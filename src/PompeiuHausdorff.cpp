#include "PompeiuHausdorff.h"
#include "BestFirstEngine.h"
#include "IDAStarEngine.h"
#include "geometric_bounds.h"
#include <chrono>
#include <cmath>
#include <stdexcept>

PompeiuHausdorff::PompeiuHausdorff(
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VA,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FA,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    const double tol,
    const double max_factor,
    const bool normalize)
    : PompeiuHausdorff(VA, FA, VB, FB, tol, max_factor, normalize, phd::search::SearchEngineType::BestFirst)
{
}

PompeiuHausdorff::PompeiuHausdorff(
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VA,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FA,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    const double tol,
    const double max_factor,
    const bool normalize,
    const phd::search::SearchEngineType engine_type)
{
    if (normalize) {
        const double x_span = VA.col(0).maxCoeff() - VA.col(0).minCoeff();
        const double y_span = VA.col(1).maxCoeff() - VA.col(1).minCoeff();
        const double z_span = VA.col(2).maxCoeff() - VA.col(2).minCoeff();
        dA = std::sqrt(x_span * x_span + y_span * y_span + z_span * z_span);
    } else {
        dA = 1.0;
    }

    // Build spatial acceleration tree for mesh B
    const double bvh_start_time = std::chrono::duration<double>(
        std::chrono::system_clock::now().time_since_epoch()).count();

    igl::AABB<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>, 3> treeB;
    treeB.init(VB, FB);

    const double bvh_end_time = std::chrono::duration<double>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    time_taken_bvh = 1000.0 * (bvh_end_time - bvh_start_time);

    // Initial distance queries for mesh A vertices onto mesh B
    Eigen::VectorXd initial_DV(VA.rows());
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> initial_C(VA.rows(), 3);
    Eigen::VectorXi initial_I(VA.rows());
    treeB.squared_distance(VB, FB, VA, initial_DV, initial_I, initial_C);
    initial_DV = initial_DV.cwiseSqrt();
    const double initial_lower = initial_DV.maxCoeff();

    // Initial cascading upper bounds calculation for all faces
    Eigen::VectorXi initial_success_bound(FA.rows());
    Eigen::VectorXd initial_upper(FA.rows());
    if (!phd::geometry::evaluate_cascading_bounds(
            VA, FA, VB, FB, initial_DV, initial_I, initial_C,
            initial_lower, initial_upper, initial_success_bound,
            nullptr, &bound_stats)) {
        throw std::runtime_error("error in upper bound function");
    }

    // Prepare search context and execute the Best-First search engine
    phd::search::SearchContext search_context(
        VA, FA, VB, FB, treeB,
        initial_DV, initial_I, initial_C,
        initial_lower, initial_upper,
        dA, tol, max_factor);

    if (engine_type == phd::search::SearchEngineType::BestFirst) {
        phd::search::BestFirstEngine engine(search_context);
        engine.search();

        lower = engine.get_lower_bound();
        upper_max = engine.get_upper_bound();
        number_of_vertices = engine.get_number_of_vertices();
        number_of_faces = engine.get_number_of_faces();
        subdivided_faces = engine.get_subdivided_faces();

        const auto& search_stats = engine.get_bound_statistics();
        bound_stats.u0_calls += search_stats.u0_calls;
        bound_stats.u1_calls += search_stats.u1_calls;
        bound_stats.u2_calls += search_stats.u2_calls;
        bound_stats.u3_calls += search_stats.u3_calls;
        bound_stats.u4_calls += search_stats.u4_calls;

        time_taken_bounds = 1000.0 * (std::chrono::duration<double>(
            std::chrono::system_clock::now().time_since_epoch()).count() - bvh_end_time);

        engine.release_state(VA_aug, C_aug, DV_aug, I_aug, FA_aug, upper_aug, Q);
    } else if (engine_type == phd::search::SearchEngineType::IDAStar) {
        phd::search::IDAStarEngine engine(search_context);
        engine.search();

        lower = engine.get_lower_bound();
        upper_max = engine.get_upper_bound();
        number_of_vertices = engine.get_number_of_vertices();
        number_of_faces = engine.get_number_of_faces();
        subdivided_faces = engine.get_subdivided_faces();

        const auto& search_stats = engine.get_bound_statistics();
        bound_stats.u0_calls += search_stats.u0_calls;
        bound_stats.u1_calls += search_stats.u1_calls;
        bound_stats.u2_calls += search_stats.u2_calls;
        bound_stats.u3_calls += search_stats.u3_calls;
        bound_stats.u4_calls += search_stats.u4_calls;

        time_taken_bounds = 1000.0 * (std::chrono::duration<double>(
            std::chrono::system_clock::now().time_since_epoch()).count() - bvh_end_time);
    } else {
        throw std::invalid_argument("Unknown search engine type requested");
    }
}
