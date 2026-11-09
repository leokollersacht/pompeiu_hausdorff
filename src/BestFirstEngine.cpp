#include "BestFirstEngine.h"
#include "geometric_bounds.h"
#include "subdivision.h"
#include <chrono>
#include <cmath>
#include <iostream>
#include <climits>

namespace phd {
namespace search {

BestFirstEngine::BestFirstEngine(const SearchContext& context)
    : context_(context)
    , lower_(context.initial_lower)
    , upper_max_(context.initial_upper.maxCoeff())
    , number_of_vertices_(context.VA.rows())
    , number_of_faces_(context.FA.rows())
    , bound_time_ms_(0.0)
{
    for (int face_idx = 0; face_idx < context_.FA.rows(); ++face_idx) {
        if (context_.initial_upper(face_idx) >= lower_) {
            Q_.emplace(context_.initial_upper(face_idx), face_idx);
        }
    }

    const int max_allowed_vertices = (context_.max_factor * context_.VA.rows() < INT_MAX)
        ? static_cast<int>(context_.max_factor * context_.VA.rows())
        : INT_MAX;
    const int max_allowed_faces = (context_.max_factor * context_.FA.rows() < INT_MAX)
        ? static_cast<int>(context_.max_factor * context_.FA.rows())
        : INT_MAX;

    VA_aug_.resize(std::min(16 * (number_of_vertices_ + 1), max_allowed_vertices), 3);
    if (VA_aug_.rows() < number_of_vertices_) {
        throw std::runtime_error("Exceeded maximum number of vertices");
    }
    VA_aug_.topRows(context_.VA.rows()) = context_.VA;

    C_aug_.resize(VA_aug_.rows(), 3);
    DV_aug_.resize(VA_aug_.rows());
    I_aug_.resize(VA_aug_.rows());

    C_aug_.topRows(context_.VA.rows()) = context_.initial_C;
    DV_aug_.head(context_.VA.rows()) = context_.initial_DV;
    I_aug_.head(context_.VA.rows()) = context_.initial_I;

    FA_aug_.resize(std::min(16 * (number_of_faces_ + 1), max_allowed_faces), 3);
    if (FA_aug_.rows() < number_of_faces_) {
        throw std::runtime_error("Exceeded maximum number of faces");
    }
    FA_aug_.topRows(context_.FA.rows()) = context_.FA;

    upper_aug_.resize(FA_aug_.rows());
    upper_aug_.head(context_.FA.rows()) = context_.initial_upper;
}

void BestFirstEngine::search()
{
    const double start_time = std::chrono::duration<double>(
        std::chrono::system_clock::now().time_since_epoch()).count();

    const int max_allowed_vertices = (context_.max_factor * context_.VA.rows() < INT_MAX)
        ? static_cast<int>(context_.max_factor * context_.VA.rows())
        : INT_MAX;
    const int max_allowed_faces = (context_.max_factor * context_.FA.rows() < INT_MAX)
        ? static_cast<int>(context_.max_factor * context_.FA.rows())
        : INT_MAX;

    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> VA_mid_query(3, 3);
    Eigen::VectorXd DV_mid(3);
    Eigen::VectorXi I_mid(3);
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> C_mid(3, 3);

    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> VA_evaluation_block(6, 3);
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> C_evaluation_block(6, 3);
    Eigen::VectorXd DV_evaluation_block(6);
    Eigen::VectorXi I_evaluation_block(6);

    Eigen::VectorXd upper_new(4);
    Eigen::VectorXi success_bound_new(4);

    while (upper_max_ - lower_ > context_.tol * context_.dA) {
        if (Q_.empty()) {
            std::cout << "\n\n\nERROR: queue got empty without reaching the given tolerance\n\n\n";
            break;
        }

        const int top_face_idx = Q_.top().second;
        Q_.pop();

        if (number_of_vertices_ + 3 > max_allowed_vertices) {
            throw std::runtime_error("Exceeded maximum number of vertices");
        }
        if (number_of_vertices_ + 3 > VA_aug_.rows()) {
            VA_aug_.conservativeResize(VA_aug_.rows() * 2, Eigen::NoChange);
            C_aug_.conservativeResize(VA_aug_.rows(), Eigen::NoChange);
            DV_aug_.conservativeResize(VA_aug_.rows());
            I_aug_.conservativeResize(VA_aug_.rows());
        }

        const Eigen::Vector3d vertex_0 = VA_aug_.row(FA_aug_(top_face_idx, 0));
        const Eigen::Vector3d vertex_1 = VA_aug_.row(FA_aug_(top_face_idx, 1));
        const Eigen::Vector3d vertex_2 = VA_aug_.row(FA_aug_(top_face_idx, 2));

        const auto subdivision_result = phd::geometry::subdivide_triangle_1_to_4(
            vertex_0, vertex_1, vertex_2);

        VA_aug_.block(number_of_vertices_, 0, 3, 3) = subdivision_result.midpoints;

        if (number_of_faces_ + 4 > max_allowed_faces) {
            throw std::runtime_error("Exceeded maximum number of faces");
        }
        if (number_of_faces_ + 4 > FA_aug_.rows()) {
            FA_aug_.conservativeResize(FA_aug_.rows() * 2, Eigen::NoChange);
            upper_aug_.conservativeResize(FA_aug_.rows());
        }

        FA_aug_.block(number_of_faces_, 0, 4, 3) = phd::geometry::construct_global_subfaces(
            FA_aug_(top_face_idx, 0),
            FA_aug_(top_face_idx, 1),
            FA_aug_(top_face_idx, 2),
            number_of_vertices_);

        VA_mid_query = VA_aug_.block(number_of_vertices_, 0, 3, 3);
        context_.treeB.squared_distance(
            context_.VB,
            context_.FB,
            VA_mid_query,
            DV_mid,
            I_mid,
            C_mid);
        DV_mid = DV_mid.cwiseSqrt();

        DV_aug_.segment(number_of_vertices_, 3) = DV_mid;
        I_aug_.segment(number_of_vertices_, 3) = I_mid;
        C_aug_.block(number_of_vertices_, 0, 3, 3) = C_mid;

        lower_ = std::max(DV_mid.maxCoeff(), lower_);

        VA_evaluation_block.row(0) = VA_aug_.row(FA_aug_(top_face_idx, 0));
        VA_evaluation_block.row(1) = VA_aug_.row(FA_aug_(top_face_idx, 1));
        VA_evaluation_block.row(2) = VA_aug_.row(FA_aug_(top_face_idx, 2));
        VA_evaluation_block.block(3, 0, 3, 3) = subdivision_result.midpoints;

        C_evaluation_block.row(0) = C_aug_.row(FA_aug_(top_face_idx, 0));
        C_evaluation_block.row(1) = C_aug_.row(FA_aug_(top_face_idx, 1));
        C_evaluation_block.row(2) = C_aug_.row(FA_aug_(top_face_idx, 2));
        C_evaluation_block.block(3, 0, 3, 3) = C_mid;

        DV_evaluation_block(0) = DV_aug_(FA_aug_(top_face_idx, 0));
        DV_evaluation_block(1) = DV_aug_(FA_aug_(top_face_idx, 1));
        DV_evaluation_block(2) = DV_aug_(FA_aug_(top_face_idx, 2));
        DV_evaluation_block.segment(3, 3) = DV_mid;

        I_evaluation_block(0) = I_aug_(FA_aug_(top_face_idx, 0));
        I_evaluation_block(1) = I_aug_(FA_aug_(top_face_idx, 1));
        I_evaluation_block(2) = I_aug_(FA_aug_(top_face_idx, 2));
        I_evaluation_block.segment(3, 3) = I_mid;

        if (!phd::geometry::evaluate_cascading_bounds(
                VA_evaluation_block,
                subdivision_result.local_faces,
                context_.VB,
                context_.FB,
                DV_evaluation_block,
                I_evaluation_block,
                C_evaluation_block,
                lower_,
                upper_new,
                success_bound_new,
                &kang_workspace_)) {
            throw std::runtime_error("error in upper bound function");
        }

        upper_aug_.segment(number_of_faces_, 4) = upper_new;
        upper_max_ = Q_.empty()
            ? upper_new.maxCoeff()
            : std::max(upper_new.maxCoeff(), Q_.top().first);

        for (int child_idx = 0; child_idx < 4; ++child_idx) {
            if (upper_new(child_idx) >= lower_) {
                Q_.emplace(upper_new(child_idx), number_of_faces_ + child_idx);
            }
        }

        number_of_vertices_ += 3;
        number_of_faces_ += 4;

        if (number_of_faces_ > max_allowed_faces - 4 || number_of_vertices_ > max_allowed_vertices - 3) {
            throw std::runtime_error("Exceeded maximum number of faces or vertices");
        }
    }

    const double end_time = std::chrono::duration<double>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    bound_time_ms_ = 1000.0 * (end_time - start_time);
}

double BestFirstEngine::get_lower_bound() const
{
    return lower_;
}

double BestFirstEngine::get_upper_bound() const
{
    return upper_max_;
}

int BestFirstEngine::get_number_of_vertices() const
{
    return number_of_vertices_;
}

int BestFirstEngine::get_number_of_faces() const
{
    return number_of_faces_;
}

double BestFirstEngine::get_bound_time_ms() const
{
    return bound_time_ms_;
}

const Eigen::MatrixXd& BestFirstEngine::get_VA_aug() const
{
    return VA_aug_;
}

const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& BestFirstEngine::get_C_aug() const
{
    return C_aug_;
}

const Eigen::VectorXd& BestFirstEngine::get_DV_aug() const
{
    return DV_aug_;
}

const Eigen::VectorXi& BestFirstEngine::get_I_aug() const
{
    return I_aug_;
}

const Eigen::MatrixXi& BestFirstEngine::get_FA_aug() const
{
    return FA_aug_;
}

const Eigen::VectorXd& BestFirstEngine::get_upper_aug() const
{
    return upper_aug_;
}

const PriorityQueue& BestFirstEngine::get_queue() const
{
    return Q_;
}

void BestFirstEngine::release_state(
    Eigen::MatrixXd& out_VA_aug,
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& out_C_aug,
    Eigen::VectorXd& out_DV_aug,
    Eigen::VectorXi& out_I_aug,
    Eigen::MatrixXi& out_FA_aug,
    Eigen::VectorXd& out_upper_aug,
    PriorityQueue& out_Q)
{
    out_VA_aug = std::move(VA_aug_);
    out_C_aug = std::move(C_aug_);
    out_DV_aug = std::move(DV_aug_);
    out_I_aug = std::move(I_aug_);
    out_FA_aug = std::move(FA_aug_);
    out_upper_aug = std::move(upper_aug_);
    out_Q = std::move(Q_);
}

} // namespace search
} // namespace phd
