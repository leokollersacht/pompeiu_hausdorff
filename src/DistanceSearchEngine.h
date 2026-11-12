#ifndef DISTANCE_SEARCH_ENGINE_H
#define DISTANCE_SEARCH_ENGINE_H

#include <Eigen/Core>
#include <igl/AABB.h>

namespace phd {
namespace search {

/// Type alias for the AABB spatial index on mesh B.
typedef igl::AABB<Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>, 3> AABBTree;

/**
 * @brief Available search engine implementations.
 */
enum class SearchEngineType {
    BestFirst,
    IDAStar
};

/**
 * @brief Search context containing problem inputs, spatial acceleration trees, and bounds.
 */
struct SearchContext {
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VA;
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FA;
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB;
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB;
    const AABBTree& treeB;
    const Eigen::VectorXd& initial_DV;
    const Eigen::VectorXi& initial_I;
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& initial_C;
    double initial_lower;
    const Eigen::VectorXd& initial_upper;
    double dA;
    double tol;
    double max_factor;

    SearchContext(
        const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VA_in,
        const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FA_in,
        const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB_in,
        const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB_in,
        const AABBTree& treeB_in,
        const Eigen::VectorXd& initial_DV_in,
        const Eigen::VectorXi& initial_I_in,
        const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& initial_C_in,
        double initial_lower_in,
        const Eigen::VectorXd& initial_upper_in,
        double dA_in,
        double tol_in,
        double max_factor_in)
        : VA(VA_in)
        , FA(FA_in)
        , VB(VB_in)
        , FB(FB_in)
        , treeB(treeB_in)
        , initial_DV(initial_DV_in)
        , initial_I(initial_I_in)
        , initial_C(initial_C_in)
        , initial_lower(initial_lower_in)
        , initial_upper(initial_upper_in)
        , dA(dA_in)
        , tol(tol_in)
        , max_factor(max_factor_in)
    {
    }
};

/**
 * @brief Abstract interface defining the contract for Hausdorff distance search routines.
 *
 * Provides a unified lifecycle for search algorithms (e.g., Best-First, IDA*),
 * isolating exploration strategies from underlying geometric bounds evaluation.
 */
class DistanceSearchEngine {
public:
    virtual ~DistanceSearchEngine() = default;

    /**
     * @brief Executes the search until the tolerance criterion is met or limits are exceeded.
     */
    virtual void search() = 0;

    /**
     * @brief Returns the computed global lower bound.
     */
    virtual double get_lower_bound() const = 0;

    /**
     * @brief Returns the computed global upper bound.
     */
    virtual double get_upper_bound() const = 0;

    /**
     * @brief Returns the total number of active vertices generated during search.
     */
    virtual int get_number_of_vertices() const = 0;

    /**
     * @brief Returns the total number of active faces generated during search.
     */
    virtual int get_number_of_faces() const = 0;

    /**
     * @brief Returns the time in milliseconds consumed by the search routine.
     */
    virtual double get_bound_time_ms() const = 0;
};

} // namespace search
} // namespace phd

#endif // DISTANCE_SEARCH_ENGINE_H
