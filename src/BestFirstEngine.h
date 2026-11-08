#ifndef BEST_FIRST_ENGINE_H
#define BEST_FIRST_ENGINE_H

#include "DistanceSearchEngine.h"
#include <queue>
#include <vector>

namespace phd {
namespace search {

/// Priority queue ordering triangles by upper bound in descending order.
typedef std::priority_queue<
    std::pair<double, int>,
    std::vector<std::pair<double, int>>,
    std::less<std::pair<double, int>>> PriorityQueue;

/**
 * @brief Search engine implementing the original Best-First Branch-and-Bound algorithm.
 *
 * Maintains a global priority queue of active triangles sorted by their upper bounds.
 * At each step, it subdivides the triangle with the highest upper bound (1-to-4),
 * projects midpoints to mesh B, updates the global lower bound, and cascades bounds.
 */
class BestFirstEngine : public DistanceSearchEngine {
public:
    /**
     * @brief Constructs the Best-First engine using the search problem context.
     *
     * @param context Immutable problem context and initial configuration.
     */
    explicit BestFirstEngine(const SearchContext& context);

    ~BestFirstEngine() override = default;

    /**
     * @brief Executes the Best-First search until tolerance is satisfied or limits are hit.
     */
    void search() override;

    double get_lower_bound() const override;
    double get_upper_bound() const override;
    int get_number_of_vertices() const override;
    int get_number_of_faces() const override;
    double get_bound_time_ms() const override;

    // Direct accessors for backward compatibility and analysis
    const Eigen::MatrixXd& get_VA_aug() const;
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& get_C_aug() const;
    const Eigen::VectorXd& get_DV_aug() const;
    const Eigen::VectorXi& get_I_aug() const;
    const Eigen::MatrixXi& get_FA_aug() const;
    const Eigen::VectorXd& get_upper_aug() const;
    const PriorityQueue& get_queue() const;

private:
    SearchContext context_;

    double lower_;
    double upper_max_;
    int number_of_vertices_;
    int number_of_faces_;
    double bound_time_ms_;

    Eigen::MatrixXd VA_aug_;
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> C_aug_;
    Eigen::VectorXd DV_aug_;
    Eigen::VectorXi I_aug_;
    Eigen::MatrixXi FA_aug_;
    Eigen::VectorXd upper_aug_;

    PriorityQueue Q_;
};

} // namespace search
} // namespace phd

#endif // BEST_FIRST_ENGINE_H
