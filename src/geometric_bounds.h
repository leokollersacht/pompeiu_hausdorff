#ifndef GEOMETRIC_BOUNDS_H
#define GEOMETRIC_BOUNDS_H

#include <Eigen/Core>
#include <cfloat>

namespace phd {
namespace geometry {

/**
 * @brief Identifies if all three vertices of a triangle project onto the same face of mesh B.
 *
 * @param projection_index_0 Face index on mesh B for vertex 0.
 * @param projection_index_1 Face index on mesh B for vertex 1.
 * @param projection_index_2 Face index on mesh B for vertex 2.
 * @return True if all three vertices project onto the identical face, false otherwise.
 */
bool has_identical_projections(
    int projection_index_0,
    int projection_index_1,
    int projection_index_2);

/**
 * @brief Computes the exact Hausdorff distance when all triangle vertices project to the same face.
 *
 * @param vertex_distances Distances from triangle vertices to mesh B (3x1 vector).
 * @return Maximum distance among the three vertices.
 */
double compute_exact_bound(const Eigen::Vector3d& vertex_distances);

/**
 * @brief Computes FirstBound (u1) based on vertex distances and opposite edge lengths.
 *
 * Formula: min_c ( d_c + max( e_{(c+1)%3}, e_{(c+2)%3} ) )
 *
 * @param vertex_distances Distances from triangle vertices to mesh B (3x1 vector).
 * @param edge_lengths Edge lengths of the triangle (3x1 vector, where e(c) is opposite vertex c).
 * @return Upper bound u1.
 */
double compute_first_bound(
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3d& edge_lengths);

/**
 * @brief Computes SecondBound (u2) based on circumradius, semiperimeter, and inradius.
 *
 * @param vertex_distances Distances from triangle vertices to mesh B (3x1 vector).
 * @param edge_lengths Edge lengths of the triangle (3x1 vector).
 * @return Upper bound u2.
 */
double compute_second_bound(
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3d& edge_lengths);

/**
 * @brief Reusable workspace buffers for Kang et al.'s continuous bound calculation (u3).
 *
 * Reusing these preallocated buffers avoids 5 dynamic heap allocations per triangle reaching u3.
 */
struct KangWorkspace {
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> VA_kang{3, 3};
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> FA_kang{1, 3};
    Eigen::VectorXd DV_kang{3};
    Eigen::VectorXi I_kang{3};
    Eigen::VectorXd u_kang{1};

    KangWorkspace() {
        FA_kang(0, 0) = 0;
        FA_kang(0, 1) = 1;
        FA_kang(0, 2) = 2;
    }
};

/**
 * @brief Statistics tracking the invocation counts of each upper bound function.
 */
struct BoundCallStatistics {
    uint64_t u0_calls = 0; // compute_exact_bound
    uint64_t u1_calls = 0; // compute_first_bound
    uint64_t u2_calls = 0; // compute_second_bound
    uint64_t u3_calls = 0; // compute_third_bound
    uint64_t u4_calls = 0; // compute_fourth_bound
};

/**
 * @brief Computes ThirdBound (u3) using Kang et al.'s continuous upper bound projection.
 *
 * @param triangle_vertices Coordinates of the 3 triangle vertices (3x3 row-major matrix).
 * @param vertex_distances Distances from triangle vertices to mesh B (3x1 vector).
 * @param projection_indices Face indices on mesh B where vertices project (3x1 vector).
 * @param VB Vertex coordinates of mesh B (#VB x 3).
 * @param FB Triangle face indices of mesh B (#FB x 3).
 * @param workspace Preallocated workspace buffer to eliminate heap allocations.
 * @return Upper bound u3.
 */
double compute_third_bound(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3i& projection_indices,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    KangWorkspace& workspace);

double compute_third_bound(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3i& projection_indices,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB);

/**
 * @brief Computes FourthBound (u4) based on the bisector plane of closest points on mesh B.
 *
 * @param triangle_vertices Coordinates of the 3 triangle vertices (3x3 row-major matrix).
 * @param closest_points Coordinates of closest points on mesh B for each vertex (3x3 row-major matrix).
 * @param edge_lengths Edge lengths of the triangle (3x1 vector).
 * @return Upper bound u4.
 */
double compute_fourth_bound(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& closest_points,
    const Eigen::Vector3d& edge_lengths);

/**
 * @brief Evaluates cascading upper bounds for a single triangle with early pruning against lower bound.
 *
 * Evaluates Exact -> FirstBound -> SecondBound -> ThirdBound -> FourthBound sequentially.
 * If any bound drops below the current global lower bound, evaluation terminates immediately.
 *
 * @param triangle_vertices Coordinates of the 3 vertices (3x3 row-major matrix).
 * @param vertex_distances Distances from vertices to mesh B (3x1 vector).
 * @param projection_indices Projected face indices on mesh B (3x1 vector).
 * @param closest_points Closest points on mesh B (3x3 row-major matrix).
 * @param VB Vertices of mesh B.
 * @param FB Faces of mesh B.
 * @param lower Current global lower bound.
 * @param[out] successful_bound_index Index of the bound that successfully pruned (0-4), or 5 if unpruned.
 * @return Calculated upper bound for the triangle.
 */
double evaluate_cascading_bound_for_triangle(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3i& projection_indices,
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& closest_points,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    double lower,
    int& successful_bound_index,
    KangWorkspace* workspace = nullptr);

double evaluate_cascading_bound_for_triangle(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3i& projection_indices,
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& closest_points,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    double lower,
    int& successful_bound_index,
    KangWorkspace* workspace,
    BoundCallStatistics* stats);

/**
 * @brief Batch evaluation of cascading upper bounds for all faces in mesh A.
 *
 * Replaces the monolithic implementation in upper_bounds while maintaining identical behavior.
 */
int evaluate_cascading_bounds(
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VA,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FA,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    const Eigen::VectorXd& DV,
    const Eigen::VectorXi& I,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& C,
    double lower,
    Eigen::VectorXd& upper_bounds_out,
    Eigen::VectorXi& success_bound_out,
    KangWorkspace* workspace = nullptr);

int evaluate_cascading_bounds(
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VA,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FA,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    const Eigen::VectorXd& DV,
    const Eigen::VectorXi& I,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& C,
    double lower,
    Eigen::VectorXd& upper_bounds_out,
    Eigen::VectorXi& success_bound_out,
    KangWorkspace* workspace,
    BoundCallStatistics* stats);

/**
 * @brief Fixed-size overload for evaluating cascading upper bounds for 4 subdivided child faces.
 *
 * Avoids any heap allocation by using fixed-size stack matrices.
 */
int evaluate_cascading_bounds(
    const Eigen::Matrix<double, 6, 3, Eigen::RowMajor>& VA,
    const Eigen::Matrix<int, 4, 3, Eigen::RowMajor>& FA,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    const Eigen::Matrix<double, 6, 1>& DV,
    const Eigen::Matrix<int, 6, 1>& I,
    const Eigen::Matrix<double, 6, 3, Eigen::RowMajor>& C,
    double lower,
    Eigen::Vector4d& upper_bounds_out,
    Eigen::Vector4i& success_bound_out,
    KangWorkspace* workspace = nullptr,
    BoundCallStatistics* stats = nullptr);

} // namespace geometry
} // namespace phd

#endif // GEOMETRIC_BOUNDS_H
