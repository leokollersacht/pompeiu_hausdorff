#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include <Eigen/Core>

namespace phd {
namespace geometry {

/**
 * @brief Result structure for a 1-to-4 midpoint triangle subdivision.
 */
struct SubdividedTriangle {
    /// Coordinates of the 3 midpoints (m01, m12, m20).
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> midpoints;

    /// Local face connectivity referring to indices 0..2 (parent vertices) and 3..5 (midpoints).
    Eigen::Matrix<int, 4, 3, Eigen::RowMajor> local_faces;

    SubdividedTriangle();
};

/**
 * @brief Performs 1-to-4 midpoint subdivision for a single triangle.
 *
 * Given parent vertices (v0, v1, v2), computes the 3 edge midpoints:
 * - m0 = (v0 + v1) / 2
 * - m1 = (v1 + v2) / 2
 * - m2 = (v2 + v0) / 2
 *
 * @param v0 3D position of first vertex.
 * @param v1 3D position of second vertex.
 * @param v2 3D position of third vertex.
 * @return SubdividedTriangle containing midpoints and local face indexing.
 */
SubdividedTriangle subdivide_triangle_1_to_4(
    const Eigen::Vector3d& v0,
    const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2);

/**
 * @brief Constructs the 4 subfaces using global mesh vertex indexing.
 *
 * @param parent_v0 Global index of vertex 0.
 * @param parent_v1 Global index of vertex 1.
 * @param parent_v2 Global index of vertex 2.
 * @param midpoint_start_index Global index allocated for the first midpoint.
 * @return 4x3 matrix containing global vertex indices for each child face.
 */
Eigen::Matrix<int, 4, 3, Eigen::RowMajor> construct_global_subfaces(
    int parent_v0,
    int parent_v1,
    int parent_v2,
    int midpoint_start_index);

} // namespace geometry
} // namespace phd

#endif // SUBDIVISION_H
