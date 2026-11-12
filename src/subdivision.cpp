#include "subdivision.h"

namespace phd {
namespace geometry {

SubdividedTriangle::SubdividedTriangle()
{
    local_faces << 0, 3, 5,
                   3, 1, 4,
                   4, 2, 5,
                   3, 4, 5;
}

SubdividedTriangle subdivide_triangle_1_to_4(
    const Eigen::Vector3d& v0,
    const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2)
{
    SubdividedTriangle result;
    result.midpoints.row(0) = v0 / 2.0 + v1 / 2.0;
    result.midpoints.row(1) = v1 / 2.0 + v2 / 2.0;
    result.midpoints.row(2) = v2 / 2.0 + v0 / 2.0;
    return result;
}

Eigen::Matrix<int, 4, 3, Eigen::RowMajor> construct_global_subfaces(
    int parent_v0,
    int parent_v1,
    int parent_v2,
    int midpoint_start_index)
{
    Eigen::Matrix<int, 4, 3, Eigen::RowMajor> subfaces;
    const int mid0 = midpoint_start_index;
    const int mid1 = midpoint_start_index + 1;
    const int mid2 = midpoint_start_index + 2;

    subfaces << parent_v0, mid0, mid2,
                parent_v1, mid1, mid0,
                parent_v2, mid2, mid1,
                mid0, mid1, mid2;

    return subfaces;
}

} // namespace geometry
} // namespace phd
