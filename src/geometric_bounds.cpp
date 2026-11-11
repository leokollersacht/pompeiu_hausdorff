#include "geometric_bounds.h"
#include "kang_upper_bound.h"
#include "bisector_of_two_points.h"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace phd {
namespace geometry {

bool has_identical_projections(
    int projection_index_0,
    int projection_index_1,
    int projection_index_2)
{
    return (projection_index_0 == projection_index_1 && projection_index_1 == projection_index_2);
}

double compute_exact_bound(const Eigen::Vector3d& vertex_distances)
{
    return std::max(std::max(vertex_distances(0), vertex_distances(1)), vertex_distances(2));
}

double compute_first_bound(
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3d& edge_lengths)
{
    double upper_bound = DBL_MAX;
    for (int vertex_idx = 0; vertex_idx < 3; ++vertex_idx) {
        const double vertex_dist = vertex_distances(vertex_idx);
        const double max_opposite_edge = std::max(
            edge_lengths((vertex_idx + 1) % 3),
            edge_lengths((vertex_idx + 2) % 3));
        upper_bound = std::min(upper_bound, vertex_dist + max_opposite_edge);
    }
    return upper_bound;
}

double compute_second_bound(
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3d& edge_lengths)
{
    const double semiperimeter = edge_lengths.array().sum() / 2.0;
    const double area = std::sqrt(
        semiperimeter *
        (semiperimeter - edge_lengths(0)) *
        (semiperimeter - edge_lengths(1)) *
        (semiperimeter - edge_lengths(2)));
    const double circumradius = (edge_lengths(0) * edge_lengths(1) * edge_lengths(2)) / (4.0 * area);
    const double inradius = area / semiperimeter;

    const double delta = (semiperimeter - inradius > 2.0 * circumradius)
        ? circumradius
        : (edge_lengths.maxCoeff() / 2.0);

    return vertex_distances.maxCoeff() + delta;
}

double compute_third_bound(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3i& projection_indices,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    KangWorkspace& workspace)
{
    workspace.VA_kang = triangle_vertices;
    workspace.DV_kang = vertex_distances;
    workspace.I_kang = projection_indices;

    kang_upper_bound(
        workspace.VA_kang,
        workspace.FA_kang,
        VB,
        FB,
        workspace.DV_kang,
        workspace.I_kang,
        workspace.u_kang);

    return workspace.u_kang(0);
}

double compute_third_bound(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3i& projection_indices,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB)
{
    KangWorkspace local_workspace;
    return compute_third_bound(
        triangle_vertices,
        vertex_distances,
        projection_indices,
        VB,
        FB,
        local_workspace);
}

double compute_fourth_bound(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& closest_points,
    const Eigen::Vector3d& edge_lengths)
{
    int smallest_edge_index = 2;
    if (edge_lengths(0) < edge_lengths(1) && edge_lengths(0) < edge_lengths(2)) {
        smallest_edge_index = 0;
    } else if (edge_lengths(1) < edge_lengths(0) && edge_lengths(1) < edge_lengths(2)) {
        smallest_edge_index = 1;
    }

    double upper_bound = DBL_MAX;
    for (int offset = 1; offset <= 2; ++offset) {
        const Eigen::Vector3d point_q1 = closest_points.row(smallest_edge_index);
        const Eigen::Vector3d point_q2 = closest_points.row((smallest_edge_index + offset) % 3);

        double plane_a, plane_b, plane_c, plane_d;
        if (!bisector_of_two_points(point_q1, point_q2, plane_a, plane_b, plane_c, plane_d)) {
            return DBL_MAX;
        }

        const Eigen::Vector3d v0 = triangle_vertices.row(0);
        const Eigen::Vector3d v1 = triangle_vertices.row(1);
        const Eigen::Vector3d v2 = triangle_vertices.row(2);

        const double dist_v0 = std::min((v0 - point_q1).norm(), (v0 - point_q2).norm());
        const double dist_v1 = std::min((v1 - point_q1).norm(), (v1 - point_q2).norm());
        const double dist_v2 = std::min((v2 - point_q1).norm(), (v2 - point_q2).norm());
        const double vertex_max_dist = std::max(std::max(dist_v0, dist_v1), dist_v2);

        double partial_bound = vertex_max_dist;

        for (int edge_idx = 0; edge_idx <= 2; ++edge_idx) {
            const Eigen::Vector3d segment_start = triangle_vertices.row((edge_idx + 1) % 3);
            const Eigen::Vector3d segment_end = triangle_vertices.row((edge_idx + 2) % 3);

            const double denominator = plane_a * (segment_end(0) - segment_start(0)) +
                                       plane_b * (segment_end(1) - segment_start(1)) +
                                       plane_c * (segment_end(2) - segment_start(2));
            const double numerator = -(plane_a * segment_start(0) +
                                       plane_b * segment_start(1) +
                                       plane_c * segment_start(2) + plane_d);

            const double parameter_t = numerator / denominator;
            if (parameter_t > 0.0 && parameter_t < 1.0) {
                const Eigen::Vector3d intersection = segment_start + parameter_t * (segment_end - segment_start);
                const double dist_intersection = std::max(
                    (intersection - point_q1).norm(),
                    (intersection - point_q2).norm());
                partial_bound = std::max(dist_intersection, partial_bound);
            }
        }

        upper_bound = std::min(upper_bound, partial_bound);
    }

    return upper_bound;
}

double evaluate_cascading_bound_for_triangle(
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& triangle_vertices,
    const Eigen::Vector3d& vertex_distances,
    const Eigen::Vector3i& projection_indices,
    const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>& closest_points,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    double lower,
    int& successful_bound_index,
    KangWorkspace* workspace)
{
    return evaluate_cascading_bound_for_triangle(
        triangle_vertices,
        vertex_distances,
        projection_indices,
        closest_points,
        VB,
        FB,
        lower,
        successful_bound_index,
        workspace,
        nullptr);
}

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
    BoundCallStatistics* stats)
{
    if (has_identical_projections(projection_indices(0), projection_indices(1), projection_indices(2))) {
        if (stats != nullptr) {
            stats->u0_calls++;
        }
        successful_bound_index = 0;
        return compute_exact_bound(vertex_distances);
    }

    Eigen::Vector3d edge_lengths;
    for (int edge_idx = 0; edge_idx < 3; ++edge_idx) {
        edge_lengths(edge_idx) = (
            triangle_vertices.row((edge_idx + 1) % 3) -
            triangle_vertices.row((edge_idx + 2) % 3)
        ).norm();
    }

    if (stats != nullptr) {
        stats->u1_calls++;
    }
    double current_bound = compute_first_bound(vertex_distances, edge_lengths);
    if (current_bound < lower) {
        successful_bound_index = 1;
        return current_bound;
    }

    if (stats != nullptr) {
        stats->u2_calls++;
    }
    const double second_bound = compute_second_bound(vertex_distances, edge_lengths);
    current_bound = std::min(current_bound, second_bound);
    if (current_bound < lower) {
        successful_bound_index = 2;
        return current_bound;
    }

    if (stats != nullptr) {
        stats->u3_calls++;
    }
    const double third_bound = (workspace != nullptr)
        ? compute_third_bound(triangle_vertices, vertex_distances, projection_indices, VB, FB, *workspace)
        : compute_third_bound(triangle_vertices, vertex_distances, projection_indices, VB, FB);
    current_bound = std::min(current_bound, third_bound);
    if (current_bound < lower) {
        successful_bound_index = 3;
        return current_bound;
    }

    if (stats != nullptr) {
        stats->u4_calls++;
    }
    const double fourth_bound = compute_fourth_bound(
        triangle_vertices,
        closest_points,
        edge_lengths);
    current_bound = std::min(current_bound, fourth_bound);
    if (current_bound < lower) {
        successful_bound_index = 4;
        return current_bound;
    }

    successful_bound_index = 5;
    return current_bound;
}

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
    KangWorkspace* workspace)
{
    return evaluate_cascading_bounds(
        VA,
        FA,
        VB,
        FB,
        DV,
        I,
        C,
        lower,
        upper_bounds_out,
        success_bound_out,
        workspace,
        nullptr);
}

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
    BoundCallStatistics* stats)
{
    if (upper_bounds_out.rows() != FA.rows()) {
        std::cerr << "geometric_bounds: Output upper bound vector dimension mismatch" << std::endl;
        return 0;
    }

    if (workspace == nullptr) {
        KangWorkspace local_workspace;
        return evaluate_cascading_bounds(
            VA, FA, VB, FB, DV, I, C, lower,
            upper_bounds_out, success_bound_out, &local_workspace, stats);
    }

    KangWorkspace& active_workspace = *workspace;

    for (int face_idx = 0; face_idx < FA.rows(); ++face_idx) {
        const int v0_idx = FA(face_idx, 0);
        const int v1_idx = FA(face_idx, 1);
        const int v2_idx = FA(face_idx, 2);

        Eigen::Matrix<double, 3, 3, Eigen::RowMajor> triangle_vertices;
        triangle_vertices.row(0) = VA.row(v0_idx);
        triangle_vertices.row(1) = VA.row(v1_idx);
        triangle_vertices.row(2) = VA.row(v2_idx);

        Eigen::Vector3d vertex_distances(DV(v0_idx), DV(v1_idx), DV(v2_idx));
        Eigen::Vector3i projection_indices(I(v0_idx), I(v1_idx), I(v2_idx));

        Eigen::Matrix<double, 3, 3, Eigen::RowMajor> closest_points;
        closest_points.row(0) = C.row(v0_idx);
        closest_points.row(1) = C.row(v1_idx);
        closest_points.row(2) = C.row(v2_idx);

        int bound_indicator = 5;
        upper_bounds_out(face_idx) = evaluate_cascading_bound_for_triangle(
            triangle_vertices,
            vertex_distances,
            projection_indices,
            closest_points,
            VB,
            FB,
            lower,
            bound_indicator,
            &active_workspace,
            stats);
        success_bound_out(face_idx) = bound_indicator;
    }

    return 1;
}

} // namespace geometry
} // namespace phd
