#include "IDAStarEngine.h"
#include "geometric_bounds.h"
#include "subdivision.h"
#include "telemetry.h"
#include <chrono>
#include <cmath>
#include <limits>
#include <vector>
#include <array>
#include <algorithm>
#include <climits>
#include <iostream>

namespace phd {
namespace search {

namespace {

struct TriangleNode {
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> V;
    Eigen::Vector3d DV;
    Eigen::Vector3i I;
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> C;
    double upper_bound;
};

struct DfsContext {
    const SearchContext& search_context;
    double& lower;
    double tol_scaled;
    int max_depth;
    int max_allowed_vertices;
    int max_allowed_faces;
    uint64_t& subdivided_faces;
    int& number_of_vertices;
    int& number_of_faces;
    phd::geometry::BoundCallStatistics& bound_stats;
    phd::geometry::KangWorkspace& workspace;
};

void dfs_branch(
    DfsContext& ctx,
    const TriangleNode& node,
    int depth,
    double current_threshold,
    double& next_threshold)
{
    // 1. Permanent pruning against current global lower bound
    if (node.upper_bound < ctx.lower) {
        return;
    }

    // 2. Iterative pruning: defer branches with upper bound at or below threshold
    if (node.upper_bound <= current_threshold) {
        if (node.upper_bound > next_threshold) {
            next_threshold = node.upper_bound;
        }
        return;
    }

    // 3. Stopping criterion: tolerance satisfied for this triangle
    if (node.upper_bound - ctx.lower <= ctx.tol_scaled) {
        if (node.upper_bound > next_threshold) {
            next_threshold = node.upper_bound;
        }
        return;
    }

    // 4. Safety limits on recursion depth and total elements
    if (depth >= ctx.max_depth) {
        throw std::runtime_error("Exceeded maximum recursion depth before satisfying tolerance");
    }
    if (ctx.number_of_vertices + 3 > ctx.max_allowed_vertices ||
        ctx.number_of_faces + 4 > ctx.max_allowed_faces) {
        throw std::runtime_error("Exceeded maximum number of faces or vertices");
    }

    // 5. Midpoint 1-to-4 subdivision
    const auto subdivision_result = phd::geometry::subdivide_triangle_1_to_4(
        node.V.row(0), node.V.row(1), node.V.row(2));
    ctx.subdivided_faces++;
    ctx.number_of_vertices += 3;
    ctx.number_of_faces += 4;

    // Query spatial distance on mesh B for midpoints
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> VA_mid_query = subdivision_result.midpoints;
    Eigen::VectorXd DV_mid(3);
    Eigen::VectorXi I_mid(3);
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> C_mid(3, 3);

    ctx.search_context.treeB.squared_distance(
        ctx.search_context.VB,
        ctx.search_context.FB,
        VA_mid_query,
        DV_mid,
        I_mid,
        C_mid);
    DV_mid = DV_mid.cwiseSqrt();

    // Monotonically increase global lower bound
    ctx.lower = std::max(ctx.lower, DV_mid.maxCoeff());

    // Prepare 6-point evaluation blocks for the 4 local subfaces
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> VA_eval(6, 3);
    Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> C_eval(6, 3);
    Eigen::VectorXd DV_eval(6);
    Eigen::VectorXi I_eval(6);

    VA_eval.row(0) = node.V.row(0);
    VA_eval.row(1) = node.V.row(1);
    VA_eval.row(2) = node.V.row(2);
    VA_eval.block(3, 0, 3, 3) = subdivision_result.midpoints;

    C_eval.row(0) = node.C.row(0);
    C_eval.row(1) = node.C.row(1);
    C_eval.row(2) = node.C.row(2);
    C_eval.block(3, 0, 3, 3) = C_mid;

    DV_eval(0) = node.DV(0);
    DV_eval(1) = node.DV(1);
    DV_eval(2) = node.DV(2);
    DV_eval.segment(3, 3) = DV_mid;

    I_eval(0) = node.I(0);
    I_eval(1) = node.I(1);
    I_eval(2) = node.I(2);
    I_eval.segment(3, 3) = I_mid;

    Eigen::VectorXd upper_children(4);
    Eigen::VectorXi success_children(4);

    if (!phd::geometry::evaluate_cascading_bounds(
            VA_eval,
            subdivision_result.local_faces,
            ctx.search_context.VB,
            ctx.search_context.FB,
            DV_eval,
            I_eval,
            C_eval,
            ctx.lower,
            upper_children,
            success_children,
            &ctx.workspace,
            &ctx.bound_stats)) {
        throw std::runtime_error("error in upper bound function during IDA* search");
    }

    // Sort children in descending order of upper bound to explore best branches first
    std::array<int, 4> order = {{0, 1, 2, 3}};
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return upper_children(a) > upper_children(b);
    });

    for (int c : order) {
        if (upper_children(c) >= ctx.lower) {
            TriangleNode child_node;
            for (int k = 0; k < 3; ++k) {
                const int idx = subdivision_result.local_faces(c, k);
                child_node.V.row(k) = VA_eval.row(idx);
                child_node.DV(k) = DV_eval(idx);
                child_node.I(k) = I_eval(idx);
                child_node.C.row(k) = C_eval.row(idx);
            }
            child_node.upper_bound = upper_children(c);

            dfs_branch(ctx, child_node, depth + 1, current_threshold, next_threshold);
        }
    }
}

} // anonymous namespace

IDAStarEngine::IDAStarEngine(const SearchContext& context)
    : context_(context)
    , lower_(context.initial_lower)
    , upper_max_(context.initial_upper.maxCoeff())
    , number_of_vertices_(context.VA.rows())
    , number_of_faces_(context.FA.rows())
    , bound_time_ms_(0.0)
{
}

void IDAStarEngine::search()
{
    const double start_time = std::chrono::duration<double>(
        std::chrono::system_clock::now().time_since_epoch()).count();

    uint64_t subdivided_faces = 0;
    phd::geometry::BoundCallStatistics search_bound_stats;
    phd::geometry::KangWorkspace kang_workspace;

    lower_ = context_.initial_lower;
    upper_max_ = context_.initial_upper.maxCoeff();

    const double tol_scaled = context_.tol * context_.dA;
    const int max_allowed_vertices = (context_.max_factor * context_.VA.rows() < INT_MAX)
        ? static_cast<int>(context_.max_factor * context_.VA.rows())
        : INT_MAX;
    const int max_allowed_faces = (context_.max_factor * context_.FA.rows() < INT_MAX)
        ? static_cast<int>(context_.max_factor * context_.FA.rows())
        : INT_MAX;
    const int max_depth = 100;

    auto construct_root_node = [&](int f) -> TriangleNode {
        TriangleNode node;
        for (int k = 0; k < 3; ++k) {
            const int v_idx = context_.FA(f, k);
            node.V.row(k) = context_.VA.row(v_idx);
            node.DV(k) = context_.initial_DV(v_idx);
            node.I(k) = context_.initial_I(v_idx);
            node.C.row(k) = context_.initial_C.row(v_idx);
        }
        node.upper_bound = context_.initial_upper(f);
        return node;
    };

    DfsContext ctx = {
        context_,
        lower_,
        tol_scaled,
        max_depth,
        max_allowed_vertices,
        max_allowed_faces,
        subdivided_faces,
        number_of_vertices_,
        number_of_faces_,
        search_bound_stats,
        kang_workspace
    };

    // Main Iterative Deepening loop
    while (upper_max_ - lower_ > tol_scaled) {
        const double gap = upper_max_ - lower_;
        const double threshold_step = std::max(tol_scaled, gap * 0.25);
        const double current_threshold = std::min(
            upper_max_ - threshold_step,
            std::nextafter(upper_max_, -std::numeric_limits<double>::infinity()));

        double next_threshold = -std::numeric_limits<double>::infinity();

        for (int f = 0; f < context_.FA.rows(); ++f) {
            if (context_.initial_upper(f) < lower_) {
                continue;
            }
            if (context_.initial_upper(f) <= current_threshold) {
                if (context_.initial_upper(f) > next_threshold) {
                    next_threshold = context_.initial_upper(f);
                }
                continue;
            }

            TriangleNode root_node = construct_root_node(f);
            dfs_branch(
                ctx,
                root_node,
                0,
                current_threshold,
                next_threshold);
        }

        // If no unexpanded branch was deferred by threshold, all leaves converged or were pruned
        if (next_threshold <= lower_ + tol_scaled) {
            upper_max_ = std::max(lower_, next_threshold);
            break;
        }

        // The maximum upper bound across all unexpanded leaves is at most next_threshold
        upper_max_ = next_threshold;

        if (upper_max_ - lower_ <= tol_scaled) {
            break;
        }
    }

    const double end_time = std::chrono::duration<double>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    bound_time_ms_ = 1000.0 * (end_time - start_time);

    phd::telemetry::record_search_metrics(subdivided_faces, search_bound_stats);
}

double IDAStarEngine::get_lower_bound() const
{
    return lower_;
}

double IDAStarEngine::get_upper_bound() const
{
    return upper_max_;
}

int IDAStarEngine::get_number_of_vertices() const
{
    return number_of_vertices_;
}

int IDAStarEngine::get_number_of_faces() const
{
    return number_of_faces_;
}

double IDAStarEngine::get_bound_time_ms() const
{
    return bound_time_ms_;
}

} // namespace search
} // namespace phd
