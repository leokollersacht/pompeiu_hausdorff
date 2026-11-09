#include "telemetry.h"

namespace phd {
namespace telemetry {

static thread_local ExecutionTelemetry g_last_execution{};

void record_initial_bounds(const geometry::BoundCallStatistics& stats)
{
    g_last_execution.bound_stats.u0_calls += stats.u0_calls;
    g_last_execution.bound_stats.u1_calls += stats.u1_calls;
    g_last_execution.bound_stats.u2_calls += stats.u2_calls;
    g_last_execution.bound_stats.u3_calls += stats.u3_calls;
    g_last_execution.bound_stats.u4_calls += stats.u4_calls;
}

void record_search_metrics(uint64_t subdivided_faces, const geometry::BoundCallStatistics& stats)
{
    g_last_execution.subdivided_faces += subdivided_faces;
    g_last_execution.bound_stats.u0_calls += stats.u0_calls;
    g_last_execution.bound_stats.u1_calls += stats.u1_calls;
    g_last_execution.bound_stats.u2_calls += stats.u2_calls;
    g_last_execution.bound_stats.u3_calls += stats.u3_calls;
    g_last_execution.bound_stats.u4_calls += stats.u4_calls;
}

void record_last_execution(uint64_t subdivided_faces, const geometry::BoundCallStatistics& stats)
{
    g_last_execution.subdivided_faces = subdivided_faces;
    g_last_execution.bound_stats = stats;
}

const ExecutionTelemetry& get_last_execution()
{
    return g_last_execution;
}

void reset_telemetry()
{
    g_last_execution = ExecutionTelemetry{};
}

} // namespace telemetry
} // namespace phd
