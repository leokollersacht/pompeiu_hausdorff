#include "telemetry.h"

namespace phd {
namespace telemetry {

static thread_local ExecutionTelemetry g_last_execution{};

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
