#ifndef PHD_TELEMETRY_H
#define PHD_TELEMETRY_H

#include "geometric_bounds.h"
#include <cstdint>

namespace phd {
namespace telemetry {

struct ExecutionTelemetry {
    uint64_t subdivided_faces = 0;
    geometry::BoundCallStatistics bound_stats;
};

void record_last_execution(uint64_t subdivided_faces, const geometry::BoundCallStatistics& stats);

const ExecutionTelemetry& get_last_execution();

void reset_telemetry();

} // namespace telemetry
} // namespace phd

#endif // PHD_TELEMETRY_H
