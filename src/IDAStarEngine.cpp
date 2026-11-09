#include "IDAStarEngine.h"

namespace phd {
namespace search {

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
    throw std::runtime_error("IDAStarEngine: Not implemented yet (scheduled for next stage).");
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
