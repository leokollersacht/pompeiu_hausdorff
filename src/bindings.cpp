#include "pompeiu_hausdorff.h"
#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/tuple.h>

namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(cascading_upper_bounds_ext, m) {
  m.def("pompeiu_hausdorff", &pompeiu_hausdorff,
      "VA"_a, "FA"_a, "VB"_a, "FB"_a, "tol"_a=1e-8, "max_factor"_a=1000000, "normalize"_a=false,
      R"(Compute lower and upper bounds on the Pompeiu-Hausdorff distance between two
meshes A and B

@param[in] VA  #VA by 3 list of vertex positions of mesh A 
@param[in] FA  #FA by 3 list of triangle indices into VA
@param[in] VB  #VB by 3 list of vertex positions of mesh B
@param[in] FB  #FB by 3 list of triangle indices into VB
@param[in] tol  tolerance value for the difference between upper and lower bounds
@param[in] max_factor  factor to define the maximum allowed number of faces and vertices in the subdivided mesh A with respect to the number of faces and vertices of the initial mesh A
@param[in] normalize  0 (false) or 1 (true) to normalize tolerance by the length of the diagonal of A's bounding box
)");
}
