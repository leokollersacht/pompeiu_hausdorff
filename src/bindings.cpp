#include "PompeiuHausdorff.h"
#include "pompeiu_hausdorff.h"
#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/pair.h>


namespace nb = nanobind;
using namespace nb::literals;

NB_MODULE(cascading_upper_bounds_ext, m) {
  typedef 
    std::priority_queue< std::pair< double, int > , std::vector< std::pair< double, int >  >,
    std::less< std::pair< double, int > > >
      PQ;
    nb::class_<PQ >(m, "PriorityQueue")
      .def(nb::init<>())
      .def("pop",  &PQ::pop)
      .def("top",  &PQ::top)
      .def("empty",&PQ::empty)
      .def("size", &PQ::size)
          ;

  nb::class_<PompeiuHausdorff>(m, "PompeiuHausdorff")
      .def(nb::init<>())
      .def(nb::init<
          const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>&,
          const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>&,
          const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>&,
          const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>&,
          double, double, bool>(),
           "VA"_a, "FA"_a, "VB"_a, "FB"_a, "tol"_a=1e-8, "max_factor"_a=1000000, "normalize"_a=true)
      .def_ro("lower", &PompeiuHausdorff::lower,"Computed lower bound of the Pompeiu-Hausdorff distance")
      .def_ro("upper_max", &PompeiuHausdorff::upper_max,"Computed upper bound of the Pompeiu-Hausdorff distance")
      .def_ro("dA", &PompeiuHausdorff::dA,"Length of the diagonal of mesh A's bounding box")
      .def_ro("time_taken_bvh", &PompeiuHausdorff::time_taken_bvh,"Time taken to build the BVH for mesh B")
      .def_ro("time_taken_bounds", &PompeiuHausdorff::time_taken_bounds,"Time taken to compute the bounds")
      .def_ro("number_of_vertices", &PompeiuHausdorff::number_of_vertices,"Current number of vertices in the subdivided mesh A")
      .def_ro("number_of_faces", &PompeiuHausdorff::number_of_faces,"Current number of faces in the subdivided mesh A")
      .def_ro("VA_aug", &PompeiuHausdorff::VA_aug,"Current memory allocation for vertices (top number_of_vertices rows of VA_aug are active)")
      .def_ro("C_aug", &PompeiuHausdorff::C_aug,"Current memory allocation for vertex positions in the subdivided mesh A")
      .def_ro("DV_aug", &PompeiuHausdorff::DV_aug,"Current memory allocation for squared distances in the subdivided mesh A")
      .def_ro("I_aug", &PompeiuHausdorff::I_aug,"Current memory allocation for indices in the subdivided mesh A")
      .def_ro("FA_aug", &PompeiuHausdorff::FA_aug,"Current memory allocation for faces (top number_of_faces rows of FA_aug are active)")
      .def_ro("upper_aug", &PompeiuHausdorff::upper_aug,"#FA_aug list of per-triangle upper bounds")
      // Even though this is read only the pop method above seems to modify it
      .def_ro("Q", &PompeiuHausdorff::Q,"Queue of triangles with upper bound greater than global lower bound")
      .def_ro("lower_va", &PompeiuHausdorff::lower_va,"Index into VA_aug for the vertex determining the current lower bound")
      .def_ro("lower_fb", &PompeiuHausdorff::lower_fb,"Index into FB for the face determining the current lower bound")
      .def_ro("lower_C", &PompeiuHausdorff::lower_C,"Position of point on (VB,FB) closest to the vertex determining the current lower bound")
      ;


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
