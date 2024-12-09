#include <Eigen/Core>
#include <tuple>

/// Compute lower and upper bounds on the Pompeiu-Hausdorff distance between two
/// meshes A and B
///
/// @param[in] VA  #VA by 3 list of vertex positions of mesh A 
/// @param[in] FA  #FA by 3 list of triangle indices into VA
/// @param[in] VB  #VB by 3 list of vertex positions of mesh B
/// @param[in] FB  #FB by 3 list of triangle indices into VB
/// @param[in] tol  tolerance value for the difference between upper and lower bounds
/// @param[in] max_factor  factor to define the maximum allowed number of faces and vertices in the subdivided mesh A with respect to the number of faces and vertices of the initial mesh A
/// @param[in] normalize  0 (false) or 1 (true) to normalize tolerance by the length of the diagonal of A's bounding box

std::tuple<
  double /* lower */,
  double /* upper */,
  double /* dA */,
  double /* time_taken_bvh */,
  double /* time_taken_bounds */>
pompeiu_hausdorff(
  const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VA,
  const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FA,
  const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB,
  const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FB,
  const double tol = 1e-8,
  const double max_factor = 1000000,
  const bool normalize = true);
