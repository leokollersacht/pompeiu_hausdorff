#include "pompeiu_hausdorff.h"
#include "PompeiuHausdorff.h"


#include <iostream>
std::tuple<
  double /* lower */,
  double /* upper_max */,
  double /* dA */,
  double /* time_taken_bvh */,
  double /* time_taken_bounds */>
pompeiu_hausdorff(
  const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VA,
  const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FA,
  const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB,
  const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FB,
  const double tol,
  const double max_factor,
  const bool normalize)
{
  std::cout<< "Computing Pompeiu-Hausdorff distance..." << std::endl;
  PompeiuHausdorff ph(VA, FA, VB, FB, tol, max_factor, normalize);
  std::cout<< "Done." << std::endl;
  return std::make_tuple(
    ph.lower, 
    ph.upper_max, 
    ph.dA, 
    ph.time_taken_bvh, 
    ph.time_taken_bounds);
}

