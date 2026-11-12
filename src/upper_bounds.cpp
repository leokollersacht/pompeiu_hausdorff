#include "upper_bounds.h"
#include "geometric_bounds.h"

int upper_bounds(
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VA,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FA,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& VB,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& FB,
    const Eigen::VectorXd& DV,
    const Eigen::VectorXi& I,
    const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>& C,
    const double& lower,
    Eigen::VectorXd& u,
    Eigen::VectorXi& success_bound)
{
    return phd::geometry::evaluate_cascading_bounds(
        VA, FA, VB, FB, DV, I, C, lower, u, success_bound);
}
