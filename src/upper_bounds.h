// Given two triangle soups (VA,FA) and (VB,FB), this function calculates cascading upper bounds for the Pompeiu-Hausdorff distance between each triangle from A to mesh B.

// Input:
// VA: #vertices(A) x 3 Eigen matrix containing x, y z coordinates of each vertex
// FA: #faces(A) x 3 Eigen matrix containing vertex indices of each face
// VB: #vertices(B) x 3 Eigen matrix containing x, y z coordinates of each vertex
// FB: #faces(B) x 3 Eigen matrix containing vertex indices of each face
// DV: #vertices(A) x 1 Eigen matrix containing distances from each vertex of A to B
// I: #vertices(A) x 1 Eigen vector containing indices of faces of B to which points from A are projected
// C: #vertices(A) x 3 Eigen matrix containing the closest points on B to the vertices of A
// lower: global lower bound (double)

// Output:
// u: #faces(A) x 1 Eigen vector containing the upper bound for the Pompeiu-Hausdorff distance from each triangle on A to mesh B
// succes_bound: #faces(A) x 1 Eigen vector containing the index of the upper bound that was successful at rejecting the triangle (= 5 if none of them were successful)

#include <stdio.h>
#include <iostream>
#include <cfloat>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/AABB.h>
#include "kang_upper_bound.h"
#include "bisector_of_two_points.h"

using namespace std;

int upper_bounds(const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VA, const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FA, const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB, const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FB, const Eigen::VectorXd & DV, const Eigen::VectorXi & I, const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & C, const double & lower, Eigen::VectorXd & u, Eigen::VectorXi & success_bound);
