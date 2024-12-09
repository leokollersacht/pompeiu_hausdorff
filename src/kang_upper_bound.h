// Given two triangle soups (VA,FA) and (VB,FB), it calculates an upper bound for the Pompeiu-Hausdorff distance from each triangle from A to mesh B, as decribed in "2018 - Kang et al. - Fast and robust Hausdorff distance computation from triangle mesh to quad mesh in near-zero cases"

// Input:
// VA: #vertices(A) x 3 Eigen matrix containing x, y z coordinates of each vertex
// FA: #faces(A) x 3 Eigen matrix containing vertex indices of each face
// VB: #vertices(B) x 3 Eigen matrix containing x, y z coordinates of each vertex
// FB: #faces(B) x 3 Eigen matrix containing vertex indices of each face
// DV: #vertices(A) x 1 Eigen matrix containing distances from each vertex of A to B
// I: #vertices(A) x 1 Eigen vector containing indices of faces from B to which points from A are projected

// Output:
// u: #faces(A) x 1 Eigen vector containing the upper bound for the Pompeiu-Hausdorff distance from each triangle on A to mesh B

#include <stdio.h>
#include <iostream>
#include <time.h>
#include <cfloat>
#include <Eigen/Core>
#include <Eigen/Dense>
#include "kang_intersect_edge_and_bisector.h"


using namespace std;

int kang_upper_bound(const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VA, const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FA, const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB, const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FB, const Eigen::VectorXd & DV, const Eigen::VectorXi & I, Eigen::VectorXd & u);
