// Given two points that define an edge and two triangles, this function calculates the intersection between the plane that bisects the planes that support the two triangles and the edge. If no intersection point is found, it returns the midpoint of the edge

// Input:
// P1: x, y, z coordinates of one of the endpoints of the edge
// P2: x, y, z coordinates of one of the endpoints of the edge
// VB: 6 x 3 Eigen matrix containing x, y z coordinates of each vertex of the two triangles. Fisrt three rows correspond to the first triangle and three last rows to the second triangle

// Output:
// P_int: x, y, z coordinates of the interesection point (or edge midpoint if the intersection is empty)

#include <stdio.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>
#include <cfloat>

using namespace std;

int kang_intersect_edge_and_bisector(const Eigen::Vector3d & P1, const Eigen::Vector3d & P2, const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB, Eigen::Vector3d & P_int);
