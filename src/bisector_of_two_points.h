// Given two points, this function calculates the coefficients a,b,c,d of the equation of the plane that bisects the two points.

// Input:
// Q1: x, y, z coordinates of the first point
// Q2: x, y, z coordinates of the second point

// Output:
// a,b,c,d: coefficients of the equation ax+by+cz+d=0 of the bisector

#include <stdio.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <math.h>
#include <cfloat>

using namespace std;

int bisector_of_two_points(const Eigen::Vector3d & Q1, const Eigen::Vector3d & Q2, double & a, double & b, double & c, double & d);
