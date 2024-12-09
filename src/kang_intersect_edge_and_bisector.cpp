// Given two points that define an edge and two triangles, this function calculates the intersection between the plane that bisects the planes that support the two triangles and the edge. If no intersection point is found, it returns the midpoint of the edge

// Input:
// P1: x, y, z coordinates of one of the endpoints of the edge
// P2: x, y, z coordinates of one of the endpoints of the edge
// VB: 6 x 3 Eigen matrix containing x, y z coordinates of each vertex of the two triangles. Fisrt three rows correspond to the first triangle and three last rows to the second triangle

// Output:
// P_int: x, y, z coordinates of the interesection point (or edge midpoint if the intersection is empty)

#include "kang_intersect_edge_and_bisector.h"

int kang_intersect_edge_and_bisector(const Eigen::Vector3d & P1, const Eigen::Vector3d & P2, const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB, Eigen::Vector3d & P_int){
    
    Eigen::Vector3d v1, v2, cross_product, coeff_abc;
    double a1, b1, c1, d1, a2, b2, c2, d2, a_plus, b_plus, c_plus, d_plus, a_minus, b_minus, c_minus, d_minus, a_bisec, b_bisec, c_bisec, d_bisec, r1, r2, x0, y0, z0, x1, y1, z1, t;
    
    // calculate coefficients of the plane that contain the first triangle of VB
    v1 = VB.row(2)-VB.row(0);
    v2 = VB.row(1)-VB.row(0);
    cross_product = v1.cross(v2);
    a1 = cross_product(0);
    b1 = cross_product(1);
    c1 = cross_product(2);
    d1 = -cross_product.dot(VB.row(2));
    
    // calculate coefficients of the plane that contain the second triangle of VB
    v1 = VB.row(5)-VB.row(3);
    v2 = VB.row(4)-VB.row(3);
    cross_product = v1.cross(v2);
    a2 = cross_product(0);
    b2 = cross_product(1);
    c2 = cross_product(2);
    d2 = -cross_product.dot(VB.row(5));
    
    // equations for the bisector:
    r1 = sqrt(pow(a1,2) + pow(b1,2) + pow(c1,2));
    r2 = sqrt(pow(a2,2) + pow(b2,2) + pow(c2,2));
    a_plus = a1/r1-a2/r2;
    b_plus = b1/r1-b2/r2;
    c_plus = c1/r1-c2/r2;
    d_plus = d1/r1-d2/r2;
    
    // Test if the bisector separates the two points
    coeff_abc(0) = a_plus;
    coeff_abc(1) = b_plus;
    coeff_abc(2) = c_plus;
    if ((coeff_abc.dot(P1)+d_plus)*(coeff_abc.dot(P2)+d_plus)<=0){
        a_bisec = a_plus; b_bisec = b_plus; c_bisec = c_plus; d_bisec = d_plus;
    } else {
        P_int = (P1+P2)/2;
        return 1;
    }

    // Calculate intersection point
    P_int(0) = DBL_MAX; P_int(1) = DBL_MAX; P_int(2) = DBL_MAX;
    x0 = P1(0); y0 = P1(1); z0 = P1(2);
    x1 = P2(0); y1 = P2(1); z1 = P2(2);
    t = -(a_bisec*x0 + b_bisec*y0 + c_bisec*z0+d_bisec)/(a_bisec*(x1-x0)+b_bisec*(y1-y0)+c_bisec*(z1-z0));
    if (t>1e-6 && t<1-1e-6){
        P_int = P1+t*(P2-P1);
    } else {
        // no intersection point was found. Use edge midpoint
        P_int = (P1+P2)/2;
    }

    return 1;
    
}
