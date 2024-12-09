// Given two triangle soups (VA,FA) and (VB,FB), it calculates an upper bound for the Pompeiu-Hausdorff distance from each triangle from A to mesh B, as decribed in "2018 - Kang et al. - Fast and robust Hausdorff distance computation from triangle mesh to quad mesh in near-zero cases".

// Input:
// VA: #vertices(A) x 3 Eigen matrix containing x, y z coordinates of each vertex
// FA: #faces(A) x 3 Eigen matrix containing vertex indices of each face
// VB: #vertices(B) x 3 Eigen matrix containing x, y z coordinates of each vertex
// FB: #faces(B) x 3 Eigen matrix containing vertex indices of each face
// DV: #vertices(A) x 1 Eigen matrix containing distances from each vertex from mesh A to mesh B
// I: #vertices(A) x 1 Eigen vector containing indices of faces from mesh B to which points from from A are projected

// Output:
// u: #faces(A) x 1 Eigen vector containing the upper bound for the Pompeiu-Hausdorff distance from each triangle on mesh A to mesh B

#include "kang_upper_bound.h"
#include <igl/point_simplex_squared_distance.h>


int kang_upper_bound(const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VA, const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FA, const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB, const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FB, const Eigen::VectorXd & DV, const Eigen::VectorXi & I, Eigen::VectorXd & u){
    
    // Variables that are going to be used in the loop
    
    // Indices of triangles from mesh B to which vertices from mesh A are projected
    Eigen::VectorXi T_idx(3), T_idx_unique(2);
    int T1, T2, T;
    
    // Vertices from mesh B (used to define planes supporting triangles from mesh B, from which bisectors are going to be calculated)
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> VB_2(6,3);
    
    // Edge-bisector intesrections
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> P_int(2,3);
    Eigen::Vector3d P_int_edge;
    
    // Barycenter and midpoints
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> B(1,3), m1(1,3), m2(1,3);
    
    // Point-triangle squared distances
    double Query2_sqrD, Query3_sqrD;

    // Partial upper bounds (max of them will be the upper bound)
    Eigen::VectorXd ha_partial(3);
    double ha, hb, hc;
    
    // Vertices from mesh A
    Eigen::Vector3d P1, P2;
    
    // Auxiliary variables to determine if vertices from a triangle from mesh A are projected to 2 adjacent triangles on mesh B
    Eigen::VectorXi adj_vertices(9);
    int num_adj_vertices, i1, i2, i, j;
    
    // loop over all triangles from mesh A
    for (int k=0; k<FA.rows();k++){
        
        // // if the three vertices project to the same triangle, return the exact Pompeiu-Hausdorff distance from the triangle to mesh B (Section 5.3.a of the paper)
        // if (I(FA(k,0)) == I(FA(k,1)) && I(FA(k,1)) == I(FA(k,2))){
        //     u(k) = max(max(DV(FA(k,0)), DV(FA(k,1))), DV(FA(k,2)));
        // }
        // else {
            
            // Count how many common vertices the triangles from mesh B to which the vertices from mesh A were projected have
            adj_vertices.segment(0,3) = FB.row(I(FA(k,0)));
            adj_vertices.segment(3,3) = FB.row(I(FA(k,1)));
            adj_vertices.segment(6,3) = FB.row(I(FA(k,2)));
            num_adj_vertices = 9;
            for (i1=8; i1>=0; i1--){
                for (i2=i1-1; i2>=0; i2--){
                    if (adj_vertices(i2)==adj_vertices(i1)){
                        num_adj_vertices--;
                        break;
                    }
                }
            }
            
            // if the triangles have 4 unique vertices, then they consist of two adjacent triangles (section 5.3.b of the paper)
            if (num_adj_vertices==4 && ( I(FA(k,0))==I(FA(k,1)) || I(FA(k,1))==I(FA(k,2)) || I(FA(k,0))==I(FA(k,2)) )){
                
                // collect indices of the triangles
                T_idx(0) = I(FA(k,0)); T_idx(1) = I(FA(k,1)); T_idx(2) = I(FA(k,2));

                // if T_idx(0)==T_idx(1), then T_idx(0) and T_idx(2) are the 2 tirangles
                if (T_idx(0)==T_idx(1)){
                    
                    // indices of the triangles
                    T1 = T_idx(0);
                    T2 = T_idx(2);
                    
                    // Select the two vertices from mesh A (to form an edge) and two triangles from B (their bisector will be intersected with the edge)
                    P1 = VA.row(FA(k,2));
                    P2 = VA.row(FA(k,0));
                    VB_2.row(0) = VB.row(FB(T_idx(2),0));
                    VB_2.row(1) = VB.row(FB(T_idx(2),1));
                    VB_2.row(2) = VB.row(FB(T_idx(2),2));
                    VB_2.row(3) = VB.row(FB(T_idx(0),0));
                    VB_2.row(4) = VB.row(FB(T_idx(0),1));
                    VB_2.row(5) = VB.row(FB(T_idx(0),2));
                    
                    // Bisector-edge intersection
                    if (!kang_intersect_edge_and_bisector(P1, P2, VB_2, P_int_edge)){
                        cout << "Error in edge bisector intersection" << endl;
                        return 0;
                    }
                    P_int(0,0) = P_int_edge(0);
                    P_int(0,1) = P_int_edge(1);
                    P_int(0,2) = P_int_edge(2);
                    
                    // Now select the other point to form the other edge with P1
                    P2 = VA.row(FA(k,1));
                    // Bisector-edge intersection
                    if (!kang_intersect_edge_and_bisector(P1, P2, VB_2, P_int_edge)){
                        cout << "Error in edge bisector intersection" << endl;
                        return 0;
                    }
                    P_int(1,0) = P_int_edge(0);
                    P_int(1,1) = P_int_edge(1);
                    P_int(1,2) = P_int_edge(2);
                    
                // if T_idx(0)==T_idx(2), then T_idx(0) and T_idx(1) are the 2 tirangles
                } else if (T_idx(0)==T_idx(2)){
                    
                    // indices of the triangles
                    T1 = T_idx(0);
                    T2 = T_idx(1);
                    
                    // Select the two vertices from mesh A (to form an edge) and two triangles from B (their bisector will be intersected with the edge)
                    P1 = VA.row(FA(k,1));
                    P2 = VA.row(FA(k,0));
                    VB_2.row(0) = VB.row(FB(T_idx(1),0));
                    VB_2.row(1) = VB.row(FB(T_idx(1),1));
                    VB_2.row(2) = VB.row(FB(T_idx(1),2));
                    VB_2.row(3) = VB.row(FB(T_idx(0),0));
                    VB_2.row(4) = VB.row(FB(T_idx(0),1));
                    VB_2.row(5) = VB.row(FB(T_idx(0),2));
                    
                    // Bisector-edge intersection
                    if (!kang_intersect_edge_and_bisector(P1, P2, VB_2, P_int_edge)){
                        cout << "Error in edge bisector intersection" << endl;
                        return 0;
                    }
                    P_int(0,0) = P_int_edge(0);
                    P_int(0,1) = P_int_edge(1);
                    P_int(0,2) = P_int_edge(2);
                    
                    // Now select the other point to form the other edge with P1
                    P2 = VA.row(FA(k,2));
                    // Bisector-edge intersection
                    if (!kang_intersect_edge_and_bisector(P1, P2, VB_2, P_int_edge)){
                        cout << "Error in edge bisector intersection" << endl;
                        return 0;
                    }
                    P_int(1,0) = P_int_edge(0);
                    P_int(1,1) = P_int_edge(1);
                    P_int(1,2) = P_int_edge(2);
                    
                // if T_idx(1)==T_idx(2), then T_idx(0) and T_idx(1) are the 2 tirangles
                } else if (T_idx(1)==T_idx(2)){
                    
                    // indices of the triangles
                    T1 = T_idx(0);
                    T2 = T_idx(1);
                    
                    // Select the two vertices from mesh A (to form an edge) and two triangles from B (their bisector will be intersected with the edge)
                    P1 = VA.row(FA(k,0));
                    P2 = VA.row(FA(k,1));
                    VB_2.row(0) = VB.row(FB(T_idx(0),0));
                    VB_2.row(1) = VB.row(FB(T_idx(0),1));
                    VB_2.row(2) = VB.row(FB(T_idx(0),2));
                    VB_2.row(3) = VB.row(FB(T_idx(1),0));
                    VB_2.row(4) = VB.row(FB(T_idx(1),1));
                    VB_2.row(5) = VB.row(FB(T_idx(1),2));
                    
                    // Bisector-edge intersection
                    if (!kang_intersect_edge_and_bisector(P1, P2, VB_2, P_int_edge)){
                        cout << "Error in edge bisector intersection" << endl;
                        return 0;
                    }
                    P_int(0,0) = P_int_edge(0);
                    P_int(0,1) = P_int_edge(1);
                    P_int(0,2) = P_int_edge(2);
                    
                    // Now select the other point to form the other edge with P1
                    P2 = VA.row(FA(k,2));
                    // Bisector-edge intersection
                    if (!kang_intersect_edge_and_bisector(P1, P2, VB_2, P_int_edge)){
                        cout << "Error in edge bisector intersection" << endl;
                        return 0;
                    }
                    P_int(1,0) = P_int_edge(0);
                    P_int(1,1) = P_int_edge(1);
                    P_int(1,2) = P_int_edge(2);
                    
                }
                
                // vertex distances
                hb = max(max(DV(FA(k,0)), DV(FA(k,1))), DV(FA(k,2)));
                
                Eigen::RowVector3d _;

                // Using the notation in the paper: distance from b_1 to r_{11}
                igl::point_simplex_squared_distance<3>(P_int.row(0).eval(),VB,FB,T1,Query2_sqrD,_);
                hb = max(hb,sqrt(Query2_sqrD));
                    
                // Using the notation in the paper: distance from b_2 to r_{21}
                igl::point_simplex_squared_distance<3>(P_int.row(1).eval(),VB,FB,T1,Query2_sqrD,_);
                hb = max(hb,sqrt(Query2_sqrD));
                
                // Using the notation in the paper: distance from b_1 to r_{12}
                igl::point_simplex_squared_distance<3>(P_int.row(0).eval(),VB,FB,T2,Query2_sqrD,_);
                hb = max(hb,sqrt(Query2_sqrD));
                
                // Using the notation in the paper: distance from b_2 to r_{22}
                igl::point_simplex_squared_distance<3>(P_int.row(1).eval(),VB,FB,T2,Query2_sqrD,_);
                
                // Return upper bound
                u(k) = max(hb,sqrt(Query2_sqrD));

            }
            
            // projection points belong to three triangle or two non-adjacent triangles (section 5.3.c of the paper)
            else {
                
                // initialize partial upper bound
                hc = 0;
                ha_partial(0) = 0;
                ha_partial(1) = 0;
                ha_partial(2) = 0;
                ha = DBL_MAX;
                
                // barycenter of the triangle from mesh A
                B = (VA.row(FA(k,0))+VA.row(FA(k,1))+VA.row(FA(k,2)))/3;
                
                Eigen::RowVector3d _;
                // loop over its vertices
                for (i=0; i<3; i++){
                    
                    // triangle on mesh B to which this vertex from mesh a is projected
                    T = I(FA(k,i));
                    
                    // using the notation in equation (22) of the paper, this is the distance from p_i to r_i
                    hc = max(hc,DV(FA(k,i)));
                    
                    // distance from c to r_{ic}
                    igl::point_simplex_squared_distance<3>(B.row(0).eval(),VB,FB,T,Query3_sqrD,_);
                    
                    // update upper bound
                    hc = max(hc,sqrt(Query3_sqrD));
                    
                    // distance from m_i to r_{ia}
                    m1 = (VA.row(FA(k,(i)%3))+VA.row(FA(k,(i+1)%3)))/2;
                    igl::point_simplex_squared_distance<3>(m1.row(0).eval(),VB,FB,T,Query3_sqrD,_);
                    
                    // update upper bound
                    hc = max(hc,sqrt(Query3_sqrD));

                    // distance from m_{i+2} to r_{ib}
                    m2 = (VA.row(FA(k,(i)%3))+VA.row(FA(k,(i+2)%3)))/2;
                    igl::point_simplex_squared_distance<3>(m2.row(0).eval(),VB,FB,T,Query3_sqrD,_);
                    
                    // update upper bound
                    hc = max(hc,sqrt(Query3_sqrD));
                    
                    // Comppute h_a as in equation (23) of the paper
                    for (j=0; j<3; j++){
                        if (j==i){
                            ha_partial(i) = max(ha_partial(i),DV(FA(k,i)));
                        } else {
                            Eigen::RowVector3d P_query(VA(FA(k,j),0),VA(FA(k,j),1),VA(FA(k,j),2));
                            igl::point_simplex_squared_distance<3>(P_query,VB,FB,T,Query3_sqrD,_);
                            ha_partial(i) = max(ha_partial(i),sqrt(Query3_sqrD));
                        }
                    }
                }
                ha = ha_partial.minCoeff();
                
                // Equation (23) of the paper
                u(k) = min(hc,ha);

            }

        // }

    }
    
    return 1;
    
}
