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

#include "upper_bounds.h"

int upper_bounds(const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VA, const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FA, const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB, const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FB, const Eigen::VectorXd & DV, const Eigen::VectorXi & I, const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & C, const double & lower, Eigen::VectorXd & u, Eigen::VectorXi & success_bound){
    
    if (u.rows()!=FA.rows()){
        cout << "upper_bounds.cpp: Upper bound vector has been passed with wrong number of entries (not the same as the number of triangles)" << endl;
        return 0;
    }
    
    // flags
    std::vector<bool> upper_bound_done(FA.rows());

    // Edge lengths (needed for u1 and u2)
    Eigen::VectorXd e(3);
    
    // uKang variables
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> VAKang(3,3);
    Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> FAKang(1,3);
    Eigen::VectorXd DVKang(3);
    Eigen::VectorXi IKang(3);
    Eigen::VectorXd uKang(1);
    
    for(int i = 0;i<FA.rows();i++){

        // upper bound initialization
        u(i) = DBL_MAX;
        // initialize flags as false
        upper_bound_done[i] = false;
        
        // Max vertex distance (exact HD) if the three vertices project to the same triangle.
        if (I(FA(i,0)) == I(FA(i,1)) && I(FA(i,1)) == I(FA(i,2))){
            u(i) = fmax(fmax(DV(FA(i,0)), DV(FA(i,1))), DV(FA(i,2)));
            upper_bound_done[i] = true;
            success_bound(i) = 0;
        }

        // u1 upper bound
        if (!upper_bound_done[i]) {
            double u1 = DBL_MAX;

            for(int c = 0;c<3;c++)
            {
                e(c) = (VA.row(FA(i,(c+1)%3)) - VA.row(FA(i,(c+2)%3))).norm();
            }

            for(int c = 0;c<3;c++)
            {
                const double dic = DV(FA(i,c));
                u1 = std::min(u1,dic+max(e((c+1)%3), e((c+2)%3)));
            }

            u(i) = u1;

            if (u(i)<lower){
                upper_bound_done[i] = true;
                success_bound(i) = 1;
            }

        }

        // u2 upper bound
        if (!upper_bound_done[i]) {

            double u2 = 0;
            // Semiperimeter
            double s = e.array().sum()/2.0;
            // Area
            double A = sqrt(s*(s-e(0))*(s-e(1))*(s-e(2)));
            // Circumradius
            double R = e(0)*e(1)*e(2)/(4.0*A);
            // Inradius
            double r = A/s;

            for(int c = 0;c<3;c++)
            {
                const double dic = DV(FA(i,c));
                u2 = std::max(u2,dic);
            }

            u2 += ( s-r > 2.*R ? R : e.maxCoeff()/2.0 );

            u(i) = std::min(u2,u(i));

            if (u(i)<lower){
                upper_bound_done[i] = true;
                success_bound(i) = 2;
            }

        }
            
        // u3 upper bound
        if (!upper_bound_done[i]) {
                
            VAKang.row(0) = VA.row(FA(i,0));
            VAKang.row(1) = VA.row(FA(i,1));
            VAKang.row(2) = VA.row(FA(i,2));
            FAKang(0,0) = 0; FAKang(0,1) = 1; FAKang(0,2) = 2;
            DVKang(0,0) = DV(FA(i,0)); DVKang(1,0) = DV(FA(i,1)); DVKang(2,0) = DV(FA(i,2));
            IKang(0) = I(FA(i,0)); IKang(1) = I(FA(i,1)); IKang(2) = I(FA(i,2));
                
            kang_upper_bound(VAKang, FAKang, VB, FB, DVKang, IKang, uKang);

                
            u(i) = std::min(uKang(0),u(i));
                
            if (u(i)<lower){
                upper_bound_done[i] = true;
                success_bound(i) = 3;
            }
                
        }

        // u4 upper bound
        if (!upper_bound_done[i]) {

            int sm_edge = -1;
            if ((e(0)<e(1))&&(e(0)<e(2))){
                sm_edge = 0;
            } else if ((e(1)<e(0))&&(e(1)<e(2))){
                sm_edge = 1;
            } else {
                sm_edge = 2;
            }

            if (sm_edge!=-1){

                double u4 = DBL_MAX;

                for (int a=1; a<=2; a++){

                    Eigen::Vector3d Q1 = C.row(FA(i,sm_edge));
                    Eigen::Vector3d Q2 = C.row(FA(i,(sm_edge+a)%3));

                    double a_bisec, b_bisec, c_bisec, d_bisec;
                    if (!bisector_of_two_points(Q1, Q2, a_bisec, b_bisec, c_bisec, d_bisec)){
                            cout << "Error in edge bisector intersection" << endl;
                            return 0;
                    }

                    Eigen::Vector3d v0 = VA.row(FA(i,0));
                    Eigen::Vector3d v1 = VA.row(FA(i,1));
                    Eigen::Vector3d v2 = VA.row(FA(i,2));
                    double u_vert = std::max(std::max(std::min((v0-Q1).norm(),(v0-Q2).norm()),std::min((v1-Q1).norm(),(v1-Q2).norm())),std::min((v2-Q1).norm(),(v2-Q2).norm()));
                        
                    double u_partial = u_vert;

                    for (int b=0; b<=2; b++){

                        Eigen::Vector3d P1 = VA.row(FA(i,(b+1)%3));
                        Eigen::Vector3d P2 = VA.row(FA(i,(b+2)%3));
                        double t_bisec = -(a_bisec*P1(0) + b_bisec*P1(1) + c_bisec*P1(2)+d_bisec)/(a_bisec*(P2(0)-P1(0))+b_bisec*(P2(1)-P1(1))+c_bisec*(P2(2)-P1(2)));

                        if (t_bisec>0 && t_bisec<1){
                            Eigen::Vector3d P_bisec = P1+t_bisec*(P2-P1);
                            u_partial = std::max(std::max((P_bisec-Q1).norm(),(P_bisec-Q2).norm()),u_partial);
                        }

                    }

                    u4 = std::min(u4,u_partial);

                }

                u(i) = std::min(u4,u(i));
                
                if (u(i)<lower){
                    upper_bound_done[i] = true;
                    success_bound(i) = 4; 
                }

            }

        }


        // None of the bounds were successful at rejecting the triangle
        if (!upper_bound_done[i]) {
            success_bound(i) = 5;
        }
        
        
    }
    
    return 1;
    
}
