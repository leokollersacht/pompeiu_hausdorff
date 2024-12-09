#include "pompeiu_hausdorff.h"

// libigl includes
#include <igl/AABB.h>

// Pompeiu-Hausdorff distance includes
#include "upper_bounds.h"


#include <chrono>

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
    // timing variables
    double t_start, t_end;
    double time_taken;
    double time_taken_bvh;
    double time_taken_bounds;

    double dA;
    if (normalize==1){
        double x_min_A = VA.col(0).minCoeff();
        double x_max_A = VA.col(0).maxCoeff();
        double y_min_A = VA.col(1).minCoeff();
        double y_max_A = VA.col(1).maxCoeff();
        double z_min_A = VA.col(2).minCoeff();
        double z_max_A = VA.col(2).maxCoeff();
        dA = sqrt(pow(x_max_A-x_min_A,2)+pow(y_max_A-y_min_A,2)+pow(z_max_A-z_min_A,2));
    } else {
        dA = 1.0;
    }


    // Put mesh B into a libigl::AABB
    t_start = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
    igl::AABB<Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor>,3> treeB;
    treeB.init(VB,FB);
    t_end = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
    time_taken_bvh = 1000*(t_end - t_start);
    // cout << "libigl::AABB build time: " << time_taken << " secs" << endl;

    // Start timing for initializations and beginning of the loop
    t_start = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();

    // Initial distance queries
    Eigen::VectorXd DV(VA.rows());
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> C(VA.rows(),3);
    Eigen::VectorXi I(VA.rows());
    treeB.squared_distance(VB,FB,VA,DV,I,C);
    DV = DV.cwiseSqrt();
    double lower = DV.maxCoeff();

    // initial upper bounds calculation
    Eigen::VectorXi success_bound(FA.rows());
    Eigen::VectorXd upper(FA.rows());
    if (!upper_bounds(VA,FA,VB,FB,DV,I,C,lower,upper,success_bound)){
        throw std::runtime_error("error in upper bound function");
    }
    double upper_max = upper.maxCoeff();

    // Enqueue triangles with upper bound greater than global lower bound
    std::priority_queue< std::pair< double, int > , std::vector< std::pair< double, int >  >,
    std::less< std::pair< double, int > > > Q;
    for (int k=0 ; k<FA.rows(); k++){
        if (upper[k]>=lower){
            Q.emplace(upper[k],k);
        }
    }

    // Variables needed for the main loop
    int max_vertices;
    if (max_factor*VA.rows()<INT_MAX){
        max_vertices = max_factor*VA.rows();
    } else {
        max_vertices = INT_MAX;
    }
    int max_faces = max_factor*FA.rows();
    if (max_factor*FA.rows()<INT_MAX){
        max_faces = max_factor*FA.rows();
    } else {
        max_faces = INT_MAX;
    }
    int number_of_vertices = VA.rows();
    int number_of_faces = FA.rows();

    // VA_aug.rows() → current number of vertices allocated
    //   C_aug
    //   DV_aug
    //   I_aug
    // FA_aug.rows() → current number of faces allocated
    //   upper_aug
    Eigen::MatrixXd VA_aug(std::min(16*(number_of_vertices+1),max_vertices),3);
    if(VA_aug.rows() < number_of_vertices)
    {
      throw std::runtime_error("Exceeded maximum number of vertices");
    }
    VA_aug.topRows(VA.rows()) = VA;
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> C_aug(VA_aug.rows(),3);
    Eigen::VectorXd DV_aug(VA_aug.rows());
    Eigen::VectorXi I_aug(VA_aug.rows());
    C_aug.topRows(VA.rows()) = C;
    DV_aug.head(VA.rows()) = DV;
    I_aug.head(VA.rows()) = I;

    Eigen::MatrixXi FA_aug(std::min(16*(number_of_faces+1),max_faces),3);
    if(FA_aug.rows() < number_of_faces)
    {
      throw std::runtime_error("Exceeded maximum number of faces");
    }
    FA_aug.topRows(FA.rows()) = FA;
    Eigen::VectorXd upper_aug(FA_aug.rows());
    upper_aug.head(FA.rows()) = upper;

    Eigen::VectorXd upper_new(4);
    Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> FA_new(4,3);
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> VA_new(3,3);
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> VA_new_2(6,3), C_new_2(6,3);
    Eigen::VectorXd DV_new_2(6);
    Eigen::VectorXi I_new_2(6);
    int f = Q.top().second;
    int iter = 0;
    Eigen::VectorXi success_bound_new(4);

    // Loop while tolerance is not reached
    while(upper_max-lower>tol*dA){

        // throw error if the queue is empty
        if (Q.size()==0){
            cout << endl << endl << endl << "ERROR: queue got empty without reaching the given tolerance" << endl << endl << endl;
            break;
        }

        // next triangle is the one on the top of the queue
        f = Q.top().second;
        // pop triangle from the top of the queue
        Q.pop();

        // new vertices (midpoint subdivision)
        if(number_of_vertices+3 > max_vertices)
        {
          throw std::runtime_error("Exceeded maximum number of vertices");
        }
        if( (number_of_vertices+3) > VA_aug.rows())
        {
          VA_aug.conservativeResize(VA_aug.rows()*2,Eigen::NoChange);
          C_aug.conservativeResize(VA_aug.rows(),Eigen::NoChange);
          DV_aug.conservativeResize(VA_aug.rows());
          I_aug.conservativeResize(VA_aug.rows());
        }
        VA_aug.row(number_of_vertices) = VA_aug.row(FA_aug(f,0))/2+VA_aug.row(FA_aug(f,1))/2;
        VA_aug.row(number_of_vertices+1) = VA_aug.row(FA_aug(f,1))/2+VA_aug.row(FA_aug(f,2))/2;
        VA_aug.row(number_of_vertices+2) = VA_aug.row(FA_aug(f,2))/2+VA_aug.row(FA_aug(f,0))/2;

        // new faces
        if(number_of_faces+4 > max_faces)
        {
          throw std::runtime_error("Exceeded maximum number of faces");
        }
        if( (number_of_faces+4) > FA_aug.rows())
        {
          FA_aug.conservativeResize(FA_aug.rows()*2,Eigen::NoChange);
          upper_aug.conservativeResize(FA_aug.rows());
        }

        FA_aug.block(number_of_faces,0,4,3) << FA_aug(f,0), number_of_vertices, number_of_vertices+2, FA_aug(f,1), number_of_vertices+1, number_of_vertices, FA_aug(f,2), number_of_vertices+2, number_of_vertices+1, number_of_vertices, number_of_vertices+1, number_of_vertices+2;

        // update lower bound
        VA_new = VA_aug.block(number_of_vertices,0,3,3);
        treeB.squared_distance(VB,FB,VA_new,DV,I,C);
        DV = DV.cwiseSqrt();
        DV_aug.segment(number_of_vertices,3) = DV;
        I_aug.segment(number_of_vertices,3) = I;
        C_aug.block(number_of_vertices,0,3,3) = C;
        lower = fmax(DV.maxCoeff(),lower);

        // calculate new upper bounds
        FA_new = FA_aug.block(number_of_faces,0,4,3);
        VA_new_2.row(0) = VA_aug.row(FA_aug(f,0));
        VA_new_2.row(1) = VA_aug.row(FA_aug(f,1));
        VA_new_2.row(2) = VA_aug.row(FA_aug(f,2));
        VA_new_2.row(3) = VA_aug.row(number_of_vertices);
        VA_new_2.row(4) = VA_aug.row(number_of_vertices+1);
        VA_new_2.row(5) = VA_aug.row(number_of_vertices+2);
        C_new_2.row(0) = C_aug.row(FA_aug(f,0));
        C_new_2.row(1) = C_aug.row(FA_aug(f,1));
        C_new_2.row(2) = C_aug.row(FA_aug(f,2));
        C_new_2.row(3) = C_aug.row(number_of_vertices);
        C_new_2.row(4) = C_aug.row(number_of_vertices+1);
        C_new_2.row(5) = C_aug.row(number_of_vertices+2);
        FA_new(0,0) = 0; FA_new(0,1) = 3; FA_new(0,2) = 5;
        FA_new(1,0) = 3; FA_new(1,1) = 1; FA_new(1,2) = 4;
        FA_new(2,0) = 4; FA_new(2,1) = 2; FA_new(2,2) = 5;
        FA_new(3,0) = 3; FA_new(3,1) = 4; FA_new(3,2) = 5;
        DV_new_2(0) = DV_aug(FA_aug(f,0));
        DV_new_2(1) = DV_aug(FA_aug(f,1));
        DV_new_2(2) = DV_aug(FA_aug(f,2));
        DV_new_2.segment(3,3) = DV;
        I_new_2(0) = I_aug(FA_aug(f,0));
        I_new_2(1) = I_aug(FA_aug(f,1));
        I_new_2(2) = I_aug(FA_aug(f,2));
        I_new_2.segment(3,3) = I;

        if (!upper_bounds(VA_new_2,FA_new,VB,FB,DV_new_2,I_new_2,C_new_2,lower,upper_new,success_bound_new)){
          throw std::runtime_error("error in upper bound function");
        }
        upper_aug.segment(number_of_faces,4) = upper_new;
        upper_max = fmax(upper_new.maxCoeff(),Q.top().first);

        // enqueue triangles with upper bound greater than current lower bound
        for (int k=0; k<4; k++){
            if (upper_new[k]>=lower){
                Q.emplace(upper_new[k],number_of_faces+k);
            }
        }

        // Update total number of vertices and faces (even if these faces don't get into the queue)
        number_of_vertices = number_of_vertices + 3;
        number_of_faces = number_of_faces + 4;

        // throw error if number of faces or vertices exceeds the maximum
        if (number_of_faces>max_faces-4 || number_of_vertices>max_vertices-3){
            throw std::runtime_error("Exceeded maximum number of faces or vertices");
        }

        // update number of itrations
        iter++;

    }

    t_end = std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
    time_taken_bounds = 1000*(t_end - t_start);

    return std::make_tuple(lower, upper_max, dA, time_taken_bvh, time_taken_bounds);
}

