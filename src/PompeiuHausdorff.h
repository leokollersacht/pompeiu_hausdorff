#include <Eigen/Core>
#include <queue>
class PompeiuHausdorff
{
  public: 
    // Member variables

    ///
    double lower;
    double upper_max;
    double dA;
    double time_taken_bvh;
    double time_taken_bounds;
    int number_of_vertices;
    int number_of_faces;
    /// Current memory allocation for vertices (top number_of_vertices rows of
    /// VA_aug are active)
    Eigen::MatrixXd VA_aug;
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> C_aug;
    Eigen::VectorXd DV_aug;
    Eigen::VectorXi I_aug;
    /// Current memory allocation for faces (top number_of_faces rows of FA_aug
    /// are active)
    Eigen::MatrixXi FA_aug;
    /// #FA_aug list of per-triangle upper bounds
    Eigen::VectorXd upper_aug;
    /// Queue of triangles with upper bound greater than global lower bound
    std::priority_queue< std::pair< double, int > , std::vector< std::pair< double, int >  >,
    std::less< std::pair< double, int > > > Q;

  // Should this be deleted?
  PompeiuHausdorff(){}
  /// @brief Class to compute the Pompeiu-Hausdorff distance between two meshes A
  /// and B
  ///
  /// @param[in] VA  #VA by 3 list of vertex positions of mesh A 
  /// @param[in] FA  #FA by 3 list of triangle indices into VA
  /// @param[in] VB  #VB by 3 list of vertex positions of mesh B
  /// @param[in] FB  #FB by 3 list of triangle indices into VB
  /// @param[in] tol  tolerance value for the difference between upper and lower bounds
  /// @param[in] max_factor  factor to define the maximum allowed number of faces and vertices in the subdivided mesh A with respect to the number of faces and vertices of the initial mesh A
  /// @param[in] normalize  0 (false) or 1 (true) to normalize tolerance by the length of the diagonal of A's bounding box
  PompeiuHausdorff(
    const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VA,
    const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FA,
    const Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> & VB,
    const Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> & FB,
    const double tol = 1e-8,
    const double max_factor = 1000000,
    const bool   normalize = true);
  // It seems this probably isn't needed after C++17
};
