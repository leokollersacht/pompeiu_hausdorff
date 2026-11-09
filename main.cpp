
// Eigen includes
#include <Eigen/Core>

// libigl includes
#include <igl/readOBJ.h>

#include "src/pompeiu_hausdorff.h"
// time include
#if ! _MSC_VER
#include <sys/time.h>
#else
#include "src/gettimeofday.h"
#endif

#include <chrono>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 6 && argc != 7) {
        cout << "Usage: " << argv[0] << " meshA.obj meshB.obj tol max_factor normalize [engine]\n"
             << "  meshA.obj, meshB.obj: triangle soups in .obj format\n"
             << "  tol: tolerance value between upper and lower bounds\n"
             << "  max_factor: maximum allowed factor for faces and vertices expansion\n"
             << "  normalize: 0 (false) or 1 (true) to normalize tolerance by bounding box diagonal\n"
             << "  engine (optional): 'best_first' (default) or 'idastar'" << endl;
        return 0;
    }


    double t_start = 
      std::chrono::duration<double>(
          std::chrono::system_clock::now().time_since_epoch()).count();
    // load meshes (vertices and faces)
    // Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> VA, VB;
    // Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> FA, FB;
    Eigen::MatrixXd VA, VB;
    Eigen::MatrixXi FA, FB;
    if (!igl::readOBJ(argv[1],VA,FA)){
        cout << "Error loading mesh A \n" << endl;
        return 0;
    }
    if (!igl::readOBJ(argv[2],VB,FB)){
        cout << "Error loading mesh B \n" << endl;
        return 0;
    }
    if (!(FA.cols()==3)){
        cout << "Error: mesh A is not a triangle mesh \n" << endl;
        return 0;
    }
    if (!(FB.cols()==3)){
        cout << "Error: mesh B is not a triangle mesh \n" << endl;
        return 0;
    }

    double t_end = 
      std::chrono::duration<double>(
          std::chrono::system_clock::now().time_since_epoch()).count();
    double time_taken = t_end - t_start;
    // cout << "Load meshes time: " << time_taken << " secs" << endl;

    // If normalized calculations are required, calculate the length of the diagonal of the bounding box
    int normalize = atoi(argv[5]);
    double max_factor = atof(argv[4]);
    double tol = atof(argv[3]);

    phd::search::SearchEngineType engine_type = phd::search::SearchEngineType::BestFirst;
    if (argc == 7) {
        std::string engine_arg = argv[6];
        if (engine_arg == "best_first" || engine_arg == "bestfirst") {
            engine_type = phd::search::SearchEngineType::BestFirst;
        } else if (engine_arg == "idastar" || engine_arg == "ida*") {
            engine_type = phd::search::SearchEngineType::IDAStar;
        } else {
            std::cerr << "Error: Unknown search engine '" << engine_arg
                      << "'. Available engines: 'best_first', 'idastar'." << std::endl;
            return 0;
        }
    }
    
    double dA;
    double lower;
    double upper_max;
    double time_taken_bvh;
    double time_taken_bounds;
    try
    {
      std::tie(lower, upper_max, dA, time_taken_bvh, time_taken_bounds) = 
        pompeiu_hausdorff(VA, FA, VB, FB, tol, max_factor, normalize, engine_type);
    }
    catch (const std::exception& e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return 0;
    }


    // Adding final statistics print for Zheng's benchmark
    cout << fixed;
    cout << setprecision(12);
    cout << "----- Results -----" << endl;
    cout << "dA = " << dA << endl;
    cout << "lower=" << lower << endl;
    cout << "upper_max=" << upper_max << endl;
    cout << "lower/dA=" << lower/dA << endl;
    cout << "upper_max/dA=" << upper_max/dA << endl;
    cout << "bvh_time(ms)=" << time_taken_bvh << endl;
    cout << "bound_time(ms)=" << time_taken_bounds << endl;
    cout << "----------------------------------------" << endl;

    return 1;
}
