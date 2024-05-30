# pompeiu_hausdorff
Code for "Cascading upper bounds for triangle soup Pompeiu-Hausdorff distance" (SGP 2024). \
Authors: Leonardo Sacht and Alec Jacobson. \
\
Note 1: This code has been tested on a Mac OS machine \
Note 2: Command "cmake .." below will fetch all the necessary libraries (Eigen, libel, CGAL, boost). The size of these libraries is approximately 2GB. \
\
-------- Compilation --------\
\
mkdir build \
cd build \
cmake .. \
make \
\
-------- Example usage --------\
\
./pompeiu_hausdorff ../meshes/107100.obj ../meshes/107100_sf.obj 1e-8 1000000 1 