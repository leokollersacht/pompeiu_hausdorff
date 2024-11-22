# pompeiu_hausdorff
Code for "Cascading upper bounds for triangle soup Pompeiu-Hausdorff distance" (SGP 2024). \
Authors: Leonardo Sacht and Alec Jacobson. \
\
Note 1: This code has been tested on a Mac OS machine \
Note 2: Command "cmake .." below will fetch all the necessary libraries (Eigen, libigl). The size of these libraries is approximately 70 MB. \
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
./pompeiu_hausdorff ../meshes/107100.obj ../meshes/107100_sf.obj 1e-8 1000000 1 \
\
-------- Input ---------- \
argv[1]: path to triangle soup A in .obj format \
argv[2]: path to triangle soup B in .obj format \
argv[3]: tolerance for the difference between upper and lower bound \
argv[4]: factor to define the maximum allowed number of faces and vertices in the subdivided mesh A with respect to the number of faces and vertices of the initial mesh A \
argv[5]: 0 (false) or 1 (true) to normalize tolerance by the length of the diagonal of A's bounding box \
-------- Output (printed) ---------- \
lower bound (absolute and relative to dA), upper bound (absolute and relative to dA), and timings
