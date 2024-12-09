# pompeiu_hausdorff
Code for "Cascading upper bounds for triangle soup Pompeiu-Hausdorff distance" (SGP 2024). \
Authors: Leonardo Sacht and Alec Jacobson. \

## C++

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

## Python

On python you can install with

```bash
pip install cascading_upper_bounds
```

Then use, for example, like this:

```python
from cascading_upper_bounds import pompeiu_hausdorff
import numpy as np
import igl

VA, FA = igl.read_triangle_mesh(f"meshes/107100.obj")
VB, FB = igl.read_triangle_mesh(f"meshes/107100_sf.obj")
tol = 1e-8
max_factor = 1000000.0
normalize = True
lower, upper_max, dA, time_taken_bvh, time_taken_bounds = pompeiu_hausdorff(VA, FA, VB, FB, tol, max_factor, normalize)

print("dA: ", dA)
print("lower bound: ", lower)
print("upper bound: ", upper_max)
print("lower bound/dA: ", lower/dA)
print("upper bound/dA: ", upper_max/dA)
print("time taken bvh(ms): ", time_taken_bvh)
print("time taken bounds(ms): ", time_taken_bounds)
```
