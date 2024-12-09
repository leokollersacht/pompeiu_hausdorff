# from the build/ dir:
#
#    cat ../tests/test.py | python  
import pytest
from cascading_upper_bounds import pompeiu_hausdorff
import numpy as np
import igl
import pathlib

def test_example():
    this_dir = pathlib.Path(__file__).parent.resolve()
    # ./pompeiu_hausdorff ../meshes/107100.obj ../meshes/107100_sf.obj 1e-8 1000000 1 \
    print(f"this_dir: {this_dir}")
    VA, FA = igl.read_triangle_mesh(f"{this_dir}/../meshes/107100.obj")
    VB, FB = igl.read_triangle_mesh(f"{this_dir}/../meshes/107100_sf.obj")
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
    
