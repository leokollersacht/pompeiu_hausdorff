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
    
    ph = PompeiuHausdorff(VA, FA, VB, FB, tol, max_factor, normalize)
    print("lower bound: ", ph.lower)
    print("upper bound: ", ph.upper_max)
    print("number of faces: ", ph.number_of_faces)
    print("# faces in queue: ", ph.Q.size())
    print("top of queue: ", ph.Q.top())
    print("lower_va: ", ph.lower_va)
    print("lower_fb: ", ph.lower_fb)
    print("lower_C: ", ph.lower_C)
    print("‖VA_aug(lower_va) - lower_C‖: ", np.linalg.norm(VA[ph.lower_va] - ph.lower_C))
    # 3 by 3 matrix where rows are triangle vertex positions
    vertices = ph.VA_aug[ph.FA_aug[ph.Q.top()[1]]]
    # area of the triangle using cross
    area = 0.5 * np.linalg.norm(np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0]))
    print("area of triangle in queue: ", area)
    barycenter = np.mean(vertices, axis=0).reshape(1, 3)
    print("barycenter: ", barycenter)
    D,I,C = igl.point_mesh_squared_distance(barycenter, VB, FB)
    print("‖VA_aug(lower_va) - barycenter‖: ",np.sqrt(D[0]))

