- 107100.obj was obtained from https://ten-thousand-models.appspot.com/detail.html?file_id=107100

- 107100_sf.obj was obtained from https://wildmeshing.github.io/tetwild/ (Dataset --> Output --> 10k surface meshes)

- cube_a.obj: Unit cube [0, 1] x [0, 1] x [0, 1] (8 vertices, 12 triangular faces).
- cube_b.obj: Unit cube translated along X by dx = 0.5, i.e. [0.5, 1.5] x [0, 1] x [0, 1].
  Theoretical Pompeiu-Hausdorff distance between cube_a and cube_b:
  - Unnormalized (normalize=0): d_H = 0.5
  - Bounding box diagonal of A: dA = sqrt(3) ~= 1.73205080757
  - Normalized (normalize=1): d_H / dA = 0.5 / sqrt(3) ~= 0.28867513459
  - Both directions d(A, B) and d(B, A) are exactly 0.5.