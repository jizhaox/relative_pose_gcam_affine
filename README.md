# README
This package implements the relative pose estimation for multi-camera systems from affine correspondences: 1AC plane solver, 2AC plane solver and 2AC vertical solver.

Source codes and Matlab mex files with demo code are provided in the package. The core solvers are written by C++. Matlab mex files are compiled using Ubuntu 16.04 + Matlab R2018b. Run `test_solver_AC.m` in folder "evaluation".

# Reference

[1] Banglei Guan, Ji Zhao, Daniel Barath, and Friedrich Fraundorfer. Minimal Solvers for Relative Pose Estimation of Multi-Camera Systems using Affine Correspondences. [**International Journal of Computer Vision**](https://rdcu.be/cYywK). DOI: https://doi.org/10.1007/s11263-022-01690-w 

[1] Banglei Guan, Ji Zhao, Daniel Barath, and Friedrich Fraundorfer. [**Minimal Cases for Computing the Generalized Relative Pose using Affine Correspondences**](https://openaccess.thecvf.com/content/ICCV2021/papers/Guan_Minimal_Cases_for_Computing_the_Generalized_Relative_Pose_Using_Affine_ICCV_2021_paper.pdf). International Conference on Computer Vision, 2021.

If you use this package in an academic work, please cite:

    @article{guan2023minimal,
      title={Minimal Solvers for Relative Pose Estimation of Multi-Camera Systems using Affine Correspondences},
      author={Guan, Banglei and Zhao, Ji and Barath, Daniel and Fraundorfer, Friedrich},
      journal={International Journal of Computer Vision},
      year={2023},
      doi={https://doi.org/10.1007/s11263-022-01690-w}
    }

    @inproceedings{guan2021minimal,
      title={Minimal Cases for Computing the Generalized Relative Pose using Affine Correspondences},
      author={Guan, Banglei and Zhao, Ji and Barath, Daniel and Fraundorfer, Friedrich},
      booktitle={International Conference on Computer Vision},
      year={2021},
      pages={6068-6077}
    }


# 1AC plane solver

A numerical solver for the relative pose estimation under planar motion using an affine correspondence. Returns a maximum of 4 solutions.
* **Solver**:  `solver_ac_pose3d_1AC_plane.mexa64`   

* **API**: `[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose3d_1AC_plane(Image1, Image2, At_33, R_cam, t_cam, Option);`

* **Input data for Demo**: 

     `Image1` (3\*1 matrix): normalized homogeneous image coordinates of an inter-camera feature point expressed in view 1.

     `Image2` (3\*1 matrix): normalized homogeneous image coordinates of an inter-camera feature point expressed in view 2.

     `At_33` (3\*3 matrix): `At_33(1:2,1:2) = A`, `At_33(3,3) = 0`, where `A` is the corresponding 2\*2 local affine transformations of the inter-camera feature point.

     `R_cam` (3\*3\*2 matrix): extrinsic rotation of two cameras expressed in the reference of the multi-camera system.

     `t_cam` (3\*2 matrix): extrinsic translation of two cameras expressed in the reference of the multi-camera system.

     `Option`: `inter` refers inter-camera feature point.

* **Output data for Demo**: 

     `Rf1tof2_recover` (3\*3\*N matrix): real number solutions of rotation matrix, N is the number of real number solutions.

     `Tf1tof2_recover` (3\*N matrix): real number solutions of translation.


# 2AC plane solver

A numerical solver for the relative pose estimation under planar motion using two affine correspondences. Returns a maximum of 4 solutions.
* **Solver**:  `solver_ac_pose3d_2AC_plane.mexa64`

* **API**: `[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose3d_2AC_plane(Image1, Image2, At_33, R_cam, t_cam, Option);`

* **Input data for Demo**: 

     `Image1` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in view 1.

     `Image2` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in view 2.

     `At_33` (3\*3\*2 matrix): `At_33(1:2,1:2,1) = A1`, `At_33(3,3,1) = 0`, `At_33(1:2,1:2,2) = A2`, `At_33(3,3,2) = 0`, where `A1` and `A2` are the corresponding 2\*2 local affine transformations of two feature points.

     `R_cam` (3\*3\*2 matrix): extrinsic rotation of two cameras expressed in the reference of the multi-camera system.

     `t_cam` (3\*2 matrix): extrinsic translation of two cameras expressed in the reference of the multi-camera system.

     `Option`: `inter` refers inter-camera feature points. `intra` refers intra-camera feature points.


* **Output data for Demo**: 

     `Rf1tof2_recover` (3\*3\*N matrix): real number solutions of rotation matrix, N is the number of real number solutions.

     `Tf1tof2_recover` (3\*N matrix): real number solutions of translation.


# 2AC vertical solver

A numerical solver for the relative pose estimation with known vertical direction using two affine correspondences. Returns a maximum of 6 solutions.
* **Solver**:  `solver_ac_pose4d.mexa64`

* **API**: `[Rf1tof2_recover, Tf1tof2_recover] = solver_ac_pose4d( Image1, Image2, At_33, R_cam, t_cam, Option);`

* **Input data for Demo**: 

     `Image1` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in view 1.

     `Image2` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in view 2.

     `At_33` (3\*3\*2 matrix): `At_33(1:2,1:2,1) = A1`, `At_33(3,3,1) = 0`, `At_33(1:2,1:2,2) = A2`, `At_33(3,3,2) = 0`, where `A1` and `A2` are the corresponding 2\*2 local affine transformations of two feature points.

     `R_cam` (3\*3\*2 matrix): extrinsic rotation of two cameras expressed in the reference of the multi-camera system.

     `t_cam` (3\*2 matrix): extrinsic translation of two cameras expressed in the reference of the multi-camera system.

     `Option`: `inter` refers inter-camera feature points. `intra` refers intra-camera feature points.


* **Output data for Demo**: 

     `Rf1tof2_recover` (3\*3\*N matrix): real number solutions of rotation matrix, N is the number of real number solutions.

     `Tf1tof2_recover` (3\*N matrix): real number solutions of translation.


# Compilation

**Dependency**: Eigen http://eigen.tuxfamily.org

You can compile the package in Linux like this:

    {path}$mkdir build
    {path}$cd build
    {path}/build$cmake ..
    {path}/build$make


# Matlab Interface and Evaluation

Compiled files using Ubuntu 16.04 + Matlab R2018b are provided. You can run the package in Matlab.

`evaluation/test_solver_AC.m` is the demo which shows how to call the APIs in Matlab.


# Demo Usage in Folder "test"

`test/src/test3d_inter_1AC_plane.cpp` , `test/src/test3d_inter_2AC_plane.cpp`, `test/src/test3d_intra_2AC_plane.cpp` , `test/src/test4d_inter_2AC_vertical.cpp` and `test/src/test4d_intra_2AC_vertical.cpp` are the demos which show how to call the APIs in C++. The input parameters of Demo are stored in folder "test/data".

    {path}/build$./test3d_inter_1AC_plane ../data/data3d_inter_1AC_plane.bin
    {path}/build$./test3d_inter_2AC_plane ../data/data3d_inter_2AC_plane.bin
    {path}/build$./test3d_intra_2AC_plane ../data/data3d_intra_2AC_plane.bin
    {path}/build$./test4d_inter_2AC_vertical ../data/data4d_inter_2AC_vertical.bin
    {path}/build$./test4d_intra_2AC_vertical ../data/data4d_intra_2AC_vertical.bin
