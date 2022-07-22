# README
This package implements the relative pose estimation for multi-camera systems from affine correspondences: 1AC plane solver, 2AC plane solver and 2AC solver.

# Reference

[1] Banglei Guan, Ji Zhao, Daniel Barath, and Friedrich Fraundorfer. [**Minimal Cases for Computing the Generalized Relative Pose using Affine Correspondences**](https://arxiv.org/pdf/2007.10700.pdf). International Conference on Computer Vision. 2021.

If you use this package in an academic work, please cite:

    @inproceedings{guan2021minimal,
      title={Minimal Cases for Computing the Generalized Relative Pose using Affine Correspondences},
      author={Guan, Banglei and Zhao, Ji and Barath, Daniel and Fraundorfer, Friedrich},
      booktitle={International Conference on Computer Vision},
      year={2021}
     }

# 1AC plane solver

A numerical solver for the relative pose estimation under planar motion using an affine correspondence. Returns a maximum of 4 solutions, including `Rf1tof2_recover` (3\*3 matrix) and `Tf1tof2_recover` (3\*1 matrix).
* **Source File**: `src/solver_ac_pose3d_1AC_plane.cpp`

* **API**: `void solver_ac_pose3d_1AC_plane(std::vector<Eigen::Matrix3d>& Rf1tof2_recover,
	std::vector<Eigen::Vector3d>& Tf1tof2_recover,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *TiAll);`
	
     Note that all input parameters are one-dimensional arrays. To call the API, the input data need to be stored in a one-dimensional array in column-major order.

* **Input for API**

     `Line_iAll` (6\*1 matrix): 6-dimensional vector of a Plucker line expressed in the multi-camera reference frame at time k.

     `Line_jAll` (6\*1 matrix): 6-dimensional vector of a Plucker line expressed in the multi-camera reference frame at time k+1.

     `P1All` (3\*1 matrix): normalized homogeneous image coordinates of a feature point expressed in the multi-camera reference frame at time k.
 
     `P2All` (3\*1 matrix): normalized homogeneous image coordinates of a feature point expressed in the multi-camera reference frame at time k+1.

     `AtempAll` (3\*3 matrix): `AtempAll = transpose(R2*[Ac 0; 0 0])`, where `Ac` is the corresponding 2\*2 local affine transformation of a feature point, `R2` is a rotation matrix of the other camera expressed in the multi-camera reference frame, and `transpose(M)` represents the transpose of the matrix `M`.

     `RiAll` (3\*3 matrix): `RiAll = R1`, where `R1` is a rotation matrix of one camera expressed in the multi-camera reference frame.

     `TiAll` (3\*2 matrix): translation vectors of two cameras expressed in the multi-camera reference frame.

* **Demo**: `test/test3d_1AC_plane.cpp`

# 2AC plane solver

A numerical solver for the relative pose estimation under planar motion using two affine correspondences. Returns a maximum of 4 solutions, including `Rf1tof2_recover` (3\*3 matrix) and `Tf1tof2_recover` (3\*1 matrix).
* **Source File**: `src/solver_ac_pose3d_2AC_plane.cpp`

* **API**: `void solver_ac_pose3d_2AC_plane(std::vector<Eigen::Matrix3d>& Rf1tof2_recover,
	std::vector<Eigen::Vector3d>& Tf1tof2_recover,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *TiAll);`
     
     Note that all input parameters are one-dimensional arrays. To call the API, the input data need to be stored in a one-dimensional array in column-major order.

* **Input for API**

     `Line_iAll` (6\*2 matrix): 6-dimensional vectors of two Plucker lines expressed in the multi-camera reference frame at time k.

     `Line_jAll` (6\*2 matrix): 6-dimensional vectors of two Plucker lines expressed in the multi-camera reference frame at time k+1.

     `P1All` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in the multi-camera reference frame at time k.
 
     `P2All` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in the multi-camera reference frame at time k+1.

     `AtempAll` (3\*3\*2 matrix): `AtempAll(:,:,1) = transpose(RiAll(:,:,1)*[Ac1 0; 0 0])`, `AtempAll(:,:,2) = transpose(RiAll(:,:,2)*[Ac2 0; 0 0])`, where `Ac1` and `Ac2` are the corresponding 2\*2 local affine transformations of two feature points, and `transpose(M)` represents the transpose of the matrix `M`.

     `RiAll` (3\*3\*2 matrix): rotation matrices of two cameras expressed in the multi-camera reference frame.

     `TiAll` (3\*2 matrix): translation vectors of two cameras expressed in the multi-camera reference frame.

* **Demo**: `test/test3d_2AC_plane.cpp`

# 2AC solver

A numerical solver for the relative pose estimation with known vertical direction using two affine correspondences. Returns a maximum of 6 solutions, including `Rf1tof2_recover` (3\*3 matrix) and `Tf1tof2_recover` (3\*1 matrix).
* **Source File**: `src/solver_ac_pose4d.cpp`

* **API**: `void solver_ac_pose4d(std::vector<Eigen::Matrix3d>& Rf1tof2_recover,
	std::vector<Eigen::Vector3d>& Tf1tof2_recover,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *Timu1All, double *Timu2All);`

     Note that all input parameters are one-dimensional arrays. To call the API, the input data need to be stored in a one-dimensional array in column-major order.

* **Input for API**

     `Line_iAll` (6\*2 matrix): 6-dimensional vectors of two Plucker lines expressed in the aligned multi-camera reference frame at time k.

     `Line_jAll` (6\*2 matrix): 6-dimensional vectors of two Plucker lines expressed in the aligned multi-camera reference frame at time k+1.

     `P1All` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in the aligned multi-camera reference frame at time k.
 
     `P2All` (3\*2 matrix): normalized homogeneous image coordinates of two feature points expressed in the aligned multi-camera reference frame at time k+1.

     `AtempAll` (3\*3\*2 matrix): `AtempAll(:,:,1) = transpose(RiAll2(:,:,1)*[Ac1 0; 0 0])`, `AtempAll(:,:,2) = Transpose(RiAll2(:,:,2)*[Ac2 0; 0 0])`, where `Ac1` and `Ac2` are the corresponding 2\*2 local affine transformations of two feature points, `RiAll2` (3\*3\*2 matrix) denotes rotation matrices of two cameras expressed in the aligned multi-camera reference frame at time k+1.

     `RiAll` (3\*3\*2 matrix): rotation matrices of two cameras expressed in the aligned multi-camera reference frame at time k.

     `Timu1All` (3\*2 matrix): translation vectors of two cameras expressed in the aligned multi-camera reference frame at time k.

     `Timu2All` (3\*2 matrix): translation vectors of two cameras expressed in the aligned multi-camera reference frame at time k+1.

* **Demo**: `test/test4d.cpp`

# Compilation

**Dependency**: Eigen http://eigen.tuxfamily.org

You can compile the package in Linux like this:

    {path}$mkdir build
    {path}$cd build
    {path}/build$cmake ..
    {path}/build$make
    {path}/build$./test3d_1AC_plane ../data/data3d_1AC_plane.bin
    {path}/build$./test3d_2AC_plane ../data/data3d_2AC_plane.bin
    {path}/build$./test4d ../data/data4d.bin

`test3d_1AC_plane.cpp` , `test3d_2AC_plane.cpp` and `test4d.cpp` are the demos which show how to call the APIs.

# Matlab Interface and Demo

Please run `test_matlab/test_solver_AC.m`