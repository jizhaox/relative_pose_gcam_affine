#ifndef _SOLVER_AC_POSE4D_H
#define _SOLVER_AC_POSE4D_H

#include <math.h>
#include <vector>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define NEAR_ZERO_THRESHOLD 1e-16

void mod_factor_order8(double *e, double *quot);

void format_convert(
        double *input_Image_1, double *input_Image_2, double *input_affine_tran,
        double *extrinsic_R_camera, double *extrinsic_T_camera,
        std::vector<Eigen::Vector3d>& Image1, std::vector<Eigen::Vector3d>& Image2, std::vector<Eigen::Matrix3d>& Ac,
        std::vector<Eigen::Matrix3d>& R_camera, std::vector<Eigen::Vector3d>& T_camera);

void f_multicamera_Ev_solver(double *C, double *s_equation1_forqy,
        double *input_Image_1, double *input_Image_2,
        double *input_Ac,  double *extrinsic_R_camera, double *extrinsic_T_camera, char *match_type);

void solver_ac_pose4d(std::vector<Eigen::Matrix3d>& Rf1tof2_recover,
        std::vector<Eigen::Vector3d>& Tf1tof2_recover,
        double *input_Image_1, double *input_Image_2,
        double *input_Ac,  double *extrinsic_R_camera, double *extrinsic_T_camera, char *match_type);

#endif