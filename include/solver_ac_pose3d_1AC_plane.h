#ifndef _SOLVER_AC_POSE3D_1AC_PLANE_H
#define _SOLVER_AC_POSE3D_1AC_PLANE_H

#include <math.h>
#include <vector>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define NEAR_ZERO_THRESHOLD 1e-16

void mod_factor_order6_1AC_plane(double *e, double *quot);

void f_multicamera_Planemotion_solver_1AC_plane(double *C, double *s_equation1_forqy,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *TiAll);

void solver_ac_pose3d_1AC_plane(std::vector<Eigen::Matrix3d>& Rf1tof2_recover,
	std::vector<Eigen::Vector3d>& Tf1tof2_recover,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *TiAll);

#endif
