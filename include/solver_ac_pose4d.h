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

void f_multicamera_Ev_solver(double *C, double *s_equation1_forqy,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *Timu1All, double *Timu2All);

void solver_ac_pose4d(std::vector<Eigen::Matrix3d>& Rf1tof2_recover,
	std::vector<Eigen::Vector3d>& Tf1tof2_recover,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *Timu1All, double *Timu2All);

#endif