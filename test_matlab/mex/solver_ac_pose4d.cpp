#include "solver_ac_pose4d.h"
#include "mex.h"


void mod_factor_order8(double *e, double *quot)
{
	quot[0] = e[0];
	quot[1] = e[1];
	quot[2] = e[2] - quot[0];
	quot[3] = e[3] - quot[1];
	quot[4] = e[4] - quot[2];
	quot[5] = e[5] - quot[3];
	quot[6] = e[6] - quot[4];
	return;
}

void f_multicamera_Ev_solver(double *C, double *s_equation1_forqy,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *Timu1All, double *Timu2All)
{
	double f1_New_C11All[2][3];
	double f1_New_C12All[2][3];
	double f1_New_C13All[2][3];
	double f1_New_C14All[2][3];

	double f2_New_C21All[2][3];
	double f2_New_C22All[2][3];
	double f2_New_C23All[2][3];
	double f2_New_C24All[2][3];

	double f3_New_C31All[2][3];
	double f3_New_C32All[2][3];
	double f3_New_C33All[2][3];
	double f3_New_C34All[2][3];

	double f_qAll[9];

	int point_num = 2;
	for (int i = 0; i < point_num; i++)
	{
		int n = 6;
		double L11 = Line_iAll[i*n];
		double L12 = Line_iAll[1 + i*n];
		double L13 = Line_iAll[2 + i*n];
		double L14 = Line_iAll[3 + i*n];
		double L15 = Line_iAll[4 + i*n];
		double L16 = Line_iAll[5 + i*n];

		double L21 = Line_jAll[i*n];
		double L22 = Line_jAll[1 + i*n];
		double L23 = Line_jAll[2 + i*n];
		double L24 = Line_jAll[3 + i*n];
		double L25 = Line_jAll[4 + i*n];
		double L26 = Line_jAll[5 + i*n];

		n = 9;
		double r1 = RiAll[i*n];
		double r2 = RiAll[3 + i*n];
		double r3 = RiAll[6 + i*n];
		double r4 = RiAll[1 + i*n];
		double r5 = RiAll[4 + i*n];
		double r6 = RiAll[7 + i*n];
		double r7 = RiAll[2 + i*n];
		double r8 = RiAll[5 + i*n];
		double r9 = RiAll[8 + i*n];

		n = 3;
		double t11 = Timu1All[i*n];
		double t12 = Timu1All[1 + i*n];
		double t13 = Timu1All[2 + i*n];

		double t21 = Timu2All[i*n];
		double t22 = Timu2All[1 + i*n];
		double t23 = Timu2All[2 + i*n];

		double p11 = P1All[i*n];
		double p12 = P1All[1 + i*n];
		double p13 = P1All[2 + i*n];

		double p21 = P2All[i*n];
		double p22 = P2All[1 + i*n];
		double p23 = P2All[2 + i*n];

		n = 9;
		double a1 = AtempAll[i*n];
		double a2 = AtempAll[3 + i*n];
		double a3 = AtempAll[6 + i*n];
		double a4 = AtempAll[1 + i*n];
		double a5 = AtempAll[4 + i*n];
		double a6 = AtempAll[7 + i*n];

		// bulid the coefficient matrix
		// f1
		f1_New_C11All[i][0] = L12*L23 + L13*L22;
		f1_New_C11All[i][1] = -2 * L11*L22;
		f1_New_C11All[i][2] = L12*L23 - L13*L22;

		f1_New_C12All[i][0] = L11*L23 - L13*L21;
		f1_New_C12All[i][1] = 2 * L11*L21 + 2 * L13*L23;
		f1_New_C12All[i][2] = L13*L21 - L11*L23;

		f1_New_C13All[i][0] = -L11*L22 - L12*L21;
		f1_New_C13All[i][1] = -2 * L13*L22;
		f1_New_C13All[i][2] = L11*L22 - L12*L21;

		f1_New_C14All[i][0] = L12*L25 - L14*L21 - L11*L24 + L15*L22 - L13*L26 - L16*L23;
		f1_New_C14All[i][1] = 2 * L11*L26 - 2 * L13*L24 + 2 * L14*L23 - 2 * L16*L21;
		f1_New_C14All[i][2] = L11*L24 + L14*L21 + L12*L25 + L15*L22 + L13*L26 + L16*L23;

		// f2
		f2_New_C21All[i][0] = -a2*p13 - a3*p12 - p23*r4 - p22*r7;
		f2_New_C21All[i][1] = 2 * a2*p11 + 2 * p22*r1;
		f2_New_C21All[i][2] = a2*p13 - a3*p12 - p23*r4 + p22*r7;

		f2_New_C22All[i][0] = a1*p13 - a3*p11 - p23*r1 + p21*r7;
		f2_New_C22All[i][1] = -2 * a1*p11 - 2 * a3*p13 - 2 * p21*r1 - 2 * p23*r7;
		f2_New_C22All[i][2] = a3*p11 - a1*p13 + p23*r1 - p21*r7;

		f2_New_C23All[i][0] = a1*p12 + a2*p11 + p22*r1 + p21*r4;
		f2_New_C23All[i][1] = 2 * a2*p13 + 2 * p22*r7;
		f2_New_C23All[i][2] = a1*p12 - a2*p11 - p22*r1 + p21*r4;

		f2_New_C24All[i][0] = a1*p13*t12 - a1*p12*t13 - a2*p11*t13 + a2*p13*t11 - a3*p11*t12 + a3*p12*t11 - a1*p12*t23 - a1*p13*t22 - a2*p11*t23 + a2*p13*t21 + a3*p11*t22 + a3*p12*t21 - p22*r1*t13 - p23*r1*t12 - p21*r4*t13 + p23*r4*t11 + p21*r7*t12 + p22*r7*t11 - p22*r1*t23 + p23*r1*t22 - p21*r4*t23 + p23*r4*t21 - p21*r7*t22 + p22*r7*t21;
		f2_New_C24All[i][1] = 2 * a1*p12*t11 - 2 * a1*p11*t12 + 2 * a3*p12*t13 - 2 * a3*p13*t12 + 2 * a1*p11*t22 - 2 * a2*p11*t21 - 2 * a2*p13*t23 + 2 * a3*p13*t22 - 2 * p21*r1*t12 + 2 * p21*r4*t11 + 2 * p23*r4*t13 - 2 * p23*r7*t12 + 2 * p21*r1*t22 - 2 * p22*r1*t21 - 2 * p22*r7*t23 + 2 * p23*r7*t22;
		f2_New_C24All[i][2] = a1*p12*t13 - a1*p13*t12 - a2*p11*t13 + a2*p13*t11 + a3*p11*t12 - a3*p12*t11 - a1*p12*t23 + a1*p13*t22 + a2*p11*t23 - a2*p13*t21 - a3*p11*t22 + a3*p12*t21 - p22*r1*t13 + p23*r1*t12 + p21*r4*t13 - p23*r4*t11 - p21*r7*t12 + p22*r7*t11 + p22*r1*t23 - p23*r1*t22 - p21*r4*t23 + p23*r4*t21 + p21*r7*t22 - p22*r7*t21;
                             
                // f3
		f3_New_C31All[i][0] = - a5*p13 - a6*p12 - p23*r5 - p22*r8;
		f3_New_C31All[i][1] = 2*a5*p11 + 2*p22*r2;
		f3_New_C31All[i][2] = a5*p13 - a6*p12 - p23*r5 + p22*r8;
        
		f3_New_C32All[i][0] = a4*p13 - a6*p11 - p23*r2 + p21*r8;
		f3_New_C32All[i][1] = - 2*a4*p11 - 2*a6*p13 - 2*p21*r2 - 2*p23*r8;
		f3_New_C32All[i][2] = a6*p11 - a4*p13 + p23*r2 - p21*r8;
        
		f3_New_C33All[i][0] = a4*p12 + a5*p11 + p22*r2 + p21*r5;
		f3_New_C33All[i][1] = 2*a5*p13 + 2*p22*r8;
		f3_New_C33All[i][2] = a4*p12 - a5*p11 - p22*r2 + p21*r5;

		f3_New_C34All[i][0] = a4*p13*t12 - a4*p12*t13 - a5*p11*t13 + a5*p13*t11 - a6*p11*t12 + a6*p12*t11 - a4*p12*t23 - a4*p13*t22 - a5*p11*t23 + a5*p13*t21 + a6*p11*t22 + a6*p12*t21 - p22*r2*t13 - p23*r2*t12 - p21*r5*t13 + p23*r5*t11 + p21*r8*t12 + p22*r8*t11 - p22*r2*t23 + p23*r2*t22 - p21*r5*t23 + p23*r5*t21 - p21*r8*t22 + p22*r8*t21;
		f3_New_C34All[i][1] = 2*a4*p12*t11 - 2*a4*p11*t12 + 2*a6*p12*t13 - 2*a6*p13*t12 + 2*a4*p11*t22 - 2*a5*p11*t21 - 2*a5*p13*t23 + 2*a6*p13*t22 - 2*p21*r2*t12 + 2*p21*r5*t11 + 2*p23*r5*t13 - 2*p23*r8*t12 + 2*p21*r2*t22 - 2*p22*r2*t21 - 2*p22*r8*t23 + 2*p23*r8*t22;
		f3_New_C34All[i][2] = a4*p12*t13 - a4*p13*t12 - a5*p11*t13 + a5*p13*t11 + a6*p11*t12 - a6*p12*t11 - a4*p12*t23 + a4*p13*t22 + a5*p11*t23 - a5*p13*t21 - a6*p11*t22 + a6*p12*t21 - p22*r2*t13 + p23*r2*t12 + p21*r5*t13 - p23*r5*t11 - p21*r8*t12 + p22*r8*t11 + p22*r2*t23 - p23*r2*t22 - p21*r5*t23 + p23*r5*t21 + p21*r8*t22 - p22*r8*t21;
	
    }
	// hidden variable resultant method
	double b2 = f1_New_C11All[0][0];
	double b1 = f1_New_C11All[0][1];
	double b0 = f1_New_C11All[0][2];

	double c2 = f1_New_C12All[0][0];
	double c1 = f1_New_C12All[0][1];
	double c0 = f1_New_C12All[0][2];

	double d2 = f1_New_C13All[0][0];
	double d1 = f1_New_C13All[0][1];
	double d0 = f1_New_C13All[0][2];

	double f2 = f1_New_C14All[0][0];
	double f1 = f1_New_C14All[0][1];
	double f0 = f1_New_C14All[0][2];

	double g2 = f2_New_C21All[0][0];
	double g1 = f2_New_C21All[0][1];
	double g0 = f2_New_C21All[0][2];

	double h2 = f2_New_C22All[0][0];
	double h1 = f2_New_C22All[0][1];
	double h0 = f2_New_C22All[0][2];

	double k2 = f2_New_C23All[0][0];
	double k1 = f2_New_C23All[0][1];
	double k0 = f2_New_C23All[0][2];

	double m2 = f2_New_C24All[0][0];
	double m1 = f2_New_C24All[0][1];
	double m0 = f2_New_C24All[0][2];

	double n2 = f3_New_C31All[0][0];
	double n1 = f3_New_C31All[0][1];
	double n0 = f3_New_C31All[0][2];

	double o2 = f3_New_C32All[0][0];
	double o1 = f3_New_C32All[0][1];
	double o0 = f3_New_C32All[0][2];

	double s2 = f3_New_C33All[0][0];
	double s1 = f3_New_C33All[0][1];
	double s0 = f3_New_C33All[0][2];

	double u2 = f3_New_C34All[0][0];
	double u1 = f3_New_C34All[0][1];
	double u0 = f3_New_C34All[0][2];

        double v2 = f1_New_C11All[1][0];
	double v1 = f1_New_C11All[1][1];
	double v0 = f1_New_C11All[1][2];

	double x2 = f1_New_C12All[1][0];
	double x1 = f1_New_C12All[1][1];
	double x0 = f1_New_C12All[1][2];

	double y2 = f1_New_C13All[1][0];
	double y1 = f1_New_C13All[1][1];
	double y0 = f1_New_C13All[1][2];

	double z2 = f1_New_C14All[1][0];
	double z1 = f1_New_C14All[1][1];
	double z0 = f1_New_C14All[1][2];

	f_qAll[0] = d2*g2*o2*z2 - d2*h2*n2*z2 - f2*g2*o2*y2 + f2*h2*n2*y2 - b2*k2*o2*z2 + c2*k2*n2*z2 - d2*m2*o2*v2 + f2*k2*o2*v2 + b2*h2*s2*z2 + b2*m2*o2*y2 - c2*g2*s2*z2 - c2*m2*n2*y2 + d2*h2*u2*v2 + d2*m2*n2*x2 - f2*h2*s2*v2 - f2*k2*n2*x2 - b2*h2*u2*y2 + c2*g2*u2*y2 - d2*g2*u2*x2 + f2*g2*s2*x2 - c2*k2*u2*v2 + c2*m2*s2*v2 + b2*k2*u2*x2 - b2*m2*s2*x2;
	f_qAll[1] = d1*g2*o2*z2 - d1*h2*n2*z2 + d2*g1*o2*z2 + d2*g2*o1*z2 + d2*g2*o2*z1 - d2*h1*n2*z2 - d2*h2*n1*z2 - d2*h2*n2*z1 - f1*g2*o2*y2 + f1*h2*n2*y2 - f2*g1*o2*y2 - f2*g2*o1*y2 - f2*g2*o2*y1 + f2*h1*n2*y2 + f2*h2*n1*y2 + f2*h2*n2*y1 - b1*k2*o2*z2 - b2*k1*o2*z2 - b2*k2*o1*z2 - b2*k2*o2*z1 + c1*k2*n2*z2 + c2*k1*n2*z2 + c2*k2*n1*z2 + c2*k2*n2*z1 - d1*m2*o2*v2 - d2*m1*o2*v2 - d2*m2*o1*v2 - d2*m2*o2*v1 + f1*k2*o2*v2 + f2*k1*o2*v2 + f2*k2*o1*v2 + f2*k2*o2*v1 + b1*h2*s2*z2 + b1*m2*o2*y2 + b2*h1*s2*z2 + b2*h2*s1*z2 + b2*h2*s2*z1 + b2*m1*o2*y2 + b2*m2*o1*y2 + b2*m2*o2*y1 - c1*g2*s2*z2 - c1*m2*n2*y2 - c2*g1*s2*z2 - c2*g2*s1*z2 - c2*g2*s2*z1 - c2*m1*n2*y2 - c2*m2*n1*y2 - c2*m2*n2*y1 + d1*h2*u2*v2 + d1*m2*n2*x2 + d2*h1*u2*v2 + d2*h2*u1*v2 + d2*h2*u2*v1 + d2*m1*n2*x2 + d2*m2*n1*x2 + d2*m2*n2*x1 - f1*h2*s2*v2 - f1*k2*n2*x2 - f2*h1*s2*v2 - f2*h2*s1*v2 - f2*h2*s2*v1 - f2*k1*n2*x2 - f2*k2*n1*x2 - f2*k2*n2*x1 - b1*h2*u2*y2 - b2*h1*u2*y2 - b2*h2*u1*y2 - b2*h2*u2*y1 + c1*g2*u2*y2 + c2*g1*u2*y2 + c2*g2*u1*y2 + c2*g2*u2*y1 - d1*g2*u2*x2 - d2*g1*u2*x2 - d2*g2*u1*x2 - d2*g2*u2*x1 + f1*g2*s2*x2 + f2*g1*s2*x2 + f2*g2*s1*x2 + f2*g2*s2*x1 - c1*k2*u2*v2 + c1*m2*s2*v2 - c2*k1*u2*v2 - c2*k2*u1*v2 - c2*k2*u2*v1 + c2*m1*s2*v2 + c2*m2*s1*v2 + c2*m2*s2*v1 + b1*k2*u2*x2 - b1*m2*s2*x2 + b2*k1*u2*x2 + b2*k2*u1*x2 + b2*k2*u2*x1 - b2*m1*s2*x2 - b2*m2*s1*x2 - b2*m2*s2*x1;
	f_qAll[2] = d0*g2*o2*z2 - d0*h2*n2*z2 + d1*g1*o2*z2 + d1*g2*o1*z2 + d1*g2*o2*z1 - d1*h1*n2*z2 - d1*h2*n1*z2 - d1*h2*n2*z1 + d2*g0*o2*z2 + d2*g1*o1*z2 + d2*g1*o2*z1 + d2*g2*o0*z2 + d2*g2*o1*z1 + d2*g2*o2*z0 - d2*h0*n2*z2 - d2*h1*n1*z2 - d2*h1*n2*z1 - d2*h2*n0*z2 - d2*h2*n1*z1 - d2*h2*n2*z0 - f0*g2*o2*y2 + f0*h2*n2*y2 - f1*g1*o2*y2 - f1*g2*o1*y2 - f1*g2*o2*y1 + f1*h1*n2*y2 + f1*h2*n1*y2 + f1*h2*n2*y1 - f2*g0*o2*y2 - f2*g1*o1*y2 - f2*g1*o2*y1 - f2*g2*o0*y2 - f2*g2*o1*y1 - f2*g2*o2*y0 + f2*h0*n2*y2 + f2*h1*n1*y2 + f2*h1*n2*y1 + f2*h2*n0*y2 + f2*h2*n1*y1 + f2*h2*n2*y0 - b0*k2*o2*z2 - b1*k1*o2*z2 - b1*k2*o1*z2 - b1*k2*o2*z1 - b2*k0*o2*z2 - b2*k1*o1*z2 - b2*k1*o2*z1 - b2*k2*o0*z2 - b2*k2*o1*z1 - b2*k2*o2*z0 + c0*k2*n2*z2 + c1*k1*n2*z2 + c1*k2*n1*z2 + c1*k2*n2*z1 + c2*k0*n2*z2 + c2*k1*n1*z2 + c2*k1*n2*z1 + c2*k2*n0*z2 + c2*k2*n1*z1 + c2*k2*n2*z0 - d0*m2*o2*v2 - d1*m1*o2*v2 - d1*m2*o1*v2 - d1*m2*o2*v1 - d2*m0*o2*v2 - d2*m1*o1*v2 - d2*m1*o2*v1 - d2*m2*o0*v2 - d2*m2*o1*v1 - d2*m2*o2*v0 + f0*k2*o2*v2 + f1*k1*o2*v2 + f1*k2*o1*v2 + f1*k2*o2*v1 + f2*k0*o2*v2 + f2*k1*o1*v2 + f2*k1*o2*v1 + f2*k2*o0*v2 + f2*k2*o1*v1 + f2*k2*o2*v0 + b0*h2*s2*z2 + b0*m2*o2*y2 + b1*h1*s2*z2 + b1*h2*s1*z2 + b1*h2*s2*z1 + b1*m1*o2*y2 + b1*m2*o1*y2 + b1*m2*o2*y1 + b2*h0*s2*z2 + b2*h1*s1*z2 + b2*h1*s2*z1 + b2*h2*s0*z2 + b2*h2*s1*z1 + b2*h2*s2*z0 + b2*m0*o2*y2 + b2*m1*o1*y2 + b2*m1*o2*y1 + b2*m2*o0*y2 + b2*m2*o1*y1 + b2*m2*o2*y0 - c0*g2*s2*z2 - c0*m2*n2*y2 - c1*g1*s2*z2 - c1*g2*s1*z2 - c1*g2*s2*z1 - c1*m1*n2*y2 - c1*m2*n1*y2 - c1*m2*n2*y1 - c2*g0*s2*z2 - c2*g1*s1*z2 - c2*g1*s2*z1 - c2*g2*s0*z2 - c2*g2*s1*z1 - c2*g2*s2*z0 - c2*m0*n2*y2 - c2*m1*n1*y2 - c2*m1*n2*y1 - c2*m2*n0*y2 - c2*m2*n1*y1 - c2*m2*n2*y0 + d0*h2*u2*v2 + d0*m2*n2*x2 + d1*h1*u2*v2 + d1*h2*u1*v2 + d1*h2*u2*v1 + d1*m1*n2*x2 + d1*m2*n1*x2 + d1*m2*n2*x1 + d2*h0*u2*v2 + d2*h1*u1*v2 + d2*h1*u2*v1 + d2*h2*u0*v2 + d2*h2*u1*v1 + d2*h2*u2*v0 + d2*m0*n2*x2 + d2*m1*n1*x2 + d2*m1*n2*x1 + d2*m2*n0*x2 + d2*m2*n1*x1 + d2*m2*n2*x0 - f0*h2*s2*v2 - f0*k2*n2*x2 - f1*h1*s2*v2 - f1*h2*s1*v2 - f1*h2*s2*v1 - f1*k1*n2*x2 - f1*k2*n1*x2 - f1*k2*n2*x1 - f2*h0*s2*v2 - f2*h1*s1*v2 - f2*h1*s2*v1 - f2*h2*s0*v2 - f2*h2*s1*v1 - f2*h2*s2*v0 - f2*k0*n2*x2 - f2*k1*n1*x2 - f2*k1*n2*x1 - f2*k2*n0*x2 - f2*k2*n1*x1 - f2*k2*n2*x0 - b0*h2*u2*y2 - b1*h1*u2*y2 - b1*h2*u1*y2 - b1*h2*u2*y1 - b2*h0*u2*y2 - b2*h1*u1*y2 - b2*h1*u2*y1 - b2*h2*u0*y2 - b2*h2*u1*y1 - b2*h2*u2*y0 + c0*g2*u2*y2 + c1*g1*u2*y2 + c1*g2*u1*y2 + c1*g2*u2*y1 + c2*g0*u2*y2 + c2*g1*u1*y2 + c2*g1*u2*y1 + c2*g2*u0*y2 + c2*g2*u1*y1 + c2*g2*u2*y0 - d0*g2*u2*x2 - d1*g1*u2*x2 - d1*g2*u1*x2 - d1*g2*u2*x1 - d2*g0*u2*x2 - d2*g1*u1*x2 - d2*g1*u2*x1 - d2*g2*u0*x2 - d2*g2*u1*x1 - d2*g2*u2*x0 + f0*g2*s2*x2 + f1*g1*s2*x2 + f1*g2*s1*x2 + f1*g2*s2*x1 + f2*g0*s2*x2 + f2*g1*s1*x2 + f2*g1*s2*x1 + f2*g2*s0*x2 + f2*g2*s1*x1 + f2*g2*s2*x0 - c0*k2*u2*v2 + c0*m2*s2*v2 - c1*k1*u2*v2 - c1*k2*u1*v2 - c1*k2*u2*v1 + c1*m1*s2*v2 + c1*m2*s1*v2 + c1*m2*s2*v1 - c2*k0*u2*v2 - c2*k1*u1*v2 - c2*k1*u2*v1 - c2*k2*u0*v2 - c2*k2*u1*v1 - c2*k2*u2*v0 + c2*m0*s2*v2 + c2*m1*s1*v2 + c2*m1*s2*v1 + c2*m2*s0*v2 + c2*m2*s1*v1 + c2*m2*s2*v0 + b0*k2*u2*x2 - b0*m2*s2*x2 + b1*k1*u2*x2 + b1*k2*u1*x2 + b1*k2*u2*x1 - b1*m1*s2*x2 - b1*m2*s1*x2 - b1*m2*s2*x1 + b2*k0*u2*x2 + b2*k1*u1*x2 + b2*k1*u2*x1 + b2*k2*u0*x2 + b2*k2*u1*x1 + b2*k2*u2*x0 - b2*m0*s2*x2 - b2*m1*s1*x2 - b2*m1*s2*x1 - b2*m2*s0*x2 - b2*m2*s1*x1 - b2*m2*s2*x0;
	f_qAll[3] = d0*g1*o2*z2 + d0*g2*o1*z2 + d0*g2*o2*z1 - d0*h1*n2*z2 - d0*h2*n1*z2 - d0*h2*n2*z1 + d1*g0*o2*z2 + d1*g1*o1*z2 + d1*g1*o2*z1 + d1*g2*o0*z2 + d1*g2*o1*z1 + d1*g2*o2*z0 - d1*h0*n2*z2 - d1*h1*n1*z2 - d1*h1*n2*z1 - d1*h2*n0*z2 - d1*h2*n1*z1 - d1*h2*n2*z0 + d2*g0*o1*z2 + d2*g0*o2*z1 + d2*g1*o0*z2 + d2*g1*o1*z1 + d2*g1*o2*z0 + d2*g2*o0*z1 + d2*g2*o1*z0 - d2*h0*n1*z2 - d2*h0*n2*z1 - d2*h1*n0*z2 - d2*h1*n1*z1 - d2*h1*n2*z0 - d2*h2*n0*z1 - d2*h2*n1*z0 - f0*g1*o2*y2 - f0*g2*o1*y2 - f0*g2*o2*y1 + f0*h1*n2*y2 + f0*h2*n1*y2 + f0*h2*n2*y1 - f1*g0*o2*y2 - f1*g1*o1*y2 - f1*g1*o2*y1 - f1*g2*o0*y2 - f1*g2*o1*y1 - f1*g2*o2*y0 + f1*h0*n2*y2 + f1*h1*n1*y2 + f1*h1*n2*y1 + f1*h2*n0*y2 + f1*h2*n1*y1 + f1*h2*n2*y0 - f2*g0*o1*y2 - f2*g0*o2*y1 - f2*g1*o0*y2 - f2*g1*o1*y1 - f2*g1*o2*y0 - f2*g2*o0*y1 - f2*g2*o1*y0 + f2*h0*n1*y2 + f2*h0*n2*y1 + f2*h1*n0*y2 + f2*h1*n1*y1 + f2*h1*n2*y0 + f2*h2*n0*y1 + f2*h2*n1*y0 - b0*k1*o2*z2 - b0*k2*o1*z2 - b0*k2*o2*z1 - b1*k0*o2*z2 - b1*k1*o1*z2 - b1*k1*o2*z1 - b1*k2*o0*z2 - b1*k2*o1*z1 - b1*k2*o2*z0 - b2*k0*o1*z2 - b2*k0*o2*z1 - b2*k1*o0*z2 - b2*k1*o1*z1 - b2*k1*o2*z0 - b2*k2*o0*z1 - b2*k2*o1*z0 + c0*k1*n2*z2 + c0*k2*n1*z2 + c0*k2*n2*z1 + c1*k0*n2*z2 + c1*k1*n1*z2 + c1*k1*n2*z1 + c1*k2*n0*z2 + c1*k2*n1*z1 + c1*k2*n2*z0 + c2*k0*n1*z2 + c2*k0*n2*z1 + c2*k1*n0*z2 + c2*k1*n1*z1 + c2*k1*n2*z0 + c2*k2*n0*z1 + c2*k2*n1*z0 - d0*m1*o2*v2 - d0*m2*o1*v2 - d0*m2*o2*v1 - d1*m0*o2*v2 - d1*m1*o1*v2 - d1*m1*o2*v1 - d1*m2*o0*v2 - d1*m2*o1*v1 - d1*m2*o2*v0 - d2*m0*o1*v2 - d2*m0*o2*v1 - d2*m1*o0*v2 - d2*m1*o1*v1 - d2*m1*o2*v0 - d2*m2*o0*v1 - d2*m2*o1*v0 + f0*k1*o2*v2 + f0*k2*o1*v2 + f0*k2*o2*v1 + f1*k0*o2*v2 + f1*k1*o1*v2 + f1*k1*o2*v1 + f1*k2*o0*v2 + f1*k2*o1*v1 + f1*k2*o2*v0 + f2*k0*o1*v2 + f2*k0*o2*v1 + f2*k1*o0*v2 + f2*k1*o1*v1 + f2*k1*o2*v0 + f2*k2*o0*v1 + f2*k2*o1*v0 + b0*h1*s2*z2 + b0*h2*s1*z2 + b0*h2*s2*z1 + b0*m1*o2*y2 + b0*m2*o1*y2 + b0*m2*o2*y1 + b1*h0*s2*z2 + b1*h1*s1*z2 + b1*h1*s2*z1 + b1*h2*s0*z2 + b1*h2*s1*z1 + b1*h2*s2*z0 + b1*m0*o2*y2 + b1*m1*o1*y2 + b1*m1*o2*y1 + b1*m2*o0*y2 + b1*m2*o1*y1 + b1*m2*o2*y0 + b2*h0*s1*z2 + b2*h0*s2*z1 + b2*h1*s0*z2 + b2*h1*s1*z1 + b2*h1*s2*z0 + b2*h2*s0*z1 + b2*h2*s1*z0 + b2*m0*o1*y2 + b2*m0*o2*y1 + b2*m1*o0*y2 + b2*m1*o1*y1 + b2*m1*o2*y0 + b2*m2*o0*y1 + b2*m2*o1*y0 - c0*g1*s2*z2 - c0*g2*s1*z2 - c0*g2*s2*z1 - c0*m1*n2*y2 - c0*m2*n1*y2 - c0*m2*n2*y1 - c1*g0*s2*z2 - c1*g1*s1*z2 - c1*g1*s2*z1 - c1*g2*s0*z2 - c1*g2*s1*z1 - c1*g2*s2*z0 - c1*m0*n2*y2 - c1*m1*n1*y2 - c1*m1*n2*y1 - c1*m2*n0*y2 - c1*m2*n1*y1 - c1*m2*n2*y0 - c2*g0*s1*z2 - c2*g0*s2*z1 - c2*g1*s0*z2 - c2*g1*s1*z1 - c2*g1*s2*z0 - c2*g2*s0*z1 - c2*g2*s1*z0 - c2*m0*n1*y2 - c2*m0*n2*y1 - c2*m1*n0*y2 - c2*m1*n1*y1 - c2*m1*n2*y0 - c2*m2*n0*y1 - c2*m2*n1*y0 + d0*h1*u2*v2 + d0*h2*u1*v2 + d0*h2*u2*v1 + d0*m1*n2*x2 + d0*m2*n1*x2 + d0*m2*n2*x1 + d1*h0*u2*v2 + d1*h1*u1*v2 + d1*h1*u2*v1 + d1*h2*u0*v2 + d1*h2*u1*v1 + d1*h2*u2*v0 + d1*m0*n2*x2 + d1*m1*n1*x2 + d1*m1*n2*x1 + d1*m2*n0*x2 + d1*m2*n1*x1 + d1*m2*n2*x0 + d2*h0*u1*v2 + d2*h0*u2*v1 + d2*h1*u0*v2 + d2*h1*u1*v1 + d2*h1*u2*v0 + d2*h2*u0*v1 + d2*h2*u1*v0 + d2*m0*n1*x2 + d2*m0*n2*x1 + d2*m1*n0*x2 + d2*m1*n1*x1 + d2*m1*n2*x0 + d2*m2*n0*x1 + d2*m2*n1*x0 - f0*h1*s2*v2 - f0*h2*s1*v2 - f0*h2*s2*v1 - f0*k1*n2*x2 - f0*k2*n1*x2 - f0*k2*n2*x1 - f1*h0*s2*v2 - f1*h1*s1*v2 - f1*h1*s2*v1 - f1*h2*s0*v2 - f1*h2*s1*v1 - f1*h2*s2*v0 - f1*k0*n2*x2 - f1*k1*n1*x2 - f1*k1*n2*x1 - f1*k2*n0*x2 - f1*k2*n1*x1 - f1*k2*n2*x0 - f2*h0*s1*v2 - f2*h0*s2*v1 - f2*h1*s0*v2 - f2*h1*s1*v1 - f2*h1*s2*v0 - f2*h2*s0*v1 - f2*h2*s1*v0 - f2*k0*n1*x2 - f2*k0*n2*x1 - f2*k1*n0*x2 - f2*k1*n1*x1 - f2*k1*n2*x0 - f2*k2*n0*x1 - f2*k2*n1*x0 - b0*h1*u2*y2 - b0*h2*u1*y2 - b0*h2*u2*y1 - b1*h0*u2*y2 - b1*h1*u1*y2 - b1*h1*u2*y1 - b1*h2*u0*y2 - b1*h2*u1*y1 - b1*h2*u2*y0 - b2*h0*u1*y2 - b2*h0*u2*y1 - b2*h1*u0*y2 - b2*h1*u1*y1 - b2*h1*u2*y0 - b2*h2*u0*y1 - b2*h2*u1*y0 + c0*g1*u2*y2 + c0*g2*u1*y2 + c0*g2*u2*y1 + c1*g0*u2*y2 + c1*g1*u1*y2 + c1*g1*u2*y1 + c1*g2*u0*y2 + c1*g2*u1*y1 + c1*g2*u2*y0 + c2*g0*u1*y2 + c2*g0*u2*y1 + c2*g1*u0*y2 + c2*g1*u1*y1 + c2*g1*u2*y0 + c2*g2*u0*y1 + c2*g2*u1*y0 - d0*g1*u2*x2 - d0*g2*u1*x2 - d0*g2*u2*x1 - d1*g0*u2*x2 - d1*g1*u1*x2 - d1*g1*u2*x1 - d1*g2*u0*x2 - d1*g2*u1*x1 - d1*g2*u2*x0 - d2*g0*u1*x2 - d2*g0*u2*x1 - d2*g1*u0*x2 - d2*g1*u1*x1 - d2*g1*u2*x0 - d2*g2*u0*x1 - d2*g2*u1*x0 + f0*g1*s2*x2 + f0*g2*s1*x2 + f0*g2*s2*x1 + f1*g0*s2*x2 + f1*g1*s1*x2 + f1*g1*s2*x1 + f1*g2*s0*x2 + f1*g2*s1*x1 + f1*g2*s2*x0 + f2*g0*s1*x2 + f2*g0*s2*x1 + f2*g1*s0*x2 + f2*g1*s1*x1 + f2*g1*s2*x0 + f2*g2*s0*x1 + f2*g2*s1*x0 - c0*k1*u2*v2 - c0*k2*u1*v2 - c0*k2*u2*v1 + c0*m1*s2*v2 + c0*m2*s1*v2 + c0*m2*s2*v1 - c1*k0*u2*v2 - c1*k1*u1*v2 - c1*k1*u2*v1 - c1*k2*u0*v2 - c1*k2*u1*v1 - c1*k2*u2*v0 + c1*m0*s2*v2 + c1*m1*s1*v2 + c1*m1*s2*v1 + c1*m2*s0*v2 + c1*m2*s1*v1 + c1*m2*s2*v0 - c2*k0*u1*v2 - c2*k0*u2*v1 - c2*k1*u0*v2 - c2*k1*u1*v1 - c2*k1*u2*v0 - c2*k2*u0*v1 - c2*k2*u1*v0 + c2*m0*s1*v2 + c2*m0*s2*v1 + c2*m1*s0*v2 + c2*m1*s1*v1 + c2*m1*s2*v0 + c2*m2*s0*v1 + c2*m2*s1*v0 + b0*k1*u2*x2 + b0*k2*u1*x2 + b0*k2*u2*x1 - b0*m1*s2*x2 - b0*m2*s1*x2 - b0*m2*s2*x1 + b1*k0*u2*x2 + b1*k1*u1*x2 + b1*k1*u2*x1 + b1*k2*u0*x2 + b1*k2*u1*x1 + b1*k2*u2*x0 - b1*m0*s2*x2 - b1*m1*s1*x2 - b1*m1*s2*x1 - b1*m2*s0*x2 - b1*m2*s1*x1 - b1*m2*s2*x0 + b2*k0*u1*x2 + b2*k0*u2*x1 + b2*k1*u0*x2 + b2*k1*u1*x1 + b2*k1*u2*x0 + b2*k2*u0*x1 + b2*k2*u1*x0 - b2*m0*s1*x2 - b2*m0*s2*x1 - b2*m1*s0*x2 - b2*m1*s1*x1 - b2*m1*s2*x0 - b2*m2*s0*x1 - b2*m2*s1*x0;
	f_qAll[4] = d0*g0*o2*z2 + d0*g1*o1*z2 + d0*g1*o2*z1 + d0*g2*o0*z2 + d0*g2*o1*z1 + d0*g2*o2*z0 - d0*h0*n2*z2 - d0*h1*n1*z2 - d0*h1*n2*z1 - d0*h2*n0*z2 - d0*h2*n1*z1 - d0*h2*n2*z0 + d1*g0*o1*z2 + d1*g0*o2*z1 + d1*g1*o0*z2 + d1*g1*o1*z1 + d1*g1*o2*z0 + d1*g2*o0*z1 + d1*g2*o1*z0 - d1*h0*n1*z2 - d1*h0*n2*z1 - d1*h1*n0*z2 - d1*h1*n1*z1 - d1*h1*n2*z0 - d1*h2*n0*z1 - d1*h2*n1*z0 + d2*g0*o0*z2 + d2*g0*o1*z1 + d2*g0*o2*z0 + d2*g1*o0*z1 + d2*g1*o1*z0 + d2*g2*o0*z0 - d2*h0*n0*z2 - d2*h0*n1*z1 - d2*h0*n2*z0 - d2*h1*n0*z1 - d2*h1*n1*z0 - d2*h2*n0*z0 - f0*g0*o2*y2 - f0*g1*o1*y2 - f0*g1*o2*y1 - f0*g2*o0*y2 - f0*g2*o1*y1 - f0*g2*o2*y0 + f0*h0*n2*y2 + f0*h1*n1*y2 + f0*h1*n2*y1 + f0*h2*n0*y2 + f0*h2*n1*y1 + f0*h2*n2*y0 - f1*g0*o1*y2 - f1*g0*o2*y1 - f1*g1*o0*y2 - f1*g1*o1*y1 - f1*g1*o2*y0 - f1*g2*o0*y1 - f1*g2*o1*y0 + f1*h0*n1*y2 + f1*h0*n2*y1 + f1*h1*n0*y2 + f1*h1*n1*y1 + f1*h1*n2*y0 + f1*h2*n0*y1 + f1*h2*n1*y0 - f2*g0*o0*y2 - f2*g0*o1*y1 - f2*g0*o2*y0 - f2*g1*o0*y1 - f2*g1*o1*y0 - f2*g2*o0*y0 + f2*h0*n0*y2 + f2*h0*n1*y1 + f2*h0*n2*y0 + f2*h1*n0*y1 + f2*h1*n1*y0 + f2*h2*n0*y0 - b0*k0*o2*z2 - b0*k1*o1*z2 - b0*k1*o2*z1 - b0*k2*o0*z2 - b0*k2*o1*z1 - b0*k2*o2*z0 - b1*k0*o1*z2 - b1*k0*o2*z1 - b1*k1*o0*z2 - b1*k1*o1*z1 - b1*k1*o2*z0 - b1*k2*o0*z1 - b1*k2*o1*z0 - b2*k0*o0*z2 - b2*k0*o1*z1 - b2*k0*o2*z0 - b2*k1*o0*z1 - b2*k1*o1*z0 - b2*k2*o0*z0 + c0*k0*n2*z2 + c0*k1*n1*z2 + c0*k1*n2*z1 + c0*k2*n0*z2 + c0*k2*n1*z1 + c0*k2*n2*z0 + c1*k0*n1*z2 + c1*k0*n2*z1 + c1*k1*n0*z2 + c1*k1*n1*z1 + c1*k1*n2*z0 + c1*k2*n0*z1 + c1*k2*n1*z0 + c2*k0*n0*z2 + c2*k0*n1*z1 + c2*k0*n2*z0 + c2*k1*n0*z1 + c2*k1*n1*z0 + c2*k2*n0*z0 - d0*m0*o2*v2 - d0*m1*o1*v2 - d0*m1*o2*v1 - d0*m2*o0*v2 - d0*m2*o1*v1 - d0*m2*o2*v0 - d1*m0*o1*v2 - d1*m0*o2*v1 - d1*m1*o0*v2 - d1*m1*o1*v1 - d1*m1*o2*v0 - d1*m2*o0*v1 - d1*m2*o1*v0 - d2*m0*o0*v2 - d2*m0*o1*v1 - d2*m0*o2*v0 - d2*m1*o0*v1 - d2*m1*o1*v0 - d2*m2*o0*v0 + f0*k0*o2*v2 + f0*k1*o1*v2 + f0*k1*o2*v1 + f0*k2*o0*v2 + f0*k2*o1*v1 + f0*k2*o2*v0 + f1*k0*o1*v2 + f1*k0*o2*v1 + f1*k1*o0*v2 + f1*k1*o1*v1 + f1*k1*o2*v0 + f1*k2*o0*v1 + f1*k2*o1*v0 + f2*k0*o0*v2 + f2*k0*o1*v1 + f2*k0*o2*v0 + f2*k1*o0*v1 + f2*k1*o1*v0 + f2*k2*o0*v0 + b0*h0*s2*z2 + b0*h1*s1*z2 + b0*h1*s2*z1 + b0*h2*s0*z2 + b0*h2*s1*z1 + b0*h2*s2*z0 + b0*m0*o2*y2 + b0*m1*o1*y2 + b0*m1*o2*y1 + b0*m2*o0*y2 + b0*m2*o1*y1 + b0*m2*o2*y0 + b1*h0*s1*z2 + b1*h0*s2*z1 + b1*h1*s0*z2 + b1*h1*s1*z1 + b1*h1*s2*z0 + b1*h2*s0*z1 + b1*h2*s1*z0 + b1*m0*o1*y2 + b1*m0*o2*y1 + b1*m1*o0*y2 + b1*m1*o1*y1 + b1*m1*o2*y0 + b1*m2*o0*y1 + b1*m2*o1*y0 + b2*h0*s0*z2 + b2*h0*s1*z1 + b2*h0*s2*z0 + b2*h1*s0*z1 + b2*h1*s1*z0 + b2*h2*s0*z0 + b2*m0*o0*y2 + b2*m0*o1*y1 + b2*m0*o2*y0 + b2*m1*o0*y1 + b2*m1*o1*y0 + b2*m2*o0*y0 - c0*g0*s2*z2 - c0*g1*s1*z2 - c0*g1*s2*z1 - c0*g2*s0*z2 - c0*g2*s1*z1 - c0*g2*s2*z0 - c0*m0*n2*y2 - c0*m1*n1*y2 - c0*m1*n2*y1 - c0*m2*n0*y2 - c0*m2*n1*y1 - c0*m2*n2*y0 - c1*g0*s1*z2 - c1*g0*s2*z1 - c1*g1*s0*z2 - c1*g1*s1*z1 - c1*g1*s2*z0 - c1*g2*s0*z1 - c1*g2*s1*z0 - c1*m0*n1*y2 - c1*m0*n2*y1 - c1*m1*n0*y2 - c1*m1*n1*y1 - c1*m1*n2*y0 - c1*m2*n0*y1 - c1*m2*n1*y0 - c2*g0*s0*z2 - c2*g0*s1*z1 - c2*g0*s2*z0 - c2*g1*s0*z1 - c2*g1*s1*z0 - c2*g2*s0*z0 - c2*m0*n0*y2 - c2*m0*n1*y1 - c2*m0*n2*y0 - c2*m1*n0*y1 - c2*m1*n1*y0 - c2*m2*n0*y0 + d0*h0*u2*v2 + d0*h1*u1*v2 + d0*h1*u2*v1 + d0*h2*u0*v2 + d0*h2*u1*v1 + d0*h2*u2*v0 + d0*m0*n2*x2 + d0*m1*n1*x2 + d0*m1*n2*x1 + d0*m2*n0*x2 + d0*m2*n1*x1 + d0*m2*n2*x0 + d1*h0*u1*v2 + d1*h0*u2*v1 + d1*h1*u0*v2 + d1*h1*u1*v1 + d1*h1*u2*v0 + d1*h2*u0*v1 + d1*h2*u1*v0 + d1*m0*n1*x2 + d1*m0*n2*x1 + d1*m1*n0*x2 + d1*m1*n1*x1 + d1*m1*n2*x0 + d1*m2*n0*x1 + d1*m2*n1*x0 + d2*h0*u0*v2 + d2*h0*u1*v1 + d2*h0*u2*v0 + d2*h1*u0*v1 + d2*h1*u1*v0 + d2*h2*u0*v0 + d2*m0*n0*x2 + d2*m0*n1*x1 + d2*m0*n2*x0 + d2*m1*n0*x1 + d2*m1*n1*x0 + d2*m2*n0*x0 - f0*h0*s2*v2 - f0*h1*s1*v2 - f0*h1*s2*v1 - f0*h2*s0*v2 - f0*h2*s1*v1 - f0*h2*s2*v0 - f0*k0*n2*x2 - f0*k1*n1*x2 - f0*k1*n2*x1 - f0*k2*n0*x2 - f0*k2*n1*x1 - f0*k2*n2*x0 - f1*h0*s1*v2 - f1*h0*s2*v1 - f1*h1*s0*v2 - f1*h1*s1*v1 - f1*h1*s2*v0 - f1*h2*s0*v1 - f1*h2*s1*v0 - f1*k0*n1*x2 - f1*k0*n2*x1 - f1*k1*n0*x2 - f1*k1*n1*x1 - f1*k1*n2*x0 - f1*k2*n0*x1 - f1*k2*n1*x0 - f2*h0*s0*v2 - f2*h0*s1*v1 - f2*h0*s2*v0 - f2*h1*s0*v1 - f2*h1*s1*v0 - f2*h2*s0*v0 - f2*k0*n0*x2 - f2*k0*n1*x1 - f2*k0*n2*x0 - f2*k1*n0*x1 - f2*k1*n1*x0 - f2*k2*n0*x0 - b0*h0*u2*y2 - b0*h1*u1*y2 - b0*h1*u2*y1 - b0*h2*u0*y2 - b0*h2*u1*y1 - b0*h2*u2*y0 - b1*h0*u1*y2 - b1*h0*u2*y1 - b1*h1*u0*y2 - b1*h1*u1*y1 - b1*h1*u2*y0 - b1*h2*u0*y1 - b1*h2*u1*y0 - b2*h0*u0*y2 - b2*h0*u1*y1 - b2*h0*u2*y0 - b2*h1*u0*y1 - b2*h1*u1*y0 - b2*h2*u0*y0 + c0*g0*u2*y2 + c0*g1*u1*y2 + c0*g1*u2*y1 + c0*g2*u0*y2 + c0*g2*u1*y1 + c0*g2*u2*y0 + c1*g0*u1*y2 + c1*g0*u2*y1 + c1*g1*u0*y2 + c1*g1*u1*y1 + c1*g1*u2*y0 + c1*g2*u0*y1 + c1*g2*u1*y0 + c2*g0*u0*y2 + c2*g0*u1*y1 + c2*g0*u2*y0 + c2*g1*u0*y1 + c2*g1*u1*y0 + c2*g2*u0*y0 - d0*g0*u2*x2 - d0*g1*u1*x2 - d0*g1*u2*x1 - d0*g2*u0*x2 - d0*g2*u1*x1 - d0*g2*u2*x0 - d1*g0*u1*x2 - d1*g0*u2*x1 - d1*g1*u0*x2 - d1*g1*u1*x1 - d1*g1*u2*x0 - d1*g2*u0*x1 - d1*g2*u1*x0 - d2*g0*u0*x2 - d2*g0*u1*x1 - d2*g0*u2*x0 - d2*g1*u0*x1 - d2*g1*u1*x0 - d2*g2*u0*x0 + f0*g0*s2*x2 + f0*g1*s1*x2 + f0*g1*s2*x1 + f0*g2*s0*x2 + f0*g2*s1*x1 + f0*g2*s2*x0 + f1*g0*s1*x2 + f1*g0*s2*x1 + f1*g1*s0*x2 + f1*g1*s1*x1 + f1*g1*s2*x0 + f1*g2*s0*x1 + f1*g2*s1*x0 + f2*g0*s0*x2 + f2*g0*s1*x1 + f2*g0*s2*x0 + f2*g1*s0*x1 + f2*g1*s1*x0 + f2*g2*s0*x0 - c0*k0*u2*v2 - c0*k1*u1*v2 - c0*k1*u2*v1 - c0*k2*u0*v2 - c0*k2*u1*v1 - c0*k2*u2*v0 + c0*m0*s2*v2 + c0*m1*s1*v2 + c0*m1*s2*v1 + c0*m2*s0*v2 + c0*m2*s1*v1 + c0*m2*s2*v0 - c1*k0*u1*v2 - c1*k0*u2*v1 - c1*k1*u0*v2 - c1*k1*u1*v1 - c1*k1*u2*v0 - c1*k2*u0*v1 - c1*k2*u1*v0 + c1*m0*s1*v2 + c1*m0*s2*v1 + c1*m1*s0*v2 + c1*m1*s1*v1 + c1*m1*s2*v0 + c1*m2*s0*v1 + c1*m2*s1*v0 - c2*k0*u0*v2 - c2*k0*u1*v1 - c2*k0*u2*v0 - c2*k1*u0*v1 - c2*k1*u1*v0 - c2*k2*u0*v0 + c2*m0*s0*v2 + c2*m0*s1*v1 + c2*m0*s2*v0 + c2*m1*s0*v1 + c2*m1*s1*v0 + c2*m2*s0*v0 + b0*k0*u2*x2 + b0*k1*u1*x2 + b0*k1*u2*x1 + b0*k2*u0*x2 + b0*k2*u1*x1 + b0*k2*u2*x0 - b0*m0*s2*x2 - b0*m1*s1*x2 - b0*m1*s2*x1 - b0*m2*s0*x2 - b0*m2*s1*x1 - b0*m2*s2*x0 + b1*k0*u1*x2 + b1*k0*u2*x1 + b1*k1*u0*x2 + b1*k1*u1*x1 + b1*k1*u2*x0 + b1*k2*u0*x1 + b1*k2*u1*x0 - b1*m0*s1*x2 - b1*m0*s2*x1 - b1*m1*s0*x2 - b1*m1*s1*x1 - b1*m1*s2*x0 - b1*m2*s0*x1 - b1*m2*s1*x0 + b2*k0*u0*x2 + b2*k0*u1*x1 + b2*k0*u2*x0 + b2*k1*u0*x1 + b2*k1*u1*x0 + b2*k2*u0*x0 - b2*m0*s0*x2 - b2*m0*s1*x1 - b2*m0*s2*x0 - b2*m1*s0*x1 - b2*m1*s1*x0 - b2*m2*s0*x0;
	f_qAll[5] = d0*g0*o1*z2 + d0*g0*o2*z1 + d0*g1*o0*z2 + d0*g1*o1*z1 + d0*g1*o2*z0 + d0*g2*o0*z1 + d0*g2*o1*z0 - d0*h0*n1*z2 - d0*h0*n2*z1 - d0*h1*n0*z2 - d0*h1*n1*z1 - d0*h1*n2*z0 - d0*h2*n0*z1 - d0*h2*n1*z0 + d1*g0*o0*z2 + d1*g0*o1*z1 + d1*g0*o2*z0 + d1*g1*o0*z1 + d1*g1*o1*z0 + d1*g2*o0*z0 - d1*h0*n0*z2 - d1*h0*n1*z1 - d1*h0*n2*z0 - d1*h1*n0*z1 - d1*h1*n1*z0 - d1*h2*n0*z0 + d2*g0*o0*z1 + d2*g0*o1*z0 + d2*g1*o0*z0 - d2*h0*n0*z1 - d2*h0*n1*z0 - d2*h1*n0*z0 - f0*g0*o1*y2 - f0*g0*o2*y1 - f0*g1*o0*y2 - f0*g1*o1*y1 - f0*g1*o2*y0 - f0*g2*o0*y1 - f0*g2*o1*y0 + f0*h0*n1*y2 + f0*h0*n2*y1 + f0*h1*n0*y2 + f0*h1*n1*y1 + f0*h1*n2*y0 + f0*h2*n0*y1 + f0*h2*n1*y0 - f1*g0*o0*y2 - f1*g0*o1*y1 - f1*g0*o2*y0 - f1*g1*o0*y1 - f1*g1*o1*y0 - f1*g2*o0*y0 + f1*h0*n0*y2 + f1*h0*n1*y1 + f1*h0*n2*y0 + f1*h1*n0*y1 + f1*h1*n1*y0 + f1*h2*n0*y0 - f2*g0*o0*y1 - f2*g0*o1*y0 - f2*g1*o0*y0 + f2*h0*n0*y1 + f2*h0*n1*y0 + f2*h1*n0*y0 - b0*k0*o1*z2 - b0*k0*o2*z1 - b0*k1*o0*z2 - b0*k1*o1*z1 - b0*k1*o2*z0 - b0*k2*o0*z1 - b0*k2*o1*z0 - b1*k0*o0*z2 - b1*k0*o1*z1 - b1*k0*o2*z0 - b1*k1*o0*z1 - b1*k1*o1*z0 - b1*k2*o0*z0 - b2*k0*o0*z1 - b2*k0*o1*z0 - b2*k1*o0*z0 + c0*k0*n1*z2 + c0*k0*n2*z1 + c0*k1*n0*z2 + c0*k1*n1*z1 + c0*k1*n2*z0 + c0*k2*n0*z1 + c0*k2*n1*z0 + c1*k0*n0*z2 + c1*k0*n1*z1 + c1*k0*n2*z0 + c1*k1*n0*z1 + c1*k1*n1*z0 + c1*k2*n0*z0 + c2*k0*n0*z1 + c2*k0*n1*z0 + c2*k1*n0*z0 - d0*m0*o1*v2 - d0*m0*o2*v1 - d0*m1*o0*v2 - d0*m1*o1*v1 - d0*m1*o2*v0 - d0*m2*o0*v1 - d0*m2*o1*v0 - d1*m0*o0*v2 - d1*m0*o1*v1 - d1*m0*o2*v0 - d1*m1*o0*v1 - d1*m1*o1*v0 - d1*m2*o0*v0 - d2*m0*o0*v1 - d2*m0*o1*v0 - d2*m1*o0*v0 + f0*k0*o1*v2 + f0*k0*o2*v1 + f0*k1*o0*v2 + f0*k1*o1*v1 + f0*k1*o2*v0 + f0*k2*o0*v1 + f0*k2*o1*v0 + f1*k0*o0*v2 + f1*k0*o1*v1 + f1*k0*o2*v0 + f1*k1*o0*v1 + f1*k1*o1*v0 + f1*k2*o0*v0 + f2*k0*o0*v1 + f2*k0*o1*v0 + f2*k1*o0*v0 + b0*h0*s1*z2 + b0*h0*s2*z1 + b0*h1*s0*z2 + b0*h1*s1*z1 + b0*h1*s2*z0 + b0*h2*s0*z1 + b0*h2*s1*z0 + b0*m0*o1*y2 + b0*m0*o2*y1 + b0*m1*o0*y2 + b0*m1*o1*y1 + b0*m1*o2*y0 + b0*m2*o0*y1 + b0*m2*o1*y0 + b1*h0*s0*z2 + b1*h0*s1*z1 + b1*h0*s2*z0 + b1*h1*s0*z1 + b1*h1*s1*z0 + b1*h2*s0*z0 + b1*m0*o0*y2 + b1*m0*o1*y1 + b1*m0*o2*y0 + b1*m1*o0*y1 + b1*m1*o1*y0 + b1*m2*o0*y0 + b2*h0*s0*z1 + b2*h0*s1*z0 + b2*h1*s0*z0 + b2*m0*o0*y1 + b2*m0*o1*y0 + b2*m1*o0*y0 - c0*g0*s1*z2 - c0*g0*s2*z1 - c0*g1*s0*z2 - c0*g1*s1*z1 - c0*g1*s2*z0 - c0*g2*s0*z1 - c0*g2*s1*z0 - c0*m0*n1*y2 - c0*m0*n2*y1 - c0*m1*n0*y2 - c0*m1*n1*y1 - c0*m1*n2*y0 - c0*m2*n0*y1 - c0*m2*n1*y0 - c1*g0*s0*z2 - c1*g0*s1*z1 - c1*g0*s2*z0 - c1*g1*s0*z1 - c1*g1*s1*z0 - c1*g2*s0*z0 - c1*m0*n0*y2 - c1*m0*n1*y1 - c1*m0*n2*y0 - c1*m1*n0*y1 - c1*m1*n1*y0 - c1*m2*n0*y0 - c2*g0*s0*z1 - c2*g0*s1*z0 - c2*g1*s0*z0 - c2*m0*n0*y1 - c2*m0*n1*y0 - c2*m1*n0*y0 + d0*h0*u1*v2 + d0*h0*u2*v1 + d0*h1*u0*v2 + d0*h1*u1*v1 + d0*h1*u2*v0 + d0*h2*u0*v1 + d0*h2*u1*v0 + d0*m0*n1*x2 + d0*m0*n2*x1 + d0*m1*n0*x2 + d0*m1*n1*x1 + d0*m1*n2*x0 + d0*m2*n0*x1 + d0*m2*n1*x0 + d1*h0*u0*v2 + d1*h0*u1*v1 + d1*h0*u2*v0 + d1*h1*u0*v1 + d1*h1*u1*v0 + d1*h2*u0*v0 + d1*m0*n0*x2 + d1*m0*n1*x1 + d1*m0*n2*x0 + d1*m1*n0*x1 + d1*m1*n1*x0 + d1*m2*n0*x0 + d2*h0*u0*v1 + d2*h0*u1*v0 + d2*h1*u0*v0 + d2*m0*n0*x1 + d2*m0*n1*x0 + d2*m1*n0*x0 - f0*h0*s1*v2 - f0*h0*s2*v1 - f0*h1*s0*v2 - f0*h1*s1*v1 - f0*h1*s2*v0 - f0*h2*s0*v1 - f0*h2*s1*v0 - f0*k0*n1*x2 - f0*k0*n2*x1 - f0*k1*n0*x2 - f0*k1*n1*x1 - f0*k1*n2*x0 - f0*k2*n0*x1 - f0*k2*n1*x0 - f1*h0*s0*v2 - f1*h0*s1*v1 - f1*h0*s2*v0 - f1*h1*s0*v1 - f1*h1*s1*v0 - f1*h2*s0*v0 - f1*k0*n0*x2 - f1*k0*n1*x1 - f1*k0*n2*x0 - f1*k1*n0*x1 - f1*k1*n1*x0 - f1*k2*n0*x0 - f2*h0*s0*v1 - f2*h0*s1*v0 - f2*h1*s0*v0 - f2*k0*n0*x1 - f2*k0*n1*x0 - f2*k1*n0*x0 - b0*h0*u1*y2 - b0*h0*u2*y1 - b0*h1*u0*y2 - b0*h1*u1*y1 - b0*h1*u2*y0 - b0*h2*u0*y1 - b0*h2*u1*y0 - b1*h0*u0*y2 - b1*h0*u1*y1 - b1*h0*u2*y0 - b1*h1*u0*y1 - b1*h1*u1*y0 - b1*h2*u0*y0 - b2*h0*u0*y1 - b2*h0*u1*y0 - b2*h1*u0*y0 + c0*g0*u1*y2 + c0*g0*u2*y1 + c0*g1*u0*y2 + c0*g1*u1*y1 + c0*g1*u2*y0 + c0*g2*u0*y1 + c0*g2*u1*y0 + c1*g0*u0*y2 + c1*g0*u1*y1 + c1*g0*u2*y0 + c1*g1*u0*y1 + c1*g1*u1*y0 + c1*g2*u0*y0 + c2*g0*u0*y1 + c2*g0*u1*y0 + c2*g1*u0*y0 - d0*g0*u1*x2 - d0*g0*u2*x1 - d0*g1*u0*x2 - d0*g1*u1*x1 - d0*g1*u2*x0 - d0*g2*u0*x1 - d0*g2*u1*x0 - d1*g0*u0*x2 - d1*g0*u1*x1 - d1*g0*u2*x0 - d1*g1*u0*x1 - d1*g1*u1*x0 - d1*g2*u0*x0 - d2*g0*u0*x1 - d2*g0*u1*x0 - d2*g1*u0*x0 + f0*g0*s1*x2 + f0*g0*s2*x1 + f0*g1*s0*x2 + f0*g1*s1*x1 + f0*g1*s2*x0 + f0*g2*s0*x1 + f0*g2*s1*x0 + f1*g0*s0*x2 + f1*g0*s1*x1 + f1*g0*s2*x0 + f1*g1*s0*x1 + f1*g1*s1*x0 + f1*g2*s0*x0 + f2*g0*s0*x1 + f2*g0*s1*x0 + f2*g1*s0*x0 - c0*k0*u1*v2 - c0*k0*u2*v1 - c0*k1*u0*v2 - c0*k1*u1*v1 - c0*k1*u2*v0 - c0*k2*u0*v1 - c0*k2*u1*v0 + c0*m0*s1*v2 + c0*m0*s2*v1 + c0*m1*s0*v2 + c0*m1*s1*v1 + c0*m1*s2*v0 + c0*m2*s0*v1 + c0*m2*s1*v0 - c1*k0*u0*v2 - c1*k0*u1*v1 - c1*k0*u2*v0 - c1*k1*u0*v1 - c1*k1*u1*v0 - c1*k2*u0*v0 + c1*m0*s0*v2 + c1*m0*s1*v1 + c1*m0*s2*v0 + c1*m1*s0*v1 + c1*m1*s1*v0 + c1*m2*s0*v0 - c2*k0*u0*v1 - c2*k0*u1*v0 - c2*k1*u0*v0 + c2*m0*s0*v1 + c2*m0*s1*v0 + c2*m1*s0*v0 + b0*k0*u1*x2 + b0*k0*u2*x1 + b0*k1*u0*x2 + b0*k1*u1*x1 + b0*k1*u2*x0 + b0*k2*u0*x1 + b0*k2*u1*x0 - b0*m0*s1*x2 - b0*m0*s2*x1 - b0*m1*s0*x2 - b0*m1*s1*x1 - b0*m1*s2*x0 - b0*m2*s0*x1 - b0*m2*s1*x0 + b1*k0*u0*x2 + b1*k0*u1*x1 + b1*k0*u2*x0 + b1*k1*u0*x1 + b1*k1*u1*x0 + b1*k2*u0*x0 - b1*m0*s0*x2 - b1*m0*s1*x1 - b1*m0*s2*x0 - b1*m1*s0*x1 - b1*m1*s1*x0 - b1*m2*s0*x0 + b2*k0*u0*x1 + b2*k0*u1*x0 + b2*k1*u0*x0 - b2*m0*s0*x1 - b2*m0*s1*x0 - b2*m1*s0*x0;
	f_qAll[6] = d0*g0*o0*z2 + d0*g0*o1*z1 + d0*g0*o2*z0 + d0*g1*o0*z1 + d0*g1*o1*z0 + d0*g2*o0*z0 - d0*h0*n0*z2 - d0*h0*n1*z1 - d0*h0*n2*z0 - d0*h1*n0*z1 - d0*h1*n1*z0 - d0*h2*n0*z0 + d1*g0*o0*z1 + d1*g0*o1*z0 + d1*g1*o0*z0 - d1*h0*n0*z1 - d1*h0*n1*z0 - d1*h1*n0*z0 + d2*g0*o0*z0 - d2*h0*n0*z0 - f0*g0*o0*y2 - f0*g0*o1*y1 - f0*g0*o2*y0 - f0*g1*o0*y1 - f0*g1*o1*y0 - f0*g2*o0*y0 + f0*h0*n0*y2 + f0*h0*n1*y1 + f0*h0*n2*y0 + f0*h1*n0*y1 + f0*h1*n1*y0 + f0*h2*n0*y0 - f1*g0*o0*y1 - f1*g0*o1*y0 - f1*g1*o0*y0 + f1*h0*n0*y1 + f1*h0*n1*y0 + f1*h1*n0*y0 - f2*g0*o0*y0 + f2*h0*n0*y0 - b0*k0*o0*z2 - b0*k0*o1*z1 - b0*k0*o2*z0 - b0*k1*o0*z1 - b0*k1*o1*z0 - b0*k2*o0*z0 - b1*k0*o0*z1 - b1*k0*o1*z0 - b1*k1*o0*z0 - b2*k0*o0*z0 + c0*k0*n0*z2 + c0*k0*n1*z1 + c0*k0*n2*z0 + c0*k1*n0*z1 + c0*k1*n1*z0 + c0*k2*n0*z0 + c1*k0*n0*z1 + c1*k0*n1*z0 + c1*k1*n0*z0 + c2*k0*n0*z0 - d0*m0*o0*v2 - d0*m0*o1*v1 - d0*m0*o2*v0 - d0*m1*o0*v1 - d0*m1*o1*v0 - d0*m2*o0*v0 - d1*m0*o0*v1 - d1*m0*o1*v0 - d1*m1*o0*v0 - d2*m0*o0*v0 + f0*k0*o0*v2 + f0*k0*o1*v1 + f0*k0*o2*v0 + f0*k1*o0*v1 + f0*k1*o1*v0 + f0*k2*o0*v0 + f1*k0*o0*v1 + f1*k0*o1*v0 + f1*k1*o0*v0 + f2*k0*o0*v0 + b0*h0*s0*z2 + b0*h0*s1*z1 + b0*h0*s2*z0 + b0*h1*s0*z1 + b0*h1*s1*z0 + b0*h2*s0*z0 + b0*m0*o0*y2 + b0*m0*o1*y1 + b0*m0*o2*y0 + b0*m1*o0*y1 + b0*m1*o1*y0 + b0*m2*o0*y0 + b1*h0*s0*z1 + b1*h0*s1*z0 + b1*h1*s0*z0 + b1*m0*o0*y1 + b1*m0*o1*y0 + b1*m1*o0*y0 + b2*h0*s0*z0 + b2*m0*o0*y0 - c0*g0*s0*z2 - c0*g0*s1*z1 - c0*g0*s2*z0 - c0*g1*s0*z1 - c0*g1*s1*z0 - c0*g2*s0*z0 - c0*m0*n0*y2 - c0*m0*n1*y1 - c0*m0*n2*y0 - c0*m1*n0*y1 - c0*m1*n1*y0 - c0*m2*n0*y0 - c1*g0*s0*z1 - c1*g0*s1*z0 - c1*g1*s0*z0 - c1*m0*n0*y1 - c1*m0*n1*y0 - c1*m1*n0*y0 - c2*g0*s0*z0 - c2*m0*n0*y0 + d0*h0*u0*v2 + d0*h0*u1*v1 + d0*h0*u2*v0 + d0*h1*u0*v1 + d0*h1*u1*v0 + d0*h2*u0*v0 + d0*m0*n0*x2 + d0*m0*n1*x1 + d0*m0*n2*x0 + d0*m1*n0*x1 + d0*m1*n1*x0 + d0*m2*n0*x0 + d1*h0*u0*v1 + d1*h0*u1*v0 + d1*h1*u0*v0 + d1*m0*n0*x1 + d1*m0*n1*x0 + d1*m1*n0*x0 + d2*h0*u0*v0 + d2*m0*n0*x0 - f0*h0*s0*v2 - f0*h0*s1*v1 - f0*h0*s2*v0 - f0*h1*s0*v1 - f0*h1*s1*v0 - f0*h2*s0*v0 - f0*k0*n0*x2 - f0*k0*n1*x1 - f0*k0*n2*x0 - f0*k1*n0*x1 - f0*k1*n1*x0 - f0*k2*n0*x0 - f1*h0*s0*v1 - f1*h0*s1*v0 - f1*h1*s0*v0 - f1*k0*n0*x1 - f1*k0*n1*x0 - f1*k1*n0*x0 - f2*h0*s0*v0 - f2*k0*n0*x0 - b0*h0*u0*y2 - b0*h0*u1*y1 - b0*h0*u2*y0 - b0*h1*u0*y1 - b0*h1*u1*y0 - b0*h2*u0*y0 - b1*h0*u0*y1 - b1*h0*u1*y0 - b1*h1*u0*y0 - b2*h0*u0*y0 + c0*g0*u0*y2 + c0*g0*u1*y1 + c0*g0*u2*y0 + c0*g1*u0*y1 + c0*g1*u1*y0 + c0*g2*u0*y0 + c1*g0*u0*y1 + c1*g0*u1*y0 + c1*g1*u0*y0 + c2*g0*u0*y0 - d0*g0*u0*x2 - d0*g0*u1*x1 - d0*g0*u2*x0 - d0*g1*u0*x1 - d0*g1*u1*x0 - d0*g2*u0*x0 - d1*g0*u0*x1 - d1*g0*u1*x0 - d1*g1*u0*x0 - d2*g0*u0*x0 + f0*g0*s0*x2 + f0*g0*s1*x1 + f0*g0*s2*x0 + f0*g1*s0*x1 + f0*g1*s1*x0 + f0*g2*s0*x0 + f1*g0*s0*x1 + f1*g0*s1*x0 + f1*g1*s0*x0 + f2*g0*s0*x0 - c0*k0*u0*v2 - c0*k0*u1*v1 - c0*k0*u2*v0 - c0*k1*u0*v1 - c0*k1*u1*v0 - c0*k2*u0*v0 + c0*m0*s0*v2 + c0*m0*s1*v1 + c0*m0*s2*v0 + c0*m1*s0*v1 + c0*m1*s1*v0 + c0*m2*s0*v0 - c1*k0*u0*v1 - c1*k0*u1*v0 - c1*k1*u0*v0 + c1*m0*s0*v1 + c1*m0*s1*v0 + c1*m1*s0*v0 - c2*k0*u0*v0 + c2*m0*s0*v0 + b0*k0*u0*x2 + b0*k0*u1*x1 + b0*k0*u2*x0 + b0*k1*u0*x1 + b0*k1*u1*x0 + b0*k2*u0*x0 - b0*m0*s0*x2 - b0*m0*s1*x1 - b0*m0*s2*x0 - b0*m1*s0*x1 - b0*m1*s1*x0 - b0*m2*s0*x0 + b1*k0*u0*x1 + b1*k0*u1*x0 + b1*k1*u0*x0 - b1*m0*s0*x1 - b1*m0*s1*x0 - b1*m1*s0*x0 + b2*k0*u0*x0 - b2*m0*s0*x0;
	f_qAll[7] = d0*g0*o0*z1 + d0*g0*o1*z0 + d0*g1*o0*z0 - d0*h0*n0*z1 - d0*h0*n1*z0 - d0*h1*n0*z0 + d1*g0*o0*z0 - d1*h0*n0*z0 - f0*g0*o0*y1 - f0*g0*o1*y0 - f0*g1*o0*y0 + f0*h0*n0*y1 + f0*h0*n1*y0 + f0*h1*n0*y0 - f1*g0*o0*y0 + f1*h0*n0*y0 - b0*k0*o0*z1 - b0*k0*o1*z0 - b0*k1*o0*z0 - b1*k0*o0*z0 + c0*k0*n0*z1 + c0*k0*n1*z0 + c0*k1*n0*z0 + c1*k0*n0*z0 - d0*m0*o0*v1 - d0*m0*o1*v0 - d0*m1*o0*v0 - d1*m0*o0*v0 + f0*k0*o0*v1 + f0*k0*o1*v0 + f0*k1*o0*v0 + f1*k0*o0*v0 + b0*h0*s0*z1 + b0*h0*s1*z0 + b0*h1*s0*z0 + b0*m0*o0*y1 + b0*m0*o1*y0 + b0*m1*o0*y0 + b1*h0*s0*z0 + b1*m0*o0*y0 - c0*g0*s0*z1 - c0*g0*s1*z0 - c0*g1*s0*z0 - c0*m0*n0*y1 - c0*m0*n1*y0 - c0*m1*n0*y0 - c1*g0*s0*z0 - c1*m0*n0*y0 + d0*h0*u0*v1 + d0*h0*u1*v0 + d0*h1*u0*v0 + d0*m0*n0*x1 + d0*m0*n1*x0 + d0*m1*n0*x0 + d1*h0*u0*v0 + d1*m0*n0*x0 - f0*h0*s0*v1 - f0*h0*s1*v0 - f0*h1*s0*v0 - f0*k0*n0*x1 - f0*k0*n1*x0 - f0*k1*n0*x0 - f1*h0*s0*v0 - f1*k0*n0*x0 - b0*h0*u0*y1 - b0*h0*u1*y0 - b0*h1*u0*y0 - b1*h0*u0*y0 + c0*g0*u0*y1 + c0*g0*u1*y0 + c0*g1*u0*y0 + c1*g0*u0*y0 - d0*g0*u0*x1 - d0*g0*u1*x0 - d0*g1*u0*x0 - d1*g0*u0*x0 + f0*g0*s0*x1 + f0*g0*s1*x0 + f0*g1*s0*x0 + f1*g0*s0*x0 - c0*k0*u0*v1 - c0*k0*u1*v0 - c0*k1*u0*v0 + c0*m0*s0*v1 + c0*m0*s1*v0 + c0*m1*s0*v0 - c1*k0*u0*v0 + c1*m0*s0*v0 + b0*k0*u0*x1 + b0*k0*u1*x0 + b0*k1*u0*x0 - b0*m0*s0*x1 - b0*m0*s1*x0 - b0*m1*s0*x0 + b1*k0*u0*x0 - b1*m0*s0*x0;
	f_qAll[8] = d0*g0*o0*z0 - d0*h0*n0*z0 - f0*g0*o0*y0 + f0*h0*n0*y0 - b0*k0*o0*z0 + c0*k0*n0*z0 - d0*m0*o0*v0 + f0*k0*o0*v0 + b0*h0*s0*z0 + b0*m0*o0*y0 - c0*g0*s0*z0 - c0*m0*n0*y0 + d0*h0*u0*v0 + d0*m0*n0*x0 - f0*h0*s0*v0 - f0*k0*n0*x0 - b0*h0*u0*y0 + c0*g0*u0*y0 - d0*g0*u0*x0 + f0*g0*s0*x0 - c0*k0*u0*v0 + c0*m0*s0*v0 + b0*k0*u0*x0 - b0*m0*s0*x0;

	mod_factor_order8(f_qAll, s_equation1_forqy);

	C[0] = b2; C[1] = b1; C[2] = b0;
	C[3] = c2; C[4] = c1; C[5] = c0;
	C[6] = d2; C[7] = d1; C[8] = d0;
	C[9] = f2; C[10] = f1; C[11] = f0;
	C[12] = g2; C[13] = g1; C[14] = g0;
	C[15] = h2; C[16] = h1; C[17] = h0;
	C[18] = k2; C[19] = k1; C[20] = k0;
	C[21] = m2; C[22] = m1; C[23] = m0;
	C[24] = n2; C[25] = n1; C[26] = n0;
	C[27] = o2; C[28] = o1; C[29] = o0;
	C[30] = s2; C[31] = s1; C[32] = s0;
	C[33] = u2; C[34] = u1; C[35] = u0;
	C[36] = v2; C[37] = v1; C[38] = v0;
	C[39] = x2; C[40] = x1; C[41] = x0;
	C[42] = y2; C[43] = y1; C[44] = y0;
	C[45] = z2; C[46] = z1; C[47] = z0;

	return;
}

void solver_ac_pose4d(std::vector<Eigen::Matrix3d>& Rf1tof2_recover,
	std::vector<Eigen::Vector3d>& Tf1tof2_recover,
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *Timu1All, double *Timu2All)
{
	Rf1tof2_recover.clear();
	Tf1tof2_recover.clear();

	double *C = new double[48];
	double *a = new double[7];

	f_multicamera_Ev_solver(C, a, Line_iAll, Line_jAll, P1All, P2All, AtempAll, RiAll, Timu1All, Timu2All);

	Eigen::Matrix<double, 16, 3> C_mat;
	int idx = 0;
	for (int i = 0; i < 16; i++){
		for (int j = 0; j < 3; j++){
			C_mat(i, j) = C[idx];
			idx++;
		}
	}

	if (abs(a[0]) < NEAR_ZERO_THRESHOLD){
		std::cout << "warning: the coefficient of highest order is 0!" << std::endl;
	}

	Eigen::Matrix<double, 6, 6> A;
	A << -a[1] / a[0], -a[2] / a[0], -a[3] / a[0], -a[4] / a[0], -a[5] / a[0], -a[6] / a[0],
		1, 0, 0, 0, 0, 0,
		0, 1, 0, 0, 0, 0,
		0, 0, 1, 0, 0, 0,
		0, 0, 0, 1, 0, 0,
		0, 0, 0, 0, 1, 0;
    
	Eigen::EigenSolver<Eigen::Matrix<double, 6, 6>> es(A);
	Eigen::ArrayXcd D = es.eigenvalues();

	for (int i = 0; i < D.rows(); i++)
	{
		double re = D[i].real();
		double im = D[i].imag();

		if (abs(im) > NEAR_ZERO_THRESHOLD)
			continue;
		
		Eigen::Vector3d q;
		q << re*re, re, 1;

		Eigen::Matrix<double, 16, 1> t = C_mat * q;
		Eigen::Matrix<double, 4, 4> Cq_recover;
		Cq_recover << 
			t(0), t(1), t(2), t(3),
			t(4), t(5), t(6), t(7),
			t(8), t(9), t(10), t(11),
			t(12), t(13), t(14), t(15);

		Eigen::Matrix<double, 4, 3> M = Cq_recover.block(0, 0, 4, 3);
		Eigen::Matrix<double, 4, 1> b = Cq_recover.block(0, 3, 4, 1);
		Eigen::Vector3d T_sol = -M.householderQr().solve(b);

		Eigen::Matrix3d R_sol;
		R_sol <<
			1 - re*re, 0, -2*re,
			0, 1 + re*re, 0,
			2*re, 0, 1 - re*re;
		R_sol = R_sol * (1.0/(1.0+re*re));

		Rf1tof2_recover.push_back(R_sol);
		Tf1tof2_recover.push_back(T_sol);
	}



	delete[] C;
	delete[] a;
	return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{

    if (nrhs < 8)
        mexErrMsgTxt("at least 8 inputs are required!");
    if (nlhs < 1)
        mexErrMsgTxt("at least 2 output is required!");

    if (mxIsEmpty(prhs[0]) || mxIsEmpty(prhs[1]) || mxIsEmpty(prhs[2]) || mxIsEmpty(prhs[3]) || mxIsEmpty(prhs[4])||
        mxIsEmpty(prhs[5])||mxIsEmpty(prhs[6])||mxIsEmpty(prhs[7]))
        mexErrMsgTxt("input parameter should not be an empty array!");

    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])
        || !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
        || !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])
        || !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])
        || !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4]))
    {
        mexErrMsgIdAndTxt("JPL:dp_intra_2ac:notDouble", "Input data must be type double.");
    }

    double *Line_iAll;
    double *Line_jAll;
    double *P1All;
    double *P2All;
    double *AtempAll;
    double *RiAll;
    double *Timu1All;
    double *Timu2All;

    Line_iAll=(double *)mxGetData(prhs[0]);
    Line_jAll=(double *)mxGetData(prhs[1]);
    P1All=(double *)mxGetData(prhs[2]);
    P2All=(double *)mxGetData(prhs[3]);
    AtempAll=(double *)mxGetData(prhs[4]);
    RiAll=(double *)mxGetData(prhs[5]);
    Timu1All=(double *)mxGetData(prhs[6]);
    Timu2All=(double *)mxGetData(prhs[7]);



    double *R_real_sols,*T_real_sols;
    std::vector<Eigen::Matrix3d> Rf1tof2_recover;
    std::vector<Eigen::Vector3d> Tf1tof2_recover;
    solver_ac_pose4d(Rf1tof2_recover,Tf1tof2_recover,Line_iAll, Line_jAll, P1All, P2All, AtempAll, RiAll, Timu1All, Timu2All);

    if (Rf1tof2_recover.size() > 0)
    {
        mwSize dims_R[3]={3,3,Rf1tof2_recover.size()};
        plhs[0]=mxCreateNumericArray(3,dims_R,mxDOUBLE_CLASS, mxREAL);
        plhs[1]=mxCreateDoubleMatrix(3, Rf1tof2_recover.size(), mxREAL);

        R_real_sols = mxGetPr(plhs[0]);
        T_real_sols = mxGetPr(plhs[1]);

        for (int i = 0; i < Rf1tof2_recover.size(); i++)
        {
            Eigen::Matrix<double,3,3> R=Rf1tof2_recover[i];
            Eigen::Matrix<double,3,1> T=Tf1tof2_recover[i];

            T_real_sols[i*3] = T(0);
            T_real_sols[i*3+1] = T(1);
            T_real_sols[i*3+2] = T(2);

            R_real_sols[i*9] = R(0,0);
            R_real_sols[i*9+1] = R(1,0);
            R_real_sols[i*9+2] = R(2,0);
            R_real_sols[i*9+3] = R(0,1);
            R_real_sols[i*9+4] = R(1,1);
            R_real_sols[i*9+5] = R(2,1);
            R_real_sols[i*9+6] = R(0,2);
            R_real_sols[i*9+7] = R(1,2);
            R_real_sols[i*9+8] = R(2,2);
        }
    }
    else
    {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(0, 0, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
}
