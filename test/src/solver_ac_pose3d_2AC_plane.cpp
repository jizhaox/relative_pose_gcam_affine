#include "solver_ac_pose3d_2AC_plane.h"

void mod_factor_order6_2AC_plane(double *e, double *quot)
{
    quot[0] = e[0];
    quot[1] = e[1];
    quot[2] = e[2] - quot[0];
    quot[3] = e[3] - quot[1];
    quot[4] = e[4] - quot[2];
    return;
}

void format_convert(
        double *input_Image_1, double *input_Image_2, double *input_affine_tran,
        double *extrinsic_R_camera, double *extrinsic_T_camera,
        std::vector<Eigen::Vector3d>& Image1, std::vector<Eigen::Vector3d>& Image2, std::vector<Eigen::Matrix3d>& Ac,
        std::vector<Eigen::Matrix3d>& R_camera, std::vector<Eigen::Vector3d>& T_camera)
{
    // data format convertion
    Eigen::Matrix3d R;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                // Matlab uses column-first order
                R(i, j) = extrinsic_R_camera[k * 9 + j * 3 + i];
            }
        }
        R_camera.push_back(R);
    }
    Eigen::Matrix3d A;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                // Matlab uses column-first order
                A(i, j) = input_affine_tran[k * 9 + j * 3 + i];
            }
        }
        Ac.push_back(A);
    }

    Eigen::Vector3d T;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            T(i) = extrinsic_T_camera[k * 3 + i];
        }
        T_camera.push_back(T);
    }
    Eigen::Vector3d X1, X2;
    for (int k = 0; k < 2; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            X1(i) = input_Image_1[k * 3 + i];
            X2(i) = input_Image_2[k * 3 + i];
        }
        Image1.push_back(X1);
        Image2.push_back(X2);
    }
    return;
}

void f_multicamera_Planemotion_solver_2AC_plane(double *C, double *s_equation1_forqy,
                                                double *input_Image_1, double *input_Image_2,
                                                double *input_Ac,  double *extrinsic_R_camera, double *extrinsic_T_camera, char *match_type)
{
    char const *inter="inter";
    char const *intra="intra";

    std::vector<Eigen::Matrix3d> R_camera, Ac;
    std::vector<Eigen::Vector3d> T_camera, Image1, Image2;

    format_convert(input_Image_1, input_Image_2, input_Ac, extrinsic_R_camera, extrinsic_T_camera,
                   Image1, Image2, Ac, R_camera, T_camera);

    double f1_New_C11All[2][3];
    double f1_New_C12All[2][3];
    double f1_New_C13All[2][3];

    double f2_New_C21All[2][3];
    double f2_New_C22All[2][3];
    double f2_New_C23All[2][3];

    double f_qAll[7];

    int point_num = 2;
    for (int i = 0; i < point_num; i++)
    {
        Eigen::Vector3d P1;
        Eigen::Vector3d P2;
        P1 = Image1[i];
        P2 = Image2[i];

        Eigen::Vector3d U1, U2;
        U1 = P1.normalized();
        U2 = P2.normalized();

        Eigen::Vector3d T1 = T_camera[i];
        Eigen::Matrix3d R1 = R_camera[i];
        Eigen::Vector3d T2;
        Eigen::Matrix3d R2,Atemp;
        if(strcmp(match_type,inter)==0)
        {
            if(i==0)
            {
                T2 = T_camera[i+1];
                R2 = R_camera[i+1];
                Atemp = Ac[i].transpose() * R2.transpose();
            } else
            {
                T2 = T_camera[i-1];
                R2 = R_camera[i-1];
                Atemp = Ac[i].transpose() * R2.transpose();
            }

        }
        else if(strcmp(match_type,intra)==0)
        {
            T2 = T_camera[i];
            R2 = R_camera[i];
            Atemp = Ac[i].transpose() * R2.transpose();
        }
        else
        {
            //mexErrMsgTxt("unsupported match type.");
            return;
        }

        Eigen::Vector3d R1U1 = R1 * U1;
        Eigen::Vector3d T1xR1U1 = T1.cross(R1U1);
        Eigen::Vector3d R2U2 = R2 * U2;
        Eigen::Vector3d T2xR2U2 = T2.cross(R2U2);

        Eigen::Matrix<double, 6, 1> Line_i, Line_j;
        Line_i.block(0, 0, 3, 1) = R1U1;
        Line_i.block(3, 0, 3, 1) = T1xR1U1;
        Line_j.block(0, 0, 3, 1) = R2U2;
        Line_j.block(3, 0, 3, 1) = T2xR2U2;

        Eigen::Matrix<double, 3, 2>Ti;
        Ti.block(0, 0, 3, 1) = T1;
        Ti.block(0, 1, 3, 1) = T2;

        Eigen::Vector3d P1All = R1 * P1;
        Eigen::Vector3d P2All = R2 * P2;


        double L11 = Line_i(0);
        double L12 = Line_i(1);
        double L13 = Line_i(2);
        double L14 = Line_i(3);
        double L15 = Line_i(4);
        double L16 = Line_i(5);

        double L21 = Line_j(0);
        double L22 = Line_j(1);
        double L23 = Line_j(2);
        double L24 = Line_j(3);
        double L25 = Line_j(4);
        double L26 = Line_j(5);


        double r1 = R1(0,0);
        double r2 = R1(0,1);
        double r3 = R1(0,2);
        double r4 = R1(1,0);
        double r5 = R1(1,1);
        double r6 = R1(1,2);
        double r7 = R1(2,0);
        double r8 = R1(2,1);
        double r9 = R1(2,2);

        double tx1 = T1(0);
        double ty1 = T1(1);
        double tz1 = T1(2);

        double tx2 = T2(0);
        double ty2 = T2(1);
        double tz2 = T2(2);

        double p11 = P1All(0);
        double p12 = P1All(1);
        double p13 = P1All(2);

        double p21 = P2All(0);
        double p22 = P2All(1);
        double p23 = P2All(2);

        double a1 = Atemp(0,0);
        double a2 = Atemp(0,1);
        double a3 = Atemp(0,2);
        double a4 = Atemp(1,0);
        double a5 = Atemp(1,1);
        double a6 = Atemp(1,2);

        // bulid the coefficient matrix
        // f1
        f1_New_C11All[i][0] = L12*L23 + L13*L22;
        f1_New_C11All[i][1] = -2 * L11*L22;
        f1_New_C11All[i][2] = L12*L23 - L13*L22;

        f1_New_C12All[i][0] = -L11*L22 - L12*L21;
        f1_New_C12All[i][1] = -2 * L13*L22;
        f1_New_C12All[i][2] = L11*L22 - L12*L21;

        f1_New_C13All[i][0] = L12*L25 - L14*L21 - L11*L24 + L15*L22 - L13*L26 - L16*L23;
        f1_New_C13All[i][1] = 2 * L11*L26 - 2 * L13*L24 + 2 * L14*L23 - 2 * L16*L21;
        f1_New_C13All[i][2] = L11*L24 + L14*L21 + L12*L25 + L15*L22 + L13*L26 + L16*L23;

        // f2
        f2_New_C21All[i][0] = -a2*p13 - a3*p12 - p23*r4 - p22*r7;
        f2_New_C21All[i][1] = 2 * a2*p11 + 2 * p22*r1;
        f2_New_C21All[i][2] = a2*p13 - a3*p12 - p23*r4 + p22*r7;

        f2_New_C22All[i][0] = a1*p12 + a2*p11 + p22*r1 + p21*r4;
        f2_New_C22All[i][1] = 2 * a2*p13 + 2 * p22*r7;
        f2_New_C22All[i][2] = a1*p12 - a2*p11 - p22*r1 + p21*r4;

        f2_New_C23All[i][0] = a2*p13*tx1 + a3*p12*tx1 + a2*p13*tx2 + a3*p12*tx2 + a1*p13*ty1 - a3*p11*ty1 - a1*p13*ty2 + a3*p11*ty2 - a1*p12*tz1 - a2*p11*tz1 - a1*p12*tz2 - a2*p11*tz2 + p23*r4*tx1 + p23*r4*tx2 + p22*r7*tx1 + p22*r7*tx2 - p23*r1*ty1 + p23*r1*ty2 + p21*r7*ty1 - p21*r7*ty2 - p22*r1*tz1 - p22*r1*tz2 - p21*r4*tz1 - p21*r4*tz2;
        f2_New_C23All[i][1] = 2*a1*p12*tx1 - 2*a2*p11*tx2 - 2*a1*p11*ty1 + 2*a1*p11*ty2 - 2*a3*p13*ty1 + 2*a3*p13*ty2 + 2*a3*p12*tz1 - 2*a2*p13*tz2 - 2*p22*r1*tx2 + 2*p21*r4*tx1 - 2*p21*r1*ty1 + 2*p21*r1*ty2 - 2*p23*r7*ty1 + 2*p23*r7*ty2 + 2*p23*r4*tz1 - 2*p22*r7*tz2;
        f2_New_C23All[i][2]=  a2*p13*tx1 - a3*p12*tx1 - a2*p13*tx2 + a3*p12*tx2 - a1*p13*ty1 + a3*p11*ty1 + a1*p13*ty2 - a3*p11*ty2 + a1*p12*tz1 - a2*p11*tz1 - a1*p12*tz2 + a2*p11*tz2 - p23*r4*tx1 + p23*r4*tx2 + p22*r7*tx1 - p22*r7*tx2 + p23*r1*ty1 - p23*r1*ty2 - p21*r7*ty1 + p21*r7*ty2 - p22*r1*tz1 + p22*r1*tz2 + p21*r4*tz1 - p21*r4*tz2;
    }

    // hidden variable resultant method
    double b3 = f1_New_C11All[0][0];
    double b4 = f1_New_C11All[0][1];
    double b5 = f1_New_C11All[0][2];

    double c3 = f1_New_C12All[0][0];
    double c4 = f1_New_C12All[0][1];
    double c5 = f1_New_C12All[0][2];


    double d3 = f1_New_C13All[0][0];
    double d4 = f1_New_C13All[0][1];
    double d5 = f1_New_C13All[0][2];


    double f3 = f2_New_C21All[0][0];
    double f4 = f2_New_C21All[0][1];
    double f5 = f2_New_C21All[0][2];


    double g3 = f2_New_C22All[0][0];
    double g4 = f2_New_C22All[0][1];
    double g5 = f2_New_C22All[0][2];

    double h3 = f2_New_C23All[0][0];
    double h4 = f2_New_C23All[0][1];
    double h5 = f2_New_C23All[0][2];

    double k3 = f1_New_C11All[1][0];
    double k4 = f1_New_C11All[1][1];
    double k5 = f1_New_C11All[1][2];



    double m3 = f1_New_C12All[1][0];
    double m4 = f1_New_C12All[1][1];
    double m5 = f1_New_C12All[1][2];

    double n3 = f1_New_C13All[1][0];
    double n4 = f1_New_C13All[1][1];
    double n5 = f1_New_C13All[1][2];

    f_qAll[0] = c3*h3*k3 - d3*g3*k3 + b3*g3*n3 - b3*h3*m3 - c3*f3*n3 + d3*f3*m3;
    f_qAll[1] = c3*h3*k4 + c3*h4*k3 + c4*h3*k3 - d3*g3*k4 - d3*g4*k3 - d4*g3*k3 + b3*g3*n4 + b3*g4*n3 - b3*h3*m4 - b3*h4*m3 + b4*g3*n3 - b4*h3*m3 - c3*f3*n4 - c3*f4*n3 - c4*f3*n3 + d3*f3*m4 + d3*f4*m3 + d4*f3*m3;
    f_qAll[2] = c3*h3*k5 + c3*h4*k4 + c3*h5*k3 + c4*h3*k4 + c4*h4*k3 + c5*h3*k3 - d3*g3*k5 - d3*g4*k4 - d3*g5*k3 - d4*g3*k4 - d4*g4*k3 - d5*g3*k3 + b3*g3*n5 + b3*g4*n4 + b3*g5*n3 - b3*h3*m5 - b3*h4*m4 - b3*h5*m3 + b4*g3*n4 + b4*g4*n3 - b4*h3*m4 - b4*h4*m3 + b5*g3*n3 - b5*h3*m3 - c3*f3*n5 - c3*f4*n4 - c3*f5*n3 - c4*f3*n4 - c4*f4*n3 - c5*f3*n3 + d3*f3*m5 + d3*f4*m4 + d3*f5*m3 + d4*f3*m4 + d4*f4*m3 + d5*f3*m3;
    f_qAll[3] = c3*h4*k5 + c3*h5*k4 + c4*h3*k5 + c4*h4*k4 + c4*h5*k3 + c5*h3*k4 + c5*h4*k3 - d3*g4*k5 - d3*g5*k4 - d4*g3*k5 - d4*g4*k4 - d4*g5*k3 - d5*g3*k4 - d5*g4*k3 + b3*g4*n5 + b3*g5*n4 - b3*h4*m5 - b3*h5*m4 + b4*g3*n5 + b4*g4*n4 + b4*g5*n3 - b4*h3*m5 - b4*h4*m4 - b4*h5*m3 + b5*g3*n4 + b5*g4*n3 - b5*h3*m4 - b5*h4*m3 - c3*f4*n5 - c3*f5*n4 - c4*f3*n5 - c4*f4*n4 - c4*f5*n3 - c5*f3*n4 - c5*f4*n3 + d3*f4*m5 + d3*f5*m4 + d4*f3*m5 + d4*f4*m4 + d4*f5*m3 + d5*f3*m4 + d5*f4*m3;
    f_qAll[4] = c3*h5*k5 + c4*h4*k5 + c4*h5*k4 + c5*h3*k5 + c5*h4*k4 + c5*h5*k3 - d3*g5*k5 - d4*g4*k5 - d4*g5*k4 - d5*g3*k5 - d5*g4*k4 - d5*g5*k3 + b3*g5*n5 - b3*h5*m5 + b4*g4*n5 + b4*g5*n4 - b4*h4*m5 - b4*h5*m4 + b5*g3*n5 + b5*g4*n4 + b5*g5*n3 - b5*h3*m5 - b5*h4*m4 - b5*h5*m3 - c3*f5*n5 - c4*f4*n5 - c4*f5*n4 - c5*f3*n5 - c5*f4*n4 - c5*f5*n3 + d3*f5*m5 + d4*f4*m5 + d4*f5*m4 + d5*f3*m5 + d5*f4*m4 + d5*f5*m3;
    f_qAll[5] = c4*h5*k5 + c5*h4*k5 + c5*h5*k4 - d4*g5*k5 - d5*g4*k5 - d5*g5*k4 + b4*g5*n5 - b4*h5*m5 + b5*g4*n5 + b5*g5*n4 - b5*h4*m5 - b5*h5*m4 - c4*f5*n5 - c5*f4*n5 - c5*f5*n4 + d4*f5*m5 + d5*f4*m5 + d5*f5*m4;
    f_qAll[6] = c5*h5*k5 - d5*g5*k5 + b5*g5*n5 - b5*h5*m5 - c5*f5*n5 + d5*f5*m5;

    mod_factor_order6_2AC_plane(f_qAll, s_equation1_forqy);

    C[0] = b3; C[1] = b4; C[2] = b5;
    C[3] = c3; C[4] = c4; C[5] = c5;
    C[6] = d3; C[7] = d4; C[8] = d5;
    C[9] = f3; C[10] = f4; C[11] = f5;
    C[12] = g3; C[13] = g4; C[14] = g5;
    C[15] = h3; C[16] = h4; C[17] = h5;
    C[18] = k3; C[19] = k4; C[20] = k5;
    C[21] = m3; C[22] = m4; C[23] = m5;
    C[24] = n3; C[25] = n4; C[26] = n5;

    return;
}

void solver_ac_pose3d_2AC_plane(std::vector<Eigen::Matrix3d>& Rf1tof2_recover,
                                std::vector<Eigen::Vector3d>& Tf1tof2_recover,
                                double *input_Image_1, double *input_Image_2,
                                double *input_Ac,  double *extrinsic_R_camera, double *extrinsic_T_camera, char *match_type)
{
    Rf1tof2_recover.clear();
    Tf1tof2_recover.clear();

    double *C = new double[27];
    double *a = new double[5];

    f_multicamera_Planemotion_solver_2AC_plane(C, a, input_Image_1, input_Image_2, input_Ac, extrinsic_R_camera, extrinsic_T_camera, match_type);

    Eigen::Matrix<double, 9, 3> C_mat;
    int idx = 0;
    for (int i = 0; i < 9; i++){
        for (int j = 0; j < 3; j++){
            C_mat(i, j) = C[idx];
            idx++;
        }
    }

    if (abs(a[0]) < NEAR_ZERO_THRESHOLD){
        std::cout << "warning: the coefficient of highest order is 0!" << std::endl;
    }

    double a3 = a[1] / a[0];
    double a2 = a[2] / a[0];
    double a1 = a[3] / a[0];
    double a0 = a[4] / a[0];

    complex<double> T1 = -a3 / 4.0;
    complex<double> T2 = a2*a2 - 3.0*a3*a1 + 12.0*a0;
    complex<double> T3 = (2.0*a2*a2*a2 - 9.0*a3*a2*a1 + 27.0*a1*a1 + 27.0*a3*a3*a0 - 72.0*a2*a0) / 2.0;
    complex<double> T4 = (-a3*a3*a3 + 4.0*a3*a2 - 8.0*a1) / 32.0;
    complex<double> T5 = (3.0*a3*a3 - 8.0*a2) / 48.0;
    complex<double> R1 = sqrt(T3*T3 - T2*T2*T2);
    complex<double> T3R1 = T3 + R1;
    complex<double> R2;
    if (T3R1.imag() == 0.0 && T3R1.real() < 0.0)
    {
        R2 = -pow(abs(T3R1), 1.0 / 3.0);
    }
    else
    {
        R2 = pow(T3R1, 1.0 / 3.0);
    }
    complex<double> R3 = (1.0 / 12.0)*(T2 / R2 + R2);
    complex<double> R4 = sqrt(T5 + R3);
    complex<double> R5 = 2.0 *T5 - R3;
    complex<double> R6 = T4 / R4;
    if ((T4 == 0.0) && (T5 == 0.0) && (abs(R3) < 1e-16))
        R6 = 1.0;

    complex<double> r[4];
    r[0] = T1 - R4 - sqrt(R5 - R6);
    r[1] = T1 - R4 + sqrt(R5 - R6);
    r[2] = T1 + R4 - sqrt(R5 + R6);
    r[3] = T1 + R4 + sqrt(R5 + R6);

    Eigen::ArrayXcd D = Eigen::ArrayXcd::Constant(4, 1, 0.0);
    D << r[0], r[1], r[2], r[3];

    for (int i = 0; i < D.rows(); i++)
    {
        double re = D[i].real();
        double im = D[i].imag();

        if (abs(im) > NEAR_ZERO_THRESHOLD)
            continue;

        Eigen::Vector3d q;
        q << re*re, re, 1;

        Eigen::Matrix<double, 9, 1> t = C_mat * q;
        Eigen::Matrix<double, 3, 3> Cq_recover;
        Cq_recover <<
                   t(0), t(1), t(2),
                t(3), t(4), t(5),
                t(6), t(7), t(8);

        Eigen::Matrix<double, 3, 2> M = Cq_recover.block(0, 0, 3, 2);
        Eigen::Matrix<double, 3, 1> b = Cq_recover.block(0, 2, 3, 1);
        Eigen::Vector2d txtz_recover = -M.householderQr().solve(b);

        Eigen::Vector3d T_sol;
        T_sol << txtz_recover(0), 0, txtz_recover(1);
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