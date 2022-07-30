#include <time.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>

#include <Eigen/Core>

#include "solver_ac_pose3d_1AC_plane.h"


using namespace std;
using namespace Eigen;

bool load_input_data3d(string fname, 
	double *Image1, double *Image2, double *At_33, double *R_cam,
	double *t_cam, double *R_gt, double *t_gt)
{
	ifstream fin;
	fin.open(fname.c_str(), ios::in | ios::binary);
	if (!fin.is_open()){
		std::cout << "can not open file!" << std::endl;
		return false;
	}

	double a;
	std::cout << "input parameters" << std::endl;
	std::cout << "Image1:" << std::endl;
	for (int i = 0; i < 3; i++){
		fin.read((char*)&a, sizeof(double));
		Image1[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "Image2:" << std::endl;
	for (int i = 0; i < 3; i++){
		fin.read((char*)&a, sizeof(double));
		Image2[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "At_33:" << std::endl;
	for (int i = 0; i < 9; i++){
		fin.read((char*)&a, sizeof(double));
		At_33[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "R_cam:" << std::endl;
	for (int i = 0; i < 18; i++){
		fin.read((char*)&a, sizeof(double));
		R_cam[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "t_cam:" << std::endl;
	for (int i = 0; i < 6; i++){
		fin.read((char*)&a, sizeof(double));
		t_cam[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "R_gt:" << std::endl;
	for (int i = 0; i < 9; i++){
		fin.read((char*)&a, sizeof(double));
		R_gt[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "t_gt:" << std::endl;
	for (int i = 0; i < 3; i++){
		fin.read((char*)&a, sizeof(double));
		t_gt[i] = a;
		std::cout << a << " ";
	}
    std::cout << std::endl;
	return true;
}


int main(int argc, char* argv[]){
	if (argc < 2)
	{
		std::cout << "USAGE: test3d_1AC <INPUT_PARAMETER_FILE>" << std::endl;
		return -1;
	}

	double *Image1 = new double[3];
	double *Image2 = new double[3];
	double *At_33 = new double[9];
	double *R_cam = new double[18];
	double *t_cam = new double[6];
	double *R_gt = new double[9];
	double *t_gt = new double[3];
	
	// load data
	load_input_data3d(argv[1], Image1, Image2, At_33, R_cam, t_cam, R_gt, t_gt);
	
	// Relative Pose Estimation Under Planar Motion: 1AC plane solver
	std::vector<Eigen::Matrix3d> Rf1tof2_recover;
	std::vector<Eigen::Vector3d> Tf1tof2_recover;
    
        solver_ac_pose3d_1AC_plane(Rf1tof2_recover, Tf1tof2_recover, Image1, Image2, At_33, R_cam, t_cam, "inter");
    
	//finished
	for (int i = 0; i < Rf1tof2_recover.size(); i++){
		std::cout << std::endl << "------------------------" << std::endl;
		std::cout << "solution #" << i << std::endl;
		std::cout << "rotation:" << std::endl;
		std::cout << Rf1tof2_recover[i] << std::endl;
		std::cout << std::endl << "translation:" << std::endl;
		std::cout << Tf1tof2_recover[i] << std::endl;
	}
    
    //Ground truth
    std::cout << std::endl << "rotation_gt:" << std::endl;
    for(int i = 0; i < 3 ;i++)
        std::cout << R_gt[i] << " " << R_gt[i+3]<< " " << R_gt[i+6]<< std::endl;
	std::cout << std::endl << "translation_gt:" << std::endl;
    for(int i = 0; i < 3 ;i++)
        std::cout << t_gt[i] << std::endl;

	delete[] Image1;
	delete[] Image2;
	delete[] At_33;
	delete[] R_cam;
	delete[] t_cam;
	delete[] R_gt;
	delete[] t_gt;

	return 0;
}
