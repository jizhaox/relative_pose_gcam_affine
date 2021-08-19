#include <time.h>
#include <fstream>
#include <iostream>

#include <Eigen/Core>

#include "solver_ac_pose3d_1AC_plane.h"


using namespace std;
using namespace Eigen;

bool load_input_data3d(string fname, 
	double *Line_iAll, double *Line_jAll, double *P1All, double *P2All,
	double *AtempAll, double *RiAll, double *TiAll)
{
	ifstream fin;
	fin.open(fname.c_str(), ios::in | ios::binary);
	if (!fin.is_open()){
		std::cout << "can not open file!" << std::endl;
		return false;
	}

	double a;
	std::cout << "input parameters" << std::endl;
	std::cout << "Line_iAll:" << std::endl;
	for (int i = 0; i < 6; i++){
		fin.read((char*)&a, sizeof(double));
		Line_iAll[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "Line_jAll:" << std::endl;
	for (int i = 0; i < 6; i++){
		fin.read((char*)&a, sizeof(double));
		Line_jAll[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "P1All:" << std::endl;
	for (int i = 0; i < 3; i++){
		fin.read((char*)&a, sizeof(double));
		P1All[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "P2All:" << std::endl;
	for (int i = 0; i < 3; i++){
		fin.read((char*)&a, sizeof(double));
		P2All[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "AtempAll:" << std::endl;
	for (int i = 0; i < 9; i++){
		fin.read((char*)&a, sizeof(double));
		AtempAll[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "RiAll:" << std::endl;
	for (int i = 0; i < 9; i++){
		fin.read((char*)&a, sizeof(double));
		RiAll[i] = a;
		std::cout << a << " ";
	}

	std::cout << std::endl << std::endl << "TiAll:" << std::endl;
	for (int i = 0; i < 6; i++){
		fin.read((char*)&a, sizeof(double));
		TiAll[i] = a;
		std::cout << a << " ";
	}

	return true;
}


int main(int argc, char* argv[]){
	if (argc < 2)
	{
		std::cout << "USAGE: test3d_1AC <INPUT_PARAMETER_FILE>" << std::endl;
		return -1;
	}

	double *Line_iAll = new double[6];
	double *Line_jAll = new double[6];
	double *P1All = new double[3];
	double *P2All = new double[3];
	double *AtempAll = new double[9];
	double *RiAll = new double[9];
	double *TiAll = new double[6];
	
	// load data
	load_input_data3d(argv[1], Line_iAll, Line_jAll, P1All, P2All, AtempAll, RiAll, TiAll);
	
	// Relative Pose Estimation Under Planar Motion: 2AC plane solver
	std::vector<Eigen::Matrix3d> Rf1tof2_recover;
	std::vector<Eigen::Vector3d> Tf1tof2_recover;
	solver_ac_pose3d_1AC_plane(Rf1tof2_recover, Tf1tof2_recover, Line_iAll, Line_jAll, P1All, P2All, AtempAll, RiAll, TiAll);
	
	//finished
	for (int i = 0; i < Rf1tof2_recover.size(); i++){
		std::cout << std::endl << "------------------------" << std::endl;
		std::cout << "solution #" << i << std::endl;
		std::cout << "rotation:" << std::endl;
		std::cout << Rf1tof2_recover[i] << std::endl;
		std::cout << std::endl << "translation:" << std::endl;
		std::cout << Tf1tof2_recover[i] << std::endl;
	}

	delete[] Line_iAll;
	delete[] Line_jAll;
	delete[] P1All;
	delete[] P2All;
	delete[] AtempAll;
	delete[] RiAll;
	delete[] TiAll;

	return 0;
}
