//User IO

#include <iostream>
#include <ostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <Windows.h>

#include "user_interface.h"


int enter_data_category()
{
	std::cout << "\t(1) Input data is regular(gridded).\n";
	std::cout << "\t(2) Input data is irregular(scattered).\n";
	
	int choice = 0;
	std::cin >> choice;

	std::cout << "\n\n";

	return choice;
}

int enter_interp_choice()
{
	std::cout << "*********************************************" << "\n";
	std::cout << "*    Ensemble Interpolation Algorithms      *" << "\n";
	std::cout << "*********************************************" << "\n";

	std::cout << "\t(1) Sibson; Natural Neighbor Interpolation (CGAL) [Sib81]" << "\n";
	std::cout << "\t(2) Sibson Continuous; Natural Neighbor Interpolation w/ Gradient (CGAL)[Sib81]" << "\n";
	std::cout << "\t(3) Farin Continuous; Natural Neighbor Interpolation w/ Gradient (CGAL) [Far90]" << "\n";
	std::cout << "\t(4) Bilinear Interpolation (NR)[Press 2007]" << "\n";
	std::cout << "\t(5) Polynomial Interpolation (NR)[Press 2007]" << "\n"; 
	std::cout << "\t(6) Bicubic Spline Interpolation (NR)[Press 2007]" << "\n";
	std::cout << "\t(7) Exit\n";
	std::cout << "Enter Choice: ";
	
	int choice = 0;
	std::cin >> choice;

	std::cout << "\n\n";

	return choice;
}

int enter_inputgrid_spacing(double& grid_dx_in, double& grid_dy_in)
{
	int status = 0;

	std::cout << "Input grid spacing [dx] : ";
	std::cin  >> grid_dx_in;
	std::cout << "Input grid spacing [dy] : ";
	std::cin  >> grid_dy_in;

	status = 1;

	return status;
}