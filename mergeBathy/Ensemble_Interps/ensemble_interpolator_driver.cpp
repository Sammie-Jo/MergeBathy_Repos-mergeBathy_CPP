//Ensemble Interpolation Driver
#include <iostream>
#include <ostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <Windows.h>


#include "numerical_recipes_interp.h"
#include "user_interface.h"
#include "computational_geometry_algo_interp.h"

int write_output_file(char*,const std::vector<double>&,const std::vector<double>&,const std::vector<double>&);
int read_input_file(char*,std::vector<double>&,std::vector<double>&,std::vector<double> &);

//Input grid
std::vector<double> grid_xvi;
std::vector<double> grid_yvi;
std::vector<double> grid_zvi;

//Output grid
std::vector<double> grid_xv;
std::vector<double> grid_yv;
std::vector<double> grid_zv;


int main(int argc, char* argv[])
{
	std::cout << "Ensemble Interpolation Driver\n\n";
	int         status = 0;
	const int	MAX_SIZE = 255;
	int			nPoints	 = 0;
	char		chAbsFilenameIn[MAX_SIZE];
	char        chAbsFilenameOut[MAX_SIZE];
	
	clock_t     time_start		= 0L;
	clock_t     time_finish		= 0L;
	double      time_duration	= 0L;
	
	std::vector<double> xv;
	std::vector<double> yv;
	std::vector<double> zv;
	    
	double grid_dx;
	double grid_dy;
	double grid_dx_in = 0.0;
	double grid_dy_in = 0.0;
	int nCategoryChoice = 0;
	int polyinterp_order= 0;
	
	//Get filename from user
	std::cout << "Enter xyz data file : ";
	std::cin.getline(chAbsFilenameIn,MAX_SIZE);// Test file -> C:\\CGAL-3.8\\examples\\Interpolation\\data\\points3

	//Read in user file
	status = read_input_file(chAbsFilenameIn,xv,yv,zv);
	if (status == -1) {	return status; }

	//Determine if data is regular or irregular
	nCategoryChoice = enter_data_category();
	if (nCategoryChoice == 1)
	{
		//User input file is regular(grid)
		enter_inputgrid_spacing(grid_dx_in,grid_dy_in);
		status = build_input_grid_locs(xv,yv,zv,grid_xvi,grid_yvi,grid_zvi,grid_dx_in,grid_dy_in);
		if (status == -1)
			return status;
	}
	else if (nCategoryChoice == 2)
	{
		//User input file is irregular(scattered)
		//No action at this time
	}
	else
	{
		std::cout << "Invalid data category...\n";
		status = -1;
		return status;
	}
	
	//Get ouput grid specifications
	std::cout << "Enter output grid spacing [dx] : ";
	std::cin  >> grid_dx;
	std::cout << "Enter output grid spacing [dy] : ";
	std::cin  >> grid_dy;

	status = build_output_grid_locs(xv,yv,grid_dx,grid_dy);
	if (status == -1) {	return status; }

	//If irregular grid must create Delaunay Triangulation(CGAL)
	if (nCategoryChoice == 2)
	{
		status = build_delaunay_triangulation(xv,yv,zv);
		if (status == -1) {	return status; }

		status = build_function_values(xv,yv,zv);
		if (status == -1) {	return status; }
	}
	
	int nInterpChoice = 0;
	while ((nInterpChoice = enter_interp_choice()) != 7)
	{
		time_start	= clock();
		switch (nInterpChoice)
		{
			case 1: //Sibson; Natural Neighbor Interpolation[Sib81];
				status = sibson_natural_neighbor_interp();
				break;
			case 2: //Sibson Continuous; Natural Neighbor Interpolation w/ Gradient Est[Sib81]
				status = sibson_natural_neighbor_cont_interp();
				break;
			case 3: //Farin Continuous; [Far90]
				status = farin_natural_neighbor_cont_interp();
				break;
			case 4: //Bilinear Interpolation; Numerical Recipes
				if (nCategoryChoice == 2)
				{
					std::cout << "Bilinear interpolation currently not supported for irregular grids...\n" ;
					status = -1;
					return status;
				}
				status =  numrecip_reggrid_bilinear();
				break;
			case 5://Higher Order (Accuracy) Polynomial; Numerical Recipes
				if (nCategoryChoice == 2)
				{
					std::cout << "Higher Order Polynomial interpolation currently not supported for irregular grids...\n" ;
					status = -1;
					return status;
				}
				//Get order of interpolation
				std::cout << "Enter order for polynomial interpolation [>=2] : ";
				std::cin  >> polyinterp_order;
				status = numrecip_reggrid_poly2d(polyinterp_order);
				break;
			case 6://Bicubic spline; Numerical Recipes
				if (nCategoryChoice == 2)
				{
					std::cout << "Bicubic spline interpolation currently not supported for irregular grids...\n" ;
					status = -1;
					return status;
				}
				status = numrecip_reggrid_bicubicspline();
				break;
			default:
				break;
		}

		if ((nInterpChoice >= 1) && (nInterpChoice <= 6))
		{
    		time_finish = clock();
			time_duration = (double)(time_finish - time_start) / (double)CLOCKS_PER_SEC;
			std::cout << "Duration  = " << time_duration << " sec. \n";

			//Export output file
			std::cout.flush();
			std::cout << "Enter output file name : ";
			std::cin >> chAbsFilenameOut;
			status = write_output_file(chAbsFilenameOut,grid_xv,grid_yv,grid_zv);
		}
	}
	
	std::cout << "Terminating...\n";
	
	return 0;
}


int read_input_file(char* chAbsFilename,std::vector<double>& xvect, std::vector<double>& yvect,std::vector<double> &zvect)
{
	int status = 0;

	try
	{
		std::cout << "Reading: " << chAbsFilename << " "; 
		std::ifstream inFile(chAbsFilename,std::ios::in);
		while (!inFile.eof())
		{
			 double x,y,z;

			 //Need to skip non-numeric lines.
			 inFile >> x >> y >> z;
			 xvect.push_back(x);
			 yvect.push_back(y);
			 zvect.push_back(z);
		}
		inFile.close();
		status = 1;
		std::cout << "(Done)\n"; 
	}
	catch(...)
	{
		status = -1;
		std::cout << "Exception in loadInputFile()\n";
	}

	return status;
}

int write_output_file(char* chAbsFilename,					\
					  const std::vector<double>& grid_xv,	\
					  const std::vector<double>& grid_yv,	\
					  const std::vector<double>& grid_zv)
{
	int status = 0;
	try
	{
		std::cout << "Writing: " << chAbsFilename << " "; 
		int xcnt   = (int)grid_xv.size();
		int ycnt   = (int)grid_yv.size();
		int zcnt   = xcnt*ycnt;
		std::cout << "Writing " << zcnt << " values to file " << chAbsFilename << " "; 

		std::ofstream outFile;
		outFile.open(chAbsFilename,std::ios_base::out);
		for (int xndx = 0; xndx < xcnt; xndx++)
		{
			for (int yndx = 0; yndx < ycnt; yndx++)
			{
				outFile << grid_xv[xndx] << " " << grid_yv[yndx] << " " << grid_zv[yndx + (ycnt*xndx)] << "\n";
			}
		}
		outFile.close();
		std::cout << "(Done)\n"; 
	}
	catch(...)
	{
		status = -1;
		std::cout << "Exception in write_output_file()\n";
	}

	status = 1;

	return status;
}