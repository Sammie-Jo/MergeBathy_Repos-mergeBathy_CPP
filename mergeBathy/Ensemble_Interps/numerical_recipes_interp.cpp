//Numerical Recipes Interpolation function
//Num. Recipes

//Standard Template Library
#include <iostream>
#include <ostream>
#include <fstream>
#include <time.h>
#include <math.h>

//Declarations
#include "numerical_recipes_interp.h"

//Implementations
#include "linbcg.h"
#include "ludcmp.h"
#include "krig.h"
#include "interp_1d.h"
#include "polcoef.h"
#include "interp_linear.h"
#include "interp_2d.h"
#include "interp_rbf.h"
#include "interp_curve.h"
#include "interp_laplace.h"

//Input grid
extern std::vector<double> grid_xvi;
extern std::vector<double> grid_yvi;
extern std::vector<double> grid_zvi;

//Output grid
extern std::vector<double> grid_xv;
extern std::vector<double> grid_yv;
extern std::vector<double> grid_zv;

//The following functions exist to handle the different representations of vectors and matrices
//between MATLAB and Numerical Recipes.
int convertMatlabVectToNumRecipMatrix(std::vector<double> &MatlabVector, MatDoub& NumRcpMatrix)
{
	int status			= 0;
	int nMatrixROWS		= NumRcpMatrix.nrows();
	int nMatrixCOLS		= NumRcpMatrix.ncols();
	int nMatrixElements = nMatrixROWS*nMatrixCOLS;
	int nVectorElements = (int)MatlabVector.size();
	
	if (nMatrixElements == nVectorElements)
	{
		for (int col_ndx = 0; col_ndx < nMatrixCOLS; col_ndx++)	
		{
			for (int row_ndx = 0; row_ndx < nMatrixROWS; row_ndx++)			
			{
				int matlab_index = row_ndx + (nMatrixROWS*col_ndx);
				NumRcpMatrix[row_ndx][col_ndx] = MatlabVector[matlab_index];
			}
		}
	}
	status = 1;

	return status;
}


std::vector<double>* convertNumRecipMatrixToMatlabVect(MatDoub &nr_matrix, int numROWS, int numCOLS)
{
	int nMatrixElements = numROWS*numCOLS;
	int nVectorElements = nMatrixElements;
	
	std::vector<double>* pMatlabVector = new std::vector<double>(nVectorElements);
		
	for (int row_ndx = 0; row_ndx < numROWS; row_ndx++)
	{
		for (int col_ndx = 0; col_ndx < numCOLS; col_ndx++)	
		{
			int matlab_index = row_ndx + (numROWS*col_ndx);
			pMatlabVector->at(matlab_index) = nr_matrix[row_ndx][col_ndx];
		}
	}
		
	return pMatlabVector;
}

//Higher order interpolation for smoothness
int numrecip_reggrid_bicubicspline()
{
	std::cout << "Executing Numerical Recipes Bicubic spline Interpolation(Single Thread)...\n";
	
	//Input grid
	int xcnti          = (int)grid_xvi.size();
	int ycnti          = (int)grid_yvi.size();

	//Output grid
	int status			= 0;
	int xcnt			= (int)grid_xv.size();
	int ycnt			= (int)grid_yv.size();
	int zcnt			= xcnt*ycnt;
	grid_zv				= std::vector<double>(zcnt);
	
	///////////////////////////////////////////////////////
	//Numerical Recipes interface 
	//ycnti = # rows; xcnti = # cols
	Int nROWSi = ycnti;
	Int nCOLSi = xcnti;
	Int nROWSo = ycnt;
	Int nCOLSo = xcnt;

	MatDoub InputMatrix(nROWSi,nCOLSi); 
	MatDoub OutputMatrix(nROWSo,nCOLSo);

	VecDoub xi(nCOLSi),yi(nROWSi);
		
	for (int xndx = 0; xndx < nCOLSi; xndx++)
	{
		xi.operator[](xndx) = grid_xvi[xndx];
	}
	
	for (int yndx = 0; yndx < nROWSi; yndx++)
	{
		yi.operator[](yndx) = grid_yvi[(nROWSi-1)-yndx];
	}
	
	status = convertMatlabVectToNumRecipMatrix(grid_zvi, InputMatrix);

	//Bicubic Spline Interpolant
	Spline2D_interp bicubicfunc(yi,xi,InputMatrix);

	for (int xndx = 0; xndx < nCOLSo; xndx++)
	{
		std::cout << "Processing col " << xndx << " of " << nCOLSo << "\n"; 
		for (int yndx = 0; yndx < nROWSo; yndx++)
		{
			OutputMatrix[yndx][xndx]		= bicubicfunc.interp(grid_yv[(nROWSo-1)-yndx],grid_xv[xndx]);
			grid_zv[yndx + (nROWSo*xndx)]	= OutputMatrix[yndx][xndx];
		}
	}
	
	status = 1;

	std::cout << "Done.\n";

	return status;

}

//Need to double check the column major vs. row major ordering between C++ and MATLAB.
int numrecip_reggrid_poly2d(int polyinterp_order)
{
	std::cout << "Executing Numerical Recipes Higher Order Polynomial Interpolation(Single Thread)...\n";
	
	//Assume uniform polynomial order for interp.
	int mOrder = polyinterp_order;
	int nOrder = polyinterp_order;

	//Input grid
	int xcnti          = (int)grid_xvi.size();
	int ycnti          = (int)grid_yvi.size();

	//Output grid
	int status			= 0;
	int xcnt			= (int)grid_xv.size();
	int ycnt			= (int)grid_yv.size();
	int zcnt			= xcnt*ycnt;
	grid_zv				= std::vector<double>(zcnt);
	
	///////////////////////////////////////////////////////
	//Numerical Recipes interface 
	//ycnti = # rows; xcnti = # cols
	Int nROWSi = ycnti;
	Int nCOLSi = xcnti;
	Int nROWSo = ycnt;
	Int nCOLSo = xcnt;

	MatDoub InputMatrix(nROWSi,nCOLSi); 
	MatDoub OutputMatrix(nROWSo,nCOLSo);

	VecDoub xi(nCOLSi),yi(nROWSi);
		
	for (int xndx = 0; xndx < nCOLSi; xndx++)
	{
		xi.operator[](xndx) = grid_xvi[xndx];
	}
	
	for (int yndx = 0; yndx < nROWSi; yndx++)
	{
		yi.operator[](yndx) = grid_yvi[(nROWSi-1)-yndx];
	}
	
	
	status = convertMatlabVectToNumRecipMatrix(grid_zvi, InputMatrix);

	//Polynomial Interpolant(Order >= 2)
	Poly2D_interp polyfunc(yi,xi,InputMatrix,mOrder,nOrder);

	for (int xndx = 0; xndx < nCOLSo; xndx++)
	{
		std::cout << "Processing col " << xndx << " of " << nCOLSo << "\n"; 
		for (int yndx = 0; yndx < nROWSo; yndx++)
		{
			OutputMatrix[yndx][xndx]		= polyfunc.interp(grid_yv[(nROWSo-1)-yndx],grid_xv[xndx]);
			grid_zv[yndx + (nROWSo*xndx)]	= OutputMatrix[yndx][xndx];
		}
	}
	
	status = 1;

	std::cout << "Done.\n";

	return status;

}

//Need to double check the column major vs. row major ordering between C++ and MATLAB.
int numrecip_reggrid_bilinear()
{
	std::cout << "Executing Numerical Recipes Bi-Linear Interpolation(Single Thread)...\n";
	

	//Input grid
	int xcnti          = (int)grid_xvi.size();
	int ycnti          = (int)grid_yvi.size();

	//Output grid
	int status			= 0;
	int xcnt			= (int)grid_xv.size();
	int ycnt			= (int)grid_yv.size();
	int zcnt			= xcnt*ycnt;
	grid_zv				= std::vector<double>(zcnt);
	
	///////////////////////////////////////////////////////
	//Numerical Recipes interface 
	//ycnti = # rows; xcnti = # cols
	Int nROWSi = ycnti;
	Int nCOLSi = xcnti;
	Int nROWSo = ycnt;
	Int nCOLSo = xcnt;

	MatDoub InputMatrix(nROWSi,nCOLSi); 

	MatDoub OutputMatrix(nROWSo,nCOLSo);

	VecDoub xi(nCOLSi),yi(nROWSi);
		
	for (int xndx = 0; xndx < nCOLSi; xndx++)
	{
		xi.operator[](xndx) = grid_xvi[xndx];
	}
	
	for (int yndx = 0; yndx < nROWSi; yndx++)
	{
		yi.operator[](yndx) = grid_yvi[(nROWSi-1)-yndx];
	}
	
	status = convertMatlabVectToNumRecipMatrix(grid_zvi, InputMatrix);
	
	//Bilinear Interpolant; not cartesian order but matrix dim order.
	Bilin_interp bifunc(yi,xi,InputMatrix);

	for (int xndx = 0; xndx < nCOLSo; xndx++)
	{
		for (int yndx = 0; yndx < nROWSo; yndx++)
		{
			std::cout << "Processing [Row,Col] = [" << yndx << ", " << xndx << "]\n"; 
			OutputMatrix[yndx][xndx]		= bifunc.interp(grid_yv[(nROWSo-1)-yndx],grid_xv[xndx]);
			grid_zv[yndx + (nROWSo*xndx)]	= OutputMatrix[yndx][xndx];
		}
	}

	status = 1;

	std::cout << "Done.\n";

	return status;
}