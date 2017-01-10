//Numerical Recipes Include files

/*
Numerical Recipes algorithms implemented using 3rd Edition by William H. Press
Chapter 3.; 2007.
*/

#ifndef NUMERICAL_RECIPES_INTERP_H
#define NUMERICAL_RECIPES_INTERP_H

#include <iostream>
#include <ostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <Windows.h>
#include "nr3.h"

//Matlab to Num. Recipes Conversion
int convertMatlabVectToNumRecipMatrix(std::vector<double>&, MatDoub&);
std::vector<double>* convertNumRecipMatrixToMatlabVect(MatDoub&, int, int);

//Regular Grid Interpolation Routines
int	numrecip_reggrid_bicubicspline();
int	numrecip_reggrid_poly2d(int);
int numrecip_reggrid_bilinear();

#endif
