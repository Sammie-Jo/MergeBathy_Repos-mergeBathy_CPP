/**********************************************************************
* CC0 License
**********************************************************************
* MergeBathy - Tool to combine one or more bathymetric data files onto a single input grid.
* Written in 2015 by Samantha J.Zambo(samantha.zambo@gmail.com) while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Todd Holland while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Nathaniel Plant while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Kevin Duvieilh while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Paul Elmore while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Will Avera while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Brian Bourgeois while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by A.Louise Perkins while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by David Lalejini while employed by the U.S.Naval Research Laboratory.
* To the extent possible under law, the author(s) and the U.S.Naval Research Laboratory have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide.This software is distributed without any warranty.
* You should have received a copy of the CC0 Public Domain Dedication along with this software.If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
**********************************************************************/
/**
* @file			standardOperations.h
* @brief		Defines standard functions used throughout mergeBathy
* @author		Kevin Duvieilh
* @date			30 June 2011
*
*/

#pragma once
#include <cmath>
#include <vector>
#include "grid.h"
#include <stdtypes.h>

/**
* Find the median of a vector. Comparable to MATLAB.
* @param a - First doubles.
* @param b - Second doubles.
* @return The lesser value between a and b.
*/
double median(vector<double> bar);

/**
* Find the minimum of two doubles.
* @param a - First doubles.
* @param b - Second doubles.
* @return The lesser value between a and b.
*/
double findMin(double a, double b);

/**
* Find the maximum of two doubles.
* @param a - First doubles.
* @param b - Second doubles.
* @return The greater value between a and b.
*/
double findMax(double a, double b);

/**
* Find the minimum of two integers.
* @param a - First integer.
* @param b - Second integer.
* @return The lesser value between a and b.
*/
int findMin(int a, int b);
double min (const double& x, const double& y);// const ;//{ if(x < y) { return x; } return y; }

double max (const double& x, const double& y);// const { if(x > y) { return x; } return y; }
/**
* Find the maximum of two integers.
* @param a - First integer.
* @param b - Second integer.
* @return The greater value between a and b.
*/
int findMax(int a, int b);

/**
* Round a double up or down.
* @param a - Double value to round.
* @return The rounded value.
*/
double roundDouble (double a);

/**
* Compute the mean from a vector.
* @param toMean - Vector of data points.
* @return The mean of the vector toMean.
*/
double mean(dvector *toMean, bool useZeros = True);

/**
* Compute the standard deviation from a vector.
* @param toStd - Vector of data points.
* @return The standard deviation of the vector toStd.
*/
double standardDeviation(dvector *toStd, bool useZeros = True);

/**
* Truncates value. (Matlab function).
* @param x - value to fix.
* @return The fixed value.
*/
double fix(double x);

/**
* Computes the unique X and Y values for the dimensions of a mesh grid.
* @param x0 - min x value for mesh boundary.
* @param y0 - min y value for mesh boundary.
* @param x1 - max x value for mesh boundary.
* @param y1 - max y value for mesh boundary.
* @param gridSpacingX - grid step size in x direction.
* @param gridSpacingY - grid step size in y direction.
* @param xt - vector of unique X dimension values.
* @param yt - vector of unique Y dimension values.
* @return The unique X and Y values for the dimensions of a mesh grid.
*/
void createMeshXYDims(double x0, double y0, double x1, double y1, double gridSpacingX, double gridSpacingY, vector<double> *xt, vector<double> *yt);

/**
* Computes the repeated X and Y values of a mesh grid.
* @param xt - vector of unique X dimension values.
* @param yt - vector of unique Y dimension values.
* @param xMeshVector - vector of repeated X grid values.
* @param yMeshVector - vector of repeated Y grid values.
* @param xMeshGrid - dgrid of repeated X grid values.
* @param yMeshGrid - dgrid of repeated Y grid values.
* @return The repeated X and Y values of a mesh grid.
*/
void createMeshGrid(vector<double> *xt, vector<double> *yt, vector<double> *xMeshVector, vector<double> *yMeshVector, dgrid *xMeshGrid, dgrid *yMeshGrid);

/**
* Reshapes the dgrid repeated values of a mesh grid. The (r*c) must be the same.
* @param g - dgrid of repeated grid values.
* @param newg - dgrid of repeated grid values.
* @return The reshaped grid of repeated values of a mesh grid.
*/
void reshapeGrid(dgrid *g, dgrid *newg);

template <typename T>
	string toStr ( T Number )
	{
		std::ostringstream ss;
		ss << Number;
		return ss.str();
	}

template <typename T>
	T toNum ( const string &Text )
	{
		std::istringstream ss(Text);
		T result;
		return ss >> result ? result : 0;
	}
