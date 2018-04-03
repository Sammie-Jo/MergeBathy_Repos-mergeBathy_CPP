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
* @file			bathyTool.h
* @purpose		Header file for single and monte-carlo data runs of mergeBathy.  Performs initial data manipulation before computation.
* @author		Kevin Duvieilh
* @date			16 June 2011
*
*/
#pragma once
#include <iostream>
#include <vector>
#include <map>
#include "scalecInterp.h"

/**
* Bathy tool that is used to call subsampleData and scalecInterpTile.  It is a scale-controlled interpolation of bathymetric (or other scalar) data.  It is only used for computed data points.
* @param inputDataX - Vector of input X data points.
* @param inputDataY - Vector of input Y data points.
* @param inputDataZ - Vector of input Depth data points.
* @param inputDataE - Vector of input Error data points.
* @param xMeshGrid - dgrid of the repeating X data points to be interpolated.  Each column contains the same X value normalized by the mean of X.
* @param yMeshGrid - dgrid of the repeating Y data points to be interpolated.  Each row contains the same Y value normalized by the mean of Y.
* @param xSingleVector - Vector of the non repeating X data points to be interpolated.
* @param ySingleVector - Vector of the non repeating Y data points to be interpolated.
* @param gridSpacingX - Spacing of the computational area in the X direction.
* @param gridSpacingY - Spacing of the computational area in the Y direction.
* @param x0 - Minimum value contained in the inputDataX vector.
* @param y0 - Minimum value contained in the inputDataY vector.
* @param meanXSingle - Mean value of inputDataX.
* @param meanYSingle - Mean value of inputDataY.
* @param kernelNames - Name of the current smoothing window interpolator.
* @param additionalOptions - A map of additional options that are required.
* @param subDataMultiplier - The multiplier value used in reducing the size of the subsampled data.
* @param dispIntermResults - Determine if intermediate output should be displayed to the command line.
* @param useDscale - Set to true for normal use, for compute offset set to false.
* @param neitol - Normalized Error Tolerance Value.
* @param xyzOut - OUTPUT_DATA that contains the interpolated depth, error, normalized error, and residual error. (Returned).
* @return Success or failure value.
*/

int bathyTool(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, dgrid *xMeshGrid, dgrid *yMeshGrid, vector<double> *xSingleVector, vector<double> *ySingleVector, double gridSpacingX, double gridSpacingY, double &x0, double &y0, double meanXSingle, double meanYSingle, string &kernelName,  map<string, int> additionalOptions, double subDataMulitplier, bool dispIntermResults, bool useDscale, const double neitol, OUTPUT_DATA *xyzOut);

/**
* Bathy tool that is used to call subsampleData and scalecInterp.  It is a scale-controlled interpolation of bathymetric (or other scalar) data.  It is only used for known plotting to specific data points defined by the user.
* @param inputDataX - Vector of input X data points.
* @param inputDataY - Vector of input Y data points.
* @param inputDataZ - Vector of input Depth data points.
* @param inputDataE - Vector of input Error data points.
* @param xInterpVector - Vector of the X points to be interpolated.
* @param yInterpVector - Vector of the Y points to be interpolated.
* @param gridSpacingX - Spacing of the computational area in the X direction.
* @param gridSpacingY - Spacing of the computational area in the Y direction.
* @param x0 - Minimum value contained in the inputDataX vector.
* @param y0 - Minimum value contained in the inputDataY vector.
* @param meanXSingle - Mean value of inputDataX.
* @param meanYSingle - Mean value of inputDataY.
* @param kernelNames - Name of the current smoothing window interpolate.
* @param additionalOptions - A map of additional options that are required.
* @param neitol - Normalized Error Tolerance Value.
* @param xyzOut - OUTPUT_DATA that contains the interpolated depth, error, normalized error, and residual error. (Returned).
* @return Success or failure value.
*/
int bathyToolPreDefined(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double &x0, double &y0, double meanXSingle, double meanYSingle, string &kernelName,  map<string, int> additionalOptions, const double neitol, OUTPUT_DATA *xyzOut);

