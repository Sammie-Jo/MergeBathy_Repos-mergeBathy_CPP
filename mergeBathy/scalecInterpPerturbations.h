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
* @file			scalecInterpPerturbations.h
* @brief		Header file for scalecInterp, scalecInterpTile, and scalecInterpPerturbations, the three main processing routines for mergeBathy.
* @author		Kevin Duvieilh
* @date			04 August 2011
*
*/

#pragma once
#include <iostream>
#include <vector>
#include <map>
#include "inFileStructs.h"
#include "outFileStructs.h"
#include "consistentWeights.h"


/**
* This is a special weighting scale defined in the matlab 
* code that performs a 1-D interpolation.  It processes the linear loess window and the 
* quadratic loess window.
* @param x - Weighting scale factor.
* @param xi - Vector of computed radial points. 
* @param yi - Vector of computed window weights. 
* @return The 1-D interpolated double value of the weighting factor applied to xi and yi.
*/
double interp1(double x, const vector<double> *xi, const vector<double> *yi);

/**
* This function is used to calculate the weights for linear and quadratic loess windows, formed from an impulse response away from edges.
* Reference: Cleveland, W. S. (1979). "ROBUST LOCALLY WEIGHTED REGRESSION AND SMOOTHING SCATTERPLOTS." Journal of the American Statistical Association 74(368): 829-836.
* @param d - Dimensions of the kernel.
* @param p - Power of the loess window. 1 is a linear window. 2 is a quadratic window.
* @param riVector - Vector of computed radial points. (Returned).
* @param aiVector - Vector of computed window weights. (Returned).
*/
void loessKernel(int dtemp, int p, vector<double> *ri, vector<double> *ai);

/**
* Computes the weighting scale, ri values, and ai values.
* This function should be called outside the main loop that uses the _Compute function.
* It was broken off from the main code for efficiency since this did not need to be recomputed every time. 
* @param subDataZ - Vector of the subsampled Depth data points.
* @param subDataE - Vector of the subsampled Error data points.
* @param perturbWeights - Vector of weights to be applied to each depth.  Index locations of the weights correspond to the index locations of the depths. (Returned).
* @param riVector - Vector of computed radial points. (Returned).
* @param aiVector - Vector of computed window weights. (Returned).
*/
//void scalecInterpPerturbations_PreCompute(const vector<double> *subDataZ, const vector<double> *subDataE, vector<double> *perturbWeights, vector<double> *riVector, vector<double> *aiVector, string *kernelName);
void scalecInterpPerturbations_PreCompute(const vector<double> *subDataZ, const vector<double> *subDataE, vector<double> *perturbWeights, PERTURBS *perturb);
/**
* This is the main interpolation function for mergeBathy.
* It computes the depth, error, normalized error, and residual error for a single computed data point.
* The _PreCompute function must be called ONCE for every chuck of data prior to calling this function.  (i.e. a tile from "scalecInterpTile" or the whole grid from "scalecInterp").
* IT DOES NOT REMOVE A POLYNOMIAL TREND, AS TREND IS ASSUMED ALREADY REMOVED FROM
* PERTURBATION DATA, DEFAULTS TO ZERO VALUE IF NO DATA
* @param subDataX - Vector of the subsampled X data points multiplied by 1/gridSpacing.
* @param subDataY - Vector of the subsampled Y data points multiplied by 1/gridSpacing.
* @param subDataZ - Vector of the subsampled Depth data points.
* @param subDataE - Vector of the subsampled Error data points.
* @param xGridValue - The interpolated X coordinate of the data to interpolate multiplied by 1/gridSpacing.
* @param yGridValue - The interpolated X coordinate of the data to interpolate multiplied by 1/gridSpacing.
* @param weights - Vector of weights to be applied to each depth.  Index locations of the weights correspond to the index locations of the depths.
* @param ri - Vector of computed radial points.
* @param ai - Vector of computed window weights.
* @param neitol - Normalized Error Tolerance Value.
* @param perturbationZ - Interpolated Depth value. (Returned).
* @param perturbationE - Interpolated Error value. (Returned).
* @param perturbationNEi - Interpolated Normalized Error value. (Returned).
* @param perturbationREi - Interpolated Residual Error value. (Returned).
*/
void scalecInterpPerturbations_Compute(const vector<double> *subDataX, const vector<double> *subDataY, const vector<double> *subDataZ, const vector<double> *subDataE, const vector<double> *subDataH, const vector<double> *subDataV, double *xGridValue, double *yGridValue, const vector<double> *weights, const double neitol, double dmin, double slopeVector, const vector<double> *subDataX0, const vector<double> *subDataY0, dgrid *Xiii, PERTURBS *perturb, bool MSE, bool PROP_UNCERT, bool KALMAN, bool KRIGING);

/**
* This is the secondary interpolation function for mergeBathy when kriging is being used.
* It computes the depth, error, normalized error, and residual error for a single computed data point so that a residual error can be established.  
* In the primary while loop "count" is incremented before "p" is computed.
* The _PreCompute function must be called ONCE for every chuck of data prior to calling this function.  (i.e. a tile from "scalecInterpTile" or the whole grid from "scalecInterp").
* IT DOES NOT REMOVE A POLYNOMIAL TREND, AS TREND IS ASSUMED ALREADY REMOVED FROM
* PERTURBATION DATA, DEFAULTS TO ZERO VALUE IF NO DATA
* @param subDataX - Vector of the subsampled X data points multiplied by 1/gridSpacing.
* @param subDataY - Vector of the subsampled Y data points multiplied by 1/gridSpacing.
* @param subDataZ - Vector of the subsampled Depth data points.
* @param subDataE - Vector of the subsampled Error data points.
* @param xGridValue - The interpolated X coordinate of the data to interpolate multiplied by 1/gridSpacing.
* @param yGridValue - The interpolated X coordinate of the data to interpolate multiplied by 1/gridSpacing.
* @param weights - Vector of weights to be applied to each depth.  Index locations of the weights correspond to the index locations of the depths.
* @param ri - Vector of computed radial points.
* @param ai - Vector of computed window weights.
* @param neitol - Normalized Error Tolerance Value.
* @param perturbationZ - Interpolated Depth value. (Returned).
* @param perturbationE - Interpolated Error value. (Returned).
* @param perturbationNEi - Interpolated Normalized Error value. (Returned).
* @param perturbationREi - Interpolated Residual Error value. (Returned).
*/
void scalecInterpPerturbations_Compute_ForKriging(const vector<double> *subDataX, const vector<double> *subDataY, const vector<double> *subDataZ, const vector<double> *subDataE, double *xGridValue, double *yGridValue, const vector<double> *weights, const vector<double> *ri, const vector<double> *ai, const double neitol, double *perturbationZ, double *perturbationE, double *perturbationNEi, double *perturbationREi, double *standardDev);














