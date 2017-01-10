/** 
* @file			kriging.h
* @brief		Header file for all Kriging routines.
* @author		Kevin Duvieilh
* @date			28 July 2011
*
*/

#ifndef KRIGING_H_
#define KRIGING_H_

#include <vector>
#include <iostream>
#include "constants.h"
#include "standardOperations.h"

using namespace std;

/**
* This function computes weighting values for use in the _PostCompute function.
* It also computes the inverse of the gammaDArray.  Inverting the grid has the highest computational and time requirements.
* This function is best for residual observation vectors smaller than 64 elements
* @param subX - Vector of the subsampled X coordinates to be computed.
* @param subY - Vector of the subsampled Y coordinates to be computed.
* @param residualObservations - Vector of the residual observed depths at the coordinates specified by subX and subY.
* @param twoGammaHatVector - Vector of the fine-scale variogram in line with theta_m. (Returned).
* @param distanceVectorBinCenters - The fit of the fine-scale variograms to spherical model. (Returned).
* @param aVectorFine - The Levenburg-Marquart fits for the spherical model of the semivariance. (Returned).
* @param invGammaDArray - The inverse of the computed solution for the interpolated residual surface using ordinary kriging as provided in Davis JC, Statistics and Data Analysis in Geology,  3rd edition, pp. 420-1, New York: Wiley, 2002. (Returned).
* @param A - The A matrix as per Calder's Eqn. (22), to account for anisotropy. (Returned).
*/
void ordinaryKrigingOfResiduals_PreCompute(vector<double> *subX, vector<double> *subY, vector<double> *residualObservations, vector<double> *twoGammaHatVector, vector<double> *distanceVectorBinCenters, vector<double> *aVectorFine, dgrid *invGammaDArray, dgrid *A);

/**
* This function computes weighting values for use in the _PostCompute function.
* It also computes the inverse of the gammaDArray.  Inverting the grid has the highest computational and time requirements.
* This function is best for residual observation vectors smaller than 64 elements
*
* THIS FUNCTION IS NOT YET IMPLEMENTED!!!
*
* @param subX - Vector of the subsampled X coordinates to be computed.
* @param subY - Vector of the subsampled Y coordinates to be computed.
* @param residualObservations - Vector of the residual observed depths at the coordinates specified by subX and subY.
* @param interpX - Vector of interpolation X coordinates.
* @param interpY - Vector of interpolation Y coordinates.
* @param spacingX - Grid spacing in meters in the X direction.
* @param spacingY - Grid spacing in meters in the Y direction.
* @param outDepthKrig - Vector of the computed depths for the given data points. (Returned).
* @param outErrorKrig - Vector of the computed errors for the given data points. (Returned).
*/
void ordinaryKrigingOfResiduals_PreComputeTile(vector<double> *subX, vector<double> *subY, vector<double> *residualObservations, vector<double> *interpX, vector<double> *interpY, double spacingX, double spacingY, vector<double> *outDepthKrig, vector<double> *outErrorKrig);

/************************************************************************************
This function "krigs" the data by placing the residuals back in.
It operates on a single gridded data point at a time over a smaller, resampled list of known data points
************************************************************************************/
/**
* This function "krigs" the data by placing the residuals back in.
* It operates on a single gridded data point at a time over a smaller, resampled list of known data points
* @param subX - Vector of the subsampled X coordinates to be computed.
* @param subY - Vector of the subsampled Y coordinates to be computed.
* @param residualObservations - Vector of the residual observed depths at the coordinates specified by subX and subY.
* @param xGridValue - The interpolated X coordinate of the data to krig.
* @param YGridValue - The interpolated Y coordinate of the data  to krig.
* @param twoGammaHatVector - Vector of the fine-scale variogram in line with theta_m. 
* @param distanceVectorBinCenters - The fit of the fine-scale variograms to spherical model. 
* @param aVectorFine - The Levenburg-Marquart fits for the spherical model of the semivariance. 
* @param invGammaDArray - The inverse of the computed solution for the interpolated residual surface using ordinary kriging as provided in Davis JC, Statistics and Data Analysis in Geology,  3rd edition, pp. 420-1, New York: Wiley, 2002. 
* @param A - The A matrix as per Calder's Eqn. (22), to account for anisotropy. 
* @param zKrigValue - The comptued depth for the given data point. (Returned).
* @param varZKrigValue - The comptued error for the given data point. (Returned).
*/
void ordinaryKrigingOfResiduals_PostCompute(vector<double> *subX, vector<double> *subY, vector<double> *residualObservations, double *xGridValue, double *yGridValue, vector<double> *twoGammaHatVector, vector<double> *distanceVectorBinCenters, vector<double> *aVectorFine, dgrid *invGammaDArray, dgrid *A, double *zKrigValue, double *varZKrigValue);

/**
* This function provides a polar 2-D emperical variogram. 
* Calder's Eqns. (15-17) and Cressie (1993), p.69.
* This was taken from the MATLAB version of the code.... It works properly.
* @param distanceVector - Vector of distances to compute over.
* @param angleVector - Vector of angles to compute over.
* @param distanceArray - Grid of the square roots of the difference squared in the x direction plus the difference squared in the y direction.
* @param angleArray - Grid of the arc tangent of the difference in the x direction by the difference in the y direction.
* @param deltaRSquaredArray - Grid of the difference of the residual observations squared.
* @param twoGammaHat - The grid of computed gamma values. The grid has angleVector number of rows and distanceVector number of columns. (Returned).
*/
void methodsOfMomentsEstimator(vector<double> *distanceVector, vector<double> *angleVector, dgrid *distanceArray, dgrid *angleArray, dgrid *deltaRSquaredArray, dgrid *twoGammaHat);

/**
* Rotate a matrix follow the rule of: y = PI * x / 180.  The output grid is:
* @n cos(y)		sin(y)
* @n -sin(y)		cos(y)
* @param x - Double value to rotate and grid.
* @return The grid rotated and compted from x.
*/
dgrid rotMtx(double x);


#endif
