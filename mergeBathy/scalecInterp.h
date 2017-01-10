/** 
* @file			scalecInterp.h
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
#include "grid.h"
#include "consistentWeights.h"
#include "scalecInterpPerturbations.h"

//************************************************************************************
//I. scalecInterp
//************************************************************************************
/**
* Processes irregular grid spacing and passes data to scalecInterpPerturbations (which does not remove any trend).
* @param subsampledData - A n by 5 vector containing the sub-sampled data.  Index 0 contains the X data normalized to the mean of X. Index 1 contains the Y data normalized to the mean of Y. Index 2 contains the Depth data. Index 3 contains the Error data.  Index 4 contains the Error data squared.
* @param xInterpVector - Vector of the X points to be interpolated.
* @param yInterpVector - Vector of the Y points to be interpolated.
* @param gridSpacingX - Spacing of the computational area in the X direction.
* @param gridSpacingY - Spacing of the computational area in the Y direction.
* @param meanXSingle - Mean value of inputDataX.
* @param meanYSingle - Mean value of inputDataY.
* @param kernelName - Name of the current smoothing window interpolator.
* @param additionalOptions - A map of additional options that are required.
* @param neitol - Normalized Error Tolerance Value.
* @param xyzOut - OUTPUT_DATA that contains the interpolated depth, error, normalized error, and residual error. (Returned).
* @return Success or failure value.
*/
int scalecInterp(vector< vector<double> > *subsampledData, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double meanXSingle, double meanYSingle, string &kernelName, map<string, int> additionalOptions, const double neitol, OUTPUT_DATA *xyzOut);

/**
* Part I: Catches regular grid output, which are passed to scalecInterpPerturbations (which does not remove any trend).  Used in when Kriging is NOT being done.
* Function was split to allow multi-Threading.
* Used only for data runs that do not krig the data.
* @param sdp - A SCALEC_DATA_POINTER for the input data structures.
* @param curIterNum - The current core number. 0 on single threaded runs.
* @param numCores - The number of total cores processing the data.  1 on single threaded runs.
* return Success or failure value.
*/
int scalecInterp_Process(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores);

/**
* Part II: Catches regular grid output, which are passed to scalecInterpPerturbations (which does not remove any trend).  Used in when Kriging is NOT being done.
* Function was split to allow multi-Threading.
* Used only for data runs that do not krig the data.
* @param sdp - A SCALEC_DATA_POINTER for the input data structures.
* @param curIterNum - The current core number. 0 on single threaded runs.
* @param numCores - The number of total cores processing the data.  1 on single threaded runs.
* return Success or failure value.
*/
int scalecInterp_Process2(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores);
int scalecInterp_Process2A(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores);
int scalecInterp_Process4A(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores);
int scalecInterp_Process5A(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores);
int scalecInterp_Process6A(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores);

/**
* Catches regular grid output and breaks into bite-size tiles, which are passed to scalecInterpPerturbations (which does not remove any trend).  Used in when Kriging is being done.
* Used only for data runs that do krig the data.
* @param sdp - A SCALEC_DATA_POINTER for the input data structures.
* @param curIterNum - The current core number. 0 on single threaded runs.
* @param numCores - The number of total cores processing the data.  1 on single threaded runs.
* return Success or failure value.
*/
int scalecInterp_ProcessKrig(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores);




//************************************************************************************
//II. scalecInterpTile
//************************************************************************************
/**
* Catches regular grid output and breaks into bite-size tiles, which are passed to scalecInterpPerturbations (which does not remove any trend).
* Used only for data runs that do not krig the data.
* @param subsampledData - A n by 5 vector containing the sub-sampled data.  Index 0 contains the X data normalized to the mean of X. Index 1 contains the Y data normalized to the mean of Y. Index 2 contains the Depth data. Index 3 contains the Error data.  Index 4 contains the Error data squared.
* @param xMeshGrid - dgrid of the repeating X data points to be interpolated.  Each column contains the same X value normalized by the mean of X.
* @param yMeshGrid - dgrid of the repeating Y data points to be interpolated.  Each row contains the same Y value normalized by the mean of Y.
* @param xSingleVector - Vector of the non repeating X data points to be interpolated.
* @param ySingleVector - Vector of the non repeating Y data points to be interpolated.
* @param gridSpacingX - Spacing of the computational area in the X direction.
* @param gridSpacingY - Spacing of the computational area in the Y direction.
* @param meanXSingle - Mean value of inputDataX.
* @param meanYSingle - Mean value of inputDataY.
* @param kernelName - Name of the current smoothing window interpolator.
* @param additionalOptions - A map of additional options that are required.
* @param dispIntermResults - Determine if intermediate output should be displayed to the command line.
* @param neitol - Normalized Error Tolerance Value.
* @param xyzOut - OUTPUT_DATA that contains the interpolated depth, error, normalized error, and residual error. (Returned).
* @return Success or failure value.
*/
int scalecInterpTile(vector< std::vector<double> > *subsampledData, dgrid *xMeshGrid, dgrid *yMeshGrid, vector<double> *xSingleVector, vector<double> *ySingleVector, double gridSpacingX, double gridSpacingY, double meanXSingle, double meanYSingle, string &kernelName, map<string, int> additionalOptions, bool dispIntermResults, const double neitol, OUTPUT_DATA *xyzOut);

/**
* Catches regular grid output and breaks into bite-size tiles, which are passed to scalecInterpPerturbations (which does not remove any trend).  Used in when Kriging is NOT being done.
* Used only for data runs that do not krig the data.
* @param stdp - A SCALEC_TILE_DATA_POINTER for the input data structures.
* @param curIterNum - The current core number. 0 on single threaded runs.
* @param numCores - The number of total cores processing the data.  1 on single threaded runs.
* return Success or failure value.
*/int scalecInterpTile_Process(SCALEC_TILE_DATA_POINTER stdp, const int curIterNum, const int numCores);
int scalecInterpTile_ProcessA(SCALEC_TILE_DATA_POINTER stdp, const int curIterNum, const int numCores);


/**
* Catches regular grid output and breaks into bite-size tiles, which are passed to scalecInterpPerturbations (which does not remove any trend).  Used in when Kriging is being done.
* Used only for data runs that do krig the data.
* @param stdp - A SCALEC_TILE_DATA_POINTER for the input data structures.
* @param curIterNum - The current core number.  0 on single threaded runs.
* @param numCores - The number of total cores processing the data.  1 on single threaded runs.
* return Success or failure value.
*/
int scalecInterpTile_ProcessKrig(SCALEC_TILE_DATA_POINTER stdp, const int curIterNum, const int numCores); 


