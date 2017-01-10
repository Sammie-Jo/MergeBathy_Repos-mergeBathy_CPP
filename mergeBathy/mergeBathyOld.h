/** 
* @file			mergeBathyOld.h
* @brief		Header file for single and monte-carlo data runs of mergeBathy.  Performs initial data manipulation before computation.
* @author		Kevin Duvieilh
* @date			04 August 2011
*
*/
#pragma once
#include <iostream>
#include <vector>
#include <map>
#include "bathyTool.h"
#include "Error_Estimator/Bathy_Grid.h"

/**
* Driver routine for merging analysis of bathymetry data. Calls either runSingle or runMonteCarlo.
* @param inputData - Vector TRUE_DATA that contains the values read from each file.  The length of the vector corresponds to the number of input files.  Each index of the vector corresponds the TRUE_DATA read from a specific input file. 
* @param refLon - Center longitude of the given computational area.
* @param refLat - Center latitude of the given computational area.
* @param rotationAngle - An angle [in degrees] of arbitrary rotation.  This allows the computation X and Y space to undergo a rotation, such that the along-shore or cross-shore direction can be oriented along a given axis.  This allows the user to more easily select unequal arbitrary smoothing length scales.
* @param gridSpacingX - Computational grid spacing in the X direction defined in meters.
* @param gridSpacingY - Computational grid spacing in the Y direction defined in meters.
* @param smoothingScaleX - Grid smoothing scale in the X direction defined in meters.
* @param smoothingScaleY - Grid smoothing scale in the Y direction defined in meters.
* @param kernelNames - Interpolation window (i.e. filter window) to use, the choices are 'quadloess', 'linloess', 'hanning' or 'hann', and 'boxcar'. Hanning is the standard smoothing window.
* @param outputFileName - A text file consisting of one line for each output point at which a valid depth was computed.  Each line consists of four columns, corresponding to Longitude, Latitude, Depth, and Uncertainty (the mean square interpolation error estimate) respectively.  Additional columns can be present depending on other additional command line options provided at runtime.
* @param additionalOptions - A map of additional options that are required.
* @param numMCRuns - The number of Monte Carlo computations to perform.
* @param MB_ZGridInput - MB_ZGRID_DATA structure defining the parameters for running MB_ZGrid if requested by the user.
* @param GMTSurfaceInput - GMT_SURFACE_DATA structure defining the parameters for running MB_ZGrid if requested by the user.
* @param forcedLocationPositions - FORCED_LOCATIONS structure defining the exact points to interpolate the data to.
* @return Success or failure value.
*/
int mergeBathy_PreCompute( std::vector<TRUE_DATA> *inputData, double refLon, double refLat, double rotationAngle, double gridSpacingX, double gridSpacingY, double smoothingScaleX, double smoothingScaleY,  std::string kernelName, std::string outputFileName, map<std::string, int> additionalOptions, int numMCRuns, MB_ZGRID_DATA *MB_ZGridInput, GMT_SURFACE_DATA *GMTSurfaceInput, ALG_SPLINE_DATA *ALGSplineInput, FORCED_LOCATIONS *forcedLocationPositions,int usagePreInterpLocsLatLon);

/**
* Computing routine for merging analysis of bathymetry data.  Only used in single data runs.
* @param inputDataX - Vector of input X data points.
* @param inputDataY - Vector of input Y data points.
* @param inputDataZ - Vector of input Depth data points.
* @param inputDataE - Vector of input Error data points.
* @param inputDataHErr - Vector of input Horizontal Error data points.
* @param inputDataVErr - Vector of input Vertical Error data points.
* @param xMeshGrid - dgrid of the repeating X data points to be interpolated.  Each column contains the same X value normalized by the mean of X.
* @param yMeshGrid - dgrid of the repeating Y data points to be interpolated.  Each row contains the same Y value normalized by the mean of Y.
* @param xSingleVector - Vector of the non repeating X data points to be interpolated.
* @param ySingleVector - Vector of the non repeating Y data points to be interpolated.
* @param xInterpVector - Vector of the X points to be interpolated.  The xMeshGrid converted to vector form.
* @param yInterpVector - Vector of the Y points to be interpolated.  The yMeshGrid converted to vector form.
* @param gridSpacingX - Computational grid spacing in the X direction defined in meters.
* @param gridSpacingY - Computational grid spacing in the Y direction defined in meters.
* @param smoothingScaleX - Grid smoothing scale in the X direction defined in meters.
* @param smoothingScaleY - Grid smoothing scale in the Y direction defined in meters.
* @param x0 - Minimum value contained in the inputDataX vector.
* @param y0 - Minimum value contained in the inputDataY vector.
* @param x1 - Maximum value contained in the inputDataX vector.
* @param y1 - Maximum value contained in the inputDataY vector.
* @param meanXSingle - Mean value of inputDataX.
* @param meanYSingle - Mean value of inputDataY.
* @param kernelNames - Name of the current smoothing window interpolator.
* @param additionalOptions - A map of additional options that are required.
* @param outputFileName - The name of the output file to be written.
* @param subDataMultiplier - The multiplier value used in reducing the size of the subsampled data.
* @param UTMNorthingRef - The UTM Northing reference point.
* @param UTMEastingRef - The UTM Easting reference point.
* @param rotAngle - An angle [in degrees] of arbitrary rotation.  This allows the computation X and Y space to undergo a rotation, such that the along-shore or cross-shore direction can be oriented along a given axis.  This allows the user to more easily select unequal arbitrary smoothing length scales.
* @param RefEllip - The reference ellipsoid used for converting from Longitude and Latitude to UTM.
* @param UTMZoneRef - The UTM Zone reference character.
* @param MB_ZGridInput - MB_ZGRID_DATA structure defining the parameters for running MB_ZGrid if requested by the user.
* @param GMTSurfaceInput - GMT_SURFACE_DATA structure defining the parameters for running MB_ZGrid if requested by the user.
* @return Success or failure value.
*/
int runSingle(std::vector<double> *inputDataX,  std::vector<double> *inputDataY, std::vector<double> *inputDataZ,  std::vector<double> *inputDataE,  std::vector<double> *inputDataHErr, std::vector<double> *inputDataVErr, dgrid *xMeshGrid, dgrid *yMeshGrid, std::vector<double> *xSingleVector,  std::vector<double> *ySingleVector,  std::vector<double> *xInterpVector, std::vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double smoothingScaleX, double smoothingScaleY, double &x0, double &y0, double &x1, double &y1, double meanXSingle, double meanYSingle, std::string &kernelName,  map<std::string, int> additionalOptions, std::string outputFileName, double subDataMulitplier, double UTMNorthingRef, double UTMEastingRef, double rotAngle, int RefEllip, char UTMZoneRef[4], MB_ZGRID_DATA *MB_ZGridInput, GMT_SURFACE_DATA *GMTSurfaceInput, ALG_SPLINE_DATA *ALGSplineInput, Bathy_Grid* bathyGrid, int USE_UTM);


int run(std::vector<double> *inputDataX,  std::vector<double> *inputDataY, std::vector<double> *inputDataZ,  std::vector<double> *inputDataE,  std::vector<double> *inputDataHErr, std::vector<double> *inputDataVErr, dgrid *xMeshGrid, dgrid *yMeshGrid, std::vector<double> *xSingleVector,  std::vector<double> *ySingleVector,  std::vector<double> *xInterpVector, std::vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double smoothingScaleX, double smoothingScaleY, double &x0, double &y0, double &x1, double &y1, double meanXSingle, double meanYSingle, std::string &kernelName,  map<std::string, int> additionalOptions, std::string outputFileName, double subDataMulitplier, double UTMNorthingRef, double UTMEastingRef, double rotAngle, int RefEllip, char UTMZoneRef[4], int numMCRuns, MB_ZGRID_DATA *MB_ZGridInput, GMT_SURFACE_DATA *GMTSurfaceInput, ALG_SPLINE_DATA *ALGSplineInput, Bathy_Grid* bathyGrid, int USE_UTM);
/**
* Computing routine for merging analysis of bathmetry data.  Only used in single data runs.
* @param inputDataX - Vector of input X data points.
* @param inputDataY - Vector of input Y data points.
* @param inputDataZ - Vector of input Depth data points.
* @param inputDataE - Vector of input Error data points.
* @param inputDataHErr - Vector of input Horizontal Error data points.
* @param inputDataVErr - Vector of input Vertical Error data points.
* @param xMeshGrid - dgrid of the repeating X data points to be interpolated.  Each column contains the same X value normalized by the mean of X.
* @param yMeshGrid - dgrid of the repeating Y data points to be interpolated.  Each row contains the same Y value normalized by the mean of Y.
* @param xSingleVector - Vector of the non repeating X data points to be interpolated.
* @param ySingleVector - Vector of the non repeating Y data points to be interpolated.
* @param xInterpVector - Vector of the X points to be interpolated.  The xMeshGrid converted to vector form.
* @param yInterpVector - Vector of the Y points to be interpolated.  The yMeshGrid converted to vector form.
* @param gridSpacingX - Computational grid spacing in the X direction defined in meters.
* @param gridSpacingY - Computational grid spacing in the Y direction defined in meters.
* @param smoothingScaleX - Grid smoothing scale in the X direction defined in meters.
* @param smoothingScaleY - Grid smoothing scale in the Y direction defined in meters.
* @param x0 - Minimum value contained in the inputDataX vector.
* @param y0 - Minimum value contained in the inputDataY vector.
* @param x1 - Maximum value contained in the inputDataX vector.
* @param y1 - Maximum value contained in the inputDataY vector.
* @param meanXSingle - Mean value of inputDataX.
* @param meanYSingle - Mean value of inputDataY.
* @param kernelNames - Name of the current smoothing window interpolator.
* @param additionalOptions - A map of additional options that are required.
* @param outputFileName - The name of the output file to be written.
* @param subDataMultiplier - The multiplier value used in reducing the size of the subsampled data.
* @param UTMNorthingRef - The UTM Northing reference point.
* @param UTMEastingRef - The UTM Easting reference point.
* @param rotAngle - An angle [in degrees] of arbitrary rotation.  This allows the computation X and Y space to undergo a rotation, such that the along-shore or cross-shore direction can be oriented along a given axis.  This allows the user to more easily select unequal arbitrary smoothing length scales.
* @param RefEllip - The reference ellipsoid used for converting from Longitude and Latitude to UTM.
* @param UTMZoneRef - The UTM Zone reference character.
* @param numMCRuns - The number of Monte Carlo computations to run.
* @param MB_ZGridInput - MB_ZGRID_DATA structure defining the parameters for running MB_ZGrid if requested by the user.
* @param GMTSurfaceInput - GMT_SURFACE_DATA structure defining the parameters for running MB_ZGrid if requested by the user.
* @return Success or failure value.
*/

int runMonteCarlo(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, dgrid *xMeshGrid, dgrid *yMeshGrid, vector<double> *xSingleVector, vector<double> *ySingleVector, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double smoothingScaleX, double smoothingScaleY, double &x0, double &y0, double &x1, double &y1, double meanXSingle, double meanYSingle, string &kernelName,  map<string, int> additionalOptions, string outputFileName, double subDataMulitplier, double UTMNorthingRef, double UTMEastingRef, double rotAngle, int RefEllip, char UTMZoneRef[4], int numMCRuns, MB_ZGRID_DATA *MB_ZGridInput, GMT_SURFACE_DATA *GMTSurfaceInput, ALG_SPLINE_DATA *ALGSplineInput, Bathy_Grid* bathyGrid);
