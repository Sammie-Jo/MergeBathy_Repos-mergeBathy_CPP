/** 
* @file			outFileStructs.h
* @brief		Define output file structures going from scalecInterp to mergeBathyOld.
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/

#pragma once
#include "Error_Estimator/Bathy_Grid.h"
#include <vector>
#include <string>
#include "grid.h"

//#pragma region Define OUTPUT_DATA
//************************************************************************
// Used for returning the interpolated values so they can be printed.
//************************************************************************
typedef struct
{
	/**
	* The interpolated depth.
	*/
	vector<double> depth;

	/**
	* The interpolated error.
	*/
	vector<double> error;

	/**
	* The interpolated normalized error.
	*/
	vector<double> nEi;

	/**
	* The interpolated residual error.
	*/
	vector<double> rEi;
	
	/**
	* The interpolated standard deviation
	*/
	vector<double> standardDev;

	/**
	* The interpolated propagated uncertainty depth.
	*/
	vector<double> depth0;

	/**
	* The interpolated propagated uncertainty error.
	*/
	vector<double> error0;

	/**
	* The interpolated propagated uncertainty standard deviation
	*/
//	vector<double> standardDev0;
		
	/**
	* The interpolated kalman depth.
	*/
	vector<double> depthK;
		
	/**
	* The interpolated kalman error.
	*/
	vector<double> errorK;

} OUTPUT_DATA;
//#pragma endregion


/*PreCompute with PERTURBS object to hold window determined ai and ri, and calculated perturbation values*/
//************************************************************************
// Used for passing necessary data values to perform 
// interpolation computations at a single point.
//************************************************************************
typedef struct {
	
	/**
	* kernelName - Name of the current smoothing window interpolator.
	*/
	string kernelName;
	
	/**
	* The interpolated depth.
	*/
	double perturbationZ;
	
	/**
	* The interpolated error.
	*/
	double perturbationE;
	
	/**
	* The interpolated standard deviation.
	*/
	double standardDev2;
	
	/**
	* The interpolated standard deviation.
	*/
	double standardDev;
		
	/**
	* Interpolated normalized error.
	*/
	double perturbationNEi;
	
	/**
	* Interpolated residual error.
	*/
	double perturbationREi;

	/**
	* The interpolated propagated uncertainty depth.
	*/
	double perturbationZ0;
	
	/**
	* The interpolated propagated uncertainty error.
	*/
	double perturbationE0;
	
	/**
	* The interpolated propagated uncertainty standard deviation.
	*/
	//double standardDev20;
	
	/**
	* The interpolated kalman depth.
	*/
	double perturbationZK;
	
	/**
	* The interpolated kalman error.
	*/
	double perturbationEK;
	


	/**
	* The interpolated depth kriged.
	*/
	double perturbationZKriged;
	
	/**
	* The interpolated error kriged.
	*/
	double perturbationEKriged;
	
	/**
	* The interpolated standard deviation.
	*/
//	double standardDev20;
	
	/**
	* The interpolated standard deviation.
	*/
//	double standardDev0;
		
	/**
	* Interpolated normalized error kriged.
	*/
	double perturbationNEiKriged;
	
	/**
	* Interpolated residual error kriged.
	*/
	double perturbationREiKriged;

	/**
	* The interpolated propagated uncertainty depth kriged.
	*/
	double perturbationZ0Kriged;
	
	/**
	* The interpolated propagated uncertainty error kriged.
	*/
	double perturbationE0Kriged;
	
	/**
	* The interpolated propagated uncertainty standard deviation.
	*/
//	double standardDev20;
	
	/**
	* The interpolated kalman depth kriged.
	*/
	double perturbationZKKriged;
	
	/**
	* The interpolated kalman error kriged.
	*/
	double perturbationEKKriged;





	/**
	* The radial distances based on window smoothing scale.
	*/
	vector<double> riVector;
	
	/**
	* The weights based on window smoothing scale.
	*/
	vector<double> aiVector;


}PERTURBS;

//#pragma region Define SCALEC_TILE_DATA, *SCALEC_TILE_DATA_POINTER
//************************************************************************
// Used for passing necessary data structures to the scalecInterpTile computation routines (Kriged or Non-Kriged).
//************************************************************************
typedef struct {
	/**
	* subsampledData - A n by 5 vector containing the sub-sampled data.  Index 0 contains the X data normalized to the mean of X. Index 1 contains the Y data normalized to the mean of Y. Index 2 contains the Depth data. Index 3 contains the Error data.  Index 4 contains the Error data squared.
	*/
    const vector< vector<double> > *subsampledData;

	/**
	* xSingleVector - Vector of the non repeating X data points to be interpolated.
	*/
	const vector<double> *xSingleVector;

	/**
	* ySingleVector - Vector of the non repeating Y data points to be interpolated.
	*/
	const vector<double> *ySingleVector;
	
	/**
	* xMeshGrid - dgrid of the repeating X data points to be interpolated.  Each column contains the same X value normalized by the mean of X.
	*/
	const dgrid *xMeshGrid;

	/**
	* yMeshGrid - dgrid of the repeating Y data points to be interpolated.  Each row contains the same Y value normalized by the mean of Y.
	*/
	const dgrid *yMeshGrid;

	/**
	* innerLoopIndexVector - Vector of vectors of the Y indices to interpolate.
	*/
	const vector< vector<int> > *innerLoopIndexVector;

	/**
	* btrend - Vector of the calculated trend surface values.
	*/
	const dvector *btrend;

	/**
	* kernelNames - Name of the current smoothing window interpolator.
	*/
	string *kernelName;

	/**
	* Lx - Modified X smoothing scale.
	*/
	const double *Lx;

	/**
	* Ly - Modified Y smoothing scale.
	*/
	const double *Ly;

	/**
	* spacingX - X smoothing scale.
	*/
	const double *spacingX;

	/**
	* spacingY - Y smoothing scale.
	*/
	const double *spacingY;

	/**
	* kx - Number of tiles in the X axis.
	*/
	const int *kx;

	/**
	* ky - Number of tiles in the Y axis.
	*/
	const int *ky;

	/**
	* nkx - Number of X points in each tile.
	*/
	const int *nkx;

	/**
	* nky - Number of Y points in each tile.
	*/
	const int *nky;

	/**
	* LMAX_x - Maximum overlap in the X axis.
	*/
	const double *LMAX_x;

	/**
	* LMAX_y - Maximum overlap in the Y axis.
	*/
	const double *LMAX_y;

	/**
	* neitol - Normalized Error Tolerance Value.
	*/
	const double *neitol;

	/**
	* outputDepth - Interpolated Depth value grid. (Returned).
	*/
	dgrid *outputDepth;

	/**
	* outputError - Interpolated Error value grid. (Returned).
	*/
	dgrid *outputError;

	/**
	* outputNEi - Interpolated Normalized Error value grid. (Returned).
	*/
	dgrid *outputNEi;

	/**
	* outputREi - Interpolated Residual Error value grid. (Returned).
	*/
	dgrid *outputREi;

	/**
	* outputStandardDev - Interpolated Standard Deviation value grid. (Returned).
	*/
	dgrid *standardDev;

	/**
	* outputDepth0 - Interpolated Depth value grid. (Returned).
	*/
	dgrid *outputDepth0;

	/**
	* outputError0 - Interpolated KALMAN Error value grid. (Returned).
	*/
	dgrid *outputError0;

	/*
	* standarDev0 - Interpolated Error value grid with propagated uncertainty. (Returned).
	*/
//	dgrid *standardDev0;

	/**
	* outputDepthKZ - Interpolated KALMAN Depth value grid. (Returned).
	*/
	dgrid *outputDepthK;

	/**
	* outputErrorK - Interpolated KALMAN Error value grid. (Returned).
	*/
	dgrid *outputErrorK;

	/**
	* x0 - Vector of input X data points with preserved spatial location.
	*/
	vector<double> *x0;
	
	/**
	* x0 - Vector of input X data points with preserved spatial location.
	*/
	vector<double> *y0;

	/**
	* yInterpLocs0 - dgrid of Y interpolation grid locations.
	*/
	dgrid *xInterpLocs0;
	
	/**
	* yInterpLocs0 - dgrid of Y interpolation grid locations.
	*/
	dgrid *yInterpLocs0;

	/**
	* new_bathyGrid - Bathy_Grid object for constructing a Delaunay Triangulation.
	*/
	Bathy_Grid *new_bathyGrid;

	/**
	* interpMethod - Interpolation method to perform by Delaunay Triangulation's estimate function.
	*/
	string interpMethod;

	///**
	//* depthInterpMethod - Depth interpolation method to perform by Delaunay Triangulation's estimate function.
	//*/
	//string depthInterpMethod;
	///**
	//* errorInterpMethod - Error interpolation method to perform by Delaunay Triangulation's estimate function.
	//*/
	//string errorInterpMethod;

	/**
	* MSE - Mean Square Error Estimator flag.
	*/
	bool MSE;			
	
	/**
	* PROP_UNCERT - Propagated Uncertainty Estimator flag.
	*/
	bool PROP_UNCERT;	
	
	/**
	* KALMAN - Kalman Estimator flag.
	*/
	bool KALMAN;

	/**
	* KRIGING - Kriging flag.
	*/bool KRIGING;


} SCALEC_TILE_DATA, *SCALEC_TILE_DATA_POINTER;
//#pragma endregion

//#pragma region Define SCALEC_DATA, *SCALEC_DATA_POINTER
//************************************************************************
// Used for passing necessary data structures to the scalecInterpTile computation routines (Kriged or Non-Kriged).
//************************************************************************
typedef struct {
	/**
	* subsampledData - A n by 5 vector containing the subsampled data.  Index 0 contains the X data normalized to the mean of X. Index 1 contains the Y data normalized to the mean of Y. Index 2 contains the Dept  h data. Index 3 contains the Error data.  Index 4 contains the Error data squared.
	*/
    const vector< vector<double> > *subsampledData;

	/**
	* xSingleVector - Vector of the non repeating X data points to be interpolated.
	*/
	const vector<double> *xInterpVector;

	/**
	* ySingleVector - Vector of the non repeating Y data points to be interpolated.
	*/
	const vector<double> *yInterpVector;

	/**
	* btrend - Vector of the calculated trend surface values.
	*/
	const dvector *btrend;

	/**
	* kernelName - Name of the current smoothing window interpolator.
	*/
	string *kernelName;
	
	/**
	* Lx - Modified X smoothing scale.
	*/
	const double *Lx;

	/**
	* Ly - Modified Y smoothing scale.
	*/
	const double *Ly;

	/**
	* minSubX - Minimum subsampled X value.
	*/
	const double *minSubX;

	/**
	* minSubY - Minimum subsampled Y value.
	*/
	const double *minSubY;

	/**
	* spacingX - X smoothing scale.
	*/
	const double *spacingX;

	/**
	* spacingY - Y smoothing scale.
	*/
	const double *spacingY;

	/**
	* kx - Number of tiles in the X axis.
	*/
	const double *kx;

	/**
	* ky - Number of tiles in the Y axis.
	*/
	const double *ky;

	/**
	* nkx - Number of X points in each tile.
	*/
	const double *numStepsX;

	/**
	* nky - Number of Y points in each tile.
	*/
	const double *numStepsY;

	/**
	* LMAX_x - Maximum overlap in the X axis.
	*/
	const double *LMAX_x;

	/**
	* LMAX_y - Maximum overlap in the Y axis.
	*/
	const double *LMAX_y;

	/**
	* neitol - Normalized Error Tolerance Value.
	*/
	const double *neitol;

	/**
	* outData - The output data structures. (Returned).
	*/
	OUTPUT_DATA *outData;

	/**
	* MSE - Mean Square Error Estimator flag.
	*/
	bool MSE;			
	
	/**
	* PROP_UNCERT - Propagated Uncertainty Estimator flag.
	*/
	bool PROP_UNCERT;	
	
	/**
	* KALMAN - Kalman Estimator flag.
	*/
	bool KALMAN;

	/**
	* KRIGING - Kriging flag.
	*/bool KRIGING;

	/**
	* residualObservationsKrigedZ - Kriging Residuals.
	*/
	vector<double> *residualObservationsKrigedZ;

	/**
	* residualObservationsKrigedZ0 - Propagated Uncertainty Kriging Residuals.
	*/
	vector<double> *residualObservationsKrigedZ0;

	/**
	* residualObservationsKrigedZK - Kalman Kriging Residuals.
	*/
	vector<double> *residualObservationsKrigedZK;

	/**
	* x_idx - Vector of input X data points.
	*/
	vector<double> *x_idx;
	
	/**
	* y_idx - Vector of input Y data points.
	*/
	vector<double> *y_idx;
	
	/**
	* x_idxKriged - Vector of input X data points.
	*/
	vector<double> *x_idxKriged;
	
	/**
	* y_idxKriged - Vector of input Y data points.
	*/
	vector<double> *y_idxKriged;
	
	/**
	* z_idx - Vector of input Z data points.
	*/
	vector<double> *z_idx;
	
	/**
	* e_idx - Vector of input E data points.
	*/
	vector<double> *e_idx;
	
	/**
	* h_idx - Vector of input H data points.
	*/
	vector<double> *h_idx;
	
	/**
	* v_idx - Vector of input V data points.
	*/
	vector<double> *v_idx;

	/**
	* x0_idx - Vector of input X data points with preserved spatial location.
	*/
	vector<double> *x0_idx;
	
	/**
	* y0_idx - Vector of input X data points with preserved spatial location..
	*/
	vector<double> *y0_idx;

	/**
	* xInterpLocs0 - Vector of X interpolation grid locations.
	*/
	vector<double> *xInterpLocs0;
	
	/**
	* yInterpLocs0 - Vector of Y interpolation grid locations.
	*/
	vector<double> *yInterpLocs0;
	
	/**
	* new_bathyGrid - Bathy_Grid object for constructing a Delaunay Triangulation.
	*/
	Bathy_Grid *new_bathyGrid;

	/**
	* slopeOut - Pointer to interpolation grid location slopes.
	*/
	vector<double> *slopeOut;

	/**
	* interpMethod - Interpolation method to perform by Delaunay Triangulation's estimate function.
	*/
	string interpMethod;

	///**
	//* depthInterpMethod - Depth interpolation method to perform by Delaunay Triangulation's estimate function.
	//*/
	//string depthInterpMethod;

	///**
	//* errorInterpMethod - Error interpolation method to perform by Delaunay Triangulation's estimate function.
	//*/
	//string errorInterpMethod;


} SCALEC_DATA, *SCALEC_DATA_POINTER;
//#pragma endregion
