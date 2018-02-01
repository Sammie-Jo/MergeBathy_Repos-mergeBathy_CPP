
#include "scalecInterp.h"
#include <fstream>
#include <time.h>
#include "regr_xzw.h"
#include "kriging.h"
#include "subSampleData.h"
#include "MB_Threads.h"
#include "Error_Estimator/Bathy_Grid.h" 
#include "Error_Estimator/GradientGrid.h"

int scalecInterpTile(vector< vector<double> > *subsampledData, dgrid *xMeshGrid, dgrid *yMeshGrid, vector<double> *xSingleVector, vector<double> *ySingleVector, double gridSpacingX, double gridSpacingY, double meanXSingle, double meanYSingle, string &kernelName, map<string, int> additionalOptions, bool dispIntermResults, const double neitol, OUTPUT_DATA *xyzOut)
{
	#pragma region --Hide Calculations
	/* Look at previous version for more confusing notes. -SJZ
	Note capitalization of variables xi, yi, Xi and Yi!

	The first part of scalecInterpTileWithAllUncertainty.m to about ln 141
	checks to see if we have gridded data.  This is not performed here.
	Therefore the variables do not match exactly.

	In this excluded MATLAB portion,
	xi is initially the entire regular grid in vector form, which is reshaped as a
	grid(r,c) and stored as Xi, and then xi becomes a vector of the regular grid col.

	In scalecInterpTile.cpp, this is the same as xMeshvector
	(or xMeshGrid+meanX stored in a vector) for the regular grid pts in xi,
	no reshaping to Xi since xMeshGrid+meanX is already in grid(r,c) form and we
	don't perform calculations in this form anyway, and then xSingleVector is the col vector.
	Furthermore, Xi0 (a new variable introduced for	Paul's v5 changes) is xMeshGrid in a vector.

	Vector XMeshVector (or xMeshGrid + meanX in a vector) - regular grid pts, (xi as a vector in matlab)
	Grid xMeshGrid + meanX - regular grid in grid(r,c) form, (Xi as a grid in matlab)
	vector xSingleVector - grid col vector (replaces xi as a vector in matlab)
	vector Xi0 - xMeshGrid in a vector, (Xi0 as a grid in matlab)
	*/

	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int i, j, k;

	//unsigned long nan2[2] = {0xffffffff, 0x7fffffff};
	//double NaN = *( double* )nan2;

	int subDataXLength = (const int)(*subsampledData)[0].size();

	double newSpacingX, newSpacingY;
	double meanX_std = 0;
	double meanY_std = 0;
	double meanErrorSquared = 0;
	double startT, stopT;

	if (dispIntermResults)
		printf("Transforming Data: .");

	// Preserve the original spatial relationships prior to interpolation scaling
	// Keep original input and output points in UTM.
	vector<double> x0	((*subsampledData)[0]);		//data
	vector<double> y0	((*subsampledData)[1]);
	//This is the same as xMeshVector (aka xInterpVector)
	dgrid xInterpLocs0	(*xMeshGrid +meanXSingle);	//interpolation grid locations
	dgrid yInterpLocs0	(*yMeshGrid +meanYSingle);

	//************************************************************************************
	//I. Remove overall trend from locations.
	//************************************************************************************
	//A. Remove trends in position array. 
	//1. First, shift data and interpolation locations to center of
	//   grid by removing the mean.
	//	 Center on output's center so subtract output's mean which we
	//	 already calculated for xSingleVector.
	// cpp's xi is cpp's subsampleData[0], but matlab's xi is cpp's xSingleVector
	(*subsampledData)[0] -= meanXSingle; //data
	(*subsampledData)[1] -= meanYSingle;

	// Interpolation grid locations; unique row and col values
	(*xSingleVector) -= meanXSingle; //interpolation grid
	(*ySingleVector) -= meanYSingle;

	//2. Scale grid by 1/std(x), where x(:, 1:2) = location of data points in meters, so
	//	 that the data point locations are not huge numbers - basically, -3 to 3).
	//i. Compute mean that is needed for the variance estimate.
	for (i = 0; i < subDataXLength; i++){
		meanX_std = meanX_std + (*subsampledData)[0][i];
		meanY_std = meanY_std + (*subsampledData)[1][i];
		meanErrorSquared += (*subsampledData)[4][i];
	}
	meanX_std /= (double)subDataXLength;
	meanY_std /= (double)subDataXLength;
	meanErrorSquared /= (double)subDataXLength;

	//ii. Compute the variance estimate
	double std_x = 0;
	double std_y = 0;
	for (i = 0; i < subDataXLength; i++){
		std_x = std_x + pow(((*subsampledData)[0][i] - meanX_std),2);
		std_y = std_y + pow(((*subsampledData)[1][i] - meanY_std),2);
	}
	std_x = std_x / (double)(subDataXLength - 1.00);
	std_y = std_y / (double)(subDataXLength - 1.00);
	std_x = sqrt(std_x);
	std_y = sqrt(std_y);

	if (dispIntermResults)
		printf(".");

	//iii. Scale the data and grid. L = diag(1./std_x);
	(*subsampledData)[0] /= std_x;	
	(*subsampledData)[1] /= std_y;
	(*xSingleVector) /= std_x;		//xi
	(*ySingleVector) /= std_y;		
	(*xMeshGrid) /= std_x;			//Xi
	(*yMeshGrid) /= std_y;			

	newSpacingX = gridSpacingX / std_x;
	newSpacingY = gridSpacingY / std_y;

	//************************************************************************************
	//II. Remove trend surface from the observations. It will be added back later.
	//************************************************************************************
	//A. Compute consistent window weights to pass to the linear regressions routine: regr_xzw.
	if (dispIntermResults)
		printf(".");

	double wtol = 0.0100;
	double constWeightsS3 = 0;
	vector<double> weights(subDataXLength, 2.00);
	consistentWeights(&(*subsampledData)[2], &(*subsampledData)[4], &wtol, &weights, &constWeightsS3);

	//B. Call regr_xzw.m to calculate a 2-D linear fit to the data set.
	//Initialize output variables from regr_xzw.m
	dvector btrend (3,0.00);
	dvector bi (3,0.00);

	//C.Call regr_xzw.m if variables have a variance; else, keep the padding with zeros
	vector<double> regrX(subDataXLength,1);
	if (dispIntermResults)
		printf(".");

	if(!(std_x <= 0 && std_y <= 0)){
		// do regression to remove a norm field
		// this is just for getting the data ready, so it is meant to be bullet proof, not statistically pure!
		regr_xzw(&regrX, &(*subsampledData)[0], &(*subsampledData)[1], &(*subsampledData)[2], &(*subsampledData)[3], &weights, &btrend, &bi);
	}
	regrX.clear();

	//D. Removes trend from the data (i.e. remove overall bias from the data).
	double zTrendValue;
	for (i = 0; i < subDataXLength; i++){
		zTrendValue = (1.0*btrend[0]) + ((*subsampledData)[0][i]*btrend[1]) + ((*subsampledData)[1][i]*btrend[2]);
		// compute deviations from trend
		(*subsampledData)[2][i] = (*subsampledData)[2][i] - zTrendValue; 
		//Multiply the sub-sampled data by Lx and Ly for use in interpPerturbations
	}

	//************************************************************************************
	//III. Compute tile sizes to break interpolation into smaller chunks.
	//************************************************************************************
	//A. Compute optimal tile size, which depends on the smoothing scale.
	//	First, determine the number of grid columns and rows
	//	and find the bounding box limits of the rescaled grid per Section III above.
	if (dispIntermResults)
		printf(".");

	double lk;
	double minXiTest = (*xSingleVector) [0];
	double maxXiTest = (*xSingleVector) [(*xSingleVector).size() - 1];
	double minYiTest = (*ySingleVector) [0];
	double maxYiTest = (*ySingleVector) [(*ySingleVector).size() - 1];

	//B. Compute (smoothing scale)/domain size ratio.
	if ((newSpacingX / (maxXiTest - minXiTest)) > (newSpacingY / (maxYiTest - minYiTest))){
		lk = (newSpacingX / (maxXiTest - minXiTest));
	}else{
		lk = (newSpacingY / (maxYiTest - minYiTest));
	}

	//C. Compute the optimal number of tiles, tiles in each dimension,
	//and approximate number of points per tile.
	double kopt = sqrt(((double)(*xSingleVector).size()) * ((double)(*ySingleVector).size()) * (1.00+lk));
	int kx  = (int)floor(sqrt(kopt));
	int ky  = (int)floor(sqrt(kopt));
	int nkx = (int)floor(((double)(*xSingleVector).size())/kx);
	int nky = (int)floor(((double)(*ySingleVector).size())/ky);

	//D. Compute efficiency, then display kopt and ropt to Command Window.
	if (dispIntermResults)
		printf("... Done Transforming Data\n");

	double ropt = (1.0/kopt) + kopt/((1.0+lk)*((double)(*xSingleVector).size())*((double)(*ySingleVector).size()));

	if (dispIntermResults)
		printf("\nNumber of computed tiles = %d, expected efficiency = %.6f\n\n", kx*ky, ropt);

	//************************************************************************************
	//IV. Begin interpolation of tiles.
	//************************************************************************************
	if (dispIntermResults)
		printf("\nENTERING INTERPOLATION ROUTINE:\n");

	/*This section uses the window size suggested in Calder's paper: ten times the largest sample spacing.*/
	double LMAX_x = 10.00*newSpacingX;
	double LMAX_y = 10.00*newSpacingY;
	int nDone = 0;

	vector< vector<int> > innerLoopIndexVector = vector< vector<int> >((int)ky);

	//A. Pre-compute the inner tiles
	//	By doing this once here we do not need to repeatedly do the same calculations in the inner loop below
	//	It takes slightly more memory but has a significant improvement in speed
	for (i = 0; i < ky; i++)
	{
		innerLoopIndexVector[i] = vector<int>(nky);
		for (j = 0; j < nky; j++)
		{
			(innerLoopIndexVector[i])[j] = j + (i * nky);
		}
	}
	if((innerLoopIndexVector[ky-1])[innerLoopIndexVector[ky-1].size() - 1] != (int)(*ySingleVector).size()-1){
		for (i = (innerLoopIndexVector[ky-1])[innerLoopIndexVector[ky-1].size() - 1] + 1; i < (int) (*ySingleVector).size(); i++){
			innerLoopIndexVector[ky-1].push_back(i);
		}
	}
#pragma endregion

	#pragma region --Scalec_Interp_Tile_Data Structures
	// create some temporary GRID objects
	dgrid outputDepth	= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);	// NaN); // Lewis Change, NaN to 0.00
	dgrid outputError	= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);	// NaN);
	dgrid outputNEi		= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 1.00);
	dgrid outputREi		= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);	// NaN);
	dgrid standardDev	= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);
	dgrid outputDepth0	= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);	// NaN);
	dgrid outputError0	= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);	// NaN);
//	dgrid standardDev0	= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);
	dgrid outputDepthK	= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);	// NaN);
	dgrid outputErrorK	= dgrid((*xMeshGrid).rows(), (*xMeshGrid).cols(), 0.00);	// NaN);

	double Lx = 1.0/newSpacingX;
	double Ly = 1.0/newSpacingY;

	// get estimators to perform
	string interpMethod = "BILINEAR";
	if(additionalOptions.find("-nnInterp")->second == 1)
		interpMethod = "NN";

	bool MSE = false;
	if(abs(additionalOptions.find("-mse")->second) == 1)
		MSE = true;

	bool PROP_UNCERT = false;
	if(abs(additionalOptions.find("-propUncert")->second) == 1)
		PROP_UNCERT = true;

	bool KALMAN = false;
	if(abs(additionalOptions.find("-kalman")->second) == 1)
		KALMAN = true;
	
	bool KRIGING = false;
	if(abs(additionalOptions["-kriging"] == 1))
		KRIGING = true;

	//B. Start with the new calling scheme here
	SCALEC_TILE_DATA scalecInterpTileData;

	scalecInterpTileData.subsampledData			= subsampledData;
	scalecInterpTileData.xSingleVector			= xSingleVector;
	scalecInterpTileData.ySingleVector			= ySingleVector;
	scalecInterpTileData.xMeshGrid				= xMeshGrid;
	scalecInterpTileData.yMeshGrid				= yMeshGrid;
	scalecInterpTileData.innerLoopIndexVector	= &innerLoopIndexVector;
	scalecInterpTileData.btrend					= &btrend;
	scalecInterpTileData.kernelName				= &kernelName;
	scalecInterpTileData.Lx						= &Lx;
	scalecInterpTileData.Ly						= &Ly;
	scalecInterpTileData.spacingX				= &newSpacingX;
	scalecInterpTileData.spacingY				= &newSpacingY;
	scalecInterpTileData.kx						= &kx;
	scalecInterpTileData.ky						= &ky;
	scalecInterpTileData.nkx					= &nkx;
	scalecInterpTileData.nky					= &nky;
	scalecInterpTileData.LMAX_x					= &LMAX_x;
	scalecInterpTileData.LMAX_y					= &LMAX_y;
	scalecInterpTileData.neitol					= &neitol;
	scalecInterpTileData.outputDepth			= &outputDepth;
	scalecInterpTileData.outputError			= &outputError;
	scalecInterpTileData.outputNEi				= &outputNEi;
	scalecInterpTileData.outputREi				= &outputREi;
	scalecInterpTileData.standardDev			= &standardDev;
	scalecInterpTileData.outputDepth0			= &outputDepth0;
	scalecInterpTileData.outputError0			= &outputError0;
//	scalecInterpTileData.standardDev0			= &standardDev0;
	scalecInterpTileData.outputDepthK			= &outputDepthK;
	scalecInterpTileData.outputErrorK			= &outputErrorK;
	scalecInterpTileData.x0						= &x0;				
	scalecInterpTileData.y0						= &y0;
	scalecInterpTileData.xInterpLocs0			= &xInterpLocs0;	
	scalecInterpTileData.yInterpLocs0			= &yInterpLocs0;
	scalecInterpTileData.interpMethod			= interpMethod;
	scalecInterpTileData.MSE					= MSE;
	scalecInterpTileData.PROP_UNCERT			= PROP_UNCERT;
	scalecInterpTileData.KALMAN					= KALMAN;
	scalecInterpTileData.KRIGING				= KRIGING;

	#pragma endregion

	startT = clock();

	//************************************************************************************
	//V. Check for multi-threading support and run the processing routines
	//************************************************************************************
	if (additionalOptions["-multiThread"] == 0)
	{
		if (additionalOptions["-kriging"] == 0)
		{
		//	int flag = scalecInterpTile_Process(&scalecInterpTileData, 0, 1);
			int flag = scalecInterpTile_ProcessA(&scalecInterpTileData, 0, 1);
		}
		else
		{ //Calling original function from mergeBathy v3.6 instead of the _Serial version for consistency.
		//	int flag = scalecInterpTile_ProcessKrig(&scalecInterpTileData, 0, 1);
			int flag = scalecInterpTile_ProcessA(&scalecInterpTileData, 0, 1);
		}
	}
	else
	{
		mbThreads mbT = mbThreads(additionalOptions["-multiThread"]);
		mbT.makeMBThread(&scalecInterpTileData);
		mbT.initMBThread_Tile(0);
//		mbT.initMBThread_Tile(additionalOptions["-kriging"]);
		mbT.joinMBThread();
		mbT.terminateMBThread();
	}

	innerLoopIndexVector.clear();
	xInterpLocs0.clear();
	yInterpLocs0.clear();
	x0.clear();
	y0.clear();

	stopT = clock();

	int nyi = (*xMeshGrid).rows();
	int nxi = (*xMeshGrid).cols();
	double t = (stopT-startT)/CLOCKS_PER_SEC;
	cout << "Interpolated "<< fix(nyi*nxi/t) << " points per second (tiled)."<< endl;

	cout<<"Store Output"<<endl;
	#pragma region --Store output
	//************************************************************************************
	//VI. Correct bad output data and uncertainty estimates from function
	//************************************************************************************
	//A. Remove NaNs from the data and replace with the trend surface
	int outputDepthSize		= outputDepth.size();
	(*xyzOut).depth			= vector<double>(outputDepthSize);
	(*xyzOut).error			= vector<double>(outputDepthSize);
	(*xyzOut).nEi			= vector<double>(outputDepthSize);
	(*xyzOut).rEi			= vector<double>(outputDepthSize);
	(*xyzOut).standardDev	= vector<double>(outputDepthSize);

	(*xyzOut).depth0		= vector<double>(outputDepthSize);
	(*xyzOut).error0		= vector<double>(outputDepthSize);
//	(*xyzOut).standardDev0	= vector<double>(outputDepthSize);
	(*xyzOut).depthK		= vector<double>(outputDepthSize);
	(*xyzOut).errorK		= vector<double>(outputDepthSize);

	double tgs0 = 1.0;
	double tgs1;
	double tgs2;

	//C. Store output from Section V in temporary column matrices.
	k = 0;
	for (i = 0; i < (const int)outputDepth.cols(); i++)
	{
		for (j = 0; j < (const int)outputDepth.rows(); j++)
		{
			(*xyzOut).depth[k]		= outputDepth(j,i);
			(*xyzOut).error[k]		= outputError(j,i);
			(*xyzOut).nEi[k]		= outputNEi(j,i);
			(*xyzOut).rEi[k]		= outputREi(j,i);
			(*xyzOut).standardDev[k] = standardDev(j,i);
			if(PROP_UNCERT)
			{
				(*xyzOut).depth0[k] = outputDepth0(j,i);
				(*xyzOut).error0[k] = outputError0(j,i);
			}
			if(KALMAN)
			{
				(*xyzOut).depthK[k] = outputDepthK(j,i);
				(*xyzOut).errorK[k] = outputErrorK(j,i);
			}

			tgs1 = (*xMeshGrid)(j,i);
			tgs2 = (*yMeshGrid)(j,i);
			
			//D. If NaN or 0 use the trend surface
			//	NaN works in Windows but Linux doesn't understand it so the value is forced to 0, that's why we need to check for it
			if ( !((*xyzOut).depth[k] >= 0 || (*xyzOut).depth[k] < 0) || ((*xyzOut).depth[k] == 0))
			{
				(*xyzOut).depth[k] = tgs0*btrend[0]+tgs1*btrend[1]+tgs2*btrend[2];
			}
			if ( !((*xyzOut).error[k] >= 0 || (*xyzOut).error[k] < 0) || ((*xyzOut).error[k] == 0))
			{
				(*xyzOut).error[k] = constWeightsS3 + meanErrorSquared;
				(*xyzOut).rEi[k] = constWeightsS3 + meanErrorSquared;
			}
			if(PROP_UNCERT)
			{
				if ( !((*xyzOut).depth0[k] >= 0 || (*xyzOut).depth0[k] < 0) || ((*xyzOut).depth0[k] == 0))
				{
					(*xyzOut).depth0[k] = tgs0*btrend[0]+tgs1*btrend[1]+tgs2*btrend[2];
				}
				if ( !((*xyzOut).error0[k] >= 0 || (*xyzOut).error0[k] < 0) || ((*xyzOut).error0[k] == 0))
				{
					(*xyzOut).error0[k] = constWeightsS3 + meanErrorSquared;
				}
			}
			if(KALMAN)
			{
				if ( !((*xyzOut).depthK[k] >= 0 || (*xyzOut).depthK[k] < 0) || ((*xyzOut).depthK[k] == 0))
				{
					(*xyzOut).depthK[k] = tgs0*btrend[0]+tgs1*btrend[1]+tgs2*btrend[2];
				}
			}

			k = k + 1;
		}
	}

	cout<<"DONE SCALECINTERPTILEDATA"<<endl;

	#pragma endregion

	btrend.clear();
	bi.clear();
	weights.clear();

	x0.clear();
	y0.clear();
	xInterpLocs0.clear();
	yInterpLocs0.clear();

	if (dispIntermResults)
	{
		printf("\nInterpolation Complete!\n");
		cout << "Time to Complete Interpolation: " << (stopT-startT)/CLOCKS_PER_SEC << " seconds" << endl << endl;
	}

	return 0;
}



int scalecInterpTile_ProcessA(SCALEC_TILE_DATA_POINTER stdp, const int curIterNum, const int numCores)
{
	//************************************************************************************
	// 0. Declare local variables and objects that are going to be used inside the loop
	//************************************************************************************
	dgrid* slopes;							// interpolation location grid slopes
	vector<double> subX0;					// scattered input data for Tri
	vector<double> subY0;
	vector<double> subXInterpLocs0_Vec;		// interpolation locations
	vector<double> subYInterpLocs0_Vec;
	dgrid Xiii(1,2);						// current interpolation location

	vector<int> outerLoopIndexVector;
	vector<int> idx;
	vector<double> subX_idy;				
	vector<double> subY_idy;
	vector<double> subZ_idy;
	vector<double> subE_idy;
	vector<double> subH_idy;
	vector<double> subV_idy;
	double xmin, xmax, ymin, ymax;
	double diffXi, diffYi;
	double dmin;

	vector <double> perturbWeights;
	const int subDataXLength = (const int)(*stdp->subsampledData)[0].size();
	int iliv_Loc;
	int oliv_Loc;
	int idySize;
	//int i, j;
	double tgs0 = 1.0;
	double tgs1, tgs2, tgs1_Compute, tgs2_Compute;
	double assnGridValue;

	PERTURBS perturb;
	perturb.kernelName = *(*stdp).kernelName;

	//Initialize Kriging vars
	vector<double> outputDepthKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	vector<double> outputErrorKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);

	vector<double> outputDepth0Krig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	vector<double> outputError0Krig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	
	vector<double> outputDepthKKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	vector<double> outputErrorKKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	double zKriged;
	double varZKriged;
	double z0Kriged;
	double varZ0Kriged;
	double zKKriged;
	double varZKKriged;
	double xGrid_indexKriged;
	double yGrid_indexKriged;
	vector<double> subX_idyKriged;
	vector<double> subY_idyKriged;
	vector<double> subX_idy0Kriged;
	vector<double> subY_idy0Kriged;
	vector<double> subX_idyKKriged;
	vector<double> subY_idyKKriged;
	vector<double> residualObservationsKrigedZ;
	vector<double> residualObservationsKrigedZ0;
	vector<double> residualObservationsKrigedZK;

	vector<double> subX_indexKriged;
	vector<double> subY_indexKriged;
	vector<double> subX_index0Kriged;
	vector<double> subY_index0Kriged;
	vector<double> subX_indexKKriged;
	vector<double> subY_indexKKriged;
	
	vector<double> xIndexKriged_Vector;
	vector<double> yIndexKriged_Vector;

	vector<double> twoGammaHatVector;
	vector<double> distanceVectorBinCenters;
	vector<double> aVectorFine;
	dgrid invGammaDArray;
	dgrid AGrid;
	double locSpacingX = (*stdp->spacingX);
	double locSpacingY = (*stdp->spacingY);
	int k=0;
	int iliv_Loc0, iliv_Loc1;
	int oliv_Loc0, oliv_Loc1;
	double tgs1_first, tgs2_first;
	double tgs1_last, tgs2_last;

	//************************************************************************************
	// I. Interpolate the data
	//************************************************************************************
	//loop x tiles, columns
	for (int outerLoop = curIterNum; outerLoop < (const int)(*stdp->kx); outerLoop+=numCores) 
	{
		outerLoopIndexVector = vector<int>((*stdp->nkx));
		//A. Get indices in XI, which gives x-locations, to interpolate on this loop.
		for (int i = 0; i < (const int)(*stdp->nkx); i++)
		{
			// indices to interpolate this time
			outerLoopIndexVector[i] = ( i + ((outerLoop) * (int)(*stdp->nkx))); 
		}
		if((outerLoop == (int)((*stdp->kx) - 1)) && (outerLoopIndexVector[outerLoopIndexVector.size() - 1] != (int)((*stdp->xSingleVector).size() - 1)))
		{
			for (int i = outerLoopIndexVector[outerLoopIndexVector.size() - 1] + 1; i < (const int)(*stdp->xSingleVector).size(); i++)
			{
				// catch the end here
				outerLoopIndexVector.push_back(i); 
			}
		}
		//B. Compute appropriate overlap between tiles horizontally and get the
		//	useful data. Recall that we scaled xi-array by 1./std(x). Then,
		//	get indices of x-coordinates in x-array for the data to be
		//	interpolated. Indices to be filtered further for useful y-coordinates.
		// find tile limits
		//xmin = (*stdp->xSingleVector)[(int)outerLoopIndexVector[0]] - (*stdp->LMAX_x); 
		//xmax = (*stdp->xSingleVector)[(int)outerLoopIndexVector[outerLoopIndexVector.size()-1]] + (*stdp->LMAX_x);
		
		/*Optimal overlap is half the distance of the tile according to Numerical Recipes.
		  This appears to run faster than LMAX but produces relatively the same results. -- SJZ*/
		double overlap = 0;
		overlap = ((*stdp->xSingleVector)[(int)outerLoopIndexVector[outerLoopIndexVector.size()-1]] - (*stdp->xSingleVector)[(int)outerLoopIndexVector[0]])/2;
		xmin = (*stdp->xSingleVector)[(int)outerLoopIndexVector[0]] - overlap; 
		xmax = (*stdp->xSingleVector)[(int)outerLoopIndexVector[outerLoopIndexVector.size()-1]] + overlap;
		
		for (int i = 0; i < subDataXLength; i++){
			if (((*stdp->subsampledData)[0][i] < xmax) && ((*stdp->subsampledData)[0][i] > xmin))
					idx.push_back(i);
		}
		#pragma region --idxSize
		if (idx.size() > 0) 
		{
			//loop y tiles, blocks in columns
			for (int innerLoop = 0; innerLoop < (const int)(*stdp->ky); innerLoop++) 
			{
				/*This section uses the window size suggested in Calder's paper: ten times the largest sample spacing.*/
				//C. Do the same as above but for the y coordinates
				/*ymin = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[0]] - (*stdp->LMAX_y);
				ymax = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[(*stdp->innerLoopIndexVector)[innerLoop].size()-1]] + (*stdp->LMAX_y);
				*/
				
				/*Optimal overlap is half the distance of the tile according to Numerical Recipes.
				  This appears to run faster than LMAX but produces relatively the same results. -- SJZ*/
				overlap = 0;
				overlap =  ( (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[(*stdp->innerLoopIndexVector)[innerLoop].size()-1]] - (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[0]] )/2;
				ymin = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[0]] - overlap;
				ymax = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[(*stdp->innerLoopIndexVector)[innerLoop].size()-1]] + overlap;
				
				//D. Get all of the original data that falls within the current tile to be interpolated
				for (int i = 0; i < (const int)idx.size(); i++)
				{
					if ( ((*stdp->subsampledData)[1][(int)idx[i]] < ymax) && ((*stdp->subsampledData)[1][(int)idx[i]] > ymin) )
					{
						subX_idy.push_back((*stdp->subsampledData)[0][idx[i]] * (*stdp->Lx));	
						subY_idy.push_back((*stdp->subsampledData)[1][idx[i]] * (*stdp->Ly));
						subZ_idy.push_back((*stdp->subsampledData)[2][idx[i]]);
						subE_idy.push_back((*stdp->subsampledData)[4][idx[i]]);
						subH_idy.push_back((*stdp->subsampledData)[5][idx[i]]);
						subV_idy.push_back((*stdp->subsampledData)[6][idx[i]]);
						
						if(stdp->KRIGING)
						{
							//subX_idyKriged.push_back((*stdp->subsampledData)[0][idx[i]]); //x for kriging
							//subY_idyKriged.push_back((*stdp->subsampledData)[1][idx[i]]);
							subX_idyKriged.push_back((*stdp->subsampledData)[0][idx[i]] * (*stdp->Lx)); //x for kriging SJZ
							subY_idyKriged.push_back((*stdp->subsampledData)[1][idx[i]] * (*stdp->Ly));
						}
						subX0.push_back((*stdp->x0)[idx[i]]);	//scattered input data sub-tile
						subY0.push_back((*stdp->y0)[idx[i]]);
					}
				}
				idySize = (const int)subX_idy.size();

				#pragma region idySize
				if (idySize > 2)
				{
					//E. Obtain a regular grid
					// get number of indices in tile. idyi, col. idxi, row.
					uint idyi = (uint)(*stdp->innerLoopIndexVector)[innerLoop].size(); 
					uint idxi = (uint)outerLoopIndexVector.size(); 

					// get starting indices of tile. i, col. j, row
					int i = ((*stdp->innerLoopIndexVector)[innerLoop])[0];
					int j = outerLoopIndexVector[0];

					// get sub-tile of grid interpolation locations.
					dgrid subXInterpLocs0;
					dgrid subYInterpLocs0;
					(*stdp->xInterpLocs0).subgrid(subXInterpLocs0, i, j, idyi, idxi);
					(*stdp->yInterpLocs0).subgrid(subYInterpLocs0, i, j, idyi, idxi);

					//F. Create Delaunay Tri from scattered input data sub-tile.
					Bathy_Grid new_bathyGrid = Bathy_Grid();
					new_bathyGrid.Construct_Tin(&subX0, &subY0, &subZ_idy, &subH_idy, &subV_idy);
					
					//G. Estimate depths at regular grid points
					InterpGrid		new_gmt			= InterpGrid(GMT);					
					vector<double>	subZInterpLocs0 = vector<double> ((subXInterpLocs0.vec()).size(), 0.00);
					new_gmt.estimate(&subXInterpLocs0.vec(), &subYInterpLocs0.vec(), &subZInterpLocs0, 1.96, 2.00, *stdp->Lx, new_bathyGrid.getTin(), stdp->interpMethod, "", "");

					//H. Find slopes at grid points
					slopes = new_gmt.getGrads()->getSlopeOut();

					//I. Find dmin
					subXInterpLocs0_Vec = subXInterpLocs0.vec();
					subYInterpLocs0_Vec = subYInterpLocs0.vec();
					double minXi		= abs(subXInterpLocs0_Vec[1] - subXInterpLocs0_Vec[0]);
					double minYi		= abs(subYInterpLocs0_Vec[1] - subYInterpLocs0_Vec[0]);
					dmin				= 0.0;
					for(int j = 1; j < (const int)subXInterpLocs0_Vec.size()-1; j++)
					{
						diffXi = abs(subXInterpLocs0_Vec[j+1] - subXInterpLocs0_Vec[j]);
						diffYi = abs(subYInterpLocs0_Vec[j+1] - subYInterpLocs0_Vec[j]);
						if(!(minXi > 0) && diffXi > 0)
							minXi = diffXi;
						else if(diffXi < minXi && diffXi > 0)
							minXi = diffXi;

						if(!(minYi > 0) && diffYi > 0)
							minYi = diffYi;
						else if(diffYi < minYi && diffYi > 0)
							minYi = diffYi;

						dmin = min(minXi, minYi);
					}

					if(dmin == 0)
						cerr << "dmin equals 0";

					//J. Pre-compute the weights, riVector, and aiVector across the whole interpolation plane
					perturbWeights = vector<double>(idySize,2);
					scalecInterpPerturbations_PreCompute(&subZ_idy, &subE_idy, &perturbWeights, &perturb);

					if(stdp->KRIGING)
					{
						#pragma region Interpolate --_Compute_ForKriging (the same as _Compute)
						iliv_Loc0 = ((*stdp->innerLoopIndexVector)[innerLoop])[0];	//first innerLoopIndexVector location
						oliv_Loc0 = outerLoopIndexVector[0];						//first outerLoopIndexVector location
						tgs1_first = (*stdp->xMeshGrid)(iliv_Loc0, oliv_Loc0);		//get first xValue for xa+yb+c
						tgs2_first = (*stdp->yMeshGrid)(iliv_Loc0, oliv_Loc0);		//get first yValue

						iliv_Loc1 = ((*stdp->innerLoopIndexVector)[innerLoop])[((*stdp->innerLoopIndexVector)[innerLoop]).size() - 1]; //first innerLoopIndexVector loc
						oliv_Loc1 = outerLoopIndexVector[outerLoopIndexVector.size()-1];	//first outerLoopIndexVector loc
						tgs1_last = (*stdp->xMeshGrid)(iliv_Loc1, oliv_Loc1);				//get last xValue for xa+yb+c
						tgs2_last = (*stdp->yMeshGrid)(iliv_Loc1, oliv_Loc1);				//get last yValue
						vector<double> slopesVec2 = vector<double>((*stdp->subsampledData)[0].size(), 1.00);

						//F. Get interpolation values at observation locations
						//We are doing input values here!
						for (int i = 0; i < (const int)idySize; i++)
						{
							Xiii(0,0) = subX0[i];				//doesn't matter cause we have a regular grid
							Xiii(0,1) = subX0[i];

							tgs1_Compute = subX_idy[i];			//current input xValue in tile
							tgs2_Compute = subY_idy[i];
							perturb.perturbationZKriged	= 0.0;
							perturb.perturbationEKriged	= 1.0;
							perturb.perturbationNEiKriged = 1.0;
							perturb.perturbationREiKriged = 1.0;
							if(stdp->PROP_UNCERT)
							{
								perturb.perturbationZ0Kriged = 0.0;
								perturb.perturbationE0Kriged = 1.0;
							}
							if(stdp->KALMAN)
							{
								perturb.perturbationZKKriged = 0.0;
								perturb.perturbationEKKriged = 1.0;
							}

							//We are doing the input data 
							//This is calls the original kriging function which does not have the new functionality and bug fixes.
							//scalecInterpPerturbations_Compute_ForKriging(&subX_idx, &subY_idx, &z_idx, &e_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, &perturb.riVector, &perturb.aiVector, (*sdp->neitol), &perturb.perturbationZKriged, &perturbationEKriged, &perturb.perturbationNEiKriged, &perturb.perturbationREiKriged, &standardDevKriged);

							//Current function
							scalecInterpPerturbations_Compute(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &subH_idy, &subV_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*stdp->neitol), dmin, (slopesVec2)[i], &subX0, &subY0, &Xiii, &perturb, stdp->MSE, stdp->PROP_UNCERT, stdp->KALMAN, stdp->KRIGING);
						
							//COMPARE THE DIFFERENCE BETWEEN SUBTILES: xMeshGrid AND subY_idyKriged AND 
							//Keep (un-scaled subtile input) kriged indices within MeshGrid (interpolation) subtile and remove kriged depth. 
							// We want to capture data in the tile that is not part of the
							// overlap.
							if ( (subY_idyKriged[i] > tgs2_first) && (subY_idyKriged[i] < tgs2_last)
							  && (subX_idyKriged[i] > tgs1_first) && (subX_idyKriged[i] < tgs1_last))
							{
								residualObservationsKrigedZ.push_back(subZ_idy[i] - perturb.perturbationZKriged);
								subX_indexKriged.push_back(subX_idyKriged[i]);
								subY_indexKriged.push_back(subY_idyKriged[i]);			

								//NEED TO CHECK TO SEE IF subX_indexKriged ARE THE SAME FOR ALL
								if(stdp->PROP_UNCERT)
									residualObservationsKrigedZ0.push_back(subZ_idy[i] - perturb.perturbationZ0Kriged);
								if(stdp->KALMAN)
									residualObservationsKrigedZK.push_back(subZ_idy[i] - perturb.perturbationZKKriged);
							}
						} //idySize
						#pragma endregion Interpolate _Compute_ForKriging

						#pragma region --ordinaryKrigingOfResiduals_PreCompute, Subtile Kriged Indices if if too large 
						if (subX_indexKriged.size() >= 15)
						{
							int KRIG_SUBTILE_FLAG = 0;
							//G. We need to subtile the kriged indices otherwise the matrix inversion takes too long -- TODO
							if(subX_indexKriged.size() > KRIGED_SIZE_THRESHOLD)
							{
								KRIG_SUBTILE_FLAG = 1;
								#pragma region --Subtile Kriged Indices
								try
								{
									k = 0;
									xIndexKriged_Vector = vector<double>(outerLoopIndexVector.size()*((*stdp->innerLoopIndexVector)[innerLoop]).size());
									yIndexKriged_Vector = vector<double>(outerLoopIndexVector.size()*((*stdp->innerLoopIndexVector)[innerLoop]).size());
									for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
									{
										for (int j = 0; j < (const int)((*stdp->innerLoopIndexVector)[innerLoop]).size(); j++)
										{
											iliv_Loc = ((*stdp->innerLoopIndexVector)[innerLoop])[j];
											oliv_Loc = outerLoopIndexVector[i];
											xIndexKriged_Vector[k] = (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);	//kriged indices subtile (interpolation locations)
											yIndexKriged_Vector[k] = (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);
											k++;

										}
									}

									//cout << "To Tile Call" << endl;
									ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKrig, &outputErrorKrig);

									if(stdp->PROP_UNCERT)
										ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepth0Krig, &outputError0Krig);
									if(stdp->KALMAN)
										ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKKrig, &outputErrorKKrig);

									xIndexKriged_Vector.clear();
									yIndexKriged_Vector.clear();
								}
								catch (exception &e) {
									cout << "An exception occurred. Exception Thrown: " << e.what() << '\n';
									cout << "Kriging Subtiles failed.  Will try kriging without subtiles.  This will take longer." << endl;
									KRIG_SUBTILE_FLAG = 0;
								}
								#pragma endregion Subtile Kriged Indices
							}
							if(!KRIG_SUBTILE_FLAG)
							{
								#pragma region --Do Not Subtile Kriged Indices
								//small enough to do all at once
								ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
							
								if(stdp->PROP_UNCERT)
									ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
								if(stdp->KALMAN)
									ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);

								k = 0;
								//H. Krig the data
								for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
								{
									for (int j = 0; j < (const int)((*stdp->innerLoopIndexVector)[innerLoop]).size(); j++)
									{
										zKriged = 0.00;
										varZKriged = 0.00;
										iliv_Loc = ((*stdp->innerLoopIndexVector)[innerLoop])[j];
										oliv_Loc = outerLoopIndexVector[i];
										xGrid_indexKriged = (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);
										yGrid_indexKriged = (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);
									
										//Compute the residual value
										ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);
									
										outputDepthKrig[k] = zKriged;
										outputErrorKrig[k] = varZKriged;

										if(stdp->PROP_UNCERT)
										{
											z0Kriged = 0.00;
											varZ0Kriged = 0.00;

											ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &z0Kriged, &varZ0Kriged);
									
											outputDepth0Krig[k] = z0Kriged;
											outputError0Krig[k] = varZ0Kriged;
										}
										if(stdp->KALMAN)
										{
											zKKriged = 0.00;
											varZKKriged = 0.00;
										
											ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKKriged, &varZKKriged);
									
											outputDepthKKrig[k] = zKKriged;
											outputErrorKKrig[k] = varZKKriged;
										}

										k++;
									} // innerLoopIndexVector
								} // outerLoopInde
								twoGammaHatVector.clear();
								distanceVectorBinCenters.clear();
								aVectorFine.clear();
								invGammaDArray.clear();
								AGrid.clear();
								#pragma endregion Do Not Subtile Kriged Indices
							} 
						} 
						#pragma endregion ordinaryKrigingOfResiduals_PreCompute
					}

					//Continue normal processing
					double outDepthKrig;
					double outErrorKrig;
					double outDepth0Krig;
					double outError0Krig;
					double outDepthKKrig;
					double outErrorKKrig;
					#pragma region Interpolate 
					k=0;
					//K. Do the interpolation over each independent point
					for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
					{
						for (int j = 0; j < (const int)(*stdp->innerLoopIndexVector)[innerLoop].size(); j++)
						{
							// current interpolation location
							Xiii(0,0) = subXInterpLocs0(j, i);
							Xiii(0,1) = subYInterpLocs0(j, i);
						
							iliv_Loc		= ((*stdp->innerLoopIndexVector)[innerLoop])[j];
							oliv_Loc		= outerLoopIndexVector[i];
							tgs1			= (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);	//get current xValue
							tgs2			= (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);	//get current yValue
							tgs1_Compute	= tgs1 * (*stdp->Lx);						//scale xValue
							tgs2_Compute	= tgs2 * (*stdp->Ly);						//scale yValue

							perturb.perturbationZ	= 0.0;
							perturb.perturbationE	= 1.0;
							perturb.perturbationNEi = 1.0;
							perturb.perturbationREi = 1.0;
							if(stdp->PROP_UNCERT)
							{
								perturb.perturbationZ0 = 0.0;
								perturb.perturbationE0 = 1.0;
							}
							if(stdp->KALMAN)
							{
								perturb.perturbationZK = 0.0;
								perturb.perturbationEK = 1.0;
							}
							if(stdp->KRIGING)
							{
								outDepthKrig  = outputDepthKrig[k];
								outErrorKrig  = outputErrorKrig[k];
								if(stdp->PROP_UNCERT)
								{
									outDepth0Krig = outputDepth0Krig[k];
									outError0Krig = outputError0Krig[k];
								}
								if(stdp->KALMAN)
								{
									outDepthKKrig = outputDepthKKrig[k];
									outErrorKKrig = outputErrorKKrig[k];
								}
							}
							else
							{
								outDepthKrig  = 0;
								outErrorKrig  = 0;
								outDepth0Krig = 0;
								outError0Krig = 0;
								outDepthKKrig = 0;
								outErrorKKrig = 0;
							}
							//L. Compute the value
							scalecInterpPerturbations_Compute(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &subH_idy, &subV_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*stdp->neitol), dmin, (*slopes)(j,i), &subX0, &subY0, &Xiii, &perturb, stdp->MSE, stdp->PROP_UNCERT, stdp->KALMAN, 0); //0 because were not kriging the residuals here... that's done above

							//M. Put trend back into this tile.
							assnGridValue = tgs0*(*stdp->btrend)[0]+tgs1*(*stdp->btrend)[1]+tgs2*(*stdp->btrend)[2];
							if((*stdp->outputNEi)(iliv_Loc,oliv_Loc)!=1)
								int stop=1;

							(*stdp->outputDepth)(iliv_Loc,oliv_Loc) = perturb.perturbationZ + assnGridValue + outDepthKrig;
							(*stdp->outputError)(iliv_Loc,oliv_Loc) = perturb.perturbationE + outErrorKrig;
							(*stdp->outputNEi)(iliv_Loc,oliv_Loc)	= perturb.perturbationNEi;
							(*stdp->outputREi)(iliv_Loc,oliv_Loc)	= perturb.perturbationREi;
							(*stdp->standardDev)(iliv_Loc,oliv_Loc) = perturb.standardDev2;
							if(stdp->PROP_UNCERT)
							{
								(*stdp->outputDepth0)(iliv_Loc,oliv_Loc) = perturb.perturbationZ0 + assnGridValue + outDepth0Krig;
								(*stdp->outputError0)(iliv_Loc,oliv_Loc) = perturb.perturbationE0 + outError0Krig;
		//						(*stdp->standardDev0)(iliv_Loc,oliv_Loc) = perturb.standardDev20;
							}
							if(stdp->KALMAN)
							{
								(*stdp->outputDepthK)(iliv_Loc,oliv_Loc) = perturb.perturbationZK + assnGridValue + outDepthKKrig;
								(*stdp->outputErrorK)(iliv_Loc,oliv_Loc) = perturb.perturbationEK + outErrorKKrig;
							}
							k++;
						}
					}
					#pragma endregion Interpolate
					//required clears!
					perturbWeights.clear();
					perturb.riVector.clear();
					perturb.aiVector.clear();
					new_bathyGrid.clear();
					new_gmt.clear();
					subXInterpLocs0.clear();
					subYInterpLocs0.clear();
					slopes = NULL;
					if(stdp->KRIGING)
					{
						residualObservationsKrigedZ.clear();
						subX_indexKriged.clear();
						subY_indexKriged.clear();
						if(stdp->PROP_UNCERT)
							residualObservationsKrigedZ0.clear();
						if(stdp->KALMAN)
							residualObservationsKrigedZK.clear();
					}
				} //idySize
				#pragma endregion idySize
				//required clears!
				subX_idy.clear();
				subY_idy.clear();
				subZ_idy.clear();
				subE_idy.clear();
				subH_idy.clear();
				subV_idy.clear();
				subX0.clear();
				subY0.clear();
				if(stdp->KRIGING)
				{
					subX_idyKriged.clear();
					subY_idyKriged.clear();
				}
			} //innerloop
		} //idxSize
		#pragma endregion idxSize
		outerLoopIndexVector.clear();
		idx.clear();
	} //outerloop
	if(stdp->KRIGING)
	{
		outputDepthKrig.clear();
		outputErrorKrig.clear();
		if(stdp->PROP_UNCERT)
		{
			outputDepth0Krig.clear();
			outputError0Krig.clear();
		}
		if(stdp->KALMAN)
		{
			outputDepthKKrig.clear();
			outputErrorKKrig.clear();
		}
	}
	return SUCCESS;
}













int scalecInterpTile_Process(SCALEC_TILE_DATA_POINTER stdp, const int curIterNum, const int numCores)
{
	//************************************************************************************
	// 0. Declare local variables and objects that are going to be used inside the loop
	//************************************************************************************
	dgrid* slopes;							// interpolation location grid slopes
	vector<double> subX0;					// scattered input data for Tri
	vector<double> subY0;
	vector<double> subXInterpLocs0_Vec;		// interpolation locations
	vector<double> subYInterpLocs0_Vec;
	dgrid Xiii(1,2);						// current interpolation location

	vector<int> outerLoopIndexVector;
	vector<int> idx;
	vector<double> subX_idy;				
	vector<double> subY_idy;
	vector<double> subZ_idy;
	vector<double> subE_idy;
	vector<double> subH_idy;
	vector<double> subV_idy;
	double xmin, xmax, ymin, ymax;
	double diffXi, diffYi;
	double dmin;

	vector <double> perturbWeights;
	const int subDataXLength = (const int)(*stdp->subsampledData)[0].size();
	int iliv_Loc;
	int oliv_Loc;
	int idySize;
	//int i, j;
	double tgs0 = 1.0;
	double tgs1, tgs2, tgs1_Compute, tgs2_Compute;
	double assnGridValue;

	PERTURBS perturb;
	perturb.kernelName = *(*stdp).kernelName;

	//************************************************************************************
	// I. Interpolate the data
	//************************************************************************************
	//loop x tiles, columns
	for (int outerLoop = curIterNum; outerLoop < (const int)(*stdp->kx); outerLoop+=numCores) 
	{
		outerLoopIndexVector = vector<int>((*stdp->nkx));
		//A. Get indices in XI, which gives x-locations, to interpolate on this loop.
		for (int i = 0; i < (const int)(*stdp->nkx); i++)
		{
			// indices to interpolate this time
			outerLoopIndexVector[i] = ( i + ((outerLoop) * (int)(*stdp->nkx))); 
		}
		if((outerLoop == (int)((*stdp->kx) - 1)) && (outerLoopIndexVector[outerLoopIndexVector.size() - 1] != (int)((*stdp->xSingleVector).size() - 1)))
		{
			for (int i = outerLoopIndexVector[outerLoopIndexVector.size() - 1] + 1; i < (const int)(*stdp->xSingleVector).size(); i++)
			{
				// catch the end here
				outerLoopIndexVector.push_back(i); 
			}
		}
		//B. Compute appropriate overlap between tiles horizontally and get the
		//	useful data. Recall that we scaled xi-array by 1./std(x). Then,
		//	get indices of x-coordinates in x-array for the data to be
		//	interpolated. Indices to be filtered further for useful y-coordinates.
		// find tile limits
		xmin = (*stdp->xSingleVector)[(int)outerLoopIndexVector[0]] - (*stdp->LMAX_x); 
		xmax = (*stdp->xSingleVector)[(int)outerLoopIndexVector[outerLoopIndexVector.size()-1]] + (*stdp->LMAX_x);

		for (int i = 0; i < subDataXLength; i++){
			if (((*stdp->subsampledData)[0][i] < xmax) && ((*stdp->subsampledData)[0][i] > xmin))
					idx.push_back(i);
		}
		#pragma region --idxSize
		if (idx.size() > 0) 
		{
			//loop y tiles, blocks in columns
			for (int innerLoop = 0; innerLoop < (const int)(*stdp->ky); innerLoop++) 
			{
				//C. Do the same as above but for the y coordinates
				ymin = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[0]] - (*stdp->LMAX_y);
				ymax = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[(*stdp->innerLoopIndexVector)[innerLoop].size()-1]] + (*stdp->LMAX_y);

				//D. Get all of the original data that falls within the current tile to be interpolated
				for (int i = 0; i < (const int)idx.size(); i++)
				{
					if ( ((*stdp->subsampledData)[1][(int)idx[i]] < ymax) && ((*stdp->subsampledData)[1][(int)idx[i]] > ymin) )
					{
						subX_idy.push_back((*stdp->subsampledData)[0][idx[i]] * (*stdp->Lx));	
						subY_idy.push_back((*stdp->subsampledData)[1][idx[i]] * (*stdp->Ly));
						subZ_idy.push_back((*stdp->subsampledData)[2][idx[i]]);
						subE_idy.push_back((*stdp->subsampledData)[4][idx[i]]);
						subH_idy.push_back((*stdp->subsampledData)[5][idx[i]]);
						subV_idy.push_back((*stdp->subsampledData)[6][idx[i]]);
						
						subX0.push_back((*stdp->x0)[idx[i]]);	//scattered input data sub-tile
						subY0.push_back((*stdp->y0)[idx[i]]);
					}
				}
				idySize = (const int)subX_idy.size();

				#pragma region idySize
				if (idySize > 2)
				{
					//E. Obtain a regular grid
					// get number of indices in tile. idyi, col. idxi, row.
					uint idyi = (uint)(*stdp->innerLoopIndexVector)[innerLoop].size(); 
					uint idxi = (uint)outerLoopIndexVector.size(); 

					// get starting indices of tile. i, col. j, row
					int i = ((*stdp->innerLoopIndexVector)[innerLoop])[0];
					int j = outerLoopIndexVector[0];

					// get sub-tile of grid interpolation locations.
					dgrid subXInterpLocs0;
					dgrid subYInterpLocs0;
					(*stdp->xInterpLocs0).subgrid(subXInterpLocs0, i, j, idyi, idxi);
					(*stdp->yInterpLocs0).subgrid(subYInterpLocs0, i, j, idyi, idxi);

					//F. Create Delaunay Tri from scattered input data sub-tile.
					Bathy_Grid new_bathyGrid = Bathy_Grid();
					new_bathyGrid.Construct_Tin(&subX0, &subY0, &subZ_idy, &subH_idy, &subV_idy);
					
					//G. Estimate depths at regular grid points
					InterpGrid		new_gmt			= InterpGrid(GMT);					
					vector<double>	subZInterpLocs0 = vector<double> ((subXInterpLocs0.vec()).size(), 0.00);
					new_gmt.estimate(&subXInterpLocs0.vec(), &subYInterpLocs0.vec(), &subZInterpLocs0, 1.96, 2.00, *stdp->Lx, new_bathyGrid.getTin(), stdp->interpMethod, "", "");

					//H. Find slopes at grid points
					slopes = new_gmt.getGrads()->getSlopeOut();

					//I. Find dmin
					subXInterpLocs0_Vec = subXInterpLocs0.vec();
					subYInterpLocs0_Vec = subYInterpLocs0.vec();
					double minXi		= abs(subXInterpLocs0_Vec[1] - subXInterpLocs0_Vec[0]);
					double minYi		= abs(subYInterpLocs0_Vec[1] - subYInterpLocs0_Vec[0]);
					dmin				= 0.0;
					for(int j = 1; j < (const int)subXInterpLocs0_Vec.size()-1; j++)
					{
						diffXi = abs(subXInterpLocs0_Vec[j+1] - subXInterpLocs0_Vec[j]);
						diffYi = abs(subYInterpLocs0_Vec[j+1] - subYInterpLocs0_Vec[j]);
						if(!(minXi > 0) && diffXi > 0)
							minXi = diffXi;
						else if(diffXi < minXi && diffXi > 0)
							minXi = diffXi;

						if(!(minYi > 0) && diffYi > 0)
							minYi = diffYi;
						else if(diffYi < minYi && diffYi > 0)
							minYi = diffYi;

						dmin = min(minXi, minYi);
					}

					if(dmin == 0)
						cerr << "dmin equals 0";

					//J. Pre-compute the weights, riVector, and aiVector across the whole interpolation plane
					perturbWeights = vector<double>(idySize,2);
					scalecInterpPerturbations_PreCompute(&subZ_idy, &subE_idy, &perturbWeights, &perturb);

					#pragma region Interpolate 
					//K. Do the interpolation over each independent point
					for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
					{
						for (int j = 0; j < (const int)(*stdp->innerLoopIndexVector)[innerLoop].size(); j++)
						{
							// current interpolation location
							Xiii(0,0) = subXInterpLocs0(j, i);
							Xiii(0,1) = subYInterpLocs0(j, i);
						
							iliv_Loc		= ((*stdp->innerLoopIndexVector)[innerLoop])[j];
							oliv_Loc		= outerLoopIndexVector[i];
							tgs1			= (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);	//get current xValue
							tgs2			= (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);	//get current yValue
							tgs1_Compute	= tgs1 * (*stdp->Lx);						//scale xValue
							tgs2_Compute	= tgs2 * (*stdp->Ly);						//scale yValue

							perturb.perturbationZ	= 0.0;
							perturb.perturbationE	= 1.0;
							perturb.perturbationNEi = 1.0;
							perturb.perturbationREi = 1.0;
							if(stdp->PROP_UNCERT)
							{
								perturb.perturbationZ0 = 0.0;
								perturb.perturbationE0 = 1.0;
							}
							if(stdp->KALMAN)
							{
								perturb.perturbationZK = 0.0;
								perturb.perturbationEK = 1.0;
							}

							bool KRIGING = 0;
							//L. Compute the value
							scalecInterpPerturbations_Compute(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &subH_idy, &subV_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*stdp->neitol), dmin, (*slopes)(j,i), &subX0, &subY0, &Xiii, &perturb, stdp->MSE, stdp->PROP_UNCERT, stdp->KALMAN, KRIGING);

							//M. Put trend back into this tile.
							assnGridValue = tgs0*(*stdp->btrend)[0]+tgs1*(*stdp->btrend)[1]+tgs2*(*stdp->btrend)[2];
							if((*stdp->outputNEi)(iliv_Loc,oliv_Loc)!=1)
								int stop=1;

							(*stdp->outputDepth)(iliv_Loc,oliv_Loc) = perturb.perturbationZ + assnGridValue;
							(*stdp->outputError)(iliv_Loc,oliv_Loc) = perturb.perturbationE;
							(*stdp->outputNEi)(iliv_Loc,oliv_Loc)	= perturb.perturbationNEi;
							(*stdp->outputREi)(iliv_Loc,oliv_Loc)	= perturb.perturbationREi;
							(*stdp->standardDev)(iliv_Loc,oliv_Loc) = perturb.standardDev2;
							if(stdp->PROP_UNCERT)
							{
								(*stdp->outputDepth0)(iliv_Loc,oliv_Loc) = perturb.perturbationZ0 + assnGridValue;
								(*stdp->outputError0)(iliv_Loc,oliv_Loc) = perturb.perturbationE0;
	//							(*stdp->standardDev0)(iliv_Loc,oliv_Loc) = perturb.standardDev20;
							}
							if(stdp->KALMAN)
							{
								(*stdp->outputDepthK)(iliv_Loc,oliv_Loc) = perturb.perturbationZK + assnGridValue;
								(*stdp->outputErrorK)(iliv_Loc,oliv_Loc) = perturb.perturbationEK;
							}
						}
					}
					#pragma endregion Interpolate
					//required clears!
					perturbWeights.clear();
					perturb.riVector.clear();
					perturb.aiVector.clear();
					new_bathyGrid.clear();
					new_gmt.clear();
					subXInterpLocs0.clear();
					subYInterpLocs0.clear();
					slopes = NULL;
				} //idySize
				#pragma endregion idySize
				//required clears!
				subX_idy.clear();
				subY_idy.clear();
				subZ_idy.clear();
				subE_idy.clear();
				subH_idy.clear();
				subV_idy.clear();
				subX0.clear();
				subY0.clear();
			} //innerloop
		} //idxSize
		#pragma endregion idxSize
		outerLoopIndexVector.clear();
		idx.clear();
	} //outerloop
	return SUCCESS;
}

//This function is an updated version of the one from mergeBathy_v3.7.1_Paul.  It was modified to call the current function Perturbations_Compute and passes KRIGING=1 when kriging and KRGING=0 when interpolating instead of having the 2 functions in older versions. This includes all changes and improvements.
int scalecInterpTile_ProcessKrig(SCALEC_TILE_DATA_POINTER stdp, const int curIterNum, const int numCores)
{
	//************************************************************************************
	// 0. Declare local variables and objects that are going to be used inside the loop
	//************************************************************************************
	vector<double> outputDepthKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	vector<double> outputErrorKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);

	vector<double> outputDepth0Krig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	vector<double> outputError0Krig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	
	vector<double> outputDepthKKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
	vector<double> outputErrorKKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);

	PERTURBS perturb;
	perturb.kernelName=*(*stdp).kernelName;
	bool KRIGING = true;

	dgrid* slopes;							// interpolation location grid slopes
	vector<double>* slopesVec;							
	vector<double> slopesVec2;							
	vector<double> subXInterpLocs0_Vec;		// interpolation locations
	vector<double> subYInterpLocs0_Vec;
	dgrid Xiii(1,2);						// current interpolation location
	double diffXi, diffYi;
	double dmin;

	//double perturbationZ;
	//double perturbationE;
	//double perturbationNEi;
	//double perturbationREi;
	//double standardDev;

	//double perturbationZKriged;
	//double perturbationEKriged;
	//double perturbationNEiKriged;
	//double perturbationREiKriged;
	//double thetaM, phi_2;
	double zKriged;
	double varZKriged;
	double z0Kriged;
	double varZ0Kriged;
	double zKKriged;
	double varZKKriged;
	double xGrid_indexKriged;
	double yGrid_indexKriged;
//	double standardDevKriged;

	vector<int> outerLoopIndexVector;
	vector<int> idx;
	vector<double> subX_idy;
	vector<double> subY_idy;
	vector<double> subZ_idy;
	vector<double> subE_idy;
	vector<double> subH_idy;
	vector<double> subV_idy;
	vector<double> subX0;			// scattered input data for Tri
	vector<double> subY0;
	vector<double> subX_idyKriged;
	vector<double> subY_idyKriged;
	vector<double> subX_idy0Kriged;
	vector<double> subY_idy0Kriged;
	vector<double> subX_idyKKriged;
	vector<double> subY_idyKKriged;
	vector<double> residualObservationsKrigedZ;
	vector<double> residualObservationsKrigedZ0;
	vector<double> residualObservationsKrigedZK;

	vector<double> subX_indexKriged;
	vector<double> subY_indexKriged;
	vector<double> subX_index0Kriged;
	vector<double> subY_index0Kriged;
	vector<double> subX_indexKKriged;
	vector<double> subY_indexKKriged;
	
	vector<double> xIndexKriged_Vector;
	vector<double> yIndexKriged_Vector;

	vector<double> twoGammaHatVector;
	vector<double> distanceVectorBinCenters;
	vector<double> aVectorFine;
	dgrid invGammaDArray;
	dgrid AGrid;


	double xmin, xmax, ymin, ymax;
//	int perturbSize;

	vector <double> perturbWeights;
	//vector <double> riVector;
	//vector <double> aiVector;

	const size_t subDataXLength = (*stdp->subsampledData)[0].size();
	int iliv_Loc, iliv_Loc0, iliv_Loc1;
	int oliv_Loc, oliv_Loc0, oliv_Loc1;
	size_t idySize;
	int /*i, j, */k;
	double tgs0 = 1.0;
	double tgs1, tgs2, tgs1_Compute, tgs2_Compute;
	double tgs1_first, tgs2_first;
	double tgs1_last, tgs2_last;

	double assnGridValue;
	double locSpacingX = (*stdp->spacingX);
	double locSpacingY = (*stdp->spacingY);

	//************************************************************************************
	// I. Interpolate the data
	//************************************************************************************
	//loop x tiles, columns
	for (int outerLoop = curIterNum; outerLoop < (const int)(*stdp->kx); outerLoop+=numCores)
	{
		outerLoopIndexVector = vector<int>((*stdp->nkx));
		//A. Get indices in XI, which gives x-locations, to interpolate on this loop.
		for (int i = 0; i < (const int)(*stdp->nkx); i++)
		{
			//indices to interpolate this time
			outerLoopIndexVector[i] = ( i + ((outerLoop) * (int)(*stdp->nkx))); 
		}
		if((outerLoop == (int)((*stdp->kx) - 1)) && (outerLoopIndexVector[outerLoopIndexVector.size() - 1] != (int)((*stdp->xSingleVector).size() - 1)))
		{
			for (int i = outerLoopIndexVector[outerLoopIndexVector.size() - 1] + 1; i < (const int)(*stdp->xSingleVector).size(); i++)
			{
				//catch the end here
				outerLoopIndexVector.push_back(i); 
			}
		}
		//B. Compute appropriate overlap between tiles horizontally and get the 
		//useful data. Recall that we scaled xi-array by 1./std(x). Then, 
		//get indices of x-coordinates in x-array for the data to be
		//interpolated. Indices to be filtered further for useful y-coordinates.
		//find tile limits
		xmin = (*stdp->xSingleVector)[(int)outerLoopIndexVector[0]] - (*stdp->LMAX_x);
		xmax = (*stdp->xSingleVector)[(int)outerLoopIndexVector[outerLoopIndexVector.size()-1]] + (*stdp->LMAX_x);

		for (int i = 0; i < (const int)subDataXLength; i++)
		{
			if (((*stdp->subsampledData)[0][i] < xmax) && ((*stdp->subsampledData)[0][i] > xmin))
			{
				idx.push_back(i);
			}
		}

		#pragma region --idxSize
		if (idx.size() > 0)
		{
			//loop y tiles, blocks in columns
			for (int innerLoop = 0; innerLoop < (const int)(*stdp->ky); innerLoop++)
			{
				//C. Do the same as above but for the y coordinates
				ymin = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[0]] - (*stdp->LMAX_y);
//				ymax = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[((*stdp->innerLoopIndexVector)[innerLoop]).size()-1]] + (*stdp->LMAX_y);
				ymax = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[(*stdp->innerLoopIndexVector)[innerLoop].size()-1]] + (*stdp->LMAX_y);
				for (int i = 0; i < (const int)idx.size(); i++)
				{
					//D. Get all of the original input data that falls within the current tile to use for interpolation
					if ( ((*stdp->subsampledData)[1][(int)idx[i]] < ymax) && ((*stdp->subsampledData)[1][(int)idx[i]] > ymin) ){
						subX_idy.push_back((*stdp->subsampledData)[0][idx[i]] * (*stdp->Lx));
						subY_idy.push_back((*stdp->subsampledData)[1][idx[i]] * (*stdp->Ly));
						subZ_idy.push_back((*stdp->subsampledData)[2][idx[i]]);
						subE_idy.push_back((*stdp->subsampledData)[4][idx[i]]);
						subH_idy.push_back((*stdp->subsampledData)[5][idx[i]]);
						subV_idy.push_back((*stdp->subsampledData)[6][idx[i]]);
						
						//subX_idyKriged.push_back((*stdp->subsampledData)[0][idx[i]]); //x for kriging
						//subY_idyKriged.push_back((*stdp->subsampledData)[1][idx[i]]);
						subX_idyKriged.push_back((*stdp->subsampledData)[0][idx[i]] * (*stdp->Lx)); //x for kriging SJZ
						subY_idyKriged.push_back((*stdp->subsampledData)[1][idx[i]] * (*stdp->Ly));

						subX0.push_back((*stdp->x0)[idx[i]]);	//scattered input data sub-tile
						subY0.push_back((*stdp->y0)[idx[i]]);
					}
				}
				idySize = subX_idy.size();

				#pragma region idySize
				if (idySize > 2)
				{
					//GradientGrid gradgrid;// = GradientGrid();
					//gradgrid.calc_GradientGrid((*stdp->subsampledData)[0], (*stdp->subsampledData)[1], (*stdp->subsampledData)[2]);
					//slopesVec = gradgrid.getSlopeOut_Vector();
					slopesVec2 = vector<double>((*stdp->subsampledData)[0].size(), 1.00);

					#pragma region Get regular grid to interpolate
					//E. Obtain a regular grid
					// get number of indices in tile. idyi, col. idxi, row.
					uint idyi = (uint)(*stdp->innerLoopIndexVector)[innerLoop].size(); 
					uint idxi = (uint)outerLoopIndexVector.size(); 

					// get starting indices of tile. i, col. j, row
					int i = ((*stdp->innerLoopIndexVector)[innerLoop])[0];
					int j = outerLoopIndexVector[0];

					// get sub-tile of grid interpolation locations.
					dgrid subXInterpLocs0;
					dgrid subYInterpLocs0;
					(*stdp->xInterpLocs0).subgrid(subXInterpLocs0, i, j, idyi, idxi);
					(*stdp->yInterpLocs0).subgrid(subYInterpLocs0, i, j, idyi, idxi);

					//F. Create Delaunay Tri from scattered input data sub-tile.
					Bathy_Grid new_bathyGrid = Bathy_Grid();
					new_bathyGrid.Construct_Tin(&subX0, &subY0, &subZ_idy, &subH_idy, &subV_idy);
					
					//G. Estimate depths at regular grid points
					InterpGrid		new_gmt			= InterpGrid(GMT);					
					vector<double>	subZInterpLocs0 = vector<double> ((subXInterpLocs0.vec()).size(), 0.00);
					new_gmt.estimate(&subXInterpLocs0.vec(), &subYInterpLocs0.vec(), &subZInterpLocs0, 1.96, 2.00, *stdp->Lx, new_bathyGrid.getTin(), stdp->interpMethod, "", "");

					//H. Find slopes at grid points
					slopes = new_gmt.getGrads()->getSlopeOut();
					
					//I. Find dmin
					subXInterpLocs0_Vec = subXInterpLocs0.vec();
					subYInterpLocs0_Vec = subYInterpLocs0.vec();
					double minXi		= abs(subXInterpLocs0_Vec[1] - subXInterpLocs0_Vec[0]);
					double minYi		= abs(subYInterpLocs0_Vec[1] - subYInterpLocs0_Vec[0]);
					dmin				= 0.0;
					for(int j = 1; j < (const int)subXInterpLocs0_Vec.size()-1; j++)
					{
						diffXi = abs(subXInterpLocs0_Vec[j+1] - subXInterpLocs0_Vec[j]);
						diffYi = abs(subYInterpLocs0_Vec[j+1] - subYInterpLocs0_Vec[j]);
						if(!(minXi > 0) && diffXi > 0)
							minXi = diffXi;
						else if(diffXi < minXi && diffXi > 0)
							minXi = diffXi;

						if(!(minYi > 0) && diffYi > 0)
							minYi = diffYi;
						else if(diffYi < minYi && diffYi > 0)
							minYi = diffYi;

						dmin = min(minXi, minYi);
					}

					if(dmin == 0)
						cerr << "dmin equals 0";

					//E. Precompute the weights, riVector, and aiVector across the whole interpolation plane
					perturbWeights = vector<double>(idySize,2);
					scalecInterpPerturbations_PreCompute(&subZ_idy, &subE_idy, &perturbWeights, &perturb);
					#pragma endregion

					#pragma region Interpolate --_Compute_ForKriging (the same as _Compute)
					iliv_Loc0 = ((*stdp->innerLoopIndexVector)[innerLoop])[0];	//first innerLoopIndexVector location
					oliv_Loc0 = outerLoopIndexVector[0];						//first outerLoopIndexVector location
					tgs1_first = (*stdp->xMeshGrid)(iliv_Loc0, oliv_Loc0);		//get first xValue for xa+yb+c
					tgs2_first = (*stdp->yMeshGrid)(iliv_Loc0, oliv_Loc0);		//get first yValue

					iliv_Loc1 = ((*stdp->innerLoopIndexVector)[innerLoop])[((*stdp->innerLoopIndexVector)[innerLoop]).size() - 1]; //first innerLoopIndexVector loc
					oliv_Loc1 = outerLoopIndexVector[outerLoopIndexVector.size()-1];	//first outerLoopIndexVector loc
					tgs1_last = (*stdp->xMeshGrid)(iliv_Loc1, oliv_Loc1);				//get last xValue for xa+yb+c
					tgs2_last = (*stdp->yMeshGrid)(iliv_Loc1, oliv_Loc1);				//get last yValue

					//F. Get interpolation values at observation locations
					//We are doing input values here!
					for (int i = 0; i < (const int)idySize; i++)
					{
						Xiii(0,0) = subX0[i];				//doesn't matter cause we have a regular grid
						Xiii(0,1) = subX0[i];

						tgs1_Compute = subX_idy[i];			//current input xValue in tile
						tgs2_Compute = subY_idy[i];
						perturb.perturbationZKriged	= 0.0;
						perturb.perturbationEKriged	= 1.0;
						perturb.perturbationNEiKriged = 1.0;
						perturb.perturbationREiKriged = 1.0;
						if(stdp->PROP_UNCERT)
						{
							perturb.perturbationZ0Kriged = 0.0;
							perturb.perturbationE0Kriged = 1.0;
						}
						if(stdp->KALMAN)
						{
							perturb.perturbationZKKriged = 0.0;
							perturb.perturbationEKKriged = 1.0;
						}
						

						//We are doing the input data 
						//This is calls the original kriging function which does not have the new functionality and bug fixes.
						//scalecInterpPerturbations_Compute_ForKriging(&subX_idx, &subY_idx, &z_idx, &e_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, &perturb.riVector, &perturb.aiVector, (*sdp->neitol), &perturb.perturbationZKriged, &perturbationEKriged, &perturb.perturbationNEiKriged, &perturb.perturbationREiKriged, &standardDevKriged);

						//Current function
						scalecInterpPerturbations_Compute(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &subH_idy, &subV_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*stdp->neitol), dmin, (slopesVec2)[i], &subX0, &subY0, &Xiii, &perturb, stdp->MSE, stdp->PROP_UNCERT, stdp->KALMAN, KRIGING=1);
						
						//COMPARE THE DIFFERENCE BETWEEN SUBTILES: xMeshGrid AND subY_idyKriged AND 
						//Keep (un-scaled subtile input) kriged indices within MeshGrid (interpolation) subtile and remove kriged depth. 
						if ( (subY_idyKriged[i] > tgs2_first) && (subY_idyKriged[i] < tgs2_last)
						  && (subX_idyKriged[i] > tgs1_first) && (subX_idyKriged[i] < tgs1_last))
						{
							residualObservationsKrigedZ.push_back(subZ_idy[i] - perturb.perturbationZKriged);
							subX_indexKriged.push_back(subX_idyKriged[i]);
							subY_indexKriged.push_back(subY_idyKriged[i]);			

							//NEED TO CHECK TO SEE IF subX_indexKriged ARE THE SAME FOR ALL
							if(stdp->PROP_UNCERT)
								residualObservationsKrigedZ.push_back(subZ_idy[i] - perturb.perturbationZ0Kriged);
							if(stdp->KALMAN)
								residualObservationsKrigedZ.push_back(subZ_idy[i] - perturb.perturbationZKKriged);
						}
					} //idySize
					#pragma endregion Interpolate _Compute_ForKriging

					#pragma region --ordinaryKrigingOfResiduals_PreCompute, Subtile Kriged Indices if if too large 
					if (subX_indexKriged.size() >= 15)
					{
						int KRIG_SUBTILE_FLAG = 0;
						//G. We need to subtile the kriged indices otherwise the matrix inversion takes too long -- TODO
						if(subX_indexKriged.size() > KRIGED_SIZE_THRESHOLD)
						{
							KRIG_SUBTILE_FLAG = 1;
							#pragma region --Subtile Kriged Indices
							try
							{
								k = 0;
								xIndexKriged_Vector = vector<double>(outerLoopIndexVector.size()*((*stdp->innerLoopIndexVector)[innerLoop]).size());
								yIndexKriged_Vector = vector<double>(outerLoopIndexVector.size()*((*stdp->innerLoopIndexVector)[innerLoop]).size());
								for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
								{
									for (int j = 0; j < (const int)((*stdp->innerLoopIndexVector)[innerLoop]).size(); j++)
									{
										iliv_Loc = ((*stdp->innerLoopIndexVector)[innerLoop])[j];
										oliv_Loc = outerLoopIndexVector[i];
										xIndexKriged_Vector[k] = (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);	//kriged indices subtile (interpolation locations)
										yIndexKriged_Vector[k] = (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);
										k++;
									}
								}

								//cout << "To Tile Call" << endl;
								ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKrig, &outputErrorKrig);

								if(stdp->PROP_UNCERT)
									ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepth0Krig, &outputError0Krig);
								if(stdp->KALMAN)
									ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKKrig, &outputErrorKKrig);

								xIndexKriged_Vector.clear();
								yIndexKriged_Vector.clear();
							}
							catch (exception &e) {
								cout << "An exception occurred. Exception Thrown: " << e.what() << '\n';
								cout << "Kriging Subtiles failed.  Will try kriging without subtiles.  This will take longer." << endl;
								KRIG_SUBTILE_FLAG = 0;
							}
							#pragma endregion Subtile Kriged Indices
						}
						if(!KRIG_SUBTILE_FLAG)
						{
							#pragma region --Do Not Subtile Kriged Indices
							//small enough to do all at once
							ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
							
							if(stdp->PROP_UNCERT)
								ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
							if(stdp->KALMAN)
								ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);

							k = 0;
							//H. Krig the data
							for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
							{
								for (int j = 0; j < (const int)((*stdp->innerLoopIndexVector)[innerLoop]).size(); j++)
								{
									zKriged = 0.00;
									varZKriged = 0.00;
									iliv_Loc = ((*stdp->innerLoopIndexVector)[innerLoop])[j];
									oliv_Loc = outerLoopIndexVector[i];
									xGrid_indexKriged = (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);
									yGrid_indexKriged = (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);
									
									//Compute the residual value
									ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);
									
									outputDepthKrig[k] = zKriged;
									outputErrorKrig[k] = varZKriged;

									if(stdp->PROP_UNCERT)
									{
										z0Kriged = 0.00;
										varZ0Kriged = 0.00;

										ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &z0Kriged, &varZ0Kriged);
									
										outputDepth0Krig[k] = z0Kriged;
										outputError0Krig[k] = varZ0Kriged;
									}
									if(stdp->KALMAN)
									{
										zKKriged = 0.00;
										varZKKriged = 0.00;
										
										ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKKriged, &varZKKriged);
									
										outputDepthKKrig[k] = zKKriged;
										outputErrorKKrig[k] = varZKKriged;
									}

									k++;
								} // innerLoopIndexVector
							} // outerLoopInde
							twoGammaHatVector.clear();
							distanceVectorBinCenters.clear();
							aVectorFine.clear();
							invGammaDArray.clear();
							AGrid.clear();
							#pragma endregion Do Not Subtile Kriged Indices
						} 
					} 
					#pragma endregion ordinaryKrigingOfResiduals_PreCompute


					#pragma region Interpolate 
					k = 0;
					//K. Do the interpolation over each independent point
					for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
					{
						for (int j = 0; j < (const int)((*stdp->innerLoopIndexVector)[innerLoop]).size(); j++)
						{
							// current regular grid interpolation location
							Xiii(0,0) = subXInterpLocs0(j, i);							//doesn't matter cause we are using a regular grid
							Xiii(0,1) = subYInterpLocs0(j, i);

							iliv_Loc = ((*stdp->innerLoopIndexVector)[innerLoop])[j];	//innerLoopIndexVector location
							oliv_Loc = outerLoopIndexVector[i];							//outerLoopIndexVector location
							tgs1 = (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);				//get current xValue for xa+yb+c
							tgs2 = (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);				//get current yValue
							tgs1_Compute = tgs1 * (*stdp->Lx);							//scale xValue
							tgs2_Compute = tgs2 * (*stdp->Ly);							//scale yValue
							if(tgs1==Xiii(0,0))
								int stop=1;
							else 
								int stop =1;

							perturb.perturbationZ	= 0.0;
							perturb.perturbationE	= 1.0;
							perturb.perturbationNEi = 1.0;
							perturb.perturbationREi = 1.0;
							if(stdp->PROP_UNCERT)
							{
								perturb.perturbationZ0 = 0.0;
								perturb.perturbationE0 = 1.0;
							}
							if(stdp->KALMAN)
							{
								perturb.perturbationZK = 0.0;
								perturb.perturbationEK = 1.0;
							}

							//I. Compute the value
							scalecInterpPerturbations_Compute(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &subH_idy, &subV_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*stdp->neitol), dmin, (*slopes)(j,i), &subX0, &subY0, &Xiii, &perturb, stdp->MSE, stdp->PROP_UNCERT, stdp->KALMAN, KRIGING=0);

							//J. Put assn grid calculation here..... do matrix * vector math.......
							//put trend back into this tile
							// 1*c + x*a + y*b
							assnGridValue = tgs0*(*stdp->btrend)[0] + tgs1*(*stdp->btrend)[1] + tgs2*(*stdp->btrend)[2];
							(*stdp->outputDepth)(iliv_Loc,oliv_Loc) = (perturb.perturbationZ + assnGridValue) + outputDepthKrig[k];
							(*stdp->outputError)(iliv_Loc,oliv_Loc) = perturb.perturbationE + outputErrorKrig[k];
							(*stdp->outputNEi)(iliv_Loc,oliv_Loc)	= perturb.perturbationNEi;
							(*stdp->outputREi)(iliv_Loc,oliv_Loc)	= perturb.perturbationREi;
							(*stdp->standardDev)(iliv_Loc,oliv_Loc) = perturb.standardDev;
							if(stdp->PROP_UNCERT)
							{
								(*stdp->outputDepth0)(iliv_Loc,oliv_Loc) = (perturb.perturbationZ0 + assnGridValue) + outputDepth0Krig[k];
								(*stdp->outputError0)(iliv_Loc,oliv_Loc) = perturb.perturbationE0 + outputError0Krig[k];
		//						(*stdp->standardDev0)(iliv_Loc,oliv_Loc) = perturb.standardDev0;
							}
							if(stdp->KALMAN)
							{
								(*stdp->outputDepthK)(iliv_Loc,oliv_Loc) = (perturb.perturbationZK + assnGridValue) + outputDepthKKrig[k];
								(*stdp->outputErrorK)(iliv_Loc,oliv_Loc) = perturb.perturbationEK + outputErrorKKrig[k];
							}

							k++;
						} // innerLoopIndexVector
					} // outerLoopIndexVector
					#pragma endregion Interpolate
					//required clears!
					perturbWeights.clear();
					perturb.riVector.clear();
					perturb.aiVector.clear();
					new_bathyGrid.clear();
					subXInterpLocs0.clear();
					subYInterpLocs0.clear();
					slopes = NULL;
					slopesVec = NULL;
//					gradgrid.clear();

					residualObservationsKrigedZ.clear();
					subX_indexKriged.clear();
					subY_indexKriged.clear();
					if(stdp->PROP_UNCERT)
						residualObservationsKrigedZ0.clear();
					if(stdp->KALMAN)
						residualObservationsKrigedZK.clear();
				}//idySize
				#pragma endregion idySize
				//required clears!
				subX_idy.clear();
				subY_idy.clear();
				subZ_idy.clear();
				subE_idy.clear();
				subH_idy.clear();
				subV_idy.clear();
				subX0.clear();
				subY0.clear();
				subX_idyKriged.clear();
				subY_idyKriged.clear();
			}//innerloop
		}//idxSize
		outerLoopIndexVector.clear();
		idx.clear();
		#pragma endregion idxSize
	}//outerloop
	outputDepthKrig.clear();
	outputErrorKrig.clear();
	if(stdp->PROP_UNCERT)
	{
		outputDepth0Krig.clear();
		outputError0Krig.clear();
	}
	if(stdp->KALMAN)
	{
		outputDepthKKrig.clear();
		outputErrorKKrig.clear();
	}
	return SUCCESS;
}





////This function is from the _Serial version and calls the old scalecInterp_Perturbations_Compute_ForKriging function.  This does not include any of the new changes and improvements. It is the function from the mergeBathy_v3.7.1_Paul version that was updated to pass standardDev in the Perturbations call.  This function calls the original Perturbations_Compute_ForKriging (which does not have any of the new functionality or bug fixes) when kriging and then calls the new Perturbations_Compute when interpolating.  To use, comment this one is and comment the other one out.
//int scalecInterpTile_ProcessKrig(SCALEC_TILE_DATA_POINTER stdp, const int curIterNum, const int numCores)
//{
//	//************************************************************************************
//	// 0. Declare local variables and objects that are going to be used inside the loop
//	//************************************************************************************
//	vector<double> outputDepthKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
//	vector<double> outputErrorKrig = vector<double>((*stdp->xMeshGrid).rows()*(*stdp->xMeshGrid).cols(),0);
//	PERTURBS perturb;
//	perturb.kernelName=*(*stdp).kernelName;
//	dgrid* slopes;							// interpolation location grid slopes
//	vector<double> subX0;					// scattered input data for Tri
//	vector<double> subY0;
//	vector<double> subXInterpLocs0_Vec;		// interpolation locations
//	vector<double> subYInterpLocs0_Vec;
//	dgrid Xiii(1,2);						// current interpolation location
//	double diffXi, diffYi;
//	double dmin;
//
//	/*double perturbationZ;
//	double perturbationE;
//	double perturbationNEi;
//	double perturbationREi;
//	double standardDev;
//	*/
//	double perturbationZKriged;
//	double perturbationEKriged;
//	double perturbationNEiKriged;
//	double perturbationREiKriged;
//	//double thetaM, phi_2;
//	double zKriged;
//	double varZKriged;
//	double xGrid_indexKriged;
//	double yGrid_indexKriged;
//	double standardDevKriged;
//
//	vector<int> outerLoopIndexVector;
//	vector<int> idx;
//	vector<double> subX_idy;
//	vector<double> subY_idy;
//	vector<double> subZ_idy;
//	vector<double> subE_idy;
//	vector<double> subX_idyKriged;
//	vector<double> subY_idyKriged;
//	vector<double> residualObservationsKrigedZ;
//	vector<double> subX_indexKriged;
//	vector<double> subY_indexKriged;
//	vector<double> xIndexKriged_Vector;
//	vector<double> yIndexKriged_Vector;
//	vector<double> subH_idy;
//	vector<double> subV_idy;
//
//	vector<double> twoGammaHatVector;
//	vector<double> distanceVectorBinCenters;
//	vector<double> aVectorFine;
//	dgrid invGammaDArray;
//	dgrid AGrid;
//
//
//	double xmin, xmax, ymin, ymax;
////	int perturbSize;
//
//	vector <double> perturbWeights;
//	vector <double> riVector;
//	vector <double> aiVector;
//
//	const int subDataXLength = (*stdp->subsampledData)[0].size();
//	int iliv_Loc;
//	int oliv_Loc;
//	int idySize;
//	int i, j, k;
//	double tgs0 = 1.0;
//	double tgs1, tgs2, tgs1_Compute, tgs2_Compute;
//	double assnGridValue;
//	double locSpacingX = (*stdp->spacingX);
//	double locSpacingY = (*stdp->spacingY);
//
//	//************************************************************************************
//	// I. Interpolate the data
//	//************************************************************************************
//	for (int outerLoop = curIterNum; outerLoop < (const int)(*stdp->kx); outerLoop+=numCores)
//	{
//		outerLoopIndexVector = vector<int>((*stdp->nkx));
//		//A. Get indices in XI, which gives x-locations, to interpolate on
//		//this loop.
//		for (i = 0; i < (const int)(*stdp->nkx); i++)
//		{
//			//indices to interpolate this time
//			outerLoopIndexVector[i] = ( i + ((outerLoop) * (int)(*stdp->nkx))); 
//		}
//		if((outerLoop == (int)((*stdp->kx) - 1)) && (outerLoopIndexVector[outerLoopIndexVector.size() - 1] != (int)((*stdp->xSingleVector).size() - 1)))
//		{
//			for (i = outerLoopIndexVector[outerLoopIndexVector.size() - 1] + 1; i < (const int)(*stdp->xSingleVector).size(); i++)
//			{
//				//catch the end here
//				outerLoopIndexVector.push_back(i); 
//			}
//		}
//		//B. Compute appropriate overlap between tiles horizontally and get the
//		//useful data. Recall that we scaled xi-array by 1./std(x). Then,
//		//get indices of x-coordinates in x-array for the data to be
//		//interpolated. Indices to be filtered further for useful y-coordinates.
//		//find tile limits
//		xmin = (*stdp->xSingleVector)[(int)outerLoopIndexVector[0]] - (*stdp->LMAX_x);
//		xmax = (*stdp->xSingleVector)[(int)outerLoopIndexVector[outerLoopIndexVector.size()-1]] + (*stdp->LMAX_x);
//
//		for (i = 0; i < subDataXLength; i++){
//			if (((*stdp->subsampledData)[0][i] < xmax) && ((*stdp->subsampledData)[0][i] > xmin)){
//				idx.push_back(i);
//			}
//		}
//		
//		#pragma region --idxSize
//		//loop y tiles, blocks in columns
//		if (idx.size() > 0)
//		{
//			for (int innerLoop = 0; innerLoop < (const int)(*stdp->ky); innerLoop++)
//			{
//				//C. Do the same as above but for the y coordinates
//				ymin = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[0]] - (*stdp->LMAX_y);
//				ymax = (*stdp->ySingleVector)[(int)((*stdp->innerLoopIndexVector)[innerLoop])[((*stdp->innerLoopIndexVector)[innerLoop]).size()-1]] + (*stdp->LMAX_y);
//
//				//D. Get all of the orignal data that falls within the current tile to be interpolated
//				for (i = 0; i < (const int)idx.size(); i++){
//					if ( ((*stdp->subsampledData)[1][(int)idx[i]] < ymax) && ((*stdp->subsampledData)[1][(int)idx[i]] > ymin) ){
//						subX_idy.push_back((*stdp->subsampledData)[0][idx[i]] * (*stdp->Lx));
//						subY_idy.push_back((*stdp->subsampledData)[1][idx[i]] * (*stdp->Ly));
//						subX_idyKriged.push_back((*stdp->subsampledData)[0][idx[i]]);
//						subY_idyKriged.push_back((*stdp->subsampledData)[1][idx[i]]);
//						subZ_idy.push_back((*stdp->subsampledData)[2][idx[i]]);
//						subE_idy.push_back((*stdp->subsampledData)[4][idx[i]]);
//
//						subH_idy.push_back((*stdp->subsampledData)[5][idx[i]]);
//						subV_idy.push_back((*stdp->subsampledData)[6][idx[i]]);
//						subX0.push_back((*stdp->x0)[idx[i]]);	//scattered input data sub-tile
//						subY0.push_back((*stdp->y0)[idx[i]]);
//
//					}
//				}
//				idySize = subX_idy.size();
//
//				#pragma region idySize
//				if (idySize > 2)
//				{
//					//E. Obtain a regular grid
//					// get number of indices in tile. idyi, col. idxi, row.
//					uint idyi = (uint)(*stdp->innerLoopIndexVector)[innerLoop].size(); 
//					uint idxi = (uint)outerLoopIndexVector.size(); 
//
//					// get starting indices of tile. i, col. j, row
//					i = ((*stdp->innerLoopIndexVector)[innerLoop])[0];
//					j = outerLoopIndexVector[0];
//
//					// get sub-tile of grid interpolation locations.
//					dgrid subXInterpLocs0;
//					dgrid subYInterpLocs0;
//					(*stdp->xInterpLocs0).subgrid(subXInterpLocs0, i, j, idyi, idxi);
//					(*stdp->yInterpLocs0).subgrid(subYInterpLocs0, i, j, idyi, idxi);
//
//					//F. Create Delaunay Tri from scattered input data sub-tile.
//					Bathy_Grid new_bathyGrid = Bathy_Grid();
//					new_bathyGrid.Construct_Tin(&subX0, &subY0, &subZ_idy, &subH_idy, &subV_idy);
//					
//					//G. Estimate depths at regular grid points
//					InterpGrid		new_gmt			= InterpGrid(GMT);					
//					vector<double>	subZInterpLocs0 = vector<double> ((subXInterpLocs0.vec()).size(), 0.00);
//					new_gmt.estimate(&subXInterpLocs0.vec(), &subYInterpLocs0.vec(), &subZInterpLocs0, 1.96, 2.00, *stdp->Lx, new_bathyGrid.getTin(), stdp->interpMethod);
//
//					//H. Find slopes at grid points
//					slopes = new_gmt.getGrads()->getSlopeOut();
//
//					//I. Find dmin
//					subXInterpLocs0_Vec = subXInterpLocs0.vec();
//					subYInterpLocs0_Vec = subYInterpLocs0.vec();
//					double minXi		= abs(subXInterpLocs0_Vec[1] - subXInterpLocs0_Vec[0]);
//					double minYi		= abs(subYInterpLocs0_Vec[1] - subYInterpLocs0_Vec[0]);
//					dmin				= 0.0;
//					for(int j = 1; j < (const int)subXInterpLocs0_Vec.size()-1; j++)
//					{
//						diffXi = abs(subXInterpLocs0_Vec[j+1] - subXInterpLocs0_Vec[j]);
//						diffYi = abs(subYInterpLocs0_Vec[j+1] - subYInterpLocs0_Vec[j]);
//						if(!(minXi > 0) && diffXi > 0)
//							minXi = diffXi;
//						else if(diffXi < minXi && diffXi > 0)
//							minXi = diffXi;
//
//						if(!(minYi > 0) && diffYi > 0)
//							minYi = diffYi;
//						else if(diffYi < minYi && diffYi > 0)
//							minYi = diffYi;
//
//						dmin = min(minXi, minYi);
//					}
//
//					if(dmin == 0)
//						cerr << "dmin equals 0";
//
//					//E. Precompute the weights, riVector, and aiVector across the whole interpolation plane
//					perturbWeights = vector<double>(idySize,2);
//					scalecInterpPerturbations_PreCompute(&subZ_idy, &subE_idy, &perturbWeights, &perturb);
//
//					#pragma region Interpolate 
//					////K. Do the interpolation over each independent point
//					//for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
//					//{
//					//	for (int j = 0; j < (const int)(*stdp->innerLoopIndexVector)[innerLoop].size(); j++)
//					//	{
//					//		// current interpolation location
//					//		Xiii(0,0) = subXInterpLocs0(j, i);
//					//		Xiii(0,1) = subYInterpLocs0(j, i);
//					//		perturb.perturbationZKriged = 0;
//					//		perturb.perturbationEKriged = 1;
//					//		perturb.perturbationNEiKriged = 1;
//					//		perturb.perturbationREiKriged = 1;
//					//		if(stdp->PROP_UNCERT)
//					//		{
//					//			perturb.perturbationZ0Kriged = 0.0;
//					//			perturb.perturbationE0Kriged = 1.0;
//					//		}
//					//		if(stdp->KALMAN)
//					//		{
//					//			perturb.perturbationZKKriged = 0.0;
//					//			perturb.perturbationEKKriged = 1.0;
//					//		}
//
//					//F. Get interpolation values at observation locations
//					for (i = 0; i < (const int)idySize; i++)
//					{
//						perturbationZKriged = 0;
//						perturbationEKriged = 1;
//						perturbationNEiKriged = 1;
//						perturbationREiKriged = 1;
//						tgs1_Compute = subX_idy[i];
//						tgs2_Compute = subY_idy[i];
//
//						scalecInterpPerturbations_Compute_ForKriging(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, &riVector, &aiVector, (*stdp->neitol), &perturbationZKriged, &perturbationEKriged, &perturbationNEiKriged, &perturbationREiKriged, &standardDevKriged);
//						//scalecInterpPerturbations_Compute_ForKriging(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, &riVector, &aiVector, (*stdp->neitol), dmin, (*slopes)(j,i), &subX0, &subY0, &Xiii, &perturb, stdp->MSE, stdp->PROP_UNCERT, stdp->KALMAN);
//
//						if ( (subY_idyKriged[i] > (*stdp->yMeshGrid)( ((*stdp->innerLoopIndexVector)[innerLoop])[0], outerLoopIndexVector[0]) )
//							&& (subY_idyKriged[i] <
//								(*stdp->yMeshGrid)( ((*stdp->innerLoopIndexVector)[innerLoop])[((*stdp->innerLoopIndexVector)[innerLoop]).size() - 1] , outerLoopIndexVector[outerLoopIndexVector.size()-1]))
//							&& (subX_idyKriged[i] >
//								(*stdp->xMeshGrid)(((*stdp->innerLoopIndexVector)[innerLoop])[0],outerLoopIndexVector[0]))
//							&& (subX_idyKriged[i] <
//								(*stdp->xMeshGrid)(((*stdp->innerLoopIndexVector)[innerLoop])[((*stdp->innerLoopIndexVector)[innerLoop]).size()-1],outerLoopIndexVector[outerLoopIndexVector.size()-1])) )
//						{
//							residualObservationsKrigedZ.push_back(subZ_idy[i] - perturb.perturbationZKriged);
//							subX_indexKriged.push_back(subX_idyKriged[i]);
//							subY_indexKriged.push_back(subY_idyKriged[i]);
//						}
//					}
//					//	}
//					//}
//					#pragma endregion Interpolate
//
//
//					if (subX_indexKriged.size() >= 15)
//					{
//						//G. We need to subtile the kriged indexes otherwise the matrix inversion takes too long -- TODO
//						if(subX_indexKriged.size() > KRIGED_SIZE_THRESHOLD)
//						{
//							k = 0;
//							xIndexKriged_Vector = vector<double>(outerLoopIndexVector.size()*((*stdp->innerLoopIndexVector)[innerLoop]).size());
//							yIndexKriged_Vector = vector<double>(outerLoopIndexVector.size()*((*stdp->innerLoopIndexVector)[innerLoop]).size());
//							for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
//							{
//								for (int j = 0; j < (const int)((*stdp->innerLoopIndexVector)[innerLoop]).size(); j++)
//								{
//									iliv_Loc = ((*stdp->innerLoopIndexVector)[innerLoop])[j];
//									oliv_Loc = outerLoopIndexVector[i];
//									xIndexKriged_Vector[k] = (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);
//									yIndexKriged_Vector[k] = (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);
//									k++;
//
//								}
//							}
//
//							//cout << "To Tile Call" << endl;
//							ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKrig, &outputErrorKrig);
//
//							xIndexKriged_Vector.clear();
//							yIndexKriged_Vector.clear();
//
//						}else
//						{//small enough to do all at once
//							ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
//
//							//H. Krig the data
//							k = 0;
//							for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
//							{
//								for (int j = 0; j < (const int)((*stdp->innerLoopIndexVector)[innerLoop]).size(); j++)
//								{
//									zKriged = 0.00;
//									varZKriged = 0.00;
//									iliv_Loc = ((*stdp->innerLoopIndexVector)[innerLoop])[j];
//									oliv_Loc = outerLoopIndexVector[i];
//									xGrid_indexKriged = (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);
//									yGrid_indexKriged = (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);
//									//Compute the residual value
//									ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);
//									outputDepthKrig[k] = zKriged;
//									outputErrorKrig[k] = varZKriged;
//									k++;
//								}
//							}
//							twoGammaHatVector.clear();
//							distanceVectorBinCenters.clear();
//							aVectorFine.clear();
//							invGammaDArray.clear();
//							AGrid.clear();
//						}
//					}
//
//					k = 0;
//					for (int i = 0; i < (const int)outerLoopIndexVector.size(); i++)
//					{
//						for (int j = 0; j < (const int)((*stdp->innerLoopIndexVector)[innerLoop]).size(); j++)
//						{
//							// current interpolation location
//							Xiii(0,0) = subXInterpLocs0(j, i);
//							Xiii(0,1) = subYInterpLocs0(j, i);
//
//							iliv_Loc = ((*stdp->innerLoopIndexVector)[innerLoop])[j];
//							oliv_Loc = outerLoopIndexVector[i];
//							tgs1 = (*stdp->xMeshGrid)(iliv_Loc,oliv_Loc);
//							tgs2 = (*stdp->yMeshGrid)(iliv_Loc,oliv_Loc);
//							tgs1_Compute = tgs1 * (*stdp->Lx);
//							tgs2_Compute = tgs2 * (*stdp->Ly);
//
//							perturb.perturbationZ = 0;
//							perturb.perturbationE = 1;
//							perturb.perturbationREi = 1;
//							if(stdp->PROP_UNCERT)
//							{
//								perturb.perturbationZ0 = 0.0;
//								perturb.perturbationE0 = 1.0;
//							}
//							if(stdp->KALMAN)
//							{
//								perturb.perturbationZK = 0.0;
//								perturb.perturbationEK = 1.0;
//							}
//							/*perturbationZ = 0;
//							perturbationE = 1;
//							perturbationREi = 1;*/
//
//							//I. Compute the value
//							//scalecInterpPerturbations_Compute(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, &riVector, &aiVector, (*stdp->neitol), &perturbationZ, &perturbationE, &perturbationNEi, &perturbationREi, &standardDev);
//							scalecInterpPerturbations_Compute(&subX_idy, &subY_idy, &subZ_idy, &subE_idy, &subH_idy, &subV_idy, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*stdp->neitol), dmin, (*slopes)(j,i), &subX0, &subY0, &Xiii, &perturb, stdp->MSE, stdp->PROP_UNCERT, stdp->KALMAN, bool KRIGING=0);
//
//							//J. Put assn grid calculation here..... do matrix * vector math.......
//							//put trend back into this tile
//							assnGridValue = tgs0*(*stdp->btrend)[0]+tgs1*(*stdp->btrend)[1]+tgs2*(*stdp->btrend)[2];
//							/*(*stdp->outputDepth)(iliv_Loc,oliv_Loc) =  (perturbationZ + assnGridValue) + outputDepthKrig[k];
//							(*stdp->outputError)(iliv_Loc,oliv_Loc) = perturbationE + outputErrorKrig[k];
//							(*stdp->outputNEi)(iliv_Loc,oliv_Loc) = perturbationNEi;
//							(*stdp->outputREi)(iliv_Loc,oliv_Loc) = perturbationREi;
//							(*stdp->standardDev)(iliv_Loc,oliv_Loc) = standardDev;
//							*/
//							(*stdp->outputDepth)(iliv_Loc,oliv_Loc) =  (perturb.perturbationZ + assnGridValue) + outputDepthKrig[k];
//							(*stdp->outputError)(iliv_Loc,oliv_Loc) = perturb.perturbationE + outputErrorKrig[k];
//							(*stdp->outputNEi)(iliv_Loc,oliv_Loc) = perturb.perturbationNEi;
//							(*stdp->outputREi)(iliv_Loc,oliv_Loc) = perturb.perturbationREi;
//							(*stdp->standardDev)(iliv_Loc,oliv_Loc) = perturb.standardDev;
//							if(stdp->PROP_UNCERT)
//							{
//								(*stdp->outputDepth0)(iliv_Loc,oliv_Loc) = perturb.perturbationZ0 + assnGridValue;
//								(*stdp->outputError0)(iliv_Loc,oliv_Loc) = perturb.perturbationE0;
//	//							(*stdp->standardDev0)(iliv_Loc,oliv_Loc) = perturb.standardDev20;
//							}
//							if(stdp->KALMAN)
//							{
//								(*stdp->outputDepthK)(iliv_Loc,oliv_Loc) = perturb.perturbationZK + assnGridValue;
//								(*stdp->outputErrorK)(iliv_Loc,oliv_Loc) = perturb.perturbationEK;
//							}
//							k++;
//						}
//					}
//					perturbWeights.clear();
//					perturb.riVector.clear();
//					perturb.aiVector.clear();
//					new_bathyGrid.clear();
//					subXInterpLocs0.clear();
//					subYInterpLocs0.clear();
//					slopes = NULL;
//
//					residualObservationsKrigedZ.clear();
//					subX_indexKriged.clear();
//					subY_indexKriged.clear();
//					//End the new stuff...
//				}//idySize
//				#pragma endregion idySize
//				subX_idy.clear();
//				subY_idy.clear();
//				subZ_idy.clear();
//				subE_idy.clear();
//				subX_idyKriged.clear();
//				subY_idyKriged.clear();
//				subH_idy.clear();
//				subV_idy.clear();
//				subX0.clear();
//				subY0.clear();
//			}//innerloop
//		}//idxSize
//		#pragma endregion idxSize
//		outerLoopIndexVector.clear();
//		idx.clear();
//	}
//	outputDepthKrig.clear();
//	outputErrorKrig.clear();
//
//	return SUCCESS;
//}
