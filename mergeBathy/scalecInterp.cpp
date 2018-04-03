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
#include "scalecInterp.h"
#include "constants.h"
#include "regr_xzw.h"
#include "kriging.h"
#include <fstream>
#include <time.h>
#include "MB_Threads.h"
#include <algorithm>

int scalecInterp(vector< vector<double> > *subsampledData, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double meanXSingle, double meanYSingle, string &kernelName, map<string, int> additionalOptions, const double neitol, OUTPUT_DATA *xyzOut)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int i;
	size_t subDataXLength	= (*subsampledData)[0].size();
	size_t outputDepthSize	= (*xInterpVector).size();
	double meanX_std		= 0.0000;
	double meanY_std		= 0.0000;
	double meanErrorSquared = 0.0000;

	double newSpacingX, newSpacingY;
	double Lx, Ly;
	double numberOfComputedTiles;
	double kx, ky;
	double numStepsX, numStepsY;
	double LMAX_x, LMAX_y;
	double startT, stopT;
	double minInterpX, maxInterpX, minInterpY, maxInterpY;
	
	//Interpolation Locations
	vector<double> xInterpVectorLocal((*xInterpVector));	//xi, Xii0
	vector<double> yInterpVectorLocal((*yInterpVector));

	//Preserve the original spatial relationships prior to interpolation scaling
	//Keep original input and output points in UTM
	//Later, a Delaunay Tri is built with x0, y0 and queried based on xInterpLocs0, yInterpLocs0
	vector<double> x0				((*subsampledData)[0]);	//x
	vector<double> y0				((*subsampledData)[1]);	
	vector<double> xInterpLocs0		(xInterpVectorLocal);	//xi
	vector<double> yInterpLocs0		(yInterpVectorLocal);	

	//************************************************************************************
	//I. Remove overall trend from locations.
	//************************************************************************************
	//A. Remove trends in position array. 
	//1. First, shift data and interpolation locations to center of
	//   grid by removing the mean.
	//	 Center on output's center so subtract output's mean which we
	//	 already calculated for xInterpVector.
	(*subsampledData)[0]	-= meanXSingle; //data
	(*subsampledData)[1]	-= meanYSingle;
	xInterpVectorLocal		-= meanXSingle;	//interpolation grid
	yInterpVectorLocal		-= meanYSingle;

	printf("Interpolating over known input grid\n");
	printf(".");

	//2. Scale grid by 1/std(x), where x(:, 1:2) = location of data points in meters, so
	//	 that the data point locations are not huge numbers - basically, -3 to 3).
	//i. Compute mean that is needed for the variance estimate.
	for (i = 0; i < (const int)subDataXLength; i++){
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
	for (i = 0; i < (const int)subDataXLength; i++){
		std_x = std_x + pow(((*subsampledData)[0][i] - meanX_std),2);
		std_y = std_y + pow(((*subsampledData)[1][i] - meanY_std),2);
	}
	std_x = std_x / (double)(subDataXLength - 1.00);
	std_y = std_y / (double)(subDataXLength - 1.00);
	std_x = sqrt(std_x);
	std_y = sqrt(std_y);

	printf(".");

	//iii. Scale the data and grid.
	(*subsampledData)[0] /= std_x;	
	(*subsampledData)[1] /= std_y;	
	xInterpVectorLocal /= std_x;	
	yInterpVectorLocal /= std_y;	

	newSpacingX = gridSpacingX / std_x;
	newSpacingY = gridSpacingY / std_y;

	Lx = 1.0/newSpacingX;
	Ly = 1.0/newSpacingY;

	//************************************************************************************
	//II. Remove trend surface from the observations. It will be added back later.
	//************************************************************************************
	//A. Compute consistent window weights to pass to the linear regressions routine: regr_xzw.
	printf(".");

	double wtol = 0.0100;
	double constWeightsS3 = 0;
	vector<double> weights(subDataXLength, 2.00);
	consistentWeights(&(*subsampledData)[2], &(*subsampledData)[4], &wtol, &weights, &constWeightsS3);

	//B. Call regr_xzw.m to calculate a 2-D linear fit to the data set.
	//Initialize output variables from regr_xzw.m
	dvector btrend	(3,0.00);
	dvector bi		(3,0.00);

	//C. Call regr_xzw if variables have a variance; else, keep the padding with zeros
	vector<double> regrX(subDataXLength,1);
	printf(".");
	if(!(std_x <= 0 && std_y <= 0)){
		// do regression to remove a norm field
		// this is just for getting the data ready, so it is meant to be bullet proof, not statistically pure!
		regr_xzw(&regrX, &(*subsampledData)[0], &(*subsampledData)[1], &(*subsampledData)[2], &(*subsampledData)[3], &weights, &btrend, &bi);
	}
	regrX.clear();

	//D. Removes trend from the data (i.e. remove overall bias from the data).
	double zTrendValue;
	for (i = 0; i < (const int)subDataXLength; i++){
		zTrendValue = (1.0*btrend[0]) + ((*subsampledData)[0][i]*btrend[1]) + ((*subsampledData)[1][i]*btrend[2]);
		// compute deviations from trend
		(*subsampledData)[2][i] = (*subsampledData)[2][i] - zTrendValue;
		//Multiply the sub-sampled data by Lx and Ly for use in interpPerturbations
	}

	startT = clock();

	//************************************************************************************
	//III. Now get to the processing by pre-computing the data
	//************************************************************************************
	//A. Now call scalecInterpPerturbations_PreCompute
	(*xyzOut).depth			= vector<double>(outputDepthSize, 0);
	(*xyzOut).error			= vector<double>(outputDepthSize, 0);
	(*xyzOut).nEi			= vector<double>(outputDepthSize, 1.00);
	(*xyzOut).rEi			= vector<double>(outputDepthSize, 0);
	(*xyzOut).standardDev	= vector<double>(outputDepthSize, 0);
	(*xyzOut).depth0		= vector<double>(outputDepthSize, 0);
	(*xyzOut).error0		= vector<double>(outputDepthSize, 0);
//	(*xyzOut).standardDev0	= vector<double>(outputDepthSize, 0);
	(*xyzOut).depthK		= vector<double>(outputDepthSize, 0);
	(*xyzOut).errorK		= vector<double>(outputDepthSize, 0);

	//B. Compute the minimum and maximum interpolation values
	minInterpX = xInterpVectorLocal[0];
	minInterpY = yInterpVectorLocal[0];
	maxInterpX = xInterpVectorLocal[0];
	maxInterpY = yInterpVectorLocal[0];
	for (i = 0; i < (const int)xInterpVectorLocal.size(); i++)
	{
		if (xInterpVectorLocal[i] < minInterpX)
			minInterpX = xInterpVectorLocal[i];
		else if (xInterpVectorLocal[i] > maxInterpX)
			maxInterpX = xInterpVectorLocal[i];
		if (yInterpVectorLocal[i] < minInterpY)
			minInterpY = yInterpVectorLocal[i];
		else if (yInterpVectorLocal[i] > maxInterpY)
			maxInterpY = yInterpVectorLocal[i];
	}
	maxInterpX = maxInterpX*1.01;
	maxInterpY = maxInterpY*1.01;

	//C. Compute the number of tiles to break the irregular data into
	numberOfComputedTiles = sqrt( ((double)xInterpVectorLocal.size()) * (1.00 + abs(newSpacingX / (maxInterpX - minInterpX))) );

	if (numberOfComputedTiles < sqrt( ((double)yInterpVectorLocal.size()) * (1.00 + abs(newSpacingY / (maxInterpY - minInterpY))) ))
		numberOfComputedTiles = sqrt( ((double)yInterpVectorLocal.size()) * (1.00 + abs(newSpacingY / (maxInterpY - minInterpY))) );

	//D. THIS IS THE TRICKY STUFF. Compute the spacing and tile numbers.
	kx = floor( sqrt( (double) numberOfComputedTiles ));	//ceil
	ky = floor( sqrt( (double) numberOfComputedTiles ));
	numStepsX = abs(maxInterpX - minInterpX) / (double)kx;
	numStepsY = abs(maxInterpY - minInterpY) / (double)ky;
	LMAX_x = 10.00*(newSpacingX);	// specify max overlap between tiles,
	LMAX_y = 10.00*(newSpacingY);	// often = 10*length scale   //numberOfComputedTiles
	#pragma endregion

	//E. Create Delaunay Tri
	vector<double> h0	(x0.size(),0.00);
	vector<double> v0	(x0.size(),0.00);
	Bathy_Grid new_bathyGrid;
	new_bathyGrid.Construct_Tin(&x0, &y0, &(*subsampledData)[2], &h0, &v0);
	
	vector<double> x_idx;
	vector<double> y_idx;
	vector<double> z_idx;
	vector<double> e_idx;
	vector<double> h_idx;
	vector<double> v_idx;
	vector<double> x0_idx;
	vector<double> y0_idx;
	vector<double> x_idxKriged;
	vector<double> y_idxKriged;
	vector<double> slopeOut	(xInterpLocs0.size(), 0.00);

	//Get sub-grid to interpolate over which will be the entire grid
	for (int ii = 0; ii < (const int)subsampledData->begin()->size(); ii++)
	{
		x_idx.push_back((*subsampledData)[0][ii] * Lx);
		y_idx.push_back((*subsampledData)[1][ii] * Ly);
		z_idx.push_back((*subsampledData)[2][ii]);
		e_idx.push_back((*subsampledData)[4][ii]);
		h_idx.push_back((*subsampledData)[5][ii]);
		v_idx.push_back((*subsampledData)[6][ii]);
		x0_idx.push_back(x0[ii]);
		y0_idx.push_back(y0[ii]);

		/*x_idxKriged.push_back((*subsampledData)[0][ii]);
		y_idxKriged.push_back((*subsampledData)[1][ii]);*/
		x_idxKriged.push_back((*subsampledData)[0][ii] * Lx); //SJZ
		y_idxKriged.push_back((*subsampledData)[1][ii] * Ly);
	}

	//Get the estimators to perform
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

	//F. Initialize the data structure
	//Start with the new calling scheme here
	SCALEC_DATA scalecInterpData;
	scalecInterpData.subsampledData = subsampledData;
	scalecInterpData.xInterpVector	= &xInterpVectorLocal;
	scalecInterpData.yInterpVector	= &yInterpVectorLocal;
	scalecInterpData.btrend			= &btrend;
	scalecInterpData.kernelName		= &kernelName;
	scalecInterpData.Lx				= &Lx;
	scalecInterpData.Ly				= &Ly;
	scalecInterpData.minSubX		= &minInterpX;
	scalecInterpData.minSubY		= &minInterpY;
	scalecInterpData.kx				= &kx;
	scalecInterpData.ky				= &ky;
	scalecInterpData.numStepsX		= &numStepsX;
	scalecInterpData.numStepsY		= &numStepsY;
	scalecInterpData.LMAX_x			= &LMAX_x;
	scalecInterpData.LMAX_y			= &LMAX_y;
	scalecInterpData.spacingX		= &newSpacingX;
	scalecInterpData.spacingY		= &newSpacingY;
	scalecInterpData.neitol			= &neitol;
	scalecInterpData.outData		= xyzOut;
	scalecInterpData.new_bathyGrid	= &new_bathyGrid;
	scalecInterpData.x_idx			= &x_idx;
	scalecInterpData.y_idx			= &y_idx;
	scalecInterpData.z_idx			= &z_idx;
	scalecInterpData.e_idx			= &e_idx;
	scalecInterpData.h_idx			= &h_idx;
	scalecInterpData.v_idx			= &v_idx;
	scalecInterpData.x0_idx			= &x0_idx;
	scalecInterpData.y0_idx			= &y0_idx;
	scalecInterpData.x_idxKriged	= &x_idxKriged;
	scalecInterpData.y_idxKriged	= &y_idxKriged;


	scalecInterpData.xInterpLocs0	= &xInterpLocs0;
	scalecInterpData.yInterpLocs0	= &yInterpLocs0;
	scalecInterpData.slopeOut		= &slopeOut;
	scalecInterpData.interpMethod	= interpMethod;
	scalecInterpData.MSE			= MSE;
	scalecInterpData.PROP_UNCERT	= PROP_UNCERT;
	scalecInterpData.KALMAN			= KALMAN;
	scalecInterpData.KRIGING		= KRIGING;
	scalecInterpData.residualObservationsKrigedZ		= new vector<double>(x_idxKriged.size(),0.00);
	scalecInterpData.residualObservationsKrigedZ0		= new vector<double>(x_idxKriged.size(),0.00);
	scalecInterpData.residualObservationsKrigedZK		= new vector<double>(x_idxKriged.size(),0.00);
	
	
	#pragma endregion

	//************************************************************************************
	//IV. Check for multi-threading support and run the processing routines
	//************************************************************************************
	if (additionalOptions["-multiThread"] == 0)
	{
		if (additionalOptions["-kriging"] == 0)
		{
			scalecInterp_Process(&scalecInterpData, 0,1);
		//	scalecInterp_Process2A(&scalecInterpData, 0,1);//handles w and w/o kriging 
			scalecInterp_Process4A(&scalecInterpData, 0,1);//handles w and w/o kriging 
			scalecInterp_Process5A(&scalecInterpData, 0,1);//handles w and w/o kriging 
			scalecInterp_Process6A(&scalecInterpData, 0,1);//handles w and w/o kriging 
			//scalecInterp_Process2(&scalecInterpData, 0,1);
		}else
		{
			//scalecInterp_ProcessKrig(&scalecInterpData, 0,1);
			scalecInterp_Process(&scalecInterpData, 0,1);
			scalecInterp_Process4A(&scalecInterpData, 0,1);//handles w and w/o kriging 
			scalecInterp_Process5A(&scalecInterpData, 0,1);//handles w and w/o kriging 
			scalecInterp_Process6A(&scalecInterpData, 0,1);//handles w and w/o kriging 
		}
	}else
	{
		mbThreads mbT = mbThreads(additionalOptions["-multiThread"]);
		mbT.makeMBThread(&scalecInterpData);
		mbT.initMBThread(0); 
//		mbT.initMBThread(additionalOptions["-kriging"]);
		mbT.joinMBThread();

		mbT.initMBThread2(0);
		//mbT.initMBThread2(additionalOptions["-kriging"]);
		mbT.joinMBThread();
		//mbT.terminateMBThread();

		scalecInterp_Process5A(&scalecInterpData, 0,1);//handles w and w/o kriging 

		mbT.initMBThread6(0);
		//mbT.initMBThread2(additionalOptions["-kriging"]);
		mbT.joinMBThread();
		mbT.terminateMBThread();
	}

	stopT = clock();
	cout << "\nTime to Complete Interpolation: " << (stopT-startT)/CLOCKS_PER_SEC << " seconds" << endl << endl;

	//Clear everything up
	weights.clear();
	regrX.clear();
	btrend.clear();
	bi.clear();

	xInterpVectorLocal.clear();
	yInterpVectorLocal.clear();
	xInterpLocs0.clear();
	yInterpLocs0.clear();
	x_idx.clear();
	y_idx.clear();
	e_idx.clear();
	h_idx.clear();
	v_idx.clear();
	x0_idx.clear();
	y0_idx.clear();
	slopeOut.clear();
	h0.clear();
	v0.clear();

	return 0;
}

//scalecInterp_Process Part I -  Function was broken up to allow for multi-threading
int scalecInterp_Process(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores)
{
	//Initialize some variables
	int centpnt;
	int ndipts;

	dgrid* slopes;
	int xSize, ySize;
	double del1, del2, delr,sumCount, nxtxi, xii0_upper, yii0_upper, xii0_lower, yii0_lower;
	vector<double> rdel, rdelTemp;
	vector<double>::iterator it;

	double tgs0 = 1.00;
	vector<double> xt;
	vector<double> yt;
	vector<double> xMeshVector;
	vector<double> yMeshVector;
	dgrid xx;
	dgrid yy;
	dgrid newxx;
	dgrid newyy;
	vector<double> zz;

	size_t Ni = sdp->xInterpLocs0->size();

	//**********************************************************************
	//I. Non gridded output assumed; Initiate Function with full dataset WEA
	//**********************************************************************
	for(int i = curIterNum; i < (const int)Ni; i += numCores)
	{
		//A. Obtain a Grid
		// #points about expansion point; should be even
		//I believe Will did this just to get the slope
		ndipts = 4;
		nxtxi  = MAX_INT;
		for (int j = 0; j < (const int)Ni; j++)
		{
			sumCount = 0;
			del1 = (*sdp->xInterpLocs0)[j] - (*sdp->xInterpLocs0)[i]; //vector difference
			del2 = (*sdp->yInterpLocs0)[j] - (*sdp->yInterpLocs0)[i]; //vector difference
			sumCount += pow(del1,2);
			sumCount += pow(del2,2);
			rdel.push_back(sqrt(sumCount));
			if(rdel[j] != 0 && rdel[j] < nxtxi)
				nxtxi = rdel[j];
		}

		// make delta for points in a temp sub-grid (rdelTemp) around
		// the point of interest
		delr = nxtxi*ndipts;

		// find the points in the temp sub-grid (rdelTemp)
//		auto ior = remove_if(rdel.begin(), rdel.end(), [=](double k){return (k < delr);});
//		rdel.resize(distance(rdel.begin(), ior));

		for(it = rdel.begin(); it < rdel.end(); it++)
		{
			if(*it >= delr)//!= NaN) //UNIX c98
				rdelTemp.push_back(*it);
		}

		// make the domain a tad bigger than the data extent
		xii0_upper = ((*sdp->xInterpLocs0)[i] + delr);
		xii0_lower = ((*sdp->xInterpLocs0)[i] - delr);
		yii0_upper = ((*sdp->yInterpLocs0)[i] + delr);
		yii0_lower = ((*sdp->yInterpLocs0)[i] - delr);
		// Calculate x and y for our new grid where we go from
		// xii0_lower to xii0_upper at a step size of nxtxi.
		createMeshXYDims(xii0_lower, yii0_lower, xii0_upper, yii0_upper, nxtxi, nxtxi, &xt, &yt);

		// get our new x,y dimensions
		xSize = (const int)xt.size();
		ySize = (const int)yt.size();
		if (0)
		{
			cout << "Dimensions of Computational Area "<< i << ": " << endl;
			cout << "\tRows: " << ySize << "\n\tCols: " << xSize << endl << endl;
		}

		// resize vectors to match our new dimensions
		xMeshVector.resize(xSize*ySize, 0.00);
		yMeshVector.resize(xSize*ySize, 0.00);
		xx.resize(ySize, xSize, 0.00);
		yy.resize(ySize, xSize, 0.00);
		// now create grid with those dimensions
		createMeshGrid(&xt, &yt, &xMeshVector, &yMeshVector, &xx, &yy);

		newxx.resize(xx.rows(), xx.cols(),0.00);
		newyy.resize(yy.rows(), yy.cols(),0.00);

		// reshape grid
		reshapeGrid(&xx, &newxx);
		reshapeGrid(&yy, &newyy);
		centpnt = (int)(fix(ndipts*2/2+1)-1);

		//B. Estimate zz value
		InterpGrid* new_gmt = new InterpGrid(GMT);
		zz.resize((newxx.vec()).size(), 0.00);
		(*new_gmt).estimate(&(newxx.vec()), &(newyy.vec()), &zz, 1.96, 2.00, *sdp->Lx, sdp->new_bathyGrid->getTin(), sdp->interpMethod, "", "");

		//C. Find slope at grid points
		slopes = (*new_gmt).getGrads()->getSlopeOut();
		(*((*sdp)).slopeOut)[i] = (*slopes)(centpnt, centpnt);//*&

		// clean up 
		delete new_gmt;
		zz.clear();
		xx.clear();
		yy.clear();
		xMeshVector.clear();
		yMeshVector.clear();
		xt.clear();
		yt.clear();
		newxx.clear();
		newyy.clear();
		rdel.clear();
		rdelTemp.clear();
	} 

	return SUCCESS;
}

//scalecInterp_Process Part II -  Function was broken up to allow for multi-threading
int scalecInterp_Process2(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores)
{
	double dmin = NaN;	//Non-gridded output assumed

	double tgs0 = 1.0, tgs1, tgs2;
	double tgs1_Compute;
	double tgs2_Compute;
	double assnGridValue;
	
	double	Ni = (double)(*sdp->xInterpLocs0).size();
	dgrid Xiii(1,2);	// current interpolation location

	PERTURBS perturb;
	perturb.kernelName = *(*sdp).kernelName;

	//3. Pre-compute the weights, riVector, and aiVector across the whole interpolation plane
	vector <double> perturbWeights(sdp->z_idx->size(),2);
	scalecInterpPerturbations_PreCompute(sdp->z_idx, sdp->e_idx, &perturbWeights, &perturb);

	for (int i = curIterNum; i < Ni; i += numCores)
	{
		tgs1			= (*sdp->xInterpVector)[i];		//get current xValue
		tgs2			= (*sdp->yInterpVector)[i];		//get current yValue
		tgs1_Compute	= tgs1 * (*sdp->Lx);			//scale xValue
		tgs2_Compute	= tgs2 * (*sdp->Ly);			//scale yValue
		Xiii(0,0)		= (*sdp->xInterpLocs0)[i];
		Xiii(0,1)		= (*sdp->yInterpLocs0)[i];

		//Initialize output fields
		perturb.perturbationZ = 0.0;
		perturb.perturbationE = 1.0;
		perturb.perturbationNEi = 1.0;
		perturb.perturbationREi = 1.0;
		if(sdp->PROP_UNCERT)
		{
			perturb.perturbationZ0 = 0.0;
			perturb.perturbationE0 = 1.0;
		}
		if(sdp->KALMAN)
		{
			perturb.perturbationZK = 0.0;
			perturb.perturbationEK = 1.0;
		}
		bool KRIGING = 0;
		//4. Compute the value
		scalecInterpPerturbations_Compute(sdp->x_idx, sdp->y_idx, sdp->z_idx, sdp->e_idx, sdp->h_idx, sdp->v_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*sdp->neitol), dmin, (*(sdp->slopeOut))[i], sdp->x0_idx, sdp->y0_idx, &Xiii, &perturb, sdp->MSE, sdp->PROP_UNCERT, sdp->KALMAN, KRIGING);

		//5. Put assn grid calculation here..... do matrix * vector math.......
		//	put trend back into this tile
		assnGridValue = tgs0*(*sdp->btrend)[0]+tgs1*(*sdp->btrend)[1]+tgs2*(*sdp->btrend)[2];
		(*sdp->outData).depth[i] = perturb.perturbationZ + assnGridValue;
		(*sdp->outData).error[i] = perturb.perturbationE;
		(*sdp->outData).nEi[i] = perturb.perturbationNEi;
		(*sdp->outData).rEi[i] = perturb.perturbationREi;
		(*sdp->outData).standardDev[i] = perturb.standardDev2;
		if(sdp->PROP_UNCERT)
		{
			(*sdp->outData).depth0[i] = perturb.perturbationZ0 + assnGridValue;
			(*sdp->outData).error0[i] = perturb.perturbationE0;
		//	(*sdp->outData).standardDev0[i] = perturb.standardDev20;
		}
		if(sdp->KALMAN)
		{
			(*sdp->outData).depthK[i] = perturb.perturbationZK + assnGridValue;
			(*sdp->outData).errorK[i] = perturb.perturbationEK;
		}
	}
	perturbWeights.clear();
	perturb.riVector.clear();
	perturb.aiVector.clear();

	return SUCCESS;
}


//This function is for both kriging and non-kriging runs. It is 2A broken down for proper threading into 4A-6A.
//Do not run 2 or 2A; replaced by 4A-6A.
//scalecInterp_Process Part II -  Function was broken up to allow for multi-threading
int scalecInterp_Process4A(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores)
{
	double dmin = NaN;	//Non-gridded output assumed

	double tgs0 = 1.0;// , tgs1, tgs2;
	double tgs1_Compute;
	double tgs2_Compute;
	//double assnGridValue;
	
	double	Ni = (double)(*sdp->xInterpLocs0).size();
	dgrid Xiii(1,2);	// current interpolation location

	PERTURBS perturb;
	perturb.kernelName = *(*sdp).kernelName;

	//3. Pre-compute the weights, riVector, and aiVector across the whole interpolation plane
	vector <double> perturbWeights(sdp->z_idx->size(),2);
	scalecInterpPerturbations_PreCompute(sdp->z_idx, sdp->e_idx, &perturbWeights, &perturb);

	if(sdp->KRIGING)
	{
		//************************************************************************************
		//I. Begin Residual Kriging.
		//************************************************************************************
		//for (int i = curIterNum; i < Ni; i += numCores)
		for (int i = curIterNum; i < (const int)(*sdp->x_idx).size(); i += numCores)
		{
			#pragma region --Find the Residuals at Observations
			//************************************************************************************
			//4. Get interpolation values at observation locations; _Compute_ForKriging
			//************************************************************************************
			Xiii(0,0)		= (*sdp->x0_idx)[i]; //query pts will be the input locations of which we know the value.
			Xiii(0,1)		= (*sdp->y0_idx)[i];
		
			tgs1_Compute = (*sdp->x_idx)[i];	//current input x value in tile shouldn't matter cause we have an irregular grid
			tgs2_Compute = (*sdp->y_idx)[i];
			perturb.perturbationZ	= 0.0;
			perturb.perturbationE	= 1.0;
			perturb.perturbationNEi = 1.0;
			perturb.perturbationREi = 1.0;
			if(sdp->PROP_UNCERT)
			{
				perturb.perturbationZ0 = 0.0;
				perturb.perturbationE0 = 1.0;
			}
			if(sdp->KALMAN)
			{
				perturb.perturbationZK = 0.0;
				perturb.perturbationEK = 1.0;
			}
			vector<double> slopesVec2 = vector<double>((*sdp->subsampledData)[0].size(), 1.00);

			//This was the __Compute_ForKriging call but was changed to use the current _Compute function.
			//A. Estimate the depth and errors at the observation location.
			scalecInterpPerturbations_Compute(sdp->x_idx, sdp->y_idx, sdp->z_idx, sdp->e_idx, sdp->h_idx, sdp->v_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*sdp->neitol), dmin, (slopesVec2)[i], sdp->x0_idx, sdp->y0_idx, &Xiii, &perturb, sdp->MSE, sdp->PROP_UNCERT, sdp->KALMAN, sdp->KRIGING);
		
			//B. Find the residual from the depth estimation and the actual value at the observation
			(*sdp->residualObservationsKrigedZ)[i] = ((*sdp->z_idx)[i] - perturb.perturbationZ);
			if(sdp->PROP_UNCERT)
				(*sdp->residualObservationsKrigedZ0)[i] = ((*sdp->z_idx)[i] - perturb.perturbationZ0);
			if(sdp->KALMAN)
				(*sdp->residualObservationsKrigedZK)[i] = ((*sdp->z_idx)[i] - perturb.perturbationZK);

			//Add the observation location to a tile in order to subtile to reduce complexity.
			//subX_indexKriged.push_back((*sdp->x_idxKriged)[i]); 
			//subY_indexKriged.push_back((*sdp->y_idxKriged)[i]);
		
			//end Get interpolation values at observation locations
			#pragma endregion
		}
	}
	perturbWeights.clear();
	perturb.riVector.clear();
	perturb.aiVector.clear();

	return SUCCESS;
}

//scalecInterp_Process Part II -  This function should only be passed numCores=1
int scalecInterp_Process5A(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores)
{
	PERTURBS perturb;
	perturb.kernelName = *(*sdp).kernelName;

	//Initialize kriging structures
	double locSpacingX = (*sdp->spacingX);
	double locSpacingY = (*sdp->spacingY);
	double zKriged;
	double varZKriged;
	double z0Kriged;
	double varZ0Kriged;
	double zKKriged;
	double varZKKriged;
	double xGrid_indexKriged;
	double yGrid_indexKriged;
	
	vector<double> xIndexKriged_Vector;
	vector<double> yIndexKriged_Vector;
	
	vector<double> twoGammaHatVector;
	vector<double> distanceVectorBinCenters;
	vector<double> aVectorFine;
	dgrid invGammaDArray;
	dgrid AGrid;

	vector<double> outputDepthKrig;
	vector<double> outputErrorKrig;
	vector<double> outputDepth0Krig;
	vector<double> outputError0Krig;
	vector<double> outputDepthKKrig;
	vector<double> outputErrorKKrig;

	//Initialize locations of which to obtain interpolation values
	xIndexKriged_Vector = vector<double>((*sdp->xInterpVector));
	yIndexKriged_Vector= vector<double>((*sdp->yInterpVector));
	
	//Initialize Krige Outputs
	outputDepthKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputErrorKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputDepth0Krig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputError0Krig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputDepthKKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputErrorKKrig = vector<double>((*sdp->xInterpVector).size(), 0);

	//3. Pre-compute the weights, riVector, and aiVector across the whole interpolation plane
	vector <double> perturbWeights(sdp->z_idx->size(),2);
	scalecInterpPerturbations_PreCompute(sdp->z_idx, sdp->e_idx, &perturbWeights, &perturb);

	if(sdp->KRIGING)
	{
		//************************************************************************************
		//I. Begin Residual Kriging.
		//************************************************************************************
		#pragma region --Krige the Residuals
		//When tile is sufficiently large enough, krige it.
		if (sdp->x_idxKriged->size() >= 15)
		{
			int KRIG_SUBTILE_FLAG = 0;
			//5. We need to subtile the kriged indexes; otherwise the matrix inversion takes too long -- TODO
			if(sdp->x_idxKriged->size() > KRIGED_SIZE_THRESHOLD)
			{
				KRIG_SUBTILE_FLAG = 1;
				#pragma region --Subtile Kriged Indices
				try
				{
					//Too large, need to subtile.
					//6. Krige the residuals
					//_PreComputeTile calls _PostCompute from within
					ordinaryKrigingOfResiduals_PreComputeTile(&*sdp->x_idxKriged, &*sdp->y_idxKriged, sdp->residualObservationsKrigedZ, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKrig, &outputErrorKrig);
					sdp->outData->depth = outputDepthKrig;
					sdp->outData->error = outputErrorKrig;

					if (sdp->PROP_UNCERT)
					{
						ordinaryKrigingOfResiduals_PreComputeTile(&*sdp->x_idxKriged, &*sdp->y_idxKriged, sdp->residualObservationsKrigedZ0, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepth0Krig, &outputError0Krig);
						sdp->outData->depth0 = outputDepth0Krig;
						sdp->outData->error0 = outputError0Krig;
					}
					if (sdp->KALMAN)
					{
						ordinaryKrigingOfResiduals_PreComputeTile(&*sdp->x_idxKriged, &*sdp->y_idxKriged, sdp->residualObservationsKrigedZK, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKKrig, &outputErrorKKrig);
						sdp->outData->depthK = outputDepthKKrig;
						sdp->outData->errorK = outputErrorKKrig;
					}
				}
				catch (exception &e) {
					cout << "An exception occurred. Exception Thrown: " << e.what() << '\n';
					cout << "Kriging Subtiles failed.  Will try kriging without subtiles.  This will take longer." << endl;
					KRIG_SUBTILE_FLAG = 0;
				}
				#pragma endregion
			}
			if(!KRIG_SUBTILE_FLAG)
			{
				#pragma region --Do Not Subtile Kriged Indices
				//small enough to do all at once
				ordinaryKrigingOfResiduals_PreCompute(sdp->x_idxKriged, sdp->y_idxKriged, sdp->residualObservationsKrigedZ, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
				if(sdp->PROP_UNCERT)
					ordinaryKrigingOfResiduals_PreCompute(sdp->x_idxKriged, sdp->y_idxKriged, sdp->residualObservationsKrigedZ0, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
				if(sdp->KALMAN)
					ordinaryKrigingOfResiduals_PreCompute(sdp->x_idxKriged, sdp->y_idxKriged, sdp->residualObservationsKrigedZK, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);

				//6. Krige the residuals
				for (int i = 0; i < (const int)(*sdp->xInterpVector).size(); i++)
				{
					xGrid_indexKriged = xIndexKriged_Vector[i]; //current output grid loc 
					yGrid_indexKriged = yIndexKriged_Vector[i];

					zKriged = 0.00;
					varZKriged = 0.00;
					
					ordinaryKrigingOfResiduals_PostCompute(sdp->x_idxKriged, sdp->y_idxKriged, sdp->residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);
					
					sdp->outData->depth[i] = zKriged;
					sdp->outData->error[i] = varZKriged;
	
					if(sdp->PROP_UNCERT)
					{
						z0Kriged = 0.00;
						varZ0Kriged = 0.00;
						ordinaryKrigingOfResiduals_PostCompute(sdp->x_idxKriged, sdp->y_idxKriged, sdp->residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &z0Kriged, &varZ0Kriged);
						
						sdp->outData->depth0[i] = z0Kriged;
						sdp->outData->error0[i] = varZ0Kriged;
					}
					if(sdp->KALMAN)
					{
						zKKriged = 0.00;
						varZKKriged = 0.00;
						ordinaryKrigingOfResiduals_PostCompute(sdp->x_idxKriged, sdp->y_idxKriged, sdp->residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKKriged, &varZKKriged);
						
						sdp->outData->depthK[i] = zKKriged;
						sdp->outData->errorK[i] = varZKKriged;
					}
				}
				twoGammaHatVector.clear();
				distanceVectorBinCenters.clear();
				aVectorFine.clear();
				invGammaDArray.clear();
				AGrid.clear();
				#pragma endregion
			}
		}
		//Done Kriging Residuals
		#pragma endregion 
	}

	perturbWeights.clear();
	perturb.riVector.clear();
	perturb.aiVector.clear();
	(*sdp->residualObservationsKrigedZ).clear();
	if(sdp->PROP_UNCERT)
	{
		(*sdp->residualObservationsKrigedZ0).clear();
	}	
	if(sdp->KALMAN)
	{
		(*sdp->residualObservationsKrigedZK).clear();
	}
	delete sdp->residualObservationsKrigedZ;
	delete sdp->residualObservationsKrigedZ0;
	delete sdp->residualObservationsKrigedZK;

	return SUCCESS;
}

int scalecInterp_Process6A(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores)
{
	double dmin = NaN;	//Non-gridded output assumed

	double tgs0 = 1.0, tgs1, tgs2;
	double tgs1_Compute;
	double tgs2_Compute;
	double assnGridValue;
	
	double	Ni = (double)(*sdp->xInterpLocs0).size();
	dgrid Xiii(1,2);	// current interpolation location

	PERTURBS perturb;
	perturb.kernelName = *(*sdp).kernelName;

	vector<double> xIndexKriged_Vector;
	vector<double> yIndexKriged_Vector;

	//Initialize locations of which to obtain interpolation values
	//xGrid_indexKriged = (*sdp->xInterpVector)[idxInterp[i]];
	//indexKriged_VectorLocs.push_back(idxInterp[i]);
	xIndexKriged_Vector = vector<double>((*sdp->xInterpVector));
	yIndexKriged_Vector= vector<double>((*sdp->yInterpVector));
	
	//3. Pre-compute the weights, riVector, and aiVector across the whole interpolation plane
	vector <double> perturbWeights(sdp->z_idx->size(),2);
	scalecInterpPerturbations_PreCompute(sdp->z_idx, sdp->e_idx, &perturbWeights, &perturb);

	for (int i = curIterNum; i < Ni; i += numCores)
	{
		//tgs1 = xIndexKriged_Vector[i];			//shouldn't use this because we have an irregular grid? same as xInterpsLocs0 and xInterpVector
		//tgs2 = yIndexKriged_Vector[i];			//get current yValue

		tgs1			= (*sdp->xInterpVector)[i];		//get current xValue
		tgs2			= (*sdp->yInterpVector)[i];		//get current yValue
		tgs1_Compute	= tgs1 * (*sdp->Lx);			//scale xValue
		tgs2_Compute	= tgs2 * (*sdp->Ly);			//scale yValue
		Xiii(0,0)		= (*sdp->xInterpLocs0)[i];
		Xiii(0,1)		= (*sdp->yInterpLocs0)[i];

		//Initialize output fields
		perturb.perturbationZ = 0.0;
		perturb.perturbationE = 1.0;
		perturb.perturbationNEi = 1.0;
		perturb.perturbationREi = 1.0;
		if(sdp->PROP_UNCERT)
		{
			perturb.perturbationZ0 = 0.0;
			perturb.perturbationE0 = 1.0;
		}
		if(sdp->KALMAN)
		{
			perturb.perturbationZK = 0.0;
			perturb.perturbationEK = 1.0;
		}
		
		//4. Compute the value
		scalecInterpPerturbations_Compute(sdp->x_idx, sdp->y_idx, sdp->z_idx, sdp->e_idx, sdp->h_idx, sdp->v_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*sdp->neitol), dmin, (*(sdp->slopeOut))[i], sdp->x0_idx, sdp->y0_idx, &Xiii, &perturb, sdp->MSE, sdp->PROP_UNCERT, sdp->KALMAN, 0);//0 because were not kriging the residuals here... that's done above

		//5. Put assn grid calculation here..... do matrix * vector math.......
		//	put trend back into this tile
		assnGridValue = tgs0*(*sdp->btrend)[0]+tgs1*(*sdp->btrend)[1]+tgs2*(*sdp->btrend)[2];		
		(*sdp->outData).depth[i] += perturb.perturbationZ + assnGridValue;// + outDepthKrig;
		(*sdp->outData).error[i] += perturb.perturbationE;// + outErrorKrig;
		(*sdp->outData).nEi[i] = perturb.perturbationNEi;
		(*sdp->outData).rEi[i] = perturb.perturbationREi;
		(*sdp->outData).standardDev[i] = perturb.standardDev2;
		if(sdp->PROP_UNCERT)
		{
			(*sdp->outData).depth0[i] += perturb.perturbationZ0 + assnGridValue;// + outDepth0Krig;
			(*sdp->outData).error0[i] += perturb.perturbationE0;// + outError0Krig;
		//	(*sdp->outData).standardDev0[i] = perturb.standardDev20;
		}
		if(sdp->KALMAN)
		{
			(*sdp->outData).depthK[i] += perturb.perturbationZK + assnGridValue;// + outDepthKKrig;
			(*sdp->outData).errorK[i] += perturb.perturbationEK;// + outErrorKKrig;
		}
	}
	perturbWeights.clear();
	perturb.riVector.clear();
	perturb.aiVector.clear();

	return SUCCESS;
}













//This is Paul's completely changed to match current scalecInterp calculations and use current scalecInterp_Perturbations
int scalecInterp_ProcessKrig(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores)
{
	double dmin = (double)NaN;	//Non-gridded output assumed
	double	Ni = (double)(*sdp->xInterpLocs0).size();
	dgrid Xiii(1,2);	// current interpolation location

	PERTURBS perturb;
	perturb.kernelName = *(*sdp).kernelName;

	//Initialize some variables
	double tgs0, tgs1, tgs2;
	double tgs1_Compute;
	double tgs2_Compute;
	double assnGridValue;

	vector<double> pw_idx;

	vector <double> perturbWeights;
	double zKriged;
	double varZKriged;
	double z0Kriged;
	double varZ0Kriged;
	double zKKriged;
	double varZKKriged;
	double xGrid_indexKriged;
	double yGrid_indexKriged;
	
	vector<double> residualObservationsKrigedZ;
	vector<double> residualObservationsKrigedZ0;
	vector<double> residualObservationsKrigedZK;
	vector<double> subX_indexKriged; //observation locations
	vector<double> subY_indexKriged;
	//vector<double> indexKriged_VectorLocs;
	vector<double> xIndexKriged_Vector;
	vector<double> yIndexKriged_Vector;

	vector<double> twoGammaHatVector;
	vector<double> distanceVectorBinCenters;
	vector<double> aVectorFine;
	dgrid invGammaDArray;
	dgrid AGrid;

	vector<double> outputDepthKrig;
	vector<double> outputErrorKrig;
	vector<double> outputDepth0Krig;
	vector<double> outputError0Krig;
	vector<double> outputDepthKKrig;
	vector<double> outputErrorKKrig;

	double locSpacingX = (*sdp->spacingX);
	double locSpacingY = (*sdp->spacingY);

	tgs0 = 1.00;


	//************************************************************************************
	//3. Precompute the weights, riVector, and aiVector across the whole interpolation plane
	//************************************************************************************
	perturbWeights = vector<double>(sdp->z_idx->size(),2);
	scalecInterpPerturbations_PreCompute(sdp->z_idx, sdp->e_idx, &perturbWeights, &perturb);

	//Initialize locations of which to obtain interpolation values
	//xGrid_indexKriged = (*sdp->xInterpVector)[idxInterp[i]];
	//indexKriged_VectorLocs.push_back(idxInterp[i]);
	xIndexKriged_Vector = vector<double>((*sdp->xInterpVector));
	yIndexKriged_Vector= vector<double>((*sdp->yInterpVector));
	//xIndexKriged_Vector = vector<double>((*sdp->xInterpLocs0));//I think it should be this one cause we're irregular
	//yIndexKriged_Vector= vector<double>((*sdp->yInterpLocs0));
	
	//Initialize Krige Outputs
	outputDepthKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputErrorKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputDepth0Krig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputError0Krig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputDepthKKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputErrorKKrig = vector<double>((*sdp->xInterpVector).size(), 0);

	//************************************************************************************
	//I. Begin Residual Kriging.
	//************************************************************************************
	//for (int i = curIterNum; i < Ni; i += numCores)
	for (int i = curIterNum; i < (const int)(*sdp->x_idx).size(); i += numCores)
	{
		#pragma region --Find the Residuals at Observations
		//************************************************************************************
		//4. Get interpolation values at observation locations; _Compute_ForKriging
		//************************************************************************************
		Xiii(0,0)		= (*sdp->x0_idx)[i]; //query pts will be the input locations of which we know the value.
		Xiii(0,1)		= (*sdp->y0_idx)[i];
		
		tgs1_Compute = (*sdp->x_idx)[i];	//current input x value in tile shouldn't matter cause we have an irregular grid
		tgs2_Compute = (*sdp->y_idx)[i];
		perturb.perturbationZ	= 0.0;
		perturb.perturbationE	= 1.0;
		perturb.perturbationNEi = 1.0;
		perturb.perturbationREi = 1.0;
		if(sdp->PROP_UNCERT)
		{
			perturb.perturbationZ0 = 0.0;
			perturb.perturbationE0 = 1.0;
		}
		if(sdp->KALMAN)
		{
			perturb.perturbationZK = 0.0;
			perturb.perturbationEK = 1.0;
		}

		bool KRIGING = 0;
		//This was the __Compute_ForKriging call but was changed to use the current _Compute function.
		//A. Estimate the depth and errors at the observation location.
		scalecInterpPerturbations_Compute(sdp->x_idx, sdp->y_idx, sdp->z_idx, sdp->e_idx, sdp->h_idx, sdp->v_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*sdp->neitol), dmin, (*(sdp->slopeOut))[i], sdp->x0_idx, sdp->y0_idx, &Xiii, &perturb, sdp->MSE, sdp->PROP_UNCERT, sdp->KALMAN, KRIGING);
		
		//B. Find the residual from the depth estimation and the actual value at the observation
		residualObservationsKrigedZ.push_back((*sdp->z_idx)[i] - perturb.perturbationZKriged);
		if(sdp->PROP_UNCERT)
			residualObservationsKrigedZ0.push_back((*sdp->z_idx)[i] - perturb.perturbationZ0Kriged);
		if(sdp->KALMAN)
			residualObservationsKrigedZK.push_back((*sdp->z_idx)[i] - perturb.perturbationZKKriged);

		//Add the observation location to a tile in order to subtile to reduce complexity.
		subX_indexKriged.push_back((*sdp->x_idxKriged)[i]); //I don't think we should use because we are irregular?
		subY_indexKriged.push_back((*sdp->y_idxKriged)[i]);
		//subX_indexKriged.push_back((*sdp->x0_idx)[i]); //don't know if this is right
		//subY_indexKriged.push_back((*sdp->y0_idx)[i]);
		//end Get interpolation values at observation locations
		#pragma endregion

		#pragma region --Krige the Residuals
		//When tile is sufficiently large enough, krige it.
		if (subX_indexKriged.size() >= 15)
		{
			int KRIG_SUBTILE_FLAG = 0;
			//5. We need to subtile the kriged indexes; otherwise the matrix inversion takes too long -- TODO
			if(subX_indexKriged.size() > KRIGED_SIZE_THRESHOLD)
			{
				KRIG_SUBTILE_FLAG = 1;
				#pragma region --Subtile Kriged Indices
				try
				{
					//Too large, need to subtile.
					//6. Krige the residuals
					//_PreComputeTile calls _PostCompute from within
					ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKrig, &outputErrorKrig);

					if(sdp->PROP_UNCERT)
						ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepth0Krig, &outputError0Krig);
					if(sdp->KALMAN)
						ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKKrig, &outputErrorKKrig);
				}
				catch (exception &e) {
					cout << "An exception occurred. Exception Thrown: " << e.what() << '\n';
					cout << "Kriging Subtiles failed.  Will try kriging without subtiles.  This will take longer." << endl;
					KRIG_SUBTILE_FLAG = 0;
				}
				#pragma endregion
			}
			if(!KRIG_SUBTILE_FLAG)
			{
				#pragma region --Do Not Subtile Kriged Indices
				//small enough to do all at once
				ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
				if(sdp->PROP_UNCERT)
					ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
				if(sdp->KALMAN)
					ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);

				//6. Krige the residuals
				for (int i = 0; i < (const int)(*sdp->xInterpVector).size(); i++)
				{
					xGrid_indexKriged = xIndexKriged_Vector[i];
					yGrid_indexKriged = yIndexKriged_Vector[i];

					zKriged = 0.00;
					varZKriged = 0.00;
					
					ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);
					outputDepthKrig[i] = zKriged;
					outputErrorKrig[i] = varZKriged;
	
					if(sdp->PROP_UNCERT)
					{
						z0Kriged = 0.00;
						varZ0Kriged = 0.00;
						ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &z0Kriged, &varZ0Kriged);
						outputDepth0Krig[i] = z0Kriged;
						outputError0Krig[i] = varZ0Kriged;
					}
					if(sdp->KALMAN)
					{
						zKKriged = 0.00;
						varZKKriged = 0.00;
						ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKKriged, &varZKKriged);
						outputDepthKKrig[i] = zKKriged;
						outputErrorKKrig[i] = varZKKriged;
					}
				}
				twoGammaHatVector.clear();
				distanceVectorBinCenters.clear();
				aVectorFine.clear();
				invGammaDArray.clear();
				AGrid.clear();
				#pragma endregion
			}
		}
		//Done Kriging Residuals
		#pragma endregion 

		//Continue normal processing
		for (int i = 0; i < (const int)(*sdp->xInterpVector).size(); i++)
		{
			tgs1 = xIndexKriged_Vector[i];			//shouldn't use this because we have an irregular grid?
			tgs2 = yIndexKriged_Vector[i];			//get current yValue
			tgs1_Compute = tgs1 * (*sdp->Lx);		//scale xValue
			tgs2_Compute = tgs2 * (*sdp->Ly);		//scale yValue

			Xiii(0,0)		= (*sdp->xInterpLocs0)[i];		//I believe this is the same as tgs1
			Xiii(0,1)		= (*sdp->yInterpLocs0)[i];		//Output grid interpolation locations

			//Initialize output fields
			perturb.perturbationZ = 0.0;
			perturb.perturbationE = 1.0;
			perturb.perturbationNEi = 1.0;
			perturb.perturbationREi = 1.0;
			if(sdp->PROP_UNCERT)
			{
				perturb.perturbationZ0 = 0.0;
				perturb.perturbationE0 = 1.0;
			}
			if(sdp->KALMAN)
			{
				perturb.perturbationZK = 0.0;
				perturb.perturbationEK = 1.0;
			}

			bool KRIGING = 0;
			//7. Compute the value
			scalecInterpPerturbations_Compute(sdp->x_idx, sdp->y_idx, sdp->z_idx, sdp->e_idx, sdp->h_idx, sdp->v_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*sdp->neitol), dmin, (*(sdp->slopeOut))[i], sdp->x0_idx, sdp->y0_idx, &Xiii, &perturb, sdp->MSE, sdp->PROP_UNCERT, sdp->KALMAN, KRIGING);

			//8. Put assn grid calculation here..... do matrix * vector math.......
			//put trend back into this tile
			assnGridValue = tgs0*(*sdp->btrend)[0]+tgs1*(*sdp->btrend)[1]+tgs2*(*sdp->btrend)[2];

			(*sdp->outData).depth[i] = perturb.perturbationZ + assnGridValue + outputDepthKrig[i];
			(*sdp->outData).error[i] = perturb.perturbationE + outputErrorKrig[i];
			(*sdp->outData).nEi[i] = perturb.perturbationNEi;
			(*sdp->outData).rEi[i] = perturb.perturbationREi;
			(*sdp->outData).standardDev[i] = perturb.standardDev2;
			if(sdp->PROP_UNCERT)
			{
				(*sdp->outData).depth0[i] = perturb.perturbationZ0 + assnGridValue + outputDepth0Krig[i];
				(*sdp->outData).error0[i] = perturb.perturbationE0 + outputError0Krig[i];
			//	(*sdp->outData).standardDev0[i] = perturb.standardDev20;
			}
			if(sdp->KALMAN)
			{
				(*sdp->outData).depthK[i] = perturb.perturbationZK + assnGridValue + outputDepthKKrig[i];
				(*sdp->outData).errorK[i] = perturb.perturbationEK + outputErrorKKrig[i];
			}
		}
		perturbWeights.clear();
		perturb.riVector.clear();
		perturb.aiVector.clear();

		residualObservationsKrigedZ.clear();
		subX_indexKriged.clear();
		subY_indexKriged.clear();

//		indexKriged_VectorLocs.clear();
		xIndexKriged_Vector.clear();
		yIndexKriged_Vector.clear();

		outputDepthKrig.clear();
		outputErrorKrig.clear();
		if(sdp->PROP_UNCERT)
		{
			outputDepth0Krig.clear();
			outputError0Krig.clear();
			residualObservationsKrigedZ0.clear();
		}
		
		if(sdp->KALMAN)
		{
			outputDepthKKrig.clear();
			outputErrorKKrig.clear();
			residualObservationsKrigedZK.clear();
		}
	}
	return SUCCESS;
}



////This is the paul's version mergeBathy_v3.7.1_Paul simply updated to use the current scalecInterp_Perturbations call.
//// To use, comment this one in and comment the other one out.
//int scalecInterp_ProcessKrig(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores)
//{
//	//Initialize some variables
//	int i, j, k;
//	int outerLoop, innerLoop;
//	double tgs0, tgs1, tgs2;
//	double tgs1_Compute;
//	double tgs2_Compute;
//	double assnGridValue;
//	double standardDev;
//
//	vector<int> idx;
//	vector<int> idxInterp;
//	vector<double> subX_idx;
//	vector<double> subY_idx;
//	//vector<double> z_idx;
//	//vector<double> e_idx;
//	vector<double> pw_idx;
//	vector<double> subZ_idy;
//	vector<double> subE_idy;
//	vector<double> subH_idy;
//	vector<double> subV_idy;
//	vector<double> subX0;
//	vector<double> subY0;
//	int idySize;
//	PERTURBS perturb;
//	perturb.kernelName=*(*sdp).kernelName;
//
//	vector <double> perturbWeights;
//
//	double perturbationZKriged;
//	double perturbationEKriged;
//	double perturbationNEiKriged;
//	double perturbationREiKriged;
//	double thetaM, phi_2;
//	double zKriged;
//	double varZKriged;
//	double xGrid_indexKriged;
//	double yGrid_indexKriged;
//	double standardDevKriged;
//
//	vector<double> subX_idyKriged;
//	vector<double> subY_idyKriged;
//	vector<double> residualObservationsKrigedZ;
//	vector<double> subX_indexKriged;
//	vector<double> subY_indexKriged;
//	vector<double> indexKriged_VectorLocs;
//	vector<double> xIndexKriged_Vector;
//	vector<double> yIndexKriged_Vector;
//
//	vector<double> twoGammaHatVector;
//	vector<double> distanceVectorBinCenters;
//	vector<double> aVectorFine;
//	dgrid invGammaDArray;
//	dgrid AGrid;
//
//	vector<double> outputDepthKrig;
//	vector<double> outputErrorKrig;
//
//	double locSpacingX = (*sdp->spacingX);
//	double locSpacingY = (*sdp->spacingY);
//
//
//	double xmin, xmax, ymin, ymax;
//	double xmin2, xmax2, ymin2, ymax2;
//	tgs0 = 1.00;
//	vector<double>* slopes;
//
//	//************************************************************************************
//	// I. Interpolate the data
//	//************************************************************************************
//	//loop x tiles, columns
//	for (outerLoop = curIterNum; outerLoop < (*sdp->kx); outerLoop+=numCores)
//	{
//		//cout << "outerLoop: " << outerLoop << " ";
//		//1. Compute appropriate overlap between tiles horizontally and get the
//		//useful data. Recall that we scaled xi-array by 1./std(x). Then,
//		//get indices of x-coordinates in x-array for the data to be
//		//interpolated. Indices to be filtered further for useful y-coordinates.
//		//find tile limits
//		xmin = ((*sdp->minSubX) + (outerLoop*(*sdp->numStepsX))) - (*sdp->LMAX_x);
//		xmax = ((*sdp->minSubX) + ((outerLoop+1.0)*(*sdp->numStepsX))) + (*sdp->LMAX_x);
//
//		xmin2 = ((*sdp->minSubX) + (outerLoop*(*sdp->numStepsX)));
//		xmax2 = ((*sdp->minSubX) + ((outerLoop+1.0)*(*sdp->numStepsX)));
//
//
//		for (i = 0; i < (const int)(*sdp->subsampledData)[0].size(); i++){
//			if (((*sdp->subsampledData)[0][i] < xmax) && ((*sdp->subsampledData)[0][i] > xmin))
//			{
//				idx.push_back(i);
//			}
//		}
//
//		for (i = 0; i < (const int)(*sdp->xInterpVector).size(); i++){
//			if (((*sdp->xInterpVector)[i] < xmax2) && ((*sdp->xInterpVector)[i] >= xmin2))
//			{
//				idxInterp.push_back(i);
//			}
//		}
//		#pragma region --idxSize
//		if (idx.size() > 0)
//		{
//			for (innerLoop = 0; innerLoop < (*sdp->ky); innerLoop++)
//			{
//				//2. Compute appropriate overlap between tiles horizontally and get the
//				//useful data. Recall that we scaled xi-array by 1./std(x). Then,
//				//get indices of x-coordinates in x-array for the data to be
//				//interpolated. Indices to be filtered further for useful y-coordinates.
//				//find tile limits
//				ymin = ((*sdp->minSubY) + (innerLoop*(*sdp->numStepsY))) - (*sdp->LMAX_y);
//				ymax = ((*sdp->minSubY) + ((innerLoop+1.0)*(*sdp->numStepsY))) + (*sdp->LMAX_y);
//				ymin2 = ((*sdp->minSubY) + (innerLoop*(*sdp->numStepsY)));
//				ymax2 = ((*sdp->minSubY) + ((innerLoop+1.0)*(*sdp->numStepsY)));
//
//				for (int i = 0; i < (const int)idx.size(); i++)
//				{
//					if (((*sdp->subsampledData)[1][idx[i]] < ymax) && ((*sdp->subsampledData)[1][idx[i]] > ymin))
//					{
//						subX_idx.push_back((*sdp->subsampledData)[0][idx[i]] * (*sdp->Lx));
//						subY_idx.push_back((*sdp->subsampledData)[1][idx[i]] * (*sdp->Ly));
//						subX_idyKriged.push_back((*sdp->subsampledData)[0][idx[i]]);
//						subY_idyKriged.push_back((*sdp->subsampledData)[1][idx[i]]);
//						//z_idx.push_back((*sdp->subsampledData)[2][idx[i]]);
//						//e_idx.push_back((*sdp->subsampledData)[4][idx[i]]);
//						subZ_idy.push_back((*sdp->subsampledData)[2][idx[i]]);
//						subE_idy.push_back((*sdp->subsampledData)[4][idx[i]]);
//						subH_idy.push_back((*sdp->subsampledData)[5][idx[i]]);
//						subV_idy.push_back((*sdp->subsampledData)[6][idx[i]]);
//						
//						subX0.push_back((*sdp->x0_idx)[idx[i]]);	//scattered input data sub-tile
//						subY0.push_back((*sdp->y0_idx)[idx[i]]);
//					}
//				}
//				idySize = (const int)subX_idx.size();
//				
//				#pragma region idySize
//				if (idySize > 2)//if ((subX_idx.size() >= 2))
//				{
//					//E. Obtain a regular grid
//					// get number of indices in tile. idyi, col. idxi, row.
//					uint idyi = (uint)(*sdp->innerLoopIndexVector)[innerLoop].size(); 
//					uint idxi = (uint)outerLoopIndexVector.size(); 
//
//					// get starting indices of tile. i, col. j, row
//					i = ((*sdp->innerLoopIndexVector)[innerLoop])[0];
//					j = outerLoopIndexVector[0];
//
//					// get sub-tile of grid interpolation locations.
//					dgrid subXInterpLocs0;
//					dgrid subYInterpLocs0;
//					(*sdp->xInterpLocs0).subgrid(subXInterpLocs0, i, j, idyi, idxi);
//					(*sdp->yInterpLocs0).subgrid(subYInterpLocs0, i, j, idyi, idxi);
//
//					//F. Create Delaunay Tri from scattered input data sub-tile.
//					Bathy_Grid new_bathyGrid = Bathy_Grid();
//					new_bathyGrid.Construct_Tin(&subX0, &subY0, &subZ_idy, &subH_idy, &subV_idy);
//					
//					//G. Estimate depths at regular grid points
//					InterpGrid		new_gmt			= InterpGrid(GMT);					
//					vector<double>	subZInterpLocs0 = vector<double> ((subXInterpLocs0.vec()).size(), 0.00);
//					new_gmt.estimate(&subXInterpLocs0.vec(), &subYInterpLocs0.vec(), &subZInterpLocs0, 1.96, 2.00, *sdp->Lx, new_bathyGrid.getTin(), sdp->interpMethod);
//
//					//H. Find slopes at grid points
//					//slopes = new_gmt.getGrads()->getSlopeOut();
//					slopes = new_gmt.getGrads()->getSlopeOut_Vector();
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
//					//3. Precompute the weights, riVector, and aiVector across the whole interpolation plane
//					perturbWeights = vector<double>(z_idx.size(),2);
//					scalecInterpPerturbations_PreCompute(&z_idx, &e_idx, &perturbWeights, &perturb);
//
//					for (int i = 0; i < (const int)idxInterp.size(); i++)
//					{
//						yGrid_indexKriged = (*sdp->yInterpVector)[idxInterp[i]];
//
//						if ((yGrid_indexKriged < ymax2 && yGrid_indexKriged >= ymin2))
//						{
//							xGrid_indexKriged = (*sdp->xInterpVector)[idxInterp[i]];
//							indexKriged_VectorLocs.push_back(idxInterp[i]);
//							xIndexKriged_Vector.push_back(xGrid_indexKriged);
//							yIndexKriged_Vector.push_back(yGrid_indexKriged);
//						}
//					}
//
//					outputDepthKrig = vector<double>(indexKriged_VectorLocs.size(), 0);
//					outputErrorKrig = vector<double>(indexKriged_VectorLocs.size(), 0);
//
//					//4. Get interpolation values at observation locations
//					for (int i = 0; i < (const int)subX_idx.size(); i++)
//					{				
//						perturbationZKriged = 0;
//						perturbationEKriged = 1;
//						perturbationNEiKriged = 1;
//						perturbationREiKriged = 1;
//						tgs1_Compute = subX_idx[i];
//						tgs2_Compute = subY_idx[i];
//
//						scalecInterpPerturbations_Compute_ForKriging(&subX_idx, &subY_idx, &z_idx, &e_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, &perturb.riVector, &perturb.aiVector, (*sdp->neitol), &perturbationZKriged, &perturbationEKriged, &perturbationNEiKriged, &perturbationREiKriged, &standardDevKriged);
//		
//						if ((subX_idyKriged[i] >= xmin2 && subX_idyKriged[i] <= xmax2) && (subY_idyKriged[i] >= ymin2 && subY_idyKriged[i] <= ymax2))
//						{
//							residualObservationsKrigedZ.push_back(z_idx[i] - perturbationZKriged);
//							subX_indexKriged.push_back(subX_idyKriged[i]);
//							subY_indexKriged.push_back(subY_idyKriged[i]);
//						}
//					}
//
//					if (subX_indexKriged.size() >= 15)
//					{
//						//5. We need to subtile the kriged indexes otherwise the matrix inversion takes too long -- TODO
//						if(subX_indexKriged.size() > KRIGED_SIZE_THRESHOLD)
//						{
//							//6. Krig the data
//
//							ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKrig, &outputErrorKrig);
//
//						}else
//						{//small enough to do all at once
//							ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
//
//							//6. Krig the data
//							for (int i = 0; i < (const int)indexKriged_VectorLocs.size(); i++)
//							{
//								xGrid_indexKriged = xIndexKriged_Vector[i];
//								yGrid_indexKriged = yIndexKriged_Vector[i];
//
//								zKriged = 0.00;
//								varZKriged = 0.00;
//
//								ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);
//								outputDepthKrig[i] = zKriged;
//								outputErrorKrig[i] = varZKriged;
//							}
//							twoGammaHatVector.clear();
//							distanceVectorBinCenters.clear();
//							aVectorFine.clear();
//							invGammaDArray.clear();
//							AGrid.clear();
//						}
//					}
//
//					for (int i = 0; i < (const int)indexKriged_VectorLocs.size(); i++)
//					{
//						tgs1 = xIndexKriged_Vector[i];
//						tgs2 = yIndexKriged_Vector[i];
//						tgs1_Compute = tgs1 * (*sdp->Lx);
//						tgs2_Compute = tgs2 * (*sdp->Ly);
//
//						perturb.perturbationZ = 0;
//						perturb.perturbationE = 1;
//						perturb.perturbationNEi = 1;
//						perturb.perturbationREi = 1
//						if(stdp->PROP_UNCERT)
//						{
//							perturb.perturbationZ0 = 0.0;
//							perturb.perturbationE0 = 1.0;
//						}
//						if(stdp->KALMAN)
//						{
//							perturb.perturbationZK = 0.0;
//							perturb.perturbationEK = 1.0;
//						}
//						bool KRIGING=0;
//						//7. Compute the value
//						scalecInterpPerturbations_Compute(&subX_idx, &subY_idx, &z_idx, &e_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*sdp->neitol), dmin, (*slopes)[i], &subX0, &subY0, &Xiii, &perturb, stdp->MSE, stdp->PROP_UNCERT, stdp->KALMAN, KRIGING);
//
//						//8. Put assn grid calculation here..... do matrix * vector math.......
//						//put trend back into this tile
//						assnGridValue = tgs0*(*sdp->btrend)[0]+tgs1*(*sdp->btrend)[1]+tgs2*(*sdp->btrend)[2];
//	
//						(*sdp->outData).depth[(uint)indexKriged_VectorLocs[i]] = perturb.perturbationZ + assnGridValue + outputDepthKrig[i];
//						(*sdp->outData).error[(uint)indexKriged_VectorLocs[i]] = perturb.perturbationE + outputErrorKrig[i];
//						(*sdp->outData).nEi[(uint)indexKriged_VectorLocs[i]] = perturb.perturbationNEi;
//						(*sdp->outData).rEi[(uint)indexKriged_VectorLocs[i]] = perturb.perturbationREi;
//						(*sdp->outData).standardDev[(uint)indexKriged_VectorLocs[i]] = standardDev;
//						if(stdp->PROP_UNCERT)
//						{
//							(*stdp->outputDepth0)(iliv_Loc,oliv_Loc) = perturb.perturbationZ0 + assnGridValue;
//							(*stdp->outputError0)(iliv_Loc,oliv_Loc) = perturb.perturbationE0;
////							(*stdp->standardDev0)(iliv_Loc,oliv_Loc) = perturb.standardDev20;
//						}
//						if(stdp->KALMAN)
//						{
//							(*stdp->outputDepthK)(iliv_Loc,oliv_Loc) = perturb.perturbationZK + assnGridValue;
//							(*stdp->outputErrorK)(iliv_Loc,oliv_Loc) = perturb.perturbationEK;
//						}
//					}
//					#pragma endregion Interpolate
//					perturbWeights.clear();
//					perturb.riVector.clear();
//					perturb.aiVector.clear();
//					
//					new_bathyGrid.clear();
//					subXInterpLocs0.clear();
//					subYInterpLocs0.clear();
//					slopes = NULL;
//
//					residualObservationsKrigedZ.clear();
//					subX_indexKriged.clear();
//					subY_indexKriged.clear();
//
//					indexKriged_VectorLocs.clear();
//					xIndexKriged_Vector.clear();
//					yIndexKriged_Vector.clear();
//
//					outputDepthKrig.clear();
//					outputErrorKrig.clear();
//				}//subX_idx
//				//idySize
//				#pragma endregion idySize
//				//required clears!
//				subX_idx.clear();
//				subY_idx.clear();
//			//	z_idx.clear();
//			//	e_idx.clear();
//				subX_idyKriged.clear();
//				subY_idyKriged.clear();
//
//				subZ_idy.clear();
//				subE_idy.clear();
//				subH_idy.clear();
//				subV_idy.clear();
//				subX0.clear();
//				subY0.clear();
//			}//innerloop
//		}//idxSize
//		#pragma endregion idxSize
//		idx.clear();
//		idxInterp.clear();
//	}//outerloop
//
//	return SUCCESS;
//}
//
//




//This function is for both kriging and non-kriging runs instead of a separate kriging function.
//Don't believe this worked properly for multi-threading kriging so it was broken down into 4A-6A.
//It was left here in case that's what it was meant to do; this function works if not threading.
//scalecInterp_Process Part II -  Function was broken up to allow for multi-threading
int scalecInterp_Process2A(SCALEC_DATA_POINTER sdp, const int curIterNum, const int numCores)
{
	double dmin = NaN;	//Non-gridded output assumed

	double tgs0 = 1.0, tgs1, tgs2;
	double tgs1_Compute;
	double tgs2_Compute;
	double assnGridValue;
	
	double	Ni = (double)(*sdp->xInterpLocs0).size();
	dgrid Xiii(1,2);	// current interpolation location

	PERTURBS perturb;
	perturb.kernelName = *(*sdp).kernelName;

	//Initialize kriging structures
	double locSpacingX = (*sdp->spacingX);
	double locSpacingY = (*sdp->spacingY);
	double zKriged;
	double varZKriged;
	double z0Kriged;
	double varZ0Kriged;
	double zKKriged;
	double varZKKriged;
	double xGrid_indexKriged;
	double yGrid_indexKriged;
	
	vector<double> residualObservationsKrigedZ;
	vector<double> residualObservationsKrigedZ0;
	vector<double> residualObservationsKrigedZK;
	vector<double> subX_indexKriged; //observation locations
	vector<double> subY_indexKriged;
	//vector<double> indexKriged_VectorLocs;
	vector<double> xIndexKriged_Vector;
	vector<double> yIndexKriged_Vector;
	
	vector<double> twoGammaHatVector;
	vector<double> distanceVectorBinCenters;
	vector<double> aVectorFine;
	dgrid invGammaDArray;
	dgrid AGrid;

	vector<double> outputDepthKrig;
	vector<double> outputErrorKrig;
	vector<double> outputDepth0Krig;
	vector<double> outputError0Krig;
	vector<double> outputDepthKKrig;
	vector<double> outputErrorKKrig;

	//Initialize locations of which to obtain interpolation values
	//xGrid_indexKriged = (*sdp->xInterpVector)[idxInterp[i]];
	//indexKriged_VectorLocs.push_back(idxInterp[i]);
	xIndexKriged_Vector = vector<double>((*sdp->xInterpVector));
	yIndexKriged_Vector= vector<double>((*sdp->yInterpVector));
	//xIndexKriged_Vector = vector<double>((*sdp->xInterpLocs0));
	//yIndexKriged_Vector= vector<double>((*sdp->yInterpLocs0));
	
	//Initialize Krige Outputs
	outputDepthKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputErrorKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputDepth0Krig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputError0Krig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputDepthKKrig = vector<double>((*sdp->xInterpVector).size(), 0);
	outputErrorKKrig = vector<double>((*sdp->xInterpVector).size(), 0);

	//3. Pre-compute the weights, riVector, and aiVector across the whole interpolation plane
	vector <double> perturbWeights(sdp->z_idx->size(),2);
	scalecInterpPerturbations_PreCompute(sdp->z_idx, sdp->e_idx, &perturbWeights, &perturb);

	if(sdp->KRIGING)
	{
		//************************************************************************************
		//I. Begin Residual Kriging.
		//************************************************************************************
		//for (int i = curIterNum; i < Ni; i += numCores)
		for (int i = curIterNum; i < (const int)(*sdp->x_idx).size(); i += numCores)
		{
			#pragma region --Find the Residuals at Observations
			//************************************************************************************
			//4. Get interpolation values at observation locations; _Compute_ForKriging
			//************************************************************************************
			Xiii(0,0)		= (*sdp->x0_idx)[i]; //query pts will be the input locations of which we know the value.
			Xiii(0,1)		= (*sdp->y0_idx)[i];
		
			tgs1_Compute = (*sdp->x_idx)[i];	//current input x value in tile shouldn't matter cause we have an irregular grid
			tgs2_Compute = (*sdp->y_idx)[i];
			perturb.perturbationZKriged	= 0.0;
			perturb.perturbationEKriged	= 1.0;
			perturb.perturbationNEiKriged = 1.0;
			perturb.perturbationREiKriged = 1.0;
			if(sdp->PROP_UNCERT)
			{
				perturb.perturbationZ0Kriged = 0.0;
				perturb.perturbationE0Kriged = 1.0;
			}
			if(sdp->KALMAN)
			{
				perturb.perturbationZKKriged = 0.0;
				perturb.perturbationEKKriged = 1.0;
			}
			vector<double> slopesVec2 = vector<double>((*sdp->subsampledData)[0].size(), 1.00);

			//This was the __Compute_ForKriging call but was changed to use the current _Compute function.
			//A. Estimate the depth and errors at the observation location.
			scalecInterpPerturbations_Compute(sdp->x_idx, sdp->y_idx, sdp->z_idx, sdp->e_idx, sdp->h_idx, sdp->v_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*sdp->neitol), dmin, (slopesVec2)[i], sdp->x0_idx, sdp->y0_idx, &Xiii, &perturb, sdp->MSE, sdp->PROP_UNCERT, sdp->KALMAN, sdp->KRIGING);
		
			//B. Find the residual from the depth estimation and the actual value at the observation
			residualObservationsKrigedZ.push_back((*sdp->z_idx)[i] - perturb.perturbationZKriged);
			if(sdp->PROP_UNCERT)
				residualObservationsKrigedZ0.push_back((*sdp->z_idx)[i] - perturb.perturbationZ0Kriged);
			if(sdp->KALMAN)
				residualObservationsKrigedZK.push_back((*sdp->z_idx)[i] - perturb.perturbationZKKriged);

			//Add the observation location to a tile in order to subtile to reduce complexity.
			subX_indexKriged.push_back((*sdp->x_idxKriged)[i]); 
			subY_indexKriged.push_back((*sdp->y_idxKriged)[i]);
			//subX_indexKriged.push_back((*sdp->x0_idx)[i]); //don't know if this is right
			//subY_indexKriged.push_back((*sdp->y0_idx)[i]);
			//end Get interpolation values at observation locations
			#pragma endregion
		}
		#pragma region --Krige the Residuals
		//When tile is sufficiently large enough, krige it.
		if (subX_indexKriged.size() >= 15)
		{
			int KRIG_SUBTILE_FLAG = 0;
			//5. We need to subtile the kriged indexes; otherwise the matrix inversion takes too long -- TODO
			if(subX_indexKriged.size() > KRIGED_SIZE_THRESHOLD)
			{
				KRIG_SUBTILE_FLAG = 1;
				#pragma region --Subtile Kriged Indices
				try
				{
					//Too large, need to subtile.
					//6. Krige the residuals
					//_PreComputeTile calls _PostCompute from within
					ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKrig, &outputErrorKrig);

					if(sdp->PROP_UNCERT)
						ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepth0Krig, &outputError0Krig);
					if(sdp->KALMAN)
						ordinaryKrigingOfResiduals_PreComputeTile(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &xIndexKriged_Vector, &yIndexKriged_Vector, locSpacingX, locSpacingY, &outputDepthKKrig, &outputErrorKKrig);
				}
				catch (exception &e) {
					cout << "An exception occurred. Exception Thrown: " << e.what() << '\n';
					cout << "Kriging Subtiles failed.  Will try kriging without subtiles.  This will take longer." << endl;
					KRIG_SUBTILE_FLAG = 0;
				}
				#pragma endregion
			}
			if(!KRIG_SUBTILE_FLAG)
			{
				#pragma region --Do Not Subtile Kriged Indices
				//small enough to do all at once
				ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
				if(sdp->PROP_UNCERT)
					ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ0, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);
				if(sdp->KALMAN)
					ordinaryKrigingOfResiduals_PreCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZK, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);

				//6. Krige the residuals
				for (int i = 0; i < (const int)(*sdp->xInterpVector).size(); i++)
				{
					xGrid_indexKriged = xIndexKriged_Vector[i]; //current output grid loc 
					yGrid_indexKriged = yIndexKriged_Vector[i];

					zKriged = 0.00;
					varZKriged = 0.00;
					
					ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);
					outputDepthKrig[i] = zKriged;
					outputErrorKrig[i] = varZKriged;
	
					if(sdp->PROP_UNCERT)
					{
						z0Kriged = 0.00;
						varZ0Kriged = 0.00;
						ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &z0Kriged, &varZ0Kriged);
						outputDepth0Krig[i] = z0Kriged;
						outputError0Krig[i] = varZ0Kriged;
					}
					if(sdp->KALMAN)
					{
						zKKriged = 0.00;
						varZKKriged = 0.00;
						ordinaryKrigingOfResiduals_PostCompute(&subX_indexKriged, &subY_indexKriged, &residualObservationsKrigedZ, &xGrid_indexKriged, &yGrid_indexKriged, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKKriged, &varZKKriged);
						outputDepthKKrig[i] = zKKriged;
						outputErrorKKrig[i] = varZKKriged;
					}
				}
				twoGammaHatVector.clear();
				distanceVectorBinCenters.clear();
				aVectorFine.clear();
				invGammaDArray.clear();
				AGrid.clear();
				#pragma endregion
			}
		}
		//Done Kriging Residuals
		#pragma endregion 
	}
	//Continue normal processing
	double outDepthKrig;
	double outErrorKrig;
	double outDepth0Krig;
	double outError0Krig;
	double outDepthKKrig;
	double outErrorKKrig;

	for (int i = curIterNum; i < Ni; i += numCores)
	{
		//tgs1 = xIndexKriged_Vector[i];			//shouldn't use this because we have an irregular grid? same as xInterpsLocs0 and xInterpVector
		//tgs2 = yIndexKriged_Vector[i];			//get current yValue

		tgs1			= (*sdp->xInterpVector)[i];		//get current xValue
		tgs2			= (*sdp->yInterpVector)[i];		//get current yValue
		tgs1_Compute	= tgs1 * (*sdp->Lx);			//scale xValue
		tgs2_Compute	= tgs2 * (*sdp->Ly);			//scale yValue
		Xiii(0,0)		= (*sdp->xInterpLocs0)[i];
		Xiii(0,1)		= (*sdp->yInterpLocs0)[i];

		//Initialize output fields
		perturb.perturbationZ = 0.0;
		perturb.perturbationE = 1.0;
		perturb.perturbationNEi = 1.0;
		perturb.perturbationREi = 1.0;
		if(sdp->PROP_UNCERT)
		{
			perturb.perturbationZ0 = 0.0;
			perturb.perturbationE0 = 1.0;
		}
		if(sdp->KALMAN)
		{
			perturb.perturbationZK = 0.0;
			perturb.perturbationEK = 1.0;
		}
		
		if(sdp->KRIGING)
		{
			outDepthKrig  = outputDepthKrig[i];
			outErrorKrig  = outputErrorKrig[i];
			if(sdp->PROP_UNCERT)
			{
				outDepth0Krig = outputDepth0Krig[i];
				outError0Krig = outputError0Krig[i];
			}
			if(sdp->KALMAN)
			{
				outDepthKKrig = outputDepthKKrig[i];
				outErrorKKrig = outputErrorKKrig[i];
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

		//4. Compute the value
		scalecInterpPerturbations_Compute(sdp->x_idx, sdp->y_idx, sdp->z_idx, sdp->e_idx, sdp->h_idx, sdp->v_idx, &tgs1_Compute, &tgs2_Compute, &perturbWeights, (*sdp->neitol), dmin, (*(sdp->slopeOut))[i], sdp->x0_idx, sdp->y0_idx, &Xiii, &perturb, sdp->MSE, sdp->PROP_UNCERT, sdp->KALMAN, 0);//0 because were not kriging the residuals here... that's done above

		//5. Put assn grid calculation here..... do matrix * vector math.......
		//	put trend back into this tile
		assnGridValue = tgs0*(*sdp->btrend)[0]+tgs1*(*sdp->btrend)[1]+tgs2*(*sdp->btrend)[2];		
		(*sdp->outData).depth[i] = perturb.perturbationZ + assnGridValue + outDepthKrig;
		(*sdp->outData).error[i] = perturb.perturbationE + outErrorKrig;
		(*sdp->outData).nEi[i] = perturb.perturbationNEi;
		(*sdp->outData).rEi[i] = perturb.perturbationREi;
		(*sdp->outData).standardDev[i] = perturb.standardDev2;
		if(sdp->PROP_UNCERT)
		{
			(*sdp->outData).depth0[i] = perturb.perturbationZ0 + assnGridValue + outDepth0Krig;
			(*sdp->outData).error0[i] = perturb.perturbationE0 + outError0Krig;
		//	(*sdp->outData).standardDev0[i] = perturb.standardDev20;
		}
		if(sdp->KALMAN)
		{
			(*sdp->outData).depthK[i] = perturb.perturbationZK + assnGridValue + outDepthKKrig;
			(*sdp->outData).errorK[i] = perturb.perturbationEK + outErrorKKrig;
		}
	}
	perturbWeights.clear();
	perturb.riVector.clear();
	perturb.aiVector.clear();
	residualObservationsKrigedZ.clear();
	
	subX_indexKriged.clear();
	subY_indexKriged.clear();

//		indexKriged_VectorLocs.clear();
	xIndexKriged_Vector.clear();
	yIndexKriged_Vector.clear();

	outputDepthKrig.clear();
	outputErrorKrig.clear();
	if(sdp->PROP_UNCERT)
	{
		outputDepth0Krig.clear();
		outputError0Krig.clear();
		residualObservationsKrigedZ0.clear();
	}	
	if(sdp->KALMAN)
	{
		outputDepthKKrig.clear();
		outputErrorKKrig.clear();
		residualObservationsKrigedZK.clear();
	}

	return SUCCESS;
}
