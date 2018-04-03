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
//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	#include "../WarningStates.h"	//Disable all Warnings!!!
#endif

extern "C"
{
	#include "./MB_ZGrid/mb_zgrid.h"
	#include "./GMT_Surface/processSurface.h"
}

//Restore warning state -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( pop )
#endif

#include "externalInterpolators.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include "constants.h"
#include "grid.h"
#include "LatLong-UTMconversion.h"
//#include "Error_Estimator/estimator.h"
#include "ALG/interpolation.h"			//ALGSpline
#include <time.h>
#include <map>
#include "Error_Estimator/Bathy_Grid.h"
#include "ALG/ap.h"
#include "standardOperations.h"
#include <math.h>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/natural_neighbor_coordinates_2.h>
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Delaunay_triangulation_2<K>             Delaunay_triangulation;
//typedef std::vector< std::pair< K::Point_2, K::FT  > >
//                                                      Point_coordinate_vector;

//UNIX
#if !defined(_MSC_VER)//(_WIN32) || !defined(_WIN64)
/* We are on Windows */
# define strcpy_s strcpy
#endif

//Input grid
vector<double> grid_xvi;
vector<double> grid_yvi;
vector<double> grid_zvi;

//Output grid
vector<double> grid_xv;
vector<double> grid_yv;
vector<double> grid_zv;

////Execute a maximum of THREAD_POOL_SIZE threads on a maximum
////of CPU_POOL_SIZE central processing units.
//#define THREAD_POOL_SIZE    8

//************************************************************************************
// SUBROUTINE I: Constructor for externalInterpolators
//************************************************************************************
externalInterpolators::externalInterpolators(double NorthingRef, double EastingRef, double rotAngle, int RefEllip, char ZoneRef[4])
{
	//Set values for use in converting Lat/Lon to UTM and back again
	UTMNorthingRef = NorthingRef;
	UTMEastingRef = EastingRef;
	rotationAngle = rotAngle;
	refEllipsoid = RefEllip;
	strcpy_s(UTMZoneRef, ZoneRef);
}

//************************************************************************************
// SUBROUTINE II: Function call for running MB ZGrid
//************************************************************************************
bool externalInterpolators::run_MB_ZGrid(vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *h, vector<double> *v, double x0, double y0, double x1, double y1, std::map<string, int> additionalOptions, double spacingX, double spacingY, double tension, string z_OutputFileName, int usage, Bathy_Grid* bathyGrid)
{
	//There was a loop to perform Monte Carlo Simulations here (see older versions) but this didn't make sense so it was removed. SJZ
	cout << endl << "***WARNING: externalInterpolators.cpp: Ensemble test cases force x0=y0=15000 and x1=y1=41000. This is hardcoded in and must be fixed or taken out!***" << endl << endl;
	if (additionalOptions.find("-ensembleTests")->second == 1)
	{
		cout << "Ensemble Tests Detected! Values for presplining set!" << endl;
		cout << "If not performing Ensemble Tests of a seamount and ridge than rename inputfilelist (dir and name) w/o 'Ensemble'" << endl;
		cout << "Old x0 = y0 =" << x0 << endl;
		cout << "Old x1 = y1 =" << y1 << endl;
		x0 = y0 = 15000;// 15800;
		x1 = y1 = 41000;// 16200;
		cout << "New x0 = y0 =" << x0 << endl;
		cout << "New x1 = y1 =" << y1 << endl;
	}

	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	ofstream outFile;
	double UTMNorthing;
	double UTMEasting;
	double newLat, newLon;
	double startTime, stopTime, compTime;

	int i, j, k;
	//Input Data
	float *mb_z = NULL;
	float *mb_xyz = NULL;

	float mb_x1 = (float)x0;
	float mb_y1 = (float)y0;
	float mb_dx = (float)spacingX;//mbzg.dx;
	float mb_dy = (float)spacingY;//mbzg.dy;

	float mb_cay = (float)tension;//mbzg.cay;
	int mb_n = (const int)(*x).size();
	int mb_nx; //interpolated area x size
	int mb_ny; //interpolated area y size
	int mb_nrng = 1000; //mbzg.nrng;

	//2.  Working arrays
	float *mb_zpij = NULL;
	int *mb_knxt = NULL;
	int *mb_imnew = NULL;

	startTime = clock();

	//************************************************************************************
	//I. Calculate the grid sizes based on the input data
	//************************************************************************************
	vector<double> xt;
	vector<double> yt;
	//A. Calculate X
	double locationValue = x0;
	int xtSize = 0;
	while (locationValue <= x1) {
		xt.push_back(locationValue);
		locationValue += spacingX;
		xtSize += 1;
	}
	//B. Calculate Y
	locationValue = y0;
	int ytSize = 0;
	while (locationValue <= y1) {
		yt.push_back(locationValue);
		locationValue += spacingY;
		ytSize += 1;
	}
	cout << "Number of computed rows: " << ytSize << endl;
	cout << "Number of computed cols: " << xtSize << endl;

	vector<double> xMeshVector = vector<double>(xtSize*ytSize);
	vector<double> yMeshVector = vector<double>(xtSize*ytSize);

	//C. Set the vectors
	int clx = 0;
	int cly = 0;
	int currentLoc = 0;
	for (int i = 0; i < ytSize; i++)
	{
		currentLoc = i;
		for (int j = 0; j < xtSize; j++)
		{
			xMeshVector[currentLoc] = xt[clx];
			yMeshVector[currentLoc] = yt[cly];
			currentLoc = currentLoc + ytSize;
			clx += 1;
		}
		clx = 0;
		cly += 1;
	}

	//D.  Now put data in the variables
	mb_nx = xtSize;
	mb_ny = ytSize;
	mb_z = (float*)calloc((mb_nx*mb_ny), sizeof(float));
	mb_xyz = (float*)calloc((3*(mb_n)), sizeof(float));
	mb_zpij = (float*)calloc(mb_n, sizeof(float));
	mb_knxt = (int*)calloc(mb_n, sizeof(int));
	mb_imnew = (int*)calloc((mb_nx+mb_ny), sizeof(int));

	//ALP for Unix
	//memset((char *)mb_z,0,mb_nx*mb_ny*sizeof(float));

	//E.  Initialize the input data structure
	j = 0;
	for (int i = 0; i < mb_n; i++){
		mb_xyz[j++] = (float)(*x).at(i);
		mb_xyz[j++] = (float)(*y).at(i);
		mb_xyz[j++] = (float)(*z).at(i);
	}

	//************************************************************************************
	//II. Call MB ZGrid
	//************************************************************************************
	mb_zgrid(mb_z, &mb_nx, &mb_ny, &mb_x1, &mb_y1,
			 &mb_dx, &mb_dy, mb_xyz, &mb_n, mb_zpij,
			 mb_knxt, mb_imnew, &mb_cay, &mb_nrng);

	dgrid zGrid_temp = dgrid(mb_ny, mb_nx);

	//A. Re align the data
	k = 0;
	for (i = 0; i < mb_ny; i++){
		for (j = 0; j < mb_nx; j++){
			zGrid_temp(i,j) = mb_z[k];
			k++;
		}
	}
	k = 0;
	vector<double> zVector(xtSize*ytSize);
	vector<double> eVector(xtSize*ytSize, 0);// this is a problem if we want to use as input
	for (i = 0; i < (const int)zGrid_temp.cols(); i++){
		for (j = 0; j < (const int)zGrid_temp.rows(); j++){
			zVector[k] = zGrid_temp(j,i);
			k += 1;
		}
	}

	stopTime = clock();
	compTime = stopTime-startTime;
	cout << "Time to Complete MB_ZGrid Interpolation: " << (stopTime-startTime)/CLOCKS_PER_SEC << endl << endl;
	if (additionalOptions.find("-inputInMeters")->second == 1)
	{
		cout << "MB_Zgrid output in (x, y) meters; no UTM conversions will be calculated." << endl << endl;
	}

	//************************************************************************************
	//III. Convert back to Lon/Lat and output the file
	//************************************************************************************
	vector<double> lonSurf = vector<double>((ytSize*xtSize));
	vector<double> latSurf = vector<double>((ytSize*xtSize));
	//Print the results of MB_ZGrid
	string ext;
	std::string outFileTemp;
	ext = (".txt");
	outFileTemp = z_OutputFileName.c_str() + ext;
	outFile.open(outFileTemp.c_str());
	outFile.precision(6);
	outFile.setf(std::ios::fixed, std::ios::floatfield);
	if (outFile.is_open())
	{
		for(int i = 0; i < (ytSize*xtSize); i++)
		{
			if (additionalOptions.find("-inputInMeters")->second == 0)
			{
				UTMEasting = (xMeshVector[i])*cos(deg2rad*(-rotationAngle))	- (yMeshVector[i])*sin(deg2rad*(-rotationAngle));
				UTMNorthing = (xMeshVector[i])*sin(deg2rad*(-rotationAngle)) + (yMeshVector[i])*cos(deg2rad*(-rotationAngle));
				
				UTMNorthing += UTMNorthingRef;
				UTMEasting += UTMEastingRef;

				UTMtoLL(refEllipsoid, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
				if(usage > 0)
					outFile << newLon << "\t" << newLat << "\t" << zVector[i] << endl;
				lonSurf[i] = newLon;
				latSurf[i] = newLat;
			}
			else
			{
				if(usage > 0)
					outFile << xMeshVector[i] << "\t" << yMeshVector[i] << "\t" << zVector[i] << endl;
				lonSurf[i] = xMeshVector[i];
				latSurf[i] = yMeshVector[i];
			}
		}
	
		//************************************************************************************
		//IV. Compute the uncertainty associated with each data point
		//************************************************************************************
		if (usage < 0 || abs(usage) == 2){ //Must compute error to  use as input!
			double minSpacing = spacingX;
			if (spacingY < spacingX)
				minSpacing = spacingY;

			InterpGrid* mbz= new InterpGrid(MBZ);
			string interpMethod = "BILINEAR";
			if(additionalOptions.find("-nnInterp")->second==1)
				interpMethod = "NN";
			mbz->estimate(&xMeshVector, &yMeshVector, &zVector,1.96, 2.00, minSpacing,bathyGrid->getTin(), interpMethod, "", "");
			if (usage == -2) //Store if to use as input for ensembling later
				bathyGrid->addToList(mbz);

			stopTime = clock();
			cout << "Time to Complete MB_ZGrid Interpolation with Errors: " << (compTime+(stopTime-startTime))/CLOCKS_PER_SEC << endl << endl;

			//A. Output the GMT Surface depth calculation and the uncertainty computed by the uncertainty estimator
			vector<double>e		= mbz->getE();
			vector<double>eS2	= mbz->getE2();
			vector<double>eS3	= mbz->getE3();
			vector<double>eS4	= mbz->getE4();
			vector<double>eS5	= mbz->getE5();
			double etemp		= standardDeviation(&e,false);
			for(int i = 0; i < (const int)xMeshVector.size(); i++)
			{
				if(e[i] == 0)
					e[i] = etemp;
				outFile << lonSurf[i] << "\t" << latSurf[i] << "\t" << zVector[i] << "\t" << e[i] <<endl;//<< "\t" << eS2[i]  << "\t" << eS3[i]  << "\t" << eS4[i]  << "\t" << eS5[i] << endl;
			}
			eVector = e;
		}	
		outFile.close();
	}else
	{
		cerr << "Failed to open MB ZGrid output file" << endl;
		return false;
	}

	//************************************************************************************
	//V. If we are going to reuse this output as input for mergeBathy then assign the variables here
	//************************************************************************************
	if ((usage == 2) || (usage == -2)) //Use as input
	{
		(*x).clear();
		(*y).clear();
		(*z).clear();
		(*e).clear();
		(*h).clear();
		(*v).clear();

		(*x) = vector<double>(xtSize*ytSize);
		(*y) = vector<double>(xtSize*ytSize);
		(*z) = vector<double>(xtSize*ytSize);
		(*e) = vector<double>(xtSize*ytSize);
		(*h) = vector<double>(xtSize*ytSize);
		(*v) = vector<double>(xtSize*ytSize);
		double temp = 0;
		for(int i = 0; i < (xtSize*ytSize); i++){
			(*x)[i] = xMeshVector[i];
			(*y)[i] = yMeshVector[i];
			(*z)[i] = zVector[i];
			(*e)[i] = eVector[i];
			
			//temp = sqrt(pow(eVector[i],2)/2);
			(*h)[i] = 0;//temp;
			(*v)[i] = eVector[i];//temp;
		}
	}

	xMeshVector.clear();
	yMeshVector.clear();
	zVector.clear();
	eVector.clear();
	xt.clear();
	yt.clear();
	lonSurf.clear();
	latSurf.clear();

	free(mb_z);
	free(mb_xyz);
	free(mb_zpij);
	free(mb_knxt);
	free(mb_imnew);

	mb_z = NULL;
	mb_xyz = NULL;
	mb_zpij = NULL;
	mb_knxt = NULL;
	mb_imnew = NULL;

	return true;
}

//************************************************************************************
// SUBROUTINE III: Function call for running GMT Surface
//************************************************************************************
bool externalInterpolators::run_Surface(vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *h, vector<double> *v, map<string, int> additionalOptions, double spacingX, double spacingY, double tension, string z_OutputFileName, double scaleFactor, double alpha, int usage, Bathy_Grid* bathyGrid)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	ofstream outFile;

	double UTMNorthing;
	double UTMEasting;
	double newLat, newLon;
	double minSpacing = spacingX;
	double startTime, stopTime, compTime;

	double *xConverted = NULL;
	double *yConverted = NULL;
	double *zConverted = NULL;

	double *xPostSurface = NULL;
	double *yPostSurface = NULL;
	double *zPostSurface = NULL;

	vector<double> xSurf;
	vector<double> ySurf;
	vector<double> lonSurf;
	vector<double> latSurf;
	vector<double> zSurf;
	vector<double> eSurf;

	double z0 = (*z)[0];
	double z1 = (*z)[0];
	double x0 = (*x)[0];
	double x1 = (*x)[0];
	double y0 = (*y)[0];
	double y1 = (*y)[0];
	int i;//, j, k;
	int postSurfaceSize;
	int returnValue;

	startTime = clock();
	if (spacingY < spacingX)
		minSpacing = spacingY;

	xConverted = (double*)calloc(((*x).size()), sizeof(double));
	yConverted = (double*)calloc(((*y).size()), sizeof(double));
	zConverted = (double*)calloc(((*z).size()), sizeof(double));
	if((xConverted==NULL || yConverted==NULL) || zConverted==NULL)
	{
		cout << "Could not allocate memory! Data set is too large. Exiting..." << endl;
		exit (1);
	}
	//A. Calculate min and max depths
	for (i = 0; i < (const int)(*x).size(); i++)
	{
		if ((*z)[i] < z0)
			z0 = (*z)[i];
		else if ((*z)[i] > z1)
			z1 = (*z)[i];

		if ((*x)[i] < x0)
			x0 = (*x)[i];
		else if ((*x)[i] > x1)
			x1 = (*x)[i];

		if ((*y)[i] < y0)
			y0 = (*y)[i];
		else if ((*y)[i] > y1)
			y1 = (*y)[i];

		xConverted[i] = (*x)[i];
		yConverted[i] = (*y)[i];
		zConverted[i] = (*z)[i];
	}

	//************************************************************************************
	//I. Call GMT Surface
	//************************************************************************************
	//There was a loop to perform Monte Carlo Simulations here (see older versions) but this didn't make sense so it was removed. SJZ
	cout << endl <<  "***WARNING: externalInterpolators.cpp: Ensemble test cases force x0=y0=15000 and x1=y1=41000. This is hardcoded in and must be fixed or taken out!***" << endl << endl;
	if (additionalOptions.find("-ensembleTests")->second == 1)
	{
		cout << "Ensemble Tests Detected! Values for presplining set!" << endl;
		cout << "If not performing Ensemble Tests of a seamount and ridge than rename inputfilelist (dir and name) w/o 'Ensemble'" << endl;
		cout << "Old x0 = y0 =" << x0 << endl;
		cout << "Old x1 = y1 =" << y1 << endl;
		x0 = y0 = 15000;// 15800;
		x1 = y1 = 41000;// 30080; //16200;
		cout << "New x0 = y0 =" << x0 << endl;
		cout << "New x1 = y1 =" << y1 << endl;
	}
	returnValue = processSurface(xConverted, yConverted, zConverted, (int)(*x).size(), x0, y0, z0, x1, y1, z1, spacingX, spacingY, tension, &xPostSurface, &yPostSurface, &zPostSurface, &postSurfaceSize);

	stopTime = clock();
	compTime = stopTime-startTime;
	cout << "Time to Complete GMT Surface Interpolation: " << (stopTime-startTime)/CLOCKS_PER_SEC << " seconds." << endl << endl;
	if (additionalOptions.find("-inputInMeters")->second == 1)
	{
		cout << "GMT_Surface output in (x, y) meters; no UTM conversions will be calculated." << endl << endl;
	}

	xSurf = vector<double>(postSurfaceSize);
	ySurf = vector<double>(postSurfaceSize);
	lonSurf = vector<double>(postSurfaceSize);
	latSurf = vector<double>(postSurfaceSize);
	zSurf = vector<double>(postSurfaceSize);
	eSurf = vector<double>(postSurfaceSize);
		
	//************************************************************************************
	//II - single run: Convert back to Lon/Lat and output lat, lon to file
	//                 if additionalOptions.find("-inputInMeters")->second == 0
	//				   Otherwise no conversion and output in meters to file.
	//************************************************************************************
	string ext;
	std::string outFileTemp;
	ext = (".txt");
	outFileTemp = z_OutputFileName.c_str() + ext;		
	outFile.open(outFileTemp.c_str());
	outFile.precision(6);
	outFile.setf(std::ios::fixed, std::ios::floatfield);
	if (outFile.is_open())
	{
		for(int i = 0; i < (postSurfaceSize); i++){
			xSurf[i] = xPostSurface[i];
			ySurf[i] = yPostSurface[i];
			zSurf[i] = zPostSurface[i];
			eSurf[i] = 0.00; // this is a problem if we want to use as input

			if (additionalOptions.find("-inputInMeters")->second == 0)
			{
				UTMEasting = (xSurf[i])*cos(deg2rad*(-rotationAngle)) - (ySurf[i])*sin(deg2rad*(-rotationAngle));
				UTMNorthing = (xSurf[i])*sin(deg2rad*(-rotationAngle)) + (ySurf[i])*cos(deg2rad*(-rotationAngle));

				UTMNorthing += UTMNorthingRef;
				UTMEasting += UTMEastingRef;

				UTMtoLL(refEllipsoid, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
				if(usage > 0)
					outFile << newLon << "\t" << newLat << "\t" << zSurf[i] << endl;
				lonSurf[i] = newLon;
				latSurf[i] = newLat;
			}
			else
			{
				if(usage > 0)
					outFile << xSurf[i] << "\t" << ySurf[i] << "\t" << zSurf[i] << endl;
				lonSurf[i] = xSurf[i];
				latSurf[i] = ySurf[i];
			}
		}

		startTime = clock();

		//************************************************************************************
		//III. Compute the uncertainty associated with each data point
		//************************************************************************************
		if (usage < 0 || abs(usage) == 2){ //Must compute error to  use as input!
			InterpGrid* gmt = new InterpGrid(GMT);
			string interpMethod = "BILINEAR";
			if(additionalOptions.find("-nnInterp")->second==1)
				interpMethod = "NN";
			gmt->estimate( &xSurf, &ySurf, &zSurf, scaleFactor, alpha, minSpacing, bathyGrid->getTin(), interpMethod, "", "");
			if (usage == -2) //Store if to use as input for ensembling later
				bathyGrid->addToList(gmt);

			stopTime = clock();
			cout << "Time to Complete GMT Surface Interpolation with Errors: " << (compTime+(stopTime-startTime))/CLOCKS_PER_SEC << endl << endl;

			//A.  Output the GMT Surface depth calculation and the uncertainty computed by the uncertainty estimator
			vector<double> e   = gmt->getE();
			vector<double> eS2 = gmt->getE2();
			vector<double> eS3 = gmt->getE3();
			vector<double> eS4 = gmt->getE4();
			vector<double> eS5 = gmt->getE5();
			double etemp = standardDeviation(&e,false);
			for(int i = 0; i < (const int)xSurf.size(); i++)
			{
				if(e[i] == 0)
					e[i] = etemp;
				if(usage < 0)
					outFile << lonSurf[i] << "\t" << latSurf[i] << "\t" << zSurf[i] << "\t" << e[i]  <<endl;//<< "\t" << eS2[i]  << "\t" << eS3[i]  << "\t" << eS4[i]  << "\t" << eS5[i] << endl;
			}
			eSurf = e;
		}
		outFile.close();
	}
	else
	{
		cerr << "Failed to open GMT Surface output file" << endl;
		return false;
	}

	//************************************************************************************
	//IV. If we are going to reuse this output as input for mergeBathy then assign the variables here
	//************************************************************************************
	if ((usage == 2) || (usage == -2)) //Use as input
	{
		(*x).clear();
		(*y).clear();
		(*z).clear();
		(*e).clear();
		(*h).clear();
		(*v).clear();

		(*x) = vector<double>(postSurfaceSize);
		(*y) = vector<double>(postSurfaceSize);
		(*z) = vector<double>(postSurfaceSize);
		(*e) = vector<double>(postSurfaceSize);
		(*h) = vector<double>(postSurfaceSize);
		(*v) = vector<double>(postSurfaceSize);
		double temp = 0;
		for(int i = 0; i < postSurfaceSize; i++)
		{
			(*x)[i] = xSurf[i];
			(*y)[i] = ySurf[i];
			(*z)[i] = zSurf[i];
			(*e)[i] = eSurf[i];

			//temp = sqrt(pow(eSurf[i],2)/2);
			(*h)[i] = 0;//temp;
			(*v)[i] = eSurf[i];//temp;
		}
	}

	//A. Clear up variables
	free(xConverted);
	free(yConverted);
	free(zConverted);

	free(xPostSurface);
	free(yPostSurface);
	free(zPostSurface);

	xConverted = NULL;
	yConverted = NULL;
	zConverted = NULL;

	xPostSurface = NULL;
	yPostSurface = NULL;
	zPostSurface = NULL;

	xSurf.clear();
	ySurf.clear();
	lonSurf.clear();
	latSurf.clear();
	zSurf.clear();
	eSurf.clear();

	return true;
}

//************************************************************************************
// Method for providing access to ALG spline functionality
//************************************************************************************
bool externalInterpolators::run_ALGSpline(vector<double> *pxin_gvect, vector<double> *pyin_gvect, vector<double> *pzin_gvect, double utm_xmin,  double utm_xmax, double utm_ymin,  double utm_ymax, double xin_size, double yin_size, vector<double> *e, vector<double> *h, vector<double> *v, map<string, int> additionalOptions, double spacingX, double spacingY, double tension, string z_OutputFileName, double scaleFactor, double alpha, int usage, Bathy_Grid* bathyGrid)
{
	cout << "THIS DOES NOT WORK! This function expects input to be in a grid already and it is not. The function was changed under the assumption it did not have to be a grid but this fails when building the spline.  See older versions for original function write which expected a grid prior to these changes.  Note neither one works since we don't have a grid beforehand.  If we use a grid computed from GMT or MBZ, then modify this function properly for gridded input and use the previous version as a reference." << endl;
		return false;

	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	ofstream outFile;

	double newLat, newLon;
	double minSpacing = spacingX;
	double startTime, stopTime, compTime=0;
	if (spacingY < spacingX)
		minSpacing = spacingY;

	vector<double> xSurf;
	vector<double> ySurf;
	vector<double> lonSurf;
	vector<double> latSurf;
	vector<double> zSurf;
	vector<double> eSurf;

	cout << "externalInterpolators::run_ALGSpline( )" << endl;
	cout << "\t Number of UTM Eastings(vector size)	= " << pxin_gvect->size() << "\n";
	cout << "\t Number of UTM Northings(vector size)= " << pyin_gvect->size() << "\n";
	cout << "\t Number of depths(vector size)		= " << pzin_gvect->size() << "\n";
	cout << "\t X spacing                           = " << spacingX << "\n";
	cout << "\t Y spacing                           = " << spacingY << "\n";
	cout << "\t Output filename                     = " << z_OutputFileName << "\n";
	cout << "\t Usage                               = " << usage << "\n";
	cout << "\t UTM Easting Ref                     = " << UTMEastingRef << "\n";
	cout << "\t UTM Northing Ref                    = " << UTMNorthingRef << "\n";
	cout << "\t UTM X Rel (min.;max.)               = [" << utm_xmin << "; " << utm_xmax << "]\n";
	cout << "\t UTM Y Rel (min.;max.)               = [" << utm_ymin << "; " << utm_ymax << "]\n\n";

	double utm_loc   = 0.0;
	double utm_xsize = 0.0;
	double utm_ysize = 0.0;
	
	//************************************************************************************
	//I. Calculate the grid sizes based on the input data
	//************************************************************************************
	#pragma region --Build UTM Grids
	//Sam saved the min y and x into pUTMs so that its grid is the same as MBZ at the least.
	//Don't know if this is right, maybe alg's grid is suppose to be smaller than mbz 
	//like how gmt is larger than mbz due to extra padding.
	//If this is changed than the ensemble function needs to change with it.
	//In order to avoid the ensemble change, all the interps need to be sorted the same.
	
	//A. Get UTM grid sizes
	utm_loc = utm_xmin;
	while (utm_loc <= utm_xmax)
	{
		utm_loc   = utm_loc + spacingX;
		utm_xsize = utm_xsize + 1;
	}

	utm_loc = utm_ymin;
	while (utm_loc <= utm_ymax)
	{
		utm_loc   = utm_loc + spacingY;
		utm_ysize = utm_ysize + 1;
	}
	//Done computing grid sizes

	double *pUTM_xt = new double [(int)utm_xsize];
	double *pUTM_yt = new double [(int)utm_ysize];
	
	//B. Build UTM grids
	utm_loc = utm_xmin;
	pUTM_xt[0] = utm_loc; //Added so grid dimensions are not smaller than MBZ
	for (int i = 1; i < utm_xsize; i++)
	{
		utm_loc   = utm_loc + spacingX;
		pUTM_xt[i] = utm_loc;
	}

	utm_loc = utm_ymin;
	pUTM_yt[0] = utm_loc; //Added so grid dimensions are not smaller than MBZ
	for (int i = 1; i < utm_ysize; i++)
	{
		utm_loc   = utm_loc + spacingY;
		pUTM_yt[i] = utm_loc;
	}
	//Done building UTM grids
	#pragma endregion

	cout << "Number of computed rows: " << utm_ysize << endl;
	cout << "Number of computed cols: " << utm_xsize << endl;

	#pragma region --Building X and Y Arrays (Input)
	int xsize = (int)pxin_gvect->size();
	int ysize = (int)pyin_gvect->size();
	int zsize = (int)pzin_gvect->size();

	int curLoc = 0;
	double *px = new double [(int)xsize];
	double *py = new double [(int)ysize];
	double *pz = new double [zsize];

	cout << "\nBuilding x array...\n";
	for (int j = 0; j < xsize; j++)
	{
		px[j] = pxin_gvect->at(j); 
	}
	cout << "\nBuilding y array...\n";
	for (int i = 0; i < ysize; i++)
	{
		py[i] = pyin_gvect->at(i);
	}
	cout << "\nBuilding z array...\n";
	for (int j = 0; j < zsize; j++)
	{
		pz[j] = pzin_gvect->at(j);
	}
	#pragma endregion
		
	cout << "Z count = " << curLoc << endl;
	cout << "Product = " << yin_size*xin_size << endl;

	//Reformat data to conform to ALG C++ interface
	alglib::real_1d_array			r1d_x;
	alglib::real_1d_array			r1d_y;
	r1d_x.setcontent((alglib::ae_int_t)xin_size,px);
	r1d_y.setcontent((alglib::ae_int_t)yin_size,py);

	alglib::real_2d_array			r2d_z;
	r2d_z.setlength((alglib::ae_int_t)yin_size,(alglib::ae_int_t)xin_size);
	r2d_z.setcontent((alglib::ae_int_t)yin_size,(alglib::ae_int_t)xin_size,pz);

	startTime = clock();
	
	//Construct both bilinear and bicubic spline objects; only keep one later.
	alglib::spline2dinterpolant s2di_bilinear;
	alglib::spline2dbuildbilinear(r1d_x, r1d_y, r2d_z, (alglib::ae_int_t)yin_size,(alglib::ae_int_t)xin_size, s2di_bilinear);
		
	if (additionalOptions.find("-inputInMeters")->second == 1)
	{
		cout << "ALG_Spline output in (x, y) meters; no UTM conversions will be calculated." << endl << endl;
	}

	//Create new bilinear spline grid
	double	bilinear_splinevalue = -999.0;
	double	bicubic_splinevalue = -999.0;
	double	latitude	= -999.0;
	double	longitude	= -999.0;
	double	UTMEasting  = 0.0;
	double	UTMNorthing = 0.0;
		
	xSurf = vector<double>((const int)(utm_xsize*utm_ysize));
	ySurf = vector<double>((const int)(utm_xsize*utm_ysize));
	zSurf = vector<double>((const int)(utm_xsize*utm_ysize));
	lonSurf = vector<double>((const int)(utm_xsize*utm_ysize));
	latSurf = vector<double>((const int)(utm_xsize*utm_ysize));
	eSurf = vector<double>((const int)(utm_xsize*utm_ysize));
	int k = 0;

	//************************************************************************************
	//II - single run: Convert back to Lon/Lat and output lat, lon to file
	//                 if additionalOptions.find("-inputInMeters")->second == 0
	//				   Otherwise no conversion and output in meters to file.
	//************************************************************************************
	string ext;
	std::string outFileTemp;
	ext = (".txt");
	outFileTemp = z_OutputFileName.c_str() + ext;
	outFile.open(outFileTemp.c_str());
	outFile.precision(6);
	outFile.setf(std::ios::fixed, std::ios::floatfield);
	if (outFile.is_open())
	{
		for (int i = 0; i < utm_ysize; i++)
		{
			//cout << "Processing UTM row = " << i << " of " << utm_ysize << "\n";
			for (int j = 0; j < utm_xsize; j++)
			{
				bilinear_splinevalue = alglib::spline2dcalc(s2di_bilinear,pUTM_xt[j],pUTM_yt[i]);
				xSurf[k] = pUTM_xt[j];
				ySurf[k] = pUTM_yt[i];
				zSurf[k] = bilinear_splinevalue;
				eSurf[k] = 0.00;// this is a problem if we want to use as input
				if (additionalOptions.find("-inputInMeters")->second == 0)
				{
					UTMEasting		= (xSurf[i])*cos(deg2rad*(-rotationAngle)) - (ySurf[i])*sin(deg2rad*(-rotationAngle));
					UTMNorthing		= (xSurf[i])*sin(deg2rad*(-rotationAngle)) + (ySurf[i])*cos(deg2rad*(-rotationAngle));

					UTMNorthing += UTMNorthingRef;
					UTMEasting += UTMEastingRef;
	
					UTMtoLL(refEllipsoid, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
					if(usage > 0)
						outFile << newLon << "\t" << newLat << "\t" << zSurf[i] << endl;
					lonSurf[i] = newLon;
					latSurf[i] = newLat;
				}
				else
				{
					if(usage > 0)
						outFile << xSurf[k] << "\t" << ySurf[k] << "\t" << zSurf[k] << endl;
					lonSurf[k] = xSurf[k];
					latSurf[k] = ySurf[k];
				}
				k++;
			}
		}
		stopTime = clock();
		compTime = stopTime-startTime;
		cout << "Time to Complete ALG Spline Interpolation: " << (stopTime-startTime)/CLOCKS_PER_SEC << " seconds." << endl << endl;

		//************************************************************************************
		//III. Compute the uncertainty associated with each data point
		//************************************************************************************
		if (usage < 0 || abs(usage) == 2){ //Must compute error to  use as input!
			InterpGrid* algSpline = new InterpGrid(ALGSpline);
			string interpMethod = "BILINEAR";
			if(additionalOptions.find("-nnInterp")->second==1)
				interpMethod = "NN";
			algSpline->estimate( &xSurf, &ySurf, &zSurf, scaleFactor, alpha, spacingX, bathyGrid->getTin(), interpMethod, "", "");
			if (usage == -2) //Store if to use as input for ensembling later
				bathyGrid->addToList(algSpline);

			stopTime = clock();
			cout << "Time to Complete ALGSpline Interpolation with Uncertainty: " << (compTime+=(stopTime-startTime))/CLOCKS_PER_SEC << endl << endl;

			//A.  Output the ALGSpline depth calculation and the uncertainty computed by the uncertainty estimator
			vector<double> e   = algSpline->getE();
			vector<double> eS2 = algSpline->getE2();
			vector<double> eS3 = algSpline->getE3();
			vector<double> eS4 = algSpline->getE4();
			vector<double> eS5 = algSpline->getE5();
			double etemp = standardDeviation(&e,false);
			for(int i = 0; i < (const int)xSurf.size(); i++)
			{
				if(e[i] == 0)
					e[i] = etemp;
				outFile << lonSurf[i] << "\t" << latSurf[i] << "\t" << zSurf[i] << "\t" << e[i]  <<endl;//<< "\t" << eS2[i]  << "\t" << eS3[i]  << "\t" << eS4[i]  << "\t" << eS5[i] << endl;
			}
			eSurf = e;
		}
		outFile.close();
	}
	else
	{
		cerr << "Failed to open ALGSpline Uncertainty output file" << endl;
		return false;
	}
		
	//************************************************************************************
	//IV. If we are going to reuse this output as input for mergeBathy then assign the variables here
	//************************************************************************************
	if ((usage == 2) || (usage == -2)) //Use as input
	{
		(*pxin_gvect).clear();
		(*pyin_gvect).clear();
		(*pzin_gvect).clear();
		(*e).clear();
		(*h).clear();
		(*v).clear();

		(*pxin_gvect) = vector<double>((const int)(utm_xsize*utm_ysize));
		(*pyin_gvect) = vector<double>((const int)(utm_xsize*utm_ysize));
		(*pzin_gvect) = vector<double>((const int)(utm_xsize*utm_ysize));
		(*e) = vector<double>((const int)(utm_xsize*utm_ysize));
		(*h) = vector<double>((const int)(utm_xsize*utm_ysize));
		(*v) = vector<double>((const int)(utm_xsize*utm_ysize));

		double temp = 0;
		for(int i = 0; i < (utm_xsize*utm_ysize); i++){
			(*pxin_gvect)[i] = xSurf[i];
			(*pyin_gvect)[i] = ySurf[i];
			(*pzin_gvect)[i] = zSurf[i];
			(*e)[i] = eSurf[i];

			//temp = sqrt(pow(eSurf[i],2)/2);
			(*h)[i] = 0;//temp;
			(*v)[i] = eSurf[i];//temp;
		}
	}

	xSurf.clear();
	ySurf.clear();
	lonSurf.clear();
	latSurf.clear();
	zSurf.clear();
	eSurf.clear();

	//Free allocated pointers
	delete [] px;
	delete [] py;
	delete [] pz;
	delete [] pUTM_xt;
	delete [] pUTM_yt;

	cout << "run_ALGSpline finished...\n" ;

	return true;
}

int externalInterpolators::write_output_file(string chAbsFilename, const vector<double>& grid_xv, const vector<double>& grid_yv, const vector<double>& grid_zv)
{
	int status = 0;
	try
	{
		cout << "Writing: " << chAbsFilename << " ";
		int xcnt   = (int)grid_xv.size();
		int ycnt   = (int)grid_yv.size();
		int zcnt   = xcnt*ycnt;
		cout << "Writing " << zcnt << " values to file " << chAbsFilename << " ";

		//UNIX ALP
		char cfilename[256];
		for (int i=0; i< (const int)chAbsFilename.length(); i++) cfilename[i]=chAbsFilename.at(i);

		ofstream outFile;
		outFile.open(cfilename,ios_base::out);
		for (int xndx = 0; xndx < xcnt; xndx++)
		{
			for (int yndx = 0; yndx < ycnt; yndx++)
			{
				outFile << grid_xv[xndx] << " " << grid_yv[yndx] << " " << grid_zv[yndx + (ycnt*xndx)] << "\n";
			}
		}
		outFile.close();
		cout << "(Done)\n";
	}
	catch(...)
	{
		status = -1;
		cout << "Exception in write_output_file()\n";
	}

	status = 1;

	return status;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/* For SIBSON and Other interpolator imports see CGAL examples and the CGAL_CONSOLE_TEST application for prior work.
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//************************************************************************************
// Method for providing access to Sibson Natural Neighbor
// functionality
//************************************************************************************
bool externalInterpolators::
	run_Sibson_Natural_Neighbor_Interp(vector<double> *pxin_gvect, vector<double>	*pyin_gvect, vector<double>	*pzin_gvect, double utm_xmin,  double utm_xmax, double utm_ymin,  double utm_ymax,double xin_size, double yin_size, vector<double> *e, vector<double> *h, vector<double> *v, map<string, int> additionalOptions, double spacingX, double spacingY, double tension, string z_OutputFileName, double scaleFactor, double alpha, int usage, Bathy_Grid* bathyGrid)
{
	cout << "Executing Sibson's Natural Neighbor Interpolation...\n";
	cout << "externalInterpolators::run_Sibson_Natural_Neighbor_Interp( )" << endl;
	cout << "\t Number of UTM eastings(vector size)	= " << pxin_gvect->size() << "\n";
	cout << "\t Number of UTM northings(vector size)= " << pyin_gvect->size() << "\n";
	cout << "\t Number of depths(vector size)		= " << pzin_gvect->size() << "\n";
	cout << "\t X spacing                           = " << spacingX << "\n";
	cout << "\t Y spacing                           = " << spacingY << "\n";
	cout << "\t Output filename                     = " << z_OutputFileName << "\n";
	cout << "\t Usage                               = " << usage << "\n";
	cout << "\t UTM Easting Ref                     = " << UTMEastingRef << "\n";
	cout << "\t UTM Northing Ref                    = " << UTMNorthingRef << "\n";
	cout << "\t UTM X Rel (min.;max.)               = [" << utm_xmin << "; " << utm_xmax << "]\n";
	cout << "\t UTM Y Rel (min.;max.)               = [" << utm_ymin << "; " << utm_ymax << "]\n";

//	double startTime, stopTime, compTime;
	double utm_loc   = 0.0;
	double utm_xsize = 0.0;
	double utm_ysize = 0.0;

	////ATTEMPT1:
	//Delaunay_triangulation dt;
	//for (int y=0 ; y<3 ; y++)
	//	for (int x=0 ; x<3 ; x++)
	//		dt.insert(K::Point_2(x,y));
	////coordinate computation
	//K::Point_2 p(1.2, 0.7);
	//Point_coordinate_vector coords;
	//CGAL::Triple< std::back_insert_iterator<Point_coordinate_vector>,	K::FT, bool> result = CGAL::natural_neighbor_coordinates_2(dt, p, std::back_inserter(coords));
	//if(!result.third){
	//	std::cout << "The coordinate computation was not successful." << std::endl;
	//	std::cout << "The point (" <<p << ") lies outside the convex hull."	<< std::endl;
	//}
	//K::FT  norm = result.second;
	//std::cout << "Coordinate computation successful." << std::endl;
	//std::cout << "Normalization factor: " <<norm << std::endl;
	//std::cout << "done" << std::endl;
	//return 0;

	////ATTEMPT 2:
	// //INTERPOLATION:
	//Point_coordinate_vector     coords;
	//std::pair<Coord_type, bool> interpolation_result;
	//Coord_type value = 0; // initialization to remove compiler warning
	//int n = points.size();
	//ITraits traits;

	//std::cout << "Interpolation at  "<<n  <<" grid points " << std::endl;
	//CGAL::Triple< std::back_insert_iterator<Point_coordinate_vector>, K::FT, bool> coordinate_result = CGAL::natural_neighbor_coordinates_2(T, points[i], std::back_inserter(coords));
	//K::FT norm = coordinate_result.second;
	////test if the computation was successful
	//assert(coordinate_result.third && norm>0);

	//interpolation_result = CGAL::sibson_c1_interpolation(coords.begin(),coords.end(),
 //      					  norm, points[i],
 //      					  CGAL::Data_access<Point_value_map>
 //      					  (values),
 //      					  CGAL::Data_access<Point_vector_map>
 //      					  (gradients), traits); break;




	#pragma region --Build UTM Grids
	//Sam saved the min y and x into pUTMs so that its grid is the same as MBZ at the least.
	//Don't know if this is right, maybe alg's grid is suppose to be smaller than mbz 
	//like how gmt is larger than mbz due to extra padding.
	//If this is changed than the ensemble function needs to change with it.
	//In order to avoid the ensemble change, all the interps need to be sorted the same.
	
	//A. Get UTM grid sizes
	utm_loc = utm_xmin;
	while (utm_loc <= utm_xmax)
	{
		utm_loc   = utm_loc + spacingX;
		utm_xsize = utm_xsize + 1;
	}

	utm_loc = utm_ymin;
	while (utm_loc <= utm_ymax)
	{
		utm_loc   = utm_loc + spacingY;
		utm_ysize = utm_ysize + 1;
	}
	//Done computing grid sizes

	double *pUTM_xt = new double [(int)utm_xsize];
	double *pUTM_yt = new double [(int)utm_ysize];
	
	//B. Build UTM grids
	utm_loc = utm_xmin;
	pUTM_xt[0] = utm_loc; //Added so grid dimensions are not smaller than MBZ
	for (int i = 0; i < utm_xsize; i++)
	{
		utm_loc   = utm_loc + spacingX;
		pUTM_xt[i] = utm_loc;
	}

	utm_loc = utm_ymin;
	pUTM_yt[0] = utm_loc; //Added so grid dimensions are not smaller than MBZ
	for (int i = 0; i < utm_ysize; i++)
	{
		utm_loc   = utm_loc + spacingY;
		pUTM_yt[i] = utm_loc;
	}
	//Done building UTM grids
	#pragma endregion

	cout << "Number of computed rows: " << utm_ysize << endl;
	cout << "Number of computed cols: " << utm_xsize << endl;

	#pragma region --Building X and Y Arrays (Input)
	int xsize = (int)pxin_gvect->size();
	int ysize = (int)pyin_gvect->size();
	int zsize = (int)pzin_gvect->size();

	int curLoc = 0;
	double *px = new double [(int)xsize];
	double *py = new double [(int)ysize];
	double *pz = new double [zsize];

	cout << "\nBuilding x array...\n";
	for (int j = 0; j < xsize; j++)
	{
		px[j] = pxin_gvect->at(j); 
	}
	cout << "\nBuilding y array...\n";
	for (int i = 0; i < ysize; i++)
	{
		py[i] = pyin_gvect->at(i);
	}
	cout << "\nBuilding z array...\n";
	for (int j = 0; j < zsize; j++)
	{
		pz[j] = pzin_gvect->at(j);
	}
	#pragma endregion


	//To pickup working were left off, uncomment from here down
	const int THREAD_POOL_SIZE=1;
	HANDLE      hThreads[THREAD_POOL_SIZE];
	int         slot = 0;
	DWORD       threadID;
	DWORD       rc;
	
	int status			= 0;
	int xcnt			= xin_size;//(int)grid_xv.size();
	int ycnt			= yin_size;//(int)grid_yv.size();
	int zcnt			= xcnt*ycnt;
	double ncols_total	= xcnt;
	vector<double> grid_zv_sub = vector<double>(ycnt);
	grid_zv				= vector<double>(zcnt);
	vector<double> col_strip_start;
	vector<double> col_strip_end;

	grid_segmentor(ncols_total,THREAD_POOL_SIZE,col_strip_start,col_strip_end);
	
	int    nTasks      = (int)col_strip_start.size();
	COLUMN_RNG_STRUCT  param_rng_struct;

	cout << "Number of threads available in pool = " << THREAD_POOL_SIZE << "\n";
	cout << "Number of strips to process         = " << col_strip_start.size() << "\n";

	for (int task_ndx = 1; task_ndx <= nTasks; task_ndx++)
	{
		double ncols_in_strip = col_strip_end.at(task_ndx-1) - col_strip_start.at(task_ndx-1) + 1;
		cout << "Column Strip(TASK) # = " << task_ndx << "\t" << col_strip_start.at(task_ndx-1)-1 << "\t" << col_strip_end.at(task_ndx-1)-1 << ";(" << ncols_in_strip << ")\n";
	
		//There are no more threads available in the pool; wait for one to become available before creating another.
		if (task_ndx > THREAD_POOL_SIZE)
		{
			cout << "\tWaiting for thread in pool to become available...\n";
			rc	 = WaitForMultipleObjects(THREAD_POOL_SIZE,hThreads,FALSE,INFINITE);
			slot = rc - WAIT_OBJECT_0;
		}

		param_rng_struct.begin = (int)col_strip_start.at(task_ndx-1)-1;
		param_rng_struct.end   = (int)col_strip_end.at(task_ndx-1)-1;

		hThreads[slot++] = CreateThread(NULL,0,compute_sibson_ntrlnbr_interp_col_thread,(LPVOID) &param_rng_struct,0,&threadID);
		cout << "\tThread in pool available; slot " << slot << "\n";
	}

	//Ensure all threads in pool are finished
	cout << "Waiting for all threads in pool to finish...\n";
	rc = WaitForMultipleObjects(THREAD_POOL_SIZE,hThreads,TRUE,INFINITE);

	//Terminate all threads in pool
	for(slot=0;slot<THREAD_POOL_SIZE;slot++)
	{
		CloseHandle(hThreads[slot]);
		cout << "Terminating thread; slot " << slot << "\n";
	}

	///////////////////
	//Original(Verified); single threaded version
	//for (int xndx = 0; xndx < xcnt; xndx++)
	//{
	//	std::cout << "Processing column :" << xndx+1 << " of " << xcnt << "\n";

		//Execute with a boost::thread
		//boost::thread thrd1(boost::bind(&compute_sibson_ntrlnbr_interp_col,xndx));
		//thrd1.join();
	
	    //Execute as a single main thread
		//status = compute_sibson_ntrlnbr_interp_col(xndx);
	//}

	
	
	//	return status;//sam

	status = write_output_file(z_OutputFileName,grid_xv,grid_yv,grid_zv);

	//for (int i = 0; i < utm_ysize; i++)
	//{
	//	cout << "Processing UTM row = " << i << " of " << utm_ysize << "\n";
	//	for (int j = 0; j < utm_xsize; j++)
	//	{
	//		bilinear_splinevalue		= alglib::spline2dcalc(s2di_bilinear,pUTM_xt[j],pUTM_yt[i]);
	//
	//		xSurf[k] = pUTM_xt[j];
	//		ySurf[k] = pUTM_yt[i]; //sam
	//		zSurf[k] = bilinear_splinevalue;
	//		eSurf[k] = 0.00;
	//		if (additionalOptions.find("-inputInMeters")->second == 0)
	//		{
	//
	//			UTMEasting		= (xSurf[k])*cos(deg2rad*(-rotationAngle)) - (ySurf[k])*sin(deg2rad*(-rotationAngle));
	//			UTMNorthing		= (xSurf[k])*sin(deg2rad*(-rotationAngle)) + (ySurf[k])*cos(deg2rad*(-rotationAngle));
	//			UTMEasting		= xSurf[k] + UTMEastingRef;//sam
	//			UTMNorthing		= ySurf[k] + UTMNorthingRef;

	//

	//			UTMtoLL(refEllipsoid, UTMNorthing, UTMEasting, UTMZoneRef, latitude, longitude);
	//			outFile << longitude << "\t" << latitude << "\t" << bilinear_splinevalue << endl;
	//			lonSurf[k] = longitude;
	//			latSurf[k] = latitude;
	//
	//		}
	//		else
	//		{
	//
	//			outFile << xSurf[k] << "\t" << ySurf[k] << "\t" << zSurf[k] << endl;
	//			lonSurf[k] = xSurf[k];
	//			latSurf[k] = ySurf[k];
	//		}
	//		k++;
	//	}
	//}
	//outFile.close();
	//
	
	////SAM
	//startTime = clock();

	////************************************************************************************
	////III. Compute the uncertainty associated with each data point
	////************************************************************************************
	//if (usage < 0){
	//	int inputDebug;

	//	InterpGrid* sibson_Natural_Neighbor_Interp = new InterpGrid(Sibson_Natural_Neighbor_Interp);
	//
	//	//sibson_Natural_Neighbor_Interp->TestDT_Example(pxin_gvect,pyin_gvect,pzin_gvect, h,
	//	//v,&xSurf, &ySurf, &zSurf, scaleFactor, alpha, spacingX);//sam

	////	algSpline->estimate( &pUTM_xt, &pUTM_yt, &zsurf, scaleFactor, alpha, spacingX, bathyGrid->getTin());
	//	sibson_Natural_Neighbor_Interp->estimate( &xSurf, &ySurf, &zSurf, scaleFactor, alpha, spacingX, bathyGrid->getTin());//sam
	//	bathyGrid->addToList(sibson_Natural_Neighbor_Interp);
	//
	//
	//
	//	stopTime = clock();
	//	cout << "Time to Complete Sibson_Natural_Neighbor_Interp Interpolation with Errors: " << (compTime+(stopTime-startTime))/CLOCKS_PER_SEC << endl << endl;

	//	//A.  Output the Sibson_Natural_Neighbor_Interp depth calculation and the uncertainty computed by the uncertainty estimator
	//	z_OutputFileName.append("_xyde.txt");
	//	//if(!sibson_Natural_Neighbor_Interp->printInterpErr(z_OutputFileName)) return false;
	//	outFile.open(z_OutputFileName.c_str());
	//	outFile.precision(6);
	//	outFile.setf(std::ios::fixed, std::ios::floatfield);
	//	if (outFile.is_open())
	//	{
	//		vector<double>e=sibson_Natural_Neighbor_Interp->getE();
	//		vector<double>eS2=sibson_Natural_Neighbor_Interp->getE2();
	//		vector<double>eS3=sibson_Natural_Neighbor_Interp->getE3();
	//		vector<double>eS4=sibson_Natural_Neighbor_Interp->getE4();
	//		vector<double>eS5=sibson_Natural_Neighbor_Interp->getE5();
	//		//algSpline->printInterpErr(outFile);
	//		for(int i = 0; i < xSurf.size(); i++){
	//			outFile << lonSurf[i] << "\t" << latSurf[i] << "\t" << zSurf[i] << "\t" << e[i]  << "\t" << eS2[i]  << "\t" << eS3[i]  << "\t" << eS4[i]  << "\t" << eS5[i] << endl;
	//		}
	//outFile.close();
	//		//exit (1);
	//	}else
	//	{
	//		cerr << "Failed to open Sibson_Natural_Neighbor_Interp uncertainty output file" << endl;
	//		//exit (1);
	//		return false;
	//	}
	//}
	////************************************************************************************
	////IV. If we are going to reuse this output as input for mergeBathy then assign the variables here
	////************************************************************************************
	//if ((usage == 2) || (usage == -2)) //Use as input
	//{
	//	(*pxin_gvect).clear();
	//	(*pyin_gvect).clear();
	//	(*pzin_gvect).clear();
	//	(*e).clear();
	//	(*pxin_gvect) = vector<double>(utm_xsize*utm_ysize);
	//	(*pyin_gvect) = vector<double>(utm_xsize*utm_ysize);
	//	(*pzin_gvect) = vector<double>(utm_xsize*utm_ysize);
	//	(*e) = vector<double>(utm_xsize*utm_ysize);

	//	for(int i = 0; i < (utm_xsize*utm_ysize); i++){
	//		(*pxin_gvect)[i] = xSurf[i];
	//		(*pyin_gvect)[i] = ySurf[i];
	//		(*pzin_gvect)[i] = zSurf[i];
	//		(*e)[i] = eSurf[i];
	//	}
	//}

	////double g=xin_size*utm_ysize;
	////free(g);
	////free(g);
	////free(g);

	////xPostSurface = NULL;
	////yPostSurface = NULL;
	////zPostSurface = NULL;

	//xSurf.clear();
	//ySurf.clear();
	//lonSurf.clear();
	//latSurf.clear();
	//zSurf.clear();
	//eSurf.clear();

	////Free allocated pointers
	//delete [] px;
	//delete [] py;
	//delete [] pz;
	//delete [] pUTM_xt;
	//delete [] pUTM_yt;

	cout << "Sibson_Natural_Neighbor_Interp finished...\n" ;

	return true;
}
void grid_segmentor(double num_columns,double num_threads,std::vector<double>& col_strip_start,std::vector<double>& col_strip_end)
{
	double prcnt_colseg = floor((num_columns/num_threads) + .5) / num_columns * (100.0);
	std::cout << "prcnt_colseg = " << prcnt_colseg << "\n";
    double min_colseg          = 10;
    double ncols_per_tile      = estimate_cols_per_tile(num_columns,prcnt_colseg,min_colseg);
	std::cout << "ncols_per_tile = " << ncols_per_tile << "\n";

    build_coltile_set(num_columns,ncols_per_tile,col_strip_start,col_strip_end);
	
	double ncolTiles  = (double)col_strip_start.size();
    
    //If last segment has less than min_colseg then delete it and add the
    //cols from this last segment to the previous segment
	double last_col_cnt = col_strip_end.at(((int)(ncolTiles-1))) - col_strip_start.at((int)(ncolTiles-1))+1; 
    if (last_col_cnt < min_colseg)
	{
        std::cout << "WARNING: Col segment count fault; merging columns with previous segment.\n";
        
        col_strip_end.at(((int)ncolTiles-2)) = col_strip_end.at(((int)ncolTiles-2)) + last_col_cnt;
        
		col_strip_start.pop_back();
		col_strip_end.pop_back();
	}   
}
DWORD WINAPI compute_sibson_ntrlnbr_interp_col_thread(LPVOID param)
{
	COLUMN_RNG_STRUCT* pRangeStruct = (COLUMN_RNG_STRUCT*)param;

	int xndx_start = pRangeStruct->begin;
	int xndx_end   = pRangeStruct->end;

	//Represents a strip from the main grid matrix
	//xndx_start = start column for this thread
	//xndx_end   = end column for this thread
	int status = 0;
	int xcnt   = (int)grid_xv.size();
	int ycnt   = (int)grid_yv.size();

	//create a ndx boundary test method; future
	for (int xndx = xndx_start; xndx <= xndx_end; xndx++)
	{
		for (int yndx = 0; yndx < ycnt; yndx++)
		{
			K::Point_2 p(grid_xv[xndx],grid_yv[yndx]);
			std::vector<std::pair<Point, Coord_type>> coords;
			Coord_type norm = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords)).second;
			Coord_type res  = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, CGAL::Data_access<std::map<Point, Coord_type, K::Less_xy_2>>(function_values));
			grid_zv[yndx + (ycnt*xndx)] = res;
		}
	}
	
	status = 1;

	return (DWORD)status;
}
int build_input_grid_locs(const std::vector<double>& xvect,    \
	                      const std::vector<double>& yvect,    \
						  const std::vector<double>& zvect,    \
						  std::vector<double>& grid_xvi,       \
						  std::vector<double>& grid_yvi,       \
						  std::vector<double>& grid_zvi,       \
						  double dx, double dy)
{
	//User has input a gridded file.
	int status = 0;

	
	try
	{
		int		zcnt = (int)zvect.size();
		std::cout << "Building input grid locations...";
		double dxmin_in = vector_min(xvect);
		double dxmax_in = vector_max(xvect);
		double dymin_in = vector_min(yvect);
		double dymax_in = vector_max(yvect);

		double  tempvalx= dxmin_in;
		int     xcnt    = 0;
		while (tempvalx <= dxmax_in)
		{
			grid_xvi.push_back(tempvalx);
			tempvalx	= tempvalx + dx;
			xcnt		= xcnt + 1;
		}
	
		double  tempvaly= dymin_in;
		int     ycnt	= 0;
		while (tempvaly <= dymax_in)
		{
			grid_yvi.push_back(tempvaly);
			tempvaly	= tempvaly + dy;
			ycnt		= ycnt + 1;
		}

		int  xycnt = xcnt*ycnt;
		grid_zvi		= std::vector<double>(xycnt);

		std::cout << "(Done)\n";
		std::cout << "\t" << xcnt  << "  X computed grid locations.\n";
		std::cout << "\t" << ycnt  << "  Y computed grid locations.\n";
		std::cout << "\t" << xycnt << " XY Product; Grid Sample Count.\n";
		std::cout << "\t" << zcnt  << "  Z input.\n";
		std::cout << "\t X min = " << dxmin_in << "\n";
		std::cout << "\t X max = " << dxmax_in << "\n";
		std::cout << "\t Y min = " << dymin_in << "\n";
		std::cout << "\t Y max = " << dymax_in << "\n";
		if (zcnt != xycnt)
		{
			std::cout << "Dataset input not a grid.\n";
			status = -1;
		}
		else
		{
			std::cout << "Dataset input is a grid.\n";
			status = 1;

			for (int yndx = 0; yndx < ycnt; yndx++)
			{
				for (int xndx = 0; xndx < xcnt; xndx++)
				{
					grid_zvi[xndx + (xcnt*yndx)] = zvect[xndx + (xcnt*yndx)];
				}
			}
		}
	}
	catch(...)
	{
		status = -1;
		std::cout << "Exception in build_input_grid_locs()\n";
	}

	return status;

}


int build_output_grid_locs(const std::vector<double>& xvect,   \
						   const std::vector<double>& yvect,   \
						   double dx, double dy) 
{
	int status = 0;

	try
	{
		std::cout << "Building output grid locations...";
		double dxmin_in = vector_min(xvect);
		double dxmax_in = vector_max(xvect);
		double dymin_in = vector_min(yvect);
		double dymax_in = vector_max(yvect);

		double  tempvalx= dxmin_in;
		int     xcnt    = 0;
		while (tempvalx <= dxmax_in)
		{
			grid_xv.push_back(tempvalx);
			tempvalx	= tempvalx + dx;
			xcnt		= xcnt + 1;
		}
	
		double  tempvaly= dymin_in;
		int     ycnt	= 0;
		while (tempvaly <= dymax_in)
		{
			grid_yv.push_back(tempvaly);
			tempvaly	= tempvaly + dy;
			ycnt		= ycnt + 1;
		}

		std::cout << "(Done)\n";
		std::cout << "\t" << xcnt << " X grid locations.\n";
		std::cout << "\t" << ycnt << " Y grid locations.\n";
		std::cout << "\t X min = " << dxmin_in << "\n";
		std::cout << "\t X max = " << dxmax_in << "\n";
		std::cout << "\t Y min = " << dymin_in << "\n";
		std::cout << "\t Y max = " << dymax_in << "\n";
	}
	catch(...)
	{
		status = -1;
		std::cout << "Exception in build_output_grid_locs()\n";
	}

	status = 1;

	return status;
}*/
//End Sibson


//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ORIGINAL FUNCTIONS FOR MONTE CARLO RUNS AND KRIGING!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

////************************************************************************************
//// SUBROUTINE II: Function call for running MB ZGrid
////************************************************************************************
//bool externalInterpolators::run_MB_ZGrid_ORIGINAL(vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *h, vector<double> *v, double x0, double y0, double x1, double y1, double spacingX, double spacingY, double tension, string z_OutputFileName, int usage)
//{
//	//************************************************************************************
//	// 0. Declare local variables and objects
//	//************************************************************************************
//	ofstream outFile;
//	double UTMNorthing;
//	double UTMEasting;
//	double newLat, newLon;
//	double startTime, stopTime, compTime;
//
//	int i, j, k;
//	//Input Data
//	float *mb_z = NULL;
//	float *mb_xyz = NULL;
//
//	float mb_x1 = (float)x0;
//	float mb_y1 = (float)y0;
//	float mb_dx = (float)spacingX;//mbzg.dx;
//	float mb_dy = (float)spacingY;//mbzg.dy;
//
//	float mb_cay = (float)tension;//mbzg.cay;
//	int mb_n = (*x).size();
//	int mb_nx; //interpolated area x size
//	int mb_ny; //interpolated area y size
//	int mb_nrng = 1000; //mbzg.nrng;
//
//	//2.  Working arrays
//	float *mb_zpij = NULL;
//	int *mb_knxt = NULL;
//	int *mb_imnew = NULL;
//
//	startTime = clock();
//
//	//************************************************************************************
//	//I. Calculate the grid sizes based on the input data
//	//************************************************************************************
//	vector<double> xt;
//	vector<double> yt;
//	//A. Calculate X
//	double locationValue = x0;
//	int xtSize = 0;
//	while (locationValue <= x1) {
//		xt.push_back(locationValue);
//		locationValue += spacingX;
//		xtSize += 1;
//	}
//	//B. Calculate Y
//	locationValue = y0;
//	int ytSize = 0;
//	while (locationValue <= y1) {
//		yt.push_back(locationValue);
//		locationValue += spacingY;
//		ytSize += 1;
//	}
//
//	vector<double> xMeshVector = vector<double>(xtSize*ytSize);
//	vector<double> yMeshVector = vector<double>(xtSize*ytSize);
//
//	//C. Set the vectors
//	int clx = 0;
//	int cly = 0;
//	int currentLoc = 0;
//	for (int i = 0; i < ytSize; i++)
//	{
//		currentLoc = i;
//		for (int j = 0; j < xtSize; j++)
//		{
//			xMeshVector[currentLoc] = xt[clx];
//			yMeshVector[currentLoc] = yt[cly];
//			currentLoc = currentLoc + ytSize;
//			clx += 1;
//		}
//		clx = 0;
//		cly += 1;
//	}
//
//	//D.  Now put data in the variables
//	mb_nx = xtSize;
//	mb_ny = ytSize;
//	mb_z = (float*)calloc((mb_nx*mb_ny), sizeof(float));
//	mb_xyz = (float*)calloc((3*(mb_n)), sizeof(float));
//	mb_zpij = (float*)calloc(mb_n, sizeof(float));
//	mb_knxt = (int*)calloc(mb_n, sizeof(int));
//	mb_imnew = (int*)calloc((mb_nx+mb_ny), sizeof(int));
//
//	memset((char *)mb_z,0,mb_nx*mb_ny*sizeof(float));
//
//	//E.  Initialize the input data structure
//	j = 0;
//	for (int i = 0; i < mb_n; i++){
//		mb_xyz[j++] = (float)(*x).at(i);
//		mb_xyz[j++] = (float)(*y).at(i);
//		mb_xyz[j++] = (float)(*z).at(i);
//	}
//
//	//************************************************************************************
//	//II. Call MB ZGrid
//	//************************************************************************************
//	mb_zgrid(mb_z, &mb_nx, &mb_ny, &mb_x1, &mb_y1, &mb_dx, &mb_dy, mb_xyz, &mb_n, mb_zpij, mb_knxt, mb_imnew, &mb_cay, &mb_nrng);
//
//	dgrid zGrid_temp = dgrid(mb_ny, mb_nx);
//
//	//A. Re align the data
//	k = 0;
//	for (i = 0; i < mb_ny; i++){
//		for (j = 0; j < mb_nx; j++){
//			zGrid_temp(i,j) = mb_z[k];
//			k++;
//		}
//	}
//	k = 0;
//	vector<double> zVector(xtSize*ytSize);
//	vector<double> eVector(xtSize*ytSize, 0);
//	for (i = 0; i < (const int)zGrid_temp.cols(); i++){
//		for (j = 0; j < (const int)zGrid_temp.rows(); j++){
//			zVector[k] = zGrid_temp(j,i);
//			k += 1;
//		}
//	}
//
//	stopTime = clock();
//	compTime = stopTime-startTime;
//	cout << "Time to Complete MB_ZGrid Interpolation: " << (stopTime-startTime)/CLOCKS_PER_SEC << endl << endl;
//
//	//************************************************************************************
//	//III. Convert back to Lon/Lat and output the file
//	//************************************************************************************
//	vector<double> lonSurf = vector<double>((ytSize*xtSize));
//	vector<double> latSurf = vector<double>((ytSize*xtSize));
//	//Print the results of MB_ZGrid
//	outFile.open(z_OutputFileName.c_str());
//	outFile.precision(6);
//	outFile.setf(std::ios::fixed, std::ios::floatfield);
//	if (outFile.is_open())
//	{
//		for(int i = 0; i < (ytSize*xtSize); i++){
//			UTMEasting = (xMeshVector[i])*cos(deg2rad*(-rotationAngle)) - (yMeshVector[i])*sin(deg2rad*(-rotationAngle));
//			UTMNorthing = (xMeshVector[i])*sin(deg2rad*(-rotationAngle)) + (yMeshVector[i])*cos(deg2rad*(-rotationAngle));
//
//			//This is saying we have no rotation angle!  
//			//Use this when we don't want to use the rotation angle or the rotation angle part isn't working correctly!
//			//This is putting everything back relative to the reference by just adding the reference value.
//			UTMNorthing = yMeshVector[i] + UTMNorthingRef;//SJZ THIS DOESN'T MAKE SENSE!!! WHY CALCULATED BEFORE
//			UTMEasting = xMeshVector[i] + UTMEastingRef;
//
//			UTMtoLL(refEllipsoid, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
//			outFile << newLon << "\t" << newLat << "\t" << zVector[i] << endl;
//			lonSurf[i] = newLon;
//			latSurf[i] = newLat;
//		}
//	}else
//	{
//		cerr << "Failed to open MB_ZGrid output file" << endl;
//		return false;
//	}
//	outFile.close();
//
//	//************************************************************************************
//	//IV. Compute the uncertainty associated with each data point
//	//************************************************************************************
//	if (usage < 0){
//		double minSpacing = spacingX;
//		if (spacingY < spacingX)
//			minSpacing = spacingY;
//		estimate(x, y, z, h, v, &xMeshVector, &yMeshVector, &zVector, &eVector, 1.96, 2.00, minSpacing);
//
//		stopTime = clock();
//		cout << "Time to Complete MB_ZGrid Interpolation with Errors: " << (compTime+(stopTime-startTime))/CLOCKS_PER_SEC << endl << endl;
//
//		//A. Output the GMT Surface depth calculation and the uncertainty computed by the uncertainty estimator
//		z_OutputFileName.append("_xyde.txt");
//		outFile.open(z_OutputFileName.c_str());
//		outFile.precision(6);
//		outFile.setf(std::ios::fixed, std::ios::floatfield);
//		if (outFile.is_open())
//		{
//			for(int i = 0; i < (const int)xMeshVector.size(); i++){
//				outFile << lonSurf[i] << "\t" << latSurf[i] << "\t" << zVector[i] << "\t" << eVector[i] << endl;
//			}
//			//exit (1);
//		}else
//		{
//			cerr << "Failed to open GMT Surface uncertainty output file" << endl;
//			//exit (1);
//			return false;
//		}
//	}
//
//	//************************************************************************************
//	//V. If we are going to reuse this output as input for mergeBathy then assign the variables here
//	//************************************************************************************
//	if ((usage == 2) || (usage == -2)) //Use as input
//	{
//		(*x).clear();
//		(*y).clear();
//		(*z).clear();
//		(*e).clear();
//		(*x) = vector<double>(xtSize*ytSize);
//		(*y) = vector<double>(xtSize*ytSize);
//		(*z) = vector<double>(xtSize*ytSize);
//		(*e) = vector<double>(xtSize*ytSize);
//
//		for(int i = 0; i < (xtSize*ytSize); i++){
//			(*x)[i] = xMeshVector[i];
//			(*y)[i] = yMeshVector[i];
//			(*z)[i] = zVector[i];
//			(*e)[i] = eVector[i];
//		}
//	}
//
//	xMeshVector.clear();
//	yMeshVector.clear();
//	zVector.clear();
//	eVector.clear();
//	xt.clear();
//	yt.clear();
//	lonSurf.clear();
//	latSurf.clear();
//
//	free(mb_z);
//	free(mb_xyz);
//	free(mb_zpij);
//	free(mb_knxt);
//	free(mb_imnew);
//
//	mb_z = NULL;
//	mb_xyz = NULL;
//	mb_zpij = NULL;
//	mb_knxt = NULL;
//	mb_imnew = NULL;
//
//	return true;
//}

////************************************************************************************
//// SUBROUTINE III: Function call for running GMT Surface
////************************************************************************************
//bool externalInterpolators::run_Surface_ORIGINAL(vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *h, vector<double> *v, double spacingX, double spacingY, double tension, string z_OutputFileName, double scaleFactor, double alpha, int usage)
//{
//	//************************************************************************************
//	// 0. Declare local variables and objects
//	//************************************************************************************
//	ofstream outFile;
//
//	double UTMNorthing;
//	double UTMEasting;
//	double newLat, newLon;
//	double minSpacing = spacingX;
//	double startTime, stopTime, compTime;
//
//	double *xConverted = NULL;
//	double *yConverted = NULL;
//	double *zConverted = NULL;
//
//	double *xPostSurface = NULL;
//	double *yPostSurface = NULL;
//	double *zPostSurface = NULL;
//
//	vector<double> xSurf;
//	vector<double> ySurf;
//	vector<double> lonSurf;
//	vector<double> latSurf;
//	vector<double> zSurf;
//	vector<double> eSurf;
//
//	double z0 = (*z)[0];
//	double z1 = (*z)[0];
//	double x0 = (*x)[0];
//	double x1 = (*x)[0];
//	double y0 = (*y)[0];
//	double y1 = (*y)[0];
//	int i;//, j, k;
//	int postSurfaceSize;
//	int returnValue;
//
//	startTime = clock();
//	if (spacingY < spacingX)
//		minSpacing = spacingY;
//
//	xConverted = (double*)calloc(((*x).size()), sizeof(double));
//	yConverted = (double*)calloc(((*y).size()), sizeof(double));
//	zConverted = (double*)calloc(((*z).size()), sizeof(double));
//
//	//A. Calculate min and max depths
//	for (i = 0; i < (const int)(*x).size(); i++)
//	{
//		if ((*z)[i] < z0)
//			z0 = (*z)[i];
//		else if ((*z)[i] > z1)
//			z1 = (*z)[i];
//
//		if ((*x)[i] < x0)
//			x0 = (*x)[i];
//		else if ((*x)[i] > x1)
//			x1 = (*x)[i];
//
//		if ((*y)[i] < y0)
//			y0 = (*y)[i];
//		else if ((*y)[i] > y1)
//			y1 = (*y)[i];
//
//		xConverted[i] = (*x)[i];
//		yConverted[i] = (*y)[i];
//		zConverted[i] = (*z)[i];
//	}
//
//	//************************************************************************************
//	//I. Call GMT Surface
//	//************************************************************************************
//	returnValue = processSurface(xConverted, yConverted, zConverted, (*x).size(), x0, y0, z0, x1, y1, z1, spacingX, spacingY, tension, &xPostSurface, &yPostSurface, &zPostSurface, &postSurfaceSize);
//
//	stopTime = clock();
//	compTime = stopTime-startTime;
//	cout << "Time to Complete GMT Surface Interpolation: " << (stopTime-startTime)/CLOCKS_PER_SEC << " seconds." << endl << endl;
//
//	xSurf = vector<double>(postSurfaceSize);
//	ySurf = vector<double>(postSurfaceSize);
//	lonSurf = vector<double>(postSurfaceSize);
//	latSurf = vector<double>(postSurfaceSize);
//	zSurf = vector<double>(postSurfaceSize);
//	eSurf = vector<double>(postSurfaceSize);
//
//	//************************************************************************************
//	//II. Convert back to Lon/Lat and output the file
//	//************************************************************************************
//	outFile.open(z_OutputFileName.c_str());
//	outFile.precision(6);
//	outFile.setf(std::ios::fixed, std::ios::floatfield);
//	if (outFile.is_open())
//	{
//		for(int i = 0; i < (postSurfaceSize); i++){
//			xSurf[i] = xPostSurface[i];
//			ySurf[i] = yPostSurface[i];
//			zSurf[i] = zPostSurface[i];
//			eSurf[i] = 0.00;
//
//			//This is a rotation matrix
//			UTMEasting = (xSurf[i])*cos(deg2rad*(-rotationAngle)) - (ySurf[i])*sin(deg2rad*(-rotationAngle));
//			UTMNorthing = (xSurf[i])*sin(deg2rad*(-rotationAngle)) + (ySurf[i])*cos(deg2rad*(-rotationAngle));
//
//			//This is saying we have no rotation angle!  
//			//Use this when we don't want to use the rotation angle or the rotation angle part isn't working correctly!
//			//This is putting everything back relative to the reference by just adding the reference value.
//			UTMNorthing = ySurf[i] + UTMNorthingRef; //SJZ THIS DOESN'T MAKE SENSE!!!! WHY CALCULATE IT BEFORE
//			UTMEasting = xSurf[i] + UTMEastingRef;
//
//			UTMtoLL(refEllipsoid, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
//			outFile << newLon << "\t" << newLat << "\t" << zSurf[i] << endl;
//			lonSurf[i] = newLon;
//			latSurf[i] = newLat;
//		}
//	}else
//	{
//		cerr << "Failed to open GMT Surface output file" << endl;
//		return false;
//	}
//	outFile.close();
//
//	startTime = clock();
//
//	//************************************************************************************
//	//III. Compute the uncertainty associated with each data point
//	//************************************************************************************
//	if (usage < 0){
////		int inputDebug;
//		estimate(x, y, z, h, v, &xSurf, &ySurf, &zSurf, &eSurf, scaleFactor, alpha, minSpacing);
//
//		stopTime = clock();
//		cout << "Time to Complete GMT Surface Interpolation with Errors: " << (compTime+(stopTime-startTime))/CLOCKS_PER_SEC << endl << endl;
//
//		//A.  Output the GMT Surface depth calculation and the uncertainty computed by the uncertainty estimator
//		z_OutputFileName.append("_xyde.txt");
//		outFile.open(z_OutputFileName.c_str());
//		outFile.precision(6);
//		outFile.setf(std::ios::fixed, std::ios::floatfield);
//		if (outFile.is_open())
//		{
//			for(int i = 0; i < (const int)xSurf.size(); i++){
//				outFile << lonSurf[i] << "\t" << latSurf[i] << "\t" << zSurf[i] << "\t" << eSurf[i] << endl;
//			}
//			//exit (1);
//		}else
//		{
//			cerr << "Failed to open GMT Surface uncertainty output file" << endl;
//			//exit (1);
//			return false;
//		}
//	}
//
//	//************************************************************************************
//	//IV. If we are going to reuse this output as input for mergeBathy then assign the variables here
//	//************************************************************************************
//	if ((usage == 2) || (usage == -2)) //Use as input
//	{
//		(*x).clear();
//		(*y).clear();
//		(*z).clear();
//		(*e).clear();
//		(*x) = vector<double>(postSurfaceSize);
//		(*y) = vector<double>(postSurfaceSize);
//		(*z) = vector<double>(postSurfaceSize);
//		(*e) = vector<double>(postSurfaceSize);
//
//		for(int i = 0; i < postSurfaceSize; i++){
//			(*x)[i] = xSurf[i];
//			(*y)[i] = ySurf[i];
//			(*z)[i] = zSurf[i];
//			(*e)[i] = eSurf[i];
//		}
//	}
//
//	//A. Clear up variables
//	free(xConverted);
//	free(yConverted);
//	free(zConverted);
//
//	free(xPostSurface);
//	free(yPostSurface);
//	free(zPostSurface);
//
//	xConverted = NULL;
//	yConverted = NULL;
//	zConverted = NULL;
//
//	xPostSurface = NULL;
//	yPostSurface = NULL;
//	zPostSurface = NULL;
//
//	xSurf.clear();
//	ySurf.clear();
//	lonSurf.clear();
//	latSurf.clear();
//	zSurf.clear();
//	eSurf.clear();
//
//	return true;
//}
