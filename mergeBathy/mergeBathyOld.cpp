
#include "mergeBathyOld.h"
#include "LatLong-UTMconversion.h"
#include "constants.h"
#include <fstream>
#include "externalInterpolators.h"
#include <algorithm>	
#include "standardOperations.h"
#include "computeOffset.h"
#include "rng.h"
#include "fileWriter.h"
#include "fileBagWriter.h"
#include <string.h> //UNIX

#include <sstream>
#include <math.h>
#include <cmath>
#include <functional> //mod
#include "Error_Estimator/Bathy_Grid.h"

#ifdef _MSC_VER
//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	//#include "../WarningStates.h"	//Disable all Warnings!!!
	#pragma warning ( disable : 4996 )	//Deprecated call for strcpy
#endif
#endif

//bool isEqual (int i) {
//  return ((i==agreeMax);
//}


//************************************************************************************
// SUBROUTINE I: Primary MergeBathy function call for
// preprocessing and grid alignment
//************************************************************************************
int mergeBathy_PreCompute(vector<TRUE_DATA> *inputData, double refLon, double refLat, double rotationAngle, double gridSpacingX, double gridSpacingY, double smoothingScaleX, double smoothingScaleY, string kernelName, string outputFileName, map<string, int> additionalOptions, int numMCRuns, MB_ZGRID_DATA *MB_ZGridInput, GMT_SURFACE_DATA *GMTSurfaceInput, ALG_SPLINE_DATA *ALGSplineInput, FORCED_LOCATIONS *forcedLocationPositions,int usagePreInterpLocsLatLon)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int returnValue = 0;
	int refEllipsoid = 23;
	int totalElements = 0;
	int currentLoc;

	double longitudeMean = 0.00;
	double maxLat, minLat;
	double UTMNorthingRef = 0, UTMEastingRef = 0, UTMNorthingRefAll = 0, UTMEastingRefAll = 0, UTMNorthing = 0, UTMEasting = 0;
	double x0, x1, y0, y1, lat0, lat1, lon0, lon1;
	int xtSize, ytSize;
	double locationValue;
	char UTMZone[4]="";//{'\0'};
	char UTMZoneRef[4]="";
	char UTMZoneRefAll[4]="";
	int USE_UTM = 0;
//	*UTMZone='\0';
	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> e;
	vector<double> hErr;
	vector<double> vErr;

	//A. Variables for computed locations
	double meanXt = 0.00;
	double meanYt = 0.00;
	vector<double> xt;
	vector<double> yt;
	vector<double> xMeshVector;
	vector<double> yMeshVector;
	dgrid xMeshGrid;
	dgrid yMeshGrid;

	OUTPUT_DATA xyzOut;
	double newX0;
	double newY0;
	double newX1;
	double newY1;

	//B. Wrap the rotation angle if greater than 360
	//Get rotation angle to 0-360 range
	//(useful if data must be transposed to
	//display in "mostly map-ish" orientation).
	rotationAngle = fmod(rotationAngle, 360);
	if(rotationAngle < 0) rotationAngle += 360;
	//else rotationAngle -= 360; // SJZ

	// Determine whether or not to use UTM grid zones
	if(((refLat == 0) & (refLon == 0)) & (rotationAngle == 0))
		USE_UTM = 1;

	double latitudeMean=0;
	double latitudeMeanAll=0;
	double longitudeMeanAll=0;
	//C. Get the mean location of the data for UTM2LL conversion in the end
	for(int count = 0; count < (int) (*inputData).size(); count++){
		totalElements += (int)(*inputData)[count].lon.size();
		longitudeMean += (*inputData)[count].longitudeSum;
		latitudeMeanAll += (*inputData)[count].latitudeSum;
	}
	longitudeMean /= double(totalElements);
	longitudeMeanAll = longitudeMean;
	latitudeMeanAll /= double(totalElements);

	x = vector<double>(totalElements);
	y = vector<double>(totalElements);
	z = vector<double>(totalElements);
	e = vector<double>(totalElements);
	hErr = vector<double>(totalElements);
	vErr = vector<double>(totalElements);

	//D. Compute offset so long as there multiple data sets and we did not specify otherwise
	additionalOptions.find("-computeOffset")->second = ((*inputData).size()>1 && additionalOptions.find("-computeOffset")->second == 1)?1:0;

	//************************************************************************************
	// I. If we have Lat/Lon input coordinates then convert them to UTM
	//************************************************************************************
	if (additionalOptions.find("-inputInMeters")->second == 0)
	{
		//Added 11/11/14 SJZ
		vector<double>::iterator minIt=(*inputData)[0].lon.begin();
		vector<double>::iterator maxIt=(*inputData)[0].lon.begin();
		vector<double>::iterator curMinIt;
		vector<double>::iterator curMaxIt;
		for(int count = 0; count < (const int)(*inputData).size(); count++)
		{
			curMinIt= min_element((*inputData)[count].lon.begin(),(*inputData)[count].lon.end());
			curMaxIt= max_element((*inputData)[count].lon.begin(),(*inputData)[count].lon.end());
			if(*curMinIt < *minIt)
				minIt = curMinIt;
			if(*curMaxIt > *maxIt)
				maxIt = curMaxIt;
		}
		double sumlon = 0;
		if(((((*minIt<0) & (*maxIt>0)) & (0-*minIt > 180+*minIt)) || (additionalOptions.find("-preInterpolatedLocations")->second == 1 
			&& (*max_element(forcedLocationPositions->forcedLonCoord.begin(), forcedLocationPositions->forcedLonCoord.end()) >=180 
			&& *min_element(forcedLocationPositions->forcedLonCoord.begin(),forcedLocationPositions->forcedLonCoord.end()) >= 0))) && abs(usagePreInterpLocsLatLon) != 2)
		{
			for(int count = 0; count < (const int)(*inputData).size(); count++)
			{
				(*inputData)[count].longitudeSum = 0;
				for(int i = 0; i < (const int)(*inputData)[count].lon.size(); i++)
				{
					if((*inputData)[count].lon[i] < 0)
						(*inputData)[count].lon[i] += 360;
					(*inputData)[count].longitudeSum += (*inputData)[count].lon[i];
				}
			}
			if(refLon < 0)
				refLon += 360;
			if(additionalOptions.find("-preInterpolatedLocations")->second == 1)
			{
				forcedLocationPositions->longitudeSum = 0;
				for(int i = 0; i < (const int)forcedLocationPositions->forcedLonCoord.size(); i++)
				{
					if(forcedLocationPositions->forcedLonCoord[i] < 0)
						forcedLocationPositions->forcedLonCoord[i] += 360;
					forcedLocationPositions->longitudeSum += forcedLocationPositions->forcedLonCoord[i];
				}
			}
		}


		#pragma region -- Input in Lat,Lon; Convert to UTM
		//A. Get the reference and the first coordinate converted

		//i. Get Reference Lat and Lon for all datasets
		//This was the original call for the UTMZoneRef used for LL2UTM conversion,
		// and UTMEastingRef and UTMNorthingRef but it did not match
		// matlab.
		//Therefore, it was kept for UTM2LL conversion as UTMZoneRefAll.
		LLtoUTM(refEllipsoid, refLat, refLon, UTMNorthingRefAll, UTMEastingRefAll, UTMZoneRefAll, refLon);

		//ii. Get individual dataset's longitude mean and reference Lat and Lon
		//Matlab code calculates UTMZoneRef and UTMZone based on
		//the longitudeMean of each individual dataset 
		//instead of the refLon, and longitudeMean of all datasets.
		longitudeMean = (*inputData)[0].longitudeSum/(*inputData)[0].lon.size();
		//Find UTMZoneRef for each dataset using its own longitudeMean.
		LLtoUTM(refEllipsoid, refLat, refLon, UTMNorthingRef, UTMEastingRef, UTMZoneRef, longitudeMean);

		//iii. Find UTMZone for first dataset from first point.
		//Convert the fist point to UTM. Get the UTM for the first data point in the dataset using the dataset's mean longitude.
		LLtoUTM(refEllipsoid, (*inputData)[0].lat[0], (*inputData)[0].lon[0], UTMNorthing, UTMEasting, UTMZone, longitudeMean);
		UTMEasting -= UTMEastingRef;
		UTMNorthing -= UTMNorthingRef;

		//B. Factor in the rotation angle to the UTM coordinates
		x0 = MAX_INT;
		x1 = MIN_INT;
		y0 = MAX_INT;
		y1 = MIN_INT;
		lon0 = (*inputData)[0].lon[0];
		lon1 = (*inputData)[0].lon[0];
		lat0 = (*inputData)[0].lat[0];
		lat1 = (*inputData)[0].lat[0];
		maxLat = lat0;
		minLat = lat0;
		int INIT_TEMPS = 1;
	
		//C. Compute each pair of UTM coordinates based
		//on the location of the reference coordinate 
		//calculated for each dataset from its longitude mean.
		currentLoc = 0;
		for(int count = 0; count < (const int)(*inputData).size(); count++)
		{
			//Find dataset's longitude mean
			longitudeMean = (*inputData)[count].longitudeSum/(*inputData)[count].lon.size();
			
			strcpy(UTMZone,"\0"); // SJZ
			latitudeMean = (*inputData)[count].latitudeSum/(*inputData)[count].lat.size(); // SJZ
			//Find the utmzone from the the dataset's mean latitude and mean longitude.
			LLtoUTM(refEllipsoid, latitudeMean, longitudeMean, UTMNorthing, UTMEasting, UTMZone, longitudeMean); // SJZ
			strcpy((*inputData)[count].utmzone, UTMZone); //utmzone for the majority of the dataset. // SJZ

			if(UTMZone==NULL)
				cerr<<"Did not expect this! The UTMZone came back null when calculating from the dataset's mean latitude and mean longitude"<<endl;
				//When this happens, we will instead calculate the utm zone for each individual data point. we did this by changing meanLong to LongTemp in LL2UTM
			strcpy(UTMZoneRef,UTMZone); // SJZ
			
			if(!USE_UTM) // SJZ
			{
				//strcpy(UTMZoneRef,UTMZone); // SJZ
				//Calculate dataset's UTMNorthingRef and UTMEastingRef 
				LLtoUTM(refEllipsoid, refLat, refLon, UTMNorthingRef, UTMEastingRef, UTMZoneRef, longitudeMean);
			}

			for(int i = 0; i < (const int)(*inputData)[count].lon.size(); i++)
			{
				//Find max and min latitude (isn't this the same as lat1 and lat0? SJZ)
				if ((*inputData)[count].lat[i] > maxLat)
					maxLat = (*inputData)[count].lat[i];
				else if ((*inputData)[count].lat[i] < minLat)
					minLat = (*inputData)[count].lat[i];
				
				LLtoUTM(refEllipsoid, (*inputData)[count].lat[i], (*inputData)[count].lon[i], UTMNorthing, UTMEasting, UTMZone, longitudeMean);
				if(!USE_UTM) // SJZ
				{
					UTMEasting -= UTMEastingRef;
					UTMNorthing -= UTMNorthingRef;

					//i. Factor in the rotation angle to the UTM coordinates
					x[currentLoc] = (UTMEasting)*cos(deg2rad*rotationAngle) - (UTMNorthing)*sin(deg2rad*rotationAngle);
					y[currentLoc] = (UTMEasting)*sin(deg2rad*rotationAngle) + (UTMNorthing)*cos(deg2rad*rotationAngle);
				}
				else
				{// SJZ
					x[currentLoc] = (UTMEasting);
					y[currentLoc] = (UTMNorthing);
				}
				//ii. Determine minimum and maximum extents of the grid for calculation later
				//if(INIT_TEMPS)
			/*	{
					x0 = x[currentLoc];
					x1 = x[currentLoc];
					y0 = y[currentLoc];
					y1 = y[currentLoc];
					INIT_TEMPS = 0;
				}*/

				if(x[currentLoc] < x0)
					x0 = x[currentLoc];
				/*else*/ if(x[currentLoc] > x1)
					x1 = x[currentLoc];

				if(y[currentLoc] < y0)
					y0 = y[currentLoc];
				/*else*/ if(y[currentLoc] > y1)
					y1 = y[currentLoc];

				if((*inputData)[count].lat[i] < lat0)
					lat0 =(*inputData)[count].lat[i];
				/*else*/ if((*inputData)[count].lat[i] > lat1)
					lat1 = (*inputData)[count].lat[i];

				if((*inputData)[count].lon[i] < lon0)
					lon0 = (*inputData)[count].lon[i];
				/*else*/ if((*inputData)[count].lon[i] > lon1)
					lon1 = (*inputData)[count].lon[i];

				z[currentLoc] = (*inputData)[count].depth[i];
				e[currentLoc] = (*inputData)[count].error[i];
				hErr[currentLoc] = (*inputData)[count].h_Error[i];
				vErr[currentLoc] = (*inputData)[count].v_Error[i];
				(*inputData)[count].x[i] = x[currentLoc];
				(*inputData)[count].y[i] = y[currentLoc];
				currentLoc+=1;
			}
			if(!USE_UTM)
				strcpy((*inputData)[count].utmzone, "999");
				

		}

		#pragma region Adjust for equator crossing or weird data set
		////D. If we have the special case where we have positive and negative latitudes (i.e. we cross the equator or have a weird data set)
		//if ((minLat < 0) && (maxLat > 0)){
		//	if(refLat >= 0){
		//		int loc = 0;
		//		if ((*inputData)[0].lat[0] < 0)
		//		{
		//			y0 = y[0] - 10000000.00;
		//			y1 = y[0] - 10000000.00;
		//		}else
		//		{
		//			y0 = y[0];
		//			y1 = y[0];
		//		}
		//		for(int count = 0; count < (const int)(*inputData).size(); count++){
		//			for(int i = 0; i < (const int)(*inputData)[count].lat.size(); i++){
		//				if ((*inputData)[count].lat[i] < 0){
		//					y[loc] = y[loc] - 10000000.00;
		//					(*inputData)[count].y[i] = y[loc];
		//					if(y[loc] < y0)
		//						y0 = y[loc];
		//					else if(y[loc] > y1)
		//						y1 = y[loc];
		//				}
		//				loc += 1;
		//			}
		//		}
		//	}else{
		//		int loc = 0;
		//		if ((*inputData)[0].lat[0] >= 0)
		//		{
		//			y0 = y[0] + 10000000.00;
		//			y1 = y[0] + 10000000.00;
		//		}else
		//		{
		//			y0 = y[0];
		//			y1 = y[0];
		//		}
		//		for(int count = 0; count < (const int)(*inputData).size(); count++){
		//			for(int i = 0; i < (const int)(*inputData)[count].lat.size(); i++){
		//				if ((*inputData)[count].lat[i] >= 0){
		//					y[loc] = y[loc] + 10000000.00;
		//					(*inputData)[count].y[i] = y[loc];
		//					if(y[loc] < y0)
		//						y0 = y[loc];
		//					else if(y[loc] > y1)
		//						y1 = y[loc];
		//				}
		//				loc += 1;
		//			}
		//		}
		//	}
		//}
		#pragma endregion
		cout << "Input Data UTM Zone: " << UTMZoneRef << endl;
		#pragma endregion
	}
	else
	{
		#pragma region -- Input already in UTM
		//E. We aren't using Lat/Lon positions, instead we read in X and Y positions that were already in UTM
		x0 = (*inputData)[0].lon[0];
		x1 = (*inputData)[0].lon[0];
		y0 = (*inputData)[0].lat[0];
		y1 = (*inputData)[0].lat[0];
		lon0 = (*inputData)[0].lon[0];
		lon1 = (*inputData)[0].lon[0];
		lat0 = (*inputData)[0].lat[0];
		lat1 = (*inputData)[0].lat[0];

		currentLoc = 0;
		for(int count = 0; count < (const int) (*inputData).size(); count++){
			for(int i = 0; i < (const int)(*inputData)[count].lon.size(); i++){
				//i. Determine minimum and maximum extents of the grid for calculation later
				if((*inputData)[count].lon[i] < x0)
				{
					x0 = (*inputData)[count].lon[i];
					lon0 = (*inputData)[count].lon[i];
				}else if((*inputData)[count].lon[i] > x1)
				{
					x1 = (*inputData)[count].lon[i];
					lon1 = (*inputData)[count].lon[i];
				}

				if((*inputData)[count].lat[i] < y0)
				{
					y0 = (*inputData)[count].lat[i];
					lat0 =(*inputData)[count].lat[i];
				}else if((*inputData)[count].lat[i] > y1)
				{
					y1 = (*inputData)[count].lat[i];
					lat1 = (*inputData)[count].lat[i];
				}

				x[currentLoc] = (*inputData)[count].lon[i];
				y[currentLoc] = (*inputData)[count].lat[i];
				z[currentLoc] = (*inputData)[count].depth[i];
				e[currentLoc] = (*inputData)[count].error[i];
				hErr[currentLoc] = (*inputData)[count].h_Error[i];
				vErr[currentLoc] = (*inputData)[count].v_Error[i];
				(*inputData)[count].x[i] = x[currentLoc];
				(*inputData)[count].y[i] = y[currentLoc];
				currentLoc+=1;
			}
		}

		if(!USE_UTM)//F. Do this for safety that way we don't end up using the Reference zone info when it has no data.
			LLtoUTM(refEllipsoid, refLat, refLon, UTMNorthingRef, UTMEastingRef, UTMZoneRef, refLon);
		#pragma endregion
	}





	
	int refellip;
//	char zoneAll;

	if(USE_UTM)
	{
		
		//strcpy(UTMZoneRefAll,'\0'); // SJZ
		//LLtoUTM(refEllipsoid, latitudeMeanAll, longitudeMeanAll, UTMNorthingRefAll, UTMEastingRefAll, UTMZoneRefAll, longitudeMeanAll);

		// Set scene center lat/lon for use with UTM gridding
		// Not necessary if local coordinate system specified)
		double clat = latitudeMeanAll;
		double clon = longitudeMeanAll;

		vector<char*> zone;
		char* zoneTemp=(*inputData)[0].utmzone;
		strcpy(UTMZoneRefAll, zoneTemp);
		int REPROJECT = 0;
		// Make sure all datasets are referenced to the same UTM zone
		for(int i = 0; i < (const int)(*inputData).size(); i++)
		{
			zone.push_back((*inputData)[i].utmzone); 
			if(strcmp(zoneTemp, zone[i])!=0)
				REPROJECT = 1;
		}
		// find mismatching zones
		if(REPROJECT)
		{
			cout << "UTM zone not consistent between datasets.  Please re-run mergeBathy with a non-zero reference position/rotation angle." << endl;
	
			// new code:
			// Changed 10/24/14 to re-project each "errant" dataset into
			// the "majority" zone.  Force the zones on the datasets. 
			// Assumes datasets are of the same order of magnitude so we
			// can simply compare UTM zones per dataset to find the
			// majority.  SJZ
	
			char zonestring[512] = "";
			refellip = 23;   // this number corresponds to WGS-84
			for(int i = 0; i < (const int)(*inputData).size(); i++)
			{
				strcat(zonestring,(*inputData)[i].utmzone);
			}
			vector<int> agree;
			size_t found;
			string zonestr=zonestring;
			int cnt, cntMax=0, idMax=0; 
			char* valueMax;
			for(int i = 0; i < (const int)(*inputData).size(); i++)
			{
				cnt = 0;
				found = zonestr.find((*inputData)[i].utmzone);
				while(found != std::string::npos)
				{
					cnt++;
					found = zonestr.find((*inputData)[i].utmzone,found+1);
				}
				if(cnt > cntMax)
				{
					cntMax = cnt;
					idMax = i;
					valueMax = (*inputData)[i].utmzone;
				}
				agree.push_back(cnt);
			}
			zoneTemp = valueMax;
			int agreeindex = idMax;
			strcpy(UTMZoneRefAll, zoneTemp);
		//	vector<int>::iterator agreeMaxIt = max_element(agree.begin(),agree.end());
			longitudeMean = (*inputData)[agreeindex].longitudeSum/(*inputData)[agreeindex].lon.size();
			x0 = MAX_INT;
			x1 = MIN_INT;
			y0 = MAX_INT;
			y1 = MIN_INT;
			double x0Temp = MAX_INT;
			double x1Temp = MIN_INT;
			double y0Temp = MAX_INT;
			double y1Temp = MIN_INT;
			int INIT_TEMPS = 1;
			int currentLoc = 0;
			for(int count = 0; count < (const int)(*inputData).size(); count++)
			{
				if (strcmp(zoneTemp,(*inputData)[count].utmzone)!=0)
				{
					cout << "reprojecting UTM coordinates into majority zone." << endl;
					strcpy((*inputData)[count].utmzone, zoneTemp);

					for(int i = 0; i < (const int)(*inputData)[count].lon.size(); i++)
					{			
						LLtoUTM(refEllipsoid, (*inputData)[count].lat[i], (*inputData)[count].lon[i], UTMNorthing, UTMEasting, zoneTemp, longitudeMean);
						x[currentLoc] = (UTMEasting);
						y[currentLoc] = (UTMNorthing);
						(*inputData)[count].x[i] = x[currentLoc];
						(*inputData)[count].y[i] = y[currentLoc];
						
						//ii. Determine minimum and maximum extents of the grid for calculation later
						//Set initial min, max
						if(INIT_TEMPS)
						{
							x0 = x[currentLoc];
							x1 = x[currentLoc];
							y0 = y[currentLoc];
							y1 = y[currentLoc];
							INIT_TEMPS=0;
						}

						//Set min, max
						if(x[currentLoc] < x0)
							x0 = x[currentLoc];
						if(x[currentLoc] > x1)
							x1 = x[currentLoc];
						if(y[currentLoc] < y0)
							y0 = y[currentLoc];
						if(y[currentLoc] > y1)
							y1 = y[currentLoc];

						currentLoc++;
					}
				}
				else
				{
					currentLoc += (int)(*inputData)[count].x.size();
					x0Temp = *min_element((*inputData)[count].x.begin(),(*inputData)[count].x.end());
					x1Temp = *max_element((*inputData)[count].x.begin(),(*inputData)[count].x.end());
					y0Temp = *min_element((*inputData)[count].y.begin(),(*inputData)[count].y.end());
					y1Temp = *max_element((*inputData)[count].y.begin(),(*inputData)[count].y.end());

					//ii. Determine minimum and maximum extents of the grid for calculation later
					//Set initial min, max
					if(INIT_TEMPS)
					{
						x0 = x0Temp;
						x1 = x1Temp;
						y0 = y0Temp;
						y1 = y1Temp;
						INIT_TEMPS = 0;
					}

					//Set min, max
					if(x0Temp < x0)
						x0 = x0Temp;
					if(x1Temp > x1)
						x1 = x1Temp;
					if(y0Temp < y0)
						y0 = y0Temp;
					if(y1Temp > y1)
						y1 = y1Temp;
				}
			}
		}
	}
	//else
	//	strcpy(UTMZoneRefAll, "999");//zoneAll = '999';

	x0 = floor(x0);
	x1 = ceil(x1);
	y0 = floor(y0);
	y1 = ceil(y1);

	//************************************************************************************
	// II. If we are computing our own locations based on a smoothing scale then do those calculations here
	//************************************************************************************
	if (additionalOptions.find("-preInterpolatedLocations")->second == 0)//calcgrid
	{
		#pragma region --No preInterpolatedLocations provided
		//Compute grid locations based on smoothing scale; No preInterpolatedLocations provided
		if (gridSpacingX == 0 && gridSpacingX == 0)
        {
			gridSpacingX = 100;
			gridSpacingY = 100;
		}
		if(smoothingScaleX == 0 && smoothingScaleY == 0)
		{
			smoothingScaleX = 100;
			smoothingScaleY = 100;
		}
		else if (smoothingScaleX < gridSpacingX || smoothingScaleY < gridSpacingY)
		{
			cout << "Grid spacing defaults to 100 making smoothing scales too small.  Make smoothing scales greater than the grid spacing." << endl;
			return ARGS_ERROR;
		}

		meanXt = 0.00;
		meanYt = 0.00;

		//A. Calculate X
		locationValue = x0 - gridSpacingX;
		xtSize = 0;
		while (locationValue <= (const double)(x1 + gridSpacingX)) {
			xt.push_back(locationValue);
			meanXt += locationValue;
			locationValue += gridSpacingX;
			xtSize += 1;
		}
		//Get the mean for later use
		meanXt = meanXt / (double)xt.size();

		//B. Calculate Y
		locationValue = y0 - gridSpacingY;
		ytSize = 0;
		while (locationValue <= (const double)(y1 + gridSpacingY)) {
			yt.push_back(locationValue);
			meanYt += locationValue;
			locationValue += gridSpacingY;
			ytSize += 1;
		}
		//Get the mean for later use
		meanYt = meanYt / (double)yt.size();

		int temp = xtSize*ytSize;
		xMeshVector = vector<double>(temp);
		yMeshVector = vector<double>(temp);
		xMeshGrid = dgrid(ytSize, xtSize);
		yMeshGrid = dgrid(ytSize, xtSize);

		//C. Reform the data to a grid that we can use for calculations
		int clx = 0;
		int cly = 0;
		currentLoc = 0;
		for (int i = 0; i < ytSize; i++)
		{
			currentLoc = i;
			for (int j = 0; j < xtSize; j++)
			{
				xMeshGrid(i,j) = xt[clx] - meanXt;
				yMeshGrid(i,j) = yt[cly] - meanYt;
				xMeshVector[currentLoc] = xt[clx];
				yMeshVector[currentLoc] = yt[cly];
				currentLoc = currentLoc + ytSize;
				clx += 1;
			}
			clx = 0;
			cly += 1;
		}
		#pragma endregion
	}
	else
	{
		#pragma region --Dont compute grid, Use pre-interpolated locs
		//D. Don't compute a grid, instead use the pre-interpolated locations
		xMeshVector = vector<double>((*forcedLocationPositions).forcedLonCoord.size());
		yMeshVector = vector<double>((*forcedLocationPositions).forcedLatCoord.size());
		
		//E. Convert the pre-interpolated lon-lat locations to meters in x and y
		minLat = (*forcedLocationPositions).forcedLatCoord[0];
		maxLat = (*forcedLocationPositions).forcedLatCoord[0];
		strcpy(UTMZone,"\0"); // SJZ
		longitudeMean = (*forcedLocationPositions).longitudeSum/(*forcedLocationPositions).forcedLonCoord.size();
		latitudeMean = (*forcedLocationPositions).latitudeSum/(*forcedLocationPositions).forcedLatCoord.size(); // SJZ
		LLtoUTM(refEllipsoid, latitudeMean, longitudeMean, UTMNorthing, UTMEasting, UTMZone, longitudeMean); // SJZ
		strcpy(UTMZoneRef,UTMZone); // SJZ

		if(!USE_UTM) // SJZ 
			LLtoUTM(refEllipsoid, refLat, refLon, UTMNorthingRef, UTMEastingRef, UTMZoneRef, longitudeMean);
		
		for(int i = 0; i < (const int)(*forcedLocationPositions).forcedLonCoord.size(); i++)
		{
			if ((*forcedLocationPositions).forcedLatCoord[i] > maxLat)
				maxLat = (*forcedLocationPositions).forcedLatCoord[i];

			if ((*forcedLocationPositions).forcedLatCoord[i] < minLat)
				minLat = (*forcedLocationPositions).forcedLatCoord[i];
			if(abs(usagePreInterpLocsLatLon) == 1)
			{
				LLtoUTM(refEllipsoid, (*forcedLocationPositions).forcedLatCoord[i], (*forcedLocationPositions).forcedLonCoord[i], UTMNorthing, UTMEasting, UTMZone, longitudeMean);
				if(!USE_UTM)// SJZ 
				{
					UTMEasting -= UTMEastingRef;
					UTMNorthing -= UTMNorthingRef;

					//i. Factor in the rotation angle to the UTM coordinates
					xMeshVector[i] = (UTMEasting)*cos(deg2rad*rotationAngle) - (UTMNorthing)*sin(deg2rad*rotationAngle);
					yMeshVector[i] = (UTMEasting)*sin(deg2rad*rotationAngle) + (UTMNorthing)*cos(deg2rad*rotationAngle);
				}
				else
				{// SJZ 
					xMeshVector[i] = (UTMEasting);
					yMeshVector[i] = (UTMNorthing);
				}
			}
			else
			{
				xMeshVector[i] = (*forcedLocationPositions).forcedLonCoord[i];
				yMeshVector[i] = (*forcedLocationPositions).forcedLatCoord[i];
			}
			strcpy((*forcedLocationPositions).utmzone, UTMZone);
			if(!USE_UTM)// SJZ 
				strcpy((*forcedLocationPositions).utmzone, "999");


			//Commented out 11/11/14 SJZ 
			//F. If we have the special case where we have positive and negative latitudes (i.e. we cross the equator or have a weird data set)
			/*if ((minLat < 0) && (maxLat > 0)){
				if(refLat >= 0){
					for(int i = 0; i < (const int)(*forcedLocationPositions).forcedLatCoord.size(); i++){
						if ((*forcedLocationPositions).forcedLatCoord[i] < 0){
							yMeshVector[i] = yMeshVector[i] - 10000000.00;
						}
					}
				}else{
					for(int i = 0; i < (const int)(*forcedLocationPositions).forcedLatCoord.size(); i++){
						if ((*forcedLocationPositions).forcedLatCoord[i] >= 0){
							yMeshVector[i] = yMeshVector[i] + 10000000.00;
						}
					}
				}
			}	*/	
		}
		
		//G. Get the mean data points for later use
		meanXt = 0.00;
		meanYt = 0.00;
		for (int i = 0; i < (const int)xMeshVector.size(); i++)
		{
			meanXt = meanXt + xMeshVector[i];
			meanYt = meanYt + yMeshVector[i];
		}
		meanXt = meanXt / (double)xMeshVector.size();
		meanYt = meanYt / (double)yMeshVector.size();
		#pragma endregion
	}

	//************************************************************************************
	// III. Compute Offset
	//************************************************************************************
	if (additionalOptions.find("-computeOffset")->second == 1)
	{
		cout << "\nComputing Offset between Data Sets" << endl;
		computeOffset(inputData, x0, y0, x1, y1, &xMeshVector, &yMeshVector, smoothingScaleX*2.00, smoothingScaleY*2.00, &z, &e, &hErr, &vErr, kernelName, additionalOptions);
		cout << "Done Computing Offset between Data Sets" << endl << endl;
	}
	//Done computing offset

	newX0 = x0;
	newY0 = y0;
	newX1 = x1;
	newY1 = y1;
	cout << "Dimensions of Computational Area: " << endl;
	cout << "\tRows: " << xMeshGrid.rows() << "\n\tCols: " << xMeshGrid.cols() << endl << endl;

	// Perform the actual gridding of the data. Do a small amount
	// of truncation (to mm precision) to allow for imprecise grids
	double dRoundoff = 1e-3;
	for(int i = 0; i < (const int)xMeshVector.size(); i++)
	{
		xMeshVector[i] = roundDouble(xMeshVector[i]/dRoundoff)*dRoundoff;
		yMeshVector[i] = roundDouble(yMeshVector[i]/dRoundoff)*dRoundoff;
	}

	Bathy_Grid* bathyGrid = new Bathy_Grid();
	#pragma region -- Gridding
	if(additionalOptions.find("-ZGrid")->second == 1 || additionalOptions.find("-GMTSurface")->second == 1)
	{
		bathyGrid->GriddingFlag = 1;
		bathyGrid->Construct_Tin(&x, &y, &z, &hErr, &vErr);
	} else bathyGrid->GriddingFlag = 0;
	#pragma endregion

	//************************************************************************************
	// IV. Call actual processing routines
	//************************************************************************************
	returnValue = run(&x, &y, &z, &e, &hErr, &vErr, &xMeshGrid, &yMeshGrid, &xt, &yt, &xMeshVector, &yMeshVector, gridSpacingX, gridSpacingY, smoothingScaleX, smoothingScaleY, newX0, newY0, newX1, newY1, meanXt, meanYt, kernelName,  additionalOptions, outputFileName, 1.00, UTMNorthingRefAll, UTMEastingRefAll, rotationAngle, refEllipsoid, UTMZoneRefAll, numMCRuns, MB_ZGridInput, GMTSurfaceInput, ALGSplineInput, bathyGrid, USE_UTM);

	//if (numMCRuns == -1)
	//{
	//	//A. Call single run
	//	returnValue = runSingle(&x, &y, &z, &e, &hErr, &vErr, &xMeshGrid, &yMeshGrid, &xt, &yt, &xMeshVector, &yMeshVector, gridSpacingX, gridSpacingY, smoothingScaleX, smoothingScaleY, newX0, newY0, newX1, newY1, meanXt, meanYt, kernelName,  additionalOptions, outputFileName, 1.00, UTMNorthingRefAll, UTMEastingRefAll, rotationAngle, refEllipsoid, UTMZoneRefAll, MB_ZGridInput, GMTSurfaceInput, ALGSplineInput, bathyGrid, USE_UTM);
	//}else
	//{
	//	//B. Call Monte Carlo run
	//	//returnValue = runMonteCarlo(&x, &y, &z, &e, &hErr, &vErr, &xMeshGrid, &yMeshGrid, &xt, &yt, &xMeshVector, &yMeshVector, gridSpacingX, gridSpacingY, smoothingScaleX, smoothingScaleY, newX0, newY0, newX1, newY1, meanXt, meanYt, kernelName,  additionalOptions, outputFileName, 1.00, UTMNorthingRef, UTMEastingRef, rotationAngle, refEllipsoid, UTMZoneRef, numMCRuns, MB_ZGridInput, GMTSurfaceInput);
	//	returnValue = runMonteCarlo(&x, &y, &z, &e, &hErr, &vErr, &xMeshGrid, &yMeshGrid, &xt, &yt, &xMeshVector, &yMeshVector, gridSpacingX, gridSpacingY, smoothingScaleX, smoothingScaleY, newX0, newY0, newX1, newY1, meanXt, meanYt, kernelName,  additionalOptions, outputFileName, 1.00, UTMNorthingRef, UTMEastingRef, rotationAngle, refEllipsoid, UTMZoneRef, numMCRuns, MB_ZGridInput, GMTSurfaceInput, ALGSplineInput, bathyGrid);
	//}

	//************************************************************************************
	// V. Clear up variables
	//************************************************************************************
	if(bathyGrid->GriddingFlag)
		delete bathyGrid;
	xMeshVector.clear();
	yMeshVector.clear();
	xt.clear();
	yt.clear();
	xMeshGrid.clear();
	yMeshGrid.clear();
	x.clear();
	y.clear();
	z.clear();
	e.clear();
	hErr.clear();
	vErr.clear();
	xyzOut.depth.clear();
	xyzOut.error.clear();
	xyzOut.nEi.clear();
	xyzOut.rEi.clear();
	//if (numMCRuns == -1)
	{
		xyzOut.depth0.clear();
		xyzOut.depthK.clear();
		xyzOut.errorK.clear();
		xyzOut.standardDev.clear();
	//	xyzOut.standardDev0.clear();
	}
	return returnValue;
}

//************************************************************************************
// SUBROUTINE II: Primary MergeBathy function call for processing single data runs
//************************************************************************************
int runSingle(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, dgrid *xMeshGrid, dgrid *yMeshGrid, vector<double> *xSingleVector, vector<double> *ySingleVector, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double smoothingScaleX, double smoothingScaleY, double &x0, double &y0, double &x1, double &y1, double meanXSingle, double meanYSingle, string &kernelName,  map<string, int> additionalOptions,  string outputFileName, double subDataMulitplier, double UTMNorthingRef, double UTMEastingRef, double rotAngle, int RefEllip, char UTMZoneRef[4], MB_ZGRID_DATA *MB_ZGridInput, GMT_SURFACE_DATA *GMTSurfaceInput, ALG_SPLINE_DATA *ALGSplineInput, Bathy_Grid* bathyGrid, int USE_UTM)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int returnValue = SUCCESS;
	OUTPUT_DATA xyzOut;

	double UTMEasting;
	double UTMNorthing;
	double newLon;
	double newLat;
	bool ensembleFlag = false;
	#pragma region -- Gridding
	if(bathyGrid->GriddingFlag)
	{
		//MergeBathy used to be able to only run one gridding algorithm at a time.  
		//These gridding algorithms spline sparse data in order to obtain more 
		//values before continuing.
		//If multiple gridding algorithms have 2/-2 USAGE values
		//set, the new up-sampled pre-splined values replace the input
		//x,y,z,e.  These values would then be input in the next gridding algorithm
		//and replaced with those new up-sampled pre-splined values, etc. before continuing on
		//to mergeBathy.
		//This is no longer the case.  Now MergeBathy will send the original input
		//to each spline routine for separate calculations and then combine them with
		//into an ensemble average afterwards.
		int sumG = 0;
		if(additionalOptions.find("-ZGrid")->second == 1)
			sumG++;
		if (additionalOptions.find("-GMTSurface")->second == 1)
			sumG++;
		if (additionalOptions.find("-ALGSpline")->second == 1)
			sumG++;
		//if (additionalOptions.find("-Sibson_Natural_Neighbor_Interp")->second == 1)
		//	sumG++;

		if(sumG > 1)
		{
			cout << "Ensembling of Pre-Spliners has not been tested and is merely a prototype!" << endl;
			cout << "It is advised that you change your arguments to select one and run again." << endl;
			//return !SUCCESS;
		}
		vector<double> copyInputX;
		vector<double> copyInputY;
		vector<double> copyInputZ;
		vector<double> copyInputE;
		vector<double> copyInputH;
		vector<double> copyInputV;
		//We intend to ensemble instead of re-splining already splined results (iteratively splining).
		if(sumG > 1)
		{
			copyInputX.assign(inputDataX->begin(),inputDataX->end());
			copyInputY.assign(inputDataY->begin(),inputDataY->end());
			copyInputZ.assign(inputDataZ->begin(),inputDataZ->end());
			copyInputE.assign(inputDataE->begin(),inputDataE->end());
			copyInputH.assign(inputDataHErr->begin(),inputDataHErr->end());
			copyInputV.assign(inputDataVErr->begin(),inputDataVErr->end());
		}
		//************************************************************************************
		// I. Check and do MB_ZGrid if necessary
		//************************************************************************************
		if (additionalOptions.find("-ZGrid")->second == 1)
		{
			#pragma region -- MB_ZGrid
			cout << "********************************************************" << endl;
			cout << "* Computing MB_ZGrid" << endl;
			cout << "********************************************************" << endl;
			externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);
			bool mbzReturn = extInterp.run_MB_ZGrid(inputDataX, inputDataY, inputDataZ, inputDataE, inputDataHErr, inputDataVErr, x0, y0, x1, y1, additionalOptions, (*MB_ZGridInput).spacingX, (*MB_ZGridInput).spacingY, (*MB_ZGridInput).tension, (*MB_ZGridInput).z_OutputFileName, (*MB_ZGridInput).usage,  bathyGrid);

			if (!mbzReturn)
			{
				cerr << "MB_ZGrid FAILED! ABORTING!!!" << endl;
				return MB_ZGRID_ERROR;
			}
			cout << "Done Computing MB_ZGrid" << endl << endl;
			if(abs((*MB_ZGridInput).usage) == 1)
				return returnValue;

			if(sumG > 1)
			{
				inputDataX->assign(copyInputX.begin(),copyInputX.end());
				inputDataY->assign(copyInputY.begin(),copyInputY.end());
				inputDataZ->assign(copyInputZ.begin(),copyInputZ.end());
				inputDataE->assign(copyInputE.begin(),copyInputE.end());
				inputDataHErr->assign(copyInputH.begin(),copyInputH.end());
				inputDataVErr->assign(copyInputV.begin(),copyInputV.end());
			}
			#pragma endregion
		}
		//Done with MB_ZGrid

		//************************************************************************************
		// II. Check and do GMT Surface if necessary
		//************************************************************************************
		if (additionalOptions.find("-GMTSurface")->second == 1)
		{
			#pragma region -- GMTSurface
			cout << "********************************************************" << endl;
			cout << "* Computing GMT Surface" << endl;
			cout << "********************************************************" << endl;
			cout << "\nComputing GMT Surface" << endl;
			externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);
			bool gmtReturn = extInterp.run_Surface(inputDataX, inputDataY, inputDataZ, inputDataE, inputDataHErr, inputDataVErr, additionalOptions, (*GMTSurfaceInput).spacingX, (*GMTSurfaceInput).spacingY, (*GMTSurfaceInput).tension, (*GMTSurfaceInput).z_OutputFileName, (*GMTSurfaceInput).scaleFactor, (*GMTSurfaceInput).alpha, (*GMTSurfaceInput).usage, bathyGrid);

			if (!gmtReturn)
			{
				cerr << "GMT Surface FAILED! ABORTING!!!" << endl;
				return GMT_SURF_ERROR;
			}
			cout << "Done Computing GMT Surface" << endl << endl;
			if(abs((*GMTSurfaceInput).usage) == 1)
				return returnValue;
			if(sumG > 1)
			{
				inputDataX->assign(copyInputX.begin(),copyInputX.end());
				inputDataY->assign(copyInputY.begin(),copyInputY.end());
				inputDataZ->assign(copyInputZ.begin(),copyInputZ.end());
				inputDataE->assign(copyInputE.begin(),copyInputE.end());
				inputDataHErr->assign(copyInputH.begin(),copyInputH.end());
				inputDataVErr->assign(copyInputV.begin(),copyInputV.end());
			}
			#pragma endregion
		}
		//Done with GMT Surface

		//************************************************************************************
		// III. Check and do ALGSpline Surface if necessary
		//************************************************************************************
		if (additionalOptions.find("-ALGSpline")->second == 1)
		{
			#pragma region -- ALGSpline
			cout << "********************************************************" << endl;
			cout << "* Computing ALG Spline Surface" << endl;
			cout << "********************************************************" << endl;
			double ySize = (*xMeshGrid).rows();
			double xSize = (*xMeshGrid).cols();
			externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);

			//Must accept a regular grid.  Problem sam.
			bool algReturn = extInterp.run_ALGSpline(inputDataX, inputDataY, inputDataZ, x0,  x1, y0, y1, xSize, ySize, inputDataE, inputDataHErr, inputDataVErr, additionalOptions, (*ALGSplineInput).spacingX, (*ALGSplineInput).spacingY, (*ALGSplineInput).tension, (*ALGSplineInput).z_OutputFileName, (*ALGSplineInput).scaleFactor, (*ALGSplineInput).alpha, (*ALGSplineInput).usage, bathyGrid);

			if (!algReturn)
			{
				cerr << "ALG Spline FAILED! ABORTING!!!" << endl;
				return ALG_SPLINE_ERROR;
			}
			cout << "Done Computing ALG Spline." << endl << endl;
			if(abs((*ALGSplineInput).usage) == 1)
				return returnValue;
			if(sumG > 1)
			{
				inputDataX->assign(copyInputX.begin(),copyInputX.end());
				inputDataY->assign(copyInputY.begin(),copyInputY.end());
				inputDataZ->assign(copyInputZ.begin(),copyInputZ.end());
				inputDataE->assign(copyInputE.begin(),copyInputE.end());
				inputDataHErr->assign(copyInputH.begin(),copyInputH.end());
				inputDataVErr->assign(copyInputV.begin(),copyInputV.end());
			}
			#pragma endregion
		}
		//Done with ALGSpline Surface

		#pragma region --SibsonNN
/*
	//	if (additionalOptions.find("-Sibson_Natural_Neighbor_Interp")->second == 1)
		{
			cout << "\nComputing run_Sibson_Natural_Neighbor_Interp" << endl;
			double ySize = (*xMeshGrid).rows();
			double xSize = (*xMeshGrid).cols();

			bool Sibson_Natural_Neighbor_InterpReturn = extInterp.run_Sibson_Natural_Neighbor_Interp(&inputDataX_Copy4, &inputDataY_Copy4, &inputDataZ_Copy4, x0,  x1, y0, y1, xSize, ySize, &inputDataE_Copy4, &inputDataHErr_Copy4, &inputDataVErr_Copy4, additionalOptions, (*ALGSplineInput).spacingX, (*ALGSplineInput).spacingY, (*ALGSplineInput).tension, (*ALGSplineInput).z_OutputFileName, (*ALGSplineInput).scaleFactor, (*ALGSplineInput).alpha, (*ALGSplineInput).usage, bathyGrid);

			if (!Sibson_Natural_Neighbor_InterpReturn)
			{
				cerr << "ALG Spline FAILED! ABORTING!!!" << endl;
			//	return Sibson_Natural_Neighbor_InterpError;
			}
			cout << "Done Computing Sibson_Natural_Neighbor_Interp." << endl << endl;
		}
		*/
		#pragma endregion
		if(bathyGrid->getGrids().size() > 1){
			bathyGrid->ensemble();
			vector<double> ensembleResultsX = bathyGrid->getEnsemble_X();
			vector<double> ensembleResultsY = bathyGrid->getEnsemble_Y();
			vector<double> ensembleResultsZ = bathyGrid->getEnsemble_Z();
			vector<double> ensembleResults  = bathyGrid->getEnsemble_U();
			vector<double> ensembleResultsH = bathyGrid->getEnsemble_H();
			vector<double> ensembleResultsV = bathyGrid->getEnsemble_V();

			inputDataX->assign(ensembleResultsX.begin(),ensembleResultsX.end());
			inputDataY->assign(ensembleResultsY.begin(),ensembleResultsY.end());
			inputDataZ->assign(ensembleResultsZ.begin(),ensembleResultsZ.end());
			inputDataE->assign(ensembleResults.begin(),ensembleResults.end());
			inputDataHErr->assign(ensembleResultsH.begin(),ensembleResultsH.end());
			inputDataVErr->assign(ensembleResultsV.begin(),ensembleResultsV.end());

		}
	}
	bathyGrid->clear();
	
	#pragma endregion

	//************************************************************************************
	// III. If we are computing the data points then run scalecInterpTile.  Otherwise go the longer route of scalecInterp.
	//************************************************************************************
	if (additionalOptions.find("-preInterpolatedLocations")->second == 0)
	{
		//A. Do normal interpolation
  		returnValue = bathyTool(inputDataX, inputDataY, inputDataZ, inputDataE, inputDataHErr, inputDataVErr, xMeshGrid, yMeshGrid, xSingleVector, ySingleVector, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, 1.0, true, true, NEITOL, &xyzOut);
	}else
	{
		//B. Do pre-interpolated location interpolation
		returnValue = bathyToolPreDefined(inputDataX, inputDataY, inputDataZ, inputDataE, inputDataHErr, inputDataVErr, xInterpVector, yInterpVector, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, NEITOL, &xyzOut);
	}

	//C. Clear up the variables if something went wrong
	if(returnValue != SUCCESS)
	{
		xyzOut.depth.clear();
		xyzOut.error.clear();
		xyzOut.nEi.clear();
		xyzOut.rEi.clear();

		if(abs(additionalOptions.find("-propUncert")->second) == 1)
		{
			xyzOut.depth0.clear();
			xyzOut.error0.clear();
		}
		if(abs(additionalOptions.find("-kalman")->second) == 1)
		{
			xyzOut.depthK.clear();
			xyzOut.errorK.clear();
		}
		return returnValue;
	}

	for(int k = 0; k < (const int)xyzOut.nEi.size(); k++)
	{
		xyzOut.error[k] = sqrt(xyzOut.error[k]);
		xyzOut.nEi[k] = sqrt(xyzOut.nEi[k]);
		xyzOut.rEi[k] = sqrt(xyzOut.rEi[k]);
		xyzOut.standardDev[k] = sqrt(xyzOut.standardDev[k]);
		if(abs(additionalOptions.find("-propUncert")->second) == 1)
		{
			xyzOut.error0[k] = sqrt(xyzOut.error0[k]);
		}
	}

	vector<double> xInterpTemp;
	vector<double> yInterpTemp;
	vector<double> nEiTemp;
	vector<double> rEiTemp;
	vector<double> depthTemp;
	vector<double> errorTemp;
	vector<double> depthKTemp;
	vector<double> errorKTemp;
	vector<double> depth0Temp;
	vector<double> error0Temp;

	//same as above but creates new vars so both raster and txt can be printed for debugging
	/*vector<double> alldepthTemp(xyzOut.depth);
	vector<double> alldepth0Temp(xyzOut.depth0);
	vector<double> alldepthKTemp(xyzOut.depthK);
	vector<double> allerrorTemp(xyzOut.error);
	vector<double> allerror0Temp(xyzOut.error0);
	vector<double> allerrorKTemp(xyzOut.errorK);
	vector<double> allnEiTemp(xyzOut.nEi);
	vector<double> allrEiTemp(xyzOut.rEi);*/

	// Post-processing
	double nescale;
	vector<double>::iterator temp;
	vector<double> bar;
	vector<double>::iterator it;
	if(additionalOptions.find("-outputRasterFile")->second == 0 && additionalOptions.find("-outputBagFile")->second == 0)
	{
		for(it = xyzOut.nEi.begin(); it < xyzOut.nEi.end(); it++)
		{
			if(*it != NaN)
				bar.push_back(*it);
		}

		 nescale = median(bar) + 0.1;

		if (nescale > 0.95)
			nescale = 0.95;

		// removes edge extrapolation effects
		nescale = max(nescale,0.75);
		temp = std::min_element(xyzOut.nEi.begin(),xyzOut.nEi.end());
		if (*temp > nescale)
			nescale = *temp + 0.1;

		//remove edge extrapolation effects
		if(additionalOptions.find("-modelflag")->second == 0)
		{
			if(additionalOptions.find("-nonegdepth")->second==0)
			{
				//remove only edge extrapolation effects 
				for(int cnt = 0; cnt < (const int)xyzOut.nEi.size(); cnt++)
				{
					if(xyzOut.nEi[cnt] > nescale )
					{ //really poor samples
						xyzOut.depth[cnt] = NaN;
						xyzOut.error[cnt] = NaN;
						xyzOut.nEi[cnt] = NaN;
						xyzOut.rEi[cnt] = NaN;
						xyzOut.depth0[cnt] = NaN;
						xyzOut.error0[cnt] = NaN;
						xyzOut.depthK[cnt] = NaN;
						xyzOut.errorK[cnt] = NaN;
					}
				}
			}
			else
			{
				//remove edge extrapolation effects and negative depths
				for(int cnt = 0; cnt < (const int)xyzOut.nEi.size(); cnt++)
				{
					if(xyzOut.nEi[cnt] > nescale || xyzOut.depth[cnt] < 0.1)//SJZ changed it from && 1/21/16
					{//This is weird... I think it should be || instead of &&
						xyzOut.depth[cnt] = NaN;
						xyzOut.error[cnt] = NaN;
						xyzOut.nEi[cnt] = NaN;
						xyzOut.rEi[cnt] = NaN;
						xyzOut.depth0[cnt] = NaN;
						xyzOut.error0[cnt] = NaN;
						xyzOut.depthK[cnt] = NaN;
						xyzOut.errorK[cnt] = NaN;
					}
				}
			}
			cout << "Points with normalized errors greater than "<< nescale << " removed." << endl;
		}
		else
			//Model Mode; Keep all points
			cout << "Nothing is bad!" << endl;

		//Keep only data points for regular output; use smaller data set
		for(int i = 0; i < (const int)(*xInterpVector).size(); i++)
		{
			if(xyzOut.depth[i] != NaN)
				xInterpTemp.push_back((*xInterpVector)[i]);
			if(xyzOut.depth[i] != NaN)
				yInterpTemp.push_back((*yInterpVector)[i]);
			if(xyzOut.depth[i] != NaN)
				depthTemp.push_back(xyzOut.depth[i]);
			if(xyzOut.error[i] != NaN)
				errorTemp.push_back(xyzOut.error[i]);
			if(xyzOut.nEi[i] != NaN)
				nEiTemp.push_back(xyzOut.nEi[i]);
			if(xyzOut.rEi[i] != NaN)
				rEiTemp.push_back(xyzOut.rEi[i]);
			if(xyzOut.depth0[i] != NaN)
				depth0Temp.push_back(xyzOut.depth0[i]);
			if(xyzOut.error0[i] != NaN)
				error0Temp.push_back(xyzOut.error0[i]);
			if(xyzOut.depthK[i] != NaN)
				depthKTemp.push_back(xyzOut.depthK[i]);
			if(xyzOut.errorK[i] != NaN)
				errorKTemp.push_back(xyzOut.errorK[i]);
		}
	}
	else 
	{	//Raster output
		xInterpTemp		= *xInterpVector;
		yInterpTemp		= *yInterpVector;
		depthTemp		= xyzOut.depth;
		depth0Temp		= xyzOut.depth0;
		depthKTemp		= xyzOut.depthK;
		errorTemp		= xyzOut.error;
		error0Temp		= xyzOut.error0;
		errorKTemp		= xyzOut.errorK;
		nEiTemp			= xyzOut.nEi;
		rEiTemp			= xyzOut.rEi;
	}

	vector<double> xMeshVector;
	vector<double> yMeshVector;
	vector<double> depthTemp2;
	vector<double> errorTemp2;
	vector<double> neiTemp2;
	vector<double> reiTemp2;
	vector<double> depth0Temp2;
	vector<double> error0Temp2;
	vector<double> depthKTemp2;
	vector<double> errorKTemp2;

//	vector<double> UTMNorthings((*xInterpVector).size());
//	vector<double> UTMEastings((*xInterpVector).size());
	vector<double> UTMNorthings(xInterpTemp.size());
	vector<double> UTMEastings(xInterpTemp.size());
	//vector<double> xInterpTemp2(xInterpTemp.size());
	//vector<double> yInterpTemp2(xInterpTemp.size());

	//SJZ clean up to try to run FOUO on 32bit**********************************
	//vector<double> xInterpTemp3;
	//vector<double> yInterpTemp3;
	//vector<double> depthTemp3;
	//vector<double> errorTemp3;
	vector<double> neiTemp3;
	//vector<double> reiTemp3;
	//vector<double> depth0Temp3;
	//vector<double> error0Temp3;
	//vector<double> depthKTemp3;
	//vector<double> errorKTemp3;
	//	
	//vector<double> xInterpTemp4;
	//vector<double> yInterpTemp4;
	//vector<double> depthTemp4;
	//vector<double> errorTemp4;
	//vector<double> neiTemp4;
	//vector<double> reiTemp4;
	//vector<double> depth0Temp4;
	//vector<double> error0Temp4;
	//vector<double> depthKTemp4;
	//vector<double> errorKTemp4;
	//end clean up
	int xSize = 0, ySize = 0;
	vector<double> xt;
	vector<double> yt;
	vector<double> x;
	vector<double> y;
	dgrid newxx;
	dgrid newyy;
	vector<double> newxxVec;
	vector<double> newyyVec;
	//************************************************************************************
	// IV. Convert from UTM to Lat/Lon so we can output the data.
	//************************************************************************************
	if ((additionalOptions.find("-inputInMeters")->second == 0) && (additionalOptions.find("-outputRasterFile")->second == 0 && additionalOptions.find("-outputBagFile")->second == 0))
	{
		for(int i = 0; i < (const int)xInterpTemp.size(); i++)
		{
			if(!USE_UTM) // SJZ
			{
				UTMEasting = (xInterpTemp[i])*cos(deg2rad*(-rotAngle)) - (yInterpTemp[i])*sin(deg2rad*(-rotAngle));
				UTMNorthing = (xInterpTemp[i])*sin(deg2rad*(-rotAngle)) + (yInterpTemp[i])*cos(deg2rad*(-rotAngle));

				UTMNorthing += UTMNorthingRef;
				UTMEasting += UTMEastingRef;
			}
			else
			{ // SJZ
				UTMEasting	= xInterpTemp[i];
				UTMNorthing = yInterpTemp[i];
			}
			UTMtoLL(RefEllip, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
			if(newLon >= 180) // SJZ
				newLon -= 360;
			yInterpTemp[i] = newLat;
			xInterpTemp[i] = newLon;
		}
		if (additionalOptions.find("-printMatlabMatch")->second == 1)
		{
			if (((rotAngle > 45) & (rotAngle <= 135)) | ((rotAngle > 225) & (rotAngle <= 315)))
				cerr << "Rotation Angle must be (<= 45 & >315) OR (> 135 & <= 225) to match Matlab output otherwise the data is transposed!" << endl;
		}
	}
	else if (additionalOptions.find("-outputRasterFile")->second == 1 || additionalOptions.find("-outputBagFile")->second == 1)
	{ 
		dgrid xx;
		dgrid yy;
		vector<double> zz;//, zz0, zzK, ee, ee0, eeK, neei, reei; // SJZ
		double minNorthing	= (double)MAX_INT;
		double minEasting	= (double)MAX_INT;
		double maxNorthing	= (double)MIN_INT;
		double maxEasting	= (double)MIN_INT;
		/*double minXInterp	= (double)MAX_INT;
		double minYInterp	= (double)MAX_INT;
		double maxXInterp	= (double)MIN_INT;
		double maxYInterp	= (double)MIN_INT;*/
		int INIT_TEMPS = 1;

		//1. Remove the central point offset from the data for raster format grids
		//A new raster grid will be created by interpolating from current data to raster grid locations
		//Convert only points with data to UTM.
		//for(int i = 0; i < (const int)xInterpVector->size(); i++)
		for(int i = 0; i < (const int)xInterpTemp.size(); i++)
		{
			if(!USE_UTM) // SJZ
			{ 
				UTMEasting = (xInterpTemp)[i]*cos(deg2rad*(-rotAngle)) - (yInterpTemp)[i]*sin(deg2rad*(-rotAngle));
				UTMNorthing = (xInterpTemp)[i]*sin(deg2rad*(-rotAngle)) + (yInterpTemp)[i]*cos(deg2rad*(-rotAngle));

				//UTMEasting = (*xInterpVector)[i]*cos(deg2rad*(-rotAngle)) - (*yInterpVector)[i]*sin(deg2rad*(-rotAngle));
				//UTMNorthing = (*xInterpVector)[i]*sin(deg2rad*(-rotAngle)) + (*yInterpVector)[i]*cos(deg2rad*(-rotAngle));
			
				UTMNorthing += UTMNorthingRef;
				UTMEasting  += UTMEastingRef;
			}
			else
			{ // SJZ
				UTMEasting	= xInterpTemp[i];
				UTMNorthing = yInterpTemp[i];
			}
			if(INIT_TEMPS)
			{
				minNorthing	= UTMNorthing;
				minEasting	= UTMEasting;
				maxNorthing	= UTMNorthing;
				maxEasting	= UTMEasting;
			/*	minXInterp	= MAX_INT;
				minYInterp	= MAX_INT;
				maxXInterp	= MIN_INT;
				maxYInterp	= MIN_INT;*/
				INIT_TEMPS = 0;
			}
			UTMNorthings[i] = UTMNorthing;	//y
			UTMEastings [i] = UTMEasting;	//x
			//yInterpTemp2[i] = UTMNorthing;
			//xInterpTemp2[i] = UTMEasting;
			//(*yInterpVector)[i] = UTMNorthing;
			//(*xInterpVector)[i] = UTMEasting; 
			UTMtoLL(RefEllip, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
			if(newLon >= 180) // SJZ
				newLon -= 360;
			yInterpTemp[i] = newLat;
			xInterpTemp[i] = newLon;

			minEasting  = min(minEasting, UTMEasting);
			maxEasting  = max(maxEasting, UTMEasting);
			minNorthing = min(minNorthing, UTMNorthing);
			maxNorthing = max(maxNorthing, UTMNorthing);
		}
		
		int Ni = (const int)UTMNorthings.size();
		string interpMethod;
		Bathy_Grid new_bathyGrid=Bathy_Grid();
				
	//	new_bathyGrid.Construct_TinRaster(&UTMEastings, &UTMNorthings, &alldepthTemp, &allerrorTemp, &allnEiTemp, &allrEiTemp, &alldepth0Temp, &allerror0Temp, &alldepthKTemp, &allerrorKTemp);
		new_bathyGrid.Construct_TinRaster(&UTMEastings, &UTMNorthings, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &depth0Temp, &error0Temp, &depthKTemp, &errorKTemp);//sjz

		//A. Obtain a Grid
		// Calculate x and y for our new grid where we go from
		createMeshXYDims(minEasting, minNorthing, maxEasting, maxNorthing, gridSpacingX, gridSpacingY, &xt, &yt);

		// get our new x,y dimensions
		xSize = (const int)xt.size();
		ySize = (const int)yt.size();
			
		// Resize vectors to match our new dimensions
		xMeshVector.resize(xSize*ySize, 0.00);
		yMeshVector.resize(xSize*ySize, 0.00);
		xx.resize(ySize, xSize, 0.00);
		yy.resize(ySize, xSize, 0.00);
		
		// Create grid with those dimensions
		createMeshGrid(&xt, &yt, &xMeshVector, &yMeshVector, &xx, &yy);
			
		// Reshape grid
		newxx.resize(xx.rows(), xx.cols(),0.00);
		newyy.resize(yy.rows(), yy.cols(),0.00);
		reshapeGrid(&xx, &newxx);
		reshapeGrid(&yy, &newyy);

		//B.Interpolate to get z, e, nmsei, rei, kz, kvar values at raster grid locations.
		InterpGrid new_gmt = InterpGrid(GMT);
		zz.resize((newxx.vec()).size(), 0.00);
		(new_gmt).estimate(&(newxx.vec()), &(newyy.vec()), &zz, 1.96, 2.00, 0, new_bathyGrid.getTin(), "bilinear", "bilinear", "none");

		depthTemp2 = *(new_gmt).getZ();
		errorTemp2 = (new_gmt).getE();
		depth0Temp2 = (new_gmt).getZ0();
		error0Temp2 = (new_gmt).getE0();
		depthKTemp2 = (new_gmt).getZK();
		errorKTemp2 = (new_gmt).getEK();
		neiTemp2 = (new_gmt).getNEI();
		reiTemp2 = (new_gmt).getREI();

		new_gmt.clear();
		new_bathyGrid.clear();
		xx.clear();
		yy.clear();
		xt.clear();
		yt.clear();
		
		newxxVec = newxx.vec();
		newyyVec = newyy.vec();
		
		//Keep only points in our range; set bad points to NaN
		for(int i = 0; i < (const int)(newxxVec).size(); i++)
		{
			if((newxxVec[i] < minEasting || newxxVec[i] > maxEasting) || (newyyVec[i] < minNorthing || newyyVec[i] > maxNorthing))
			{
				depthTemp2[i] = NaN;
				errorTemp2[i] = NaN;
				neiTemp2[i] = NaN;
				reiTemp2[i] = NaN;
				if(abs(additionalOptions.find("-propUncert")->second) == 1)
				{
					depth0Temp2[i] = NaN;
					error0Temp2[i] = NaN;
				}
				if(abs(additionalOptions.find("-kalman")->second) == 1)
				{
					depthKTemp2[i] = NaN;
					errorKTemp2[i] = NaN;
				}
			}
			else
			{
				/*xInterpTemp3.push_back(newxxVec[i]); 
				yInterpTemp3.push_back(newyyVec[i]); 
				depthTemp3.push_back(depthTemp2[i]); 
				errorTemp3.push_back(errorTemp2[i]);*/
				neiTemp3.push_back(neiTemp2[i]);
				/*reiTemp3.push_back(reiTemp2[i]); 
				if(abs(additionalOptions.find("-propUncert")->second) == 1)
				{
					depth0Temp3.push_back(depth0Temp2[i]);
					error0Temp3.push_back(error0Temp2[i]);
				}
				if(abs(additionalOptions.find("-kalman")->second) == 1)
				{
					depthKTemp3.push_back(depthKTemp2[i]);
					errorKTemp3.push_back(errorKTemp2[i]);
				}*/
			}
		}

		bar.clear();
		// Post-processing
		for(it = neiTemp3.begin(); it < neiTemp3.end(); it++)
		{
			if(*it != NaN)
				bar.push_back(*it);
		}

		nescale = median(bar) + 0.1;

		if (nescale > 0.95)
			nescale = 0.95;

		// removes edge extrapolation effects
		nescale = max(nescale,0.75);
		temp = std::min_element(neiTemp3.begin(),neiTemp3.end());
		if (*temp > nescale)
			nescale = *temp + 0.1;
	
		//Remove edge extrapolation effects
		if(additionalOptions.find("-modelflag")->second == 0)
		{
			if(additionalOptions.find("-nonegdepth")->second==0)
			{
				//remove only edge extrapolation effects
				for(int cnt = 0; cnt < (const int)neiTemp2.size(); cnt++)
				{
					if(neiTemp2[cnt] > nescale || neiTemp2[cnt] < 0 )
					{ //really poor samples
						depthTemp2[cnt] = NaN;
						errorTemp2[cnt] = NaN;
						//Keep all nei
						//neiTemp2[cnt] = NaN;
						reiTemp2[cnt] = NaN;
						if(abs(additionalOptions.find("-propUncert")->second) == 1)
						{
							depth0Temp2[cnt] = NaN;
							error0Temp2[cnt] = NaN;
						}
						if(abs(additionalOptions.find("-kalman")->second) == 1)
						{
							depthKTemp2[cnt] = NaN;
							errorKTemp2[cnt] = NaN;
						}
					}
				}
			}
			else
			{
				//remove edge extrapolation effects and negative depths
				for(int cnt = 0; cnt < (const int)neiTemp2.size(); cnt++)
				{ 
					if((neiTemp2[cnt] > nescale || neiTemp2[cnt] < 0) || depthTemp2[cnt] < 0.1)//SJZ changed it from && 1/21/16
					{
						depthTemp2[cnt] = NaN;
						errorTemp2[cnt] = NaN;
						//Keep all nei
						//neiTemp2[cnt] = NaN;
						reiTemp2[cnt] = NaN;
						if(abs(additionalOptions.find("-propUncert")->second) == 1)
						{
							depth0Temp2[cnt] = NaN;
							error0Temp2[cnt] = NaN;
						}
						if(abs(additionalOptions.find("-kalman")->second) == 1)
						{
							depthKTemp2[cnt] = NaN;
							errorKTemp2[cnt] = NaN;
						}
					}
				}
			}
			cout << "Points with normalized errors greater than "<< nescale << " removed." << endl;
		}
		else //Model Mode; Keep all points!
			cout << "Nothing is bad!" << endl;

		//for(int i = 0; i < (const int)(newxxVec).size(); i++) 
		//{
		//	if(depthTemp2[i] != NaN)
		//		xInterpTemp4.push_back((newxxVec)[i]);
		//	if(depthTemp2[i] != NaN)
		//		yInterpTemp4.push_back((newyyVec)[i]);
		//	if(depthTemp2[i] != NaN)
		//		depthTemp4.push_back(depthTemp2[i]);
		//	if(errorTemp2[i] != NaN)
		//		errorTemp4.push_back(errorTemp2[i]);
		//	//Keep all nei
		///*	if(neiTemp2[i] != NaN)
		//		neiTemp4.push_back(neiTemp2[i]);*/
		//	if(reiTemp2[i] != NaN)
		//		reiTemp4.push_back(reiTemp2[i]);
		//	if(abs(additionalOptions.find("-propUncert")->second) == 1)
		//	{
		//		if(depth0Temp2[i] != NaN)
		//			depth0Temp4.push_back(depth0Temp2[i]);
		//		if(error0Temp2[i] != NaN)
		//			error0Temp4.push_back(error0Temp2[i]);
		//	}
		//	if(abs(additionalOptions.find("-kalman")->second) == 1)
		//	{
		//		if(depthKTemp2[i] != NaN)
		//			depthKTemp4.push_back(depthKTemp2[i]);
		//		if(errorKTemp2[i] != NaN)
		//			errorKTemp4.push_back(errorKTemp2[i]);
		//	}
		//}
	}

	//************************************************************************************
	// V. Write the data to an output file and clean up.
	//************************************************************************************
	//Find last underscore
	string fileName = outputFileName;
	string fname;
	string f;

	if(additionalOptions.find("-appendFilename")->second == 1)
	{
		//size_t found = fileName.back();
		size_t found = *fileName.rbegin(); // UNIX
		
		//Print results to file and attach the window name used to filename.
		if(kernelName.compare("hanning") == 0 || kernelName.compare("hann") == 0)
			fname = fileName.substr(0,found+1).append("hann");
		else if (kernelName.compare("boxcar") == 0)
			fname = fileName.substr(0,found+1).append("boxcar");
		else if (kernelName.compare("quadloess") == 0)
			fname = fileName.substr(0,found+1).append("quadloess");
		else if (kernelName.compare("loess") == 0)
			fname = fileName.substr(0,found+1).append("loess");

		if (additionalOptions.find("-nnInterp")->second == 1)
			fname = fname.append("_NN");
	}
	else 
		fname = fileName;

	if (additionalOptions.find("-outputRasterFile")->second == 0 && additionalOptions.find("-outputBagFile")->second == 0)
	{
		//A. Write the data to a flat file
		if((additionalOptions.find("-mse")->second == 1 || additionalOptions.find("-propUncert")->second == 1) || additionalOptions.find("-kalman")->second == 1)
		{
			f = fname;
			writeFiles(f, &xInterpTemp, &yInterpTemp, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions, &depth0Temp, &error0Temp, &depthKTemp, &errorKTemp);
		}
		//i. Print MSE with Kalman results tagged on the end.
		if(additionalOptions.find("-printMSEwK")->second == 1)
		{
			f = fname;
			f = f.append("MSEwK");
			writeFile(f, &xInterpTemp, &yInterpTemp, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions, &depthKTemp, &errorKTemp);
		}
		//ii. Print matlab matched output;
		if (additionalOptions.find("-printMatlabMatch")->second == 1)
		{
			f = fname;
			f = f.append("Matlab");
			writeFileMatch(f, &xInterpTemp, &yInterpTemp, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions,&depthKTemp, &errorKTemp);
		}
		//iii. Each estimator prints to its own output file.
	/*	if(additionalOptions.find("-mse")->second == 1)
		{
			f = fname;
			writeFile(f, &xInterpTemp, &yInterpTemp, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions);
		}
		if(additionalOptions.find("-propUncert")->second == 1)
		{
			f = fname;
			f = f.append("P");
			writeFile(f, &xInterpTemp, &yInterpTemp, &depth0Temp, &error0Temp, &nEiTemp, &rEiTemp, &additionalOptions);
		}
		if(additionalOptions.find("-kalman")->second == 1)
		{
			f = fname;
			f = f.append("K");
			writeFile(f, &xInterpTemp, &yInterpTemp, &depthKTemp, &errorKTemp, &nEiTemp, &rEiTemp, &additionalOptions);
		}*/
	}
	else if (additionalOptions.find("-outputRasterFile")->second == 0)
	{
		//B. Write the data to a bag file
		xSize = newxx.rows();
		ySize = newxx.cols();

		f = fname;
		writeBagFile(f, &newxxVec, &newyyVec, &depthTemp2, &errorTemp2, &neiTemp2, &reiTemp2,(int)gridSpacingX, (int)gridSpacingY, xSize, ySize, &additionalOptions, UTMZoneRef, &depth0Temp2, &error0Temp2, &depthKTemp2, &errorKTemp2);

		//Debug by printing bag to raster
	//	f = fname;
	//	readBagFile(f, &newxxVec, &newyyVec, &depthTemp2, &errorTemp2, &neiTemp2, &reiTemp2, (int)gridSpacingX, (int)gridSpacingY, xSize, ySize, &additionalOptions);
	}
	else
	{
		//B. Write the data to a raster file
		xSize = newxx.rows();
		ySize = newxx.cols();

		////Print MSE raster files 
		////write mse file for debugging raster--delete when done
		//f = fname;
		//f = f.append("mse");
		////original mse values no nans
		////writeFile(f, &xInterpTemp2, &yInterpTemp2, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions);
		//writeFile(f, &UTMEastings, &UTMNorthings, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions);
		//
		////f = fname;
		////f = f.append("mseALL");
		////writeFile(f, &UTMEastings, &UTMNorthings, &alldepthTemp, &allerrorTemp, &allnEiTemp, &allrEiTemp, &additionalOptions);
		//	
		////new raster values with before its been naned
		//f = fname;
		//f = f.append("mse3");
		//writeFile(f, &xInterpTemp3, &yInterpTemp3,  &depthTemp3, &errorTemp3, &neiTemp3, &reiTemp3, &additionalOptions);

		f = fname;
		writeRasterFile(f, &newxxVec, &newyyVec, &depthTemp2, &errorTemp2, &neiTemp2, &reiTemp2, gridSpacingX, gridSpacingY, xSize, ySize, &additionalOptions, UTMZoneRef, &depth0Temp2, &error0Temp2, &depthKTemp2, &errorKTemp2);
	}

	cout << "Done Creating Output File" << endl;

	//C. Clear up more variables
	xyzOut.depth.clear();
	xyzOut.error.clear();
	xyzOut.nEi.clear();
	xyzOut.rEi.clear();
	if(abs(additionalOptions.find("-propUncert")->second) == 1)
	{
		xyzOut.depth0.clear();
		xyzOut.error0.clear();
	}
	if(abs(additionalOptions.find("-kalman")->second) == 1 )
	{
		xyzOut.depthK.clear();
		xyzOut.errorK.clear();
	}
	return returnValue;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ORIGINAL FUNCTIONS FOR MONTE CARLO RUNS!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//************************************************************************************
// SUBROUTINE III: Primary MergeBathy function call for processing Monte Carlo data runs
//************************************************************************************
int runMonteCarlo(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, dgrid *xMeshGrid, dgrid *yMeshGrid, vector<double> *xSingleVector, vector<double> *ySingleVector, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double smoothingScaleX, double smoothingScaleY, double &x0, double &y0, double &x1, double &y1, double meanXSingle, double meanYSingle, string &kernelName,  map<string, int> additionalOptions, string outputFileName, double subDataMulitplier, double UTMNorthingRef, double UTMEastingRef, double rotAngle, int RefEllip, char UTMZoneRef[4], int numMCRuns, MB_ZGRID_DATA *MB_ZGridInput, GMT_SURFACE_DATA *GMTSurfaceInput, ALG_SPLINE_DATA *ALGSplineInput, Bathy_Grid* bathyGrid)
{
	//MergeBathy used to be able to only run one gridding algorithm at a time.  
	//These gridding algorithms spline sparse data in order to obtain more 
	//values before continuing.
	//If multiple gridding algorithms have 2/-2 USAGE values
	//set, the new up-sampled pre-splined values replace the input
	//x,y,z,e.  These values would then be input in the next gridding algorithm
	//and replaced with those new up-sampled pre-splined values, etc. before continuing on
	//to mergeBathy.
	//This is no longer the case.  Now MergeBathy will send the original input
	//to each spline routine for separate calculations and then combine them with
	//into an ensemble average afterwards.
	int sumG = additionalOptions.find("-ZGrid")->second + additionalOptions.find("-GMTSurface")->second
		+ additionalOptions.find("-ALGSpline")->second;// + additionalOptions.find("-Sibson_Natural_Neighbor_Interp")->second;
	
	if(sumG > 1)
	{
		cout << "Ensembling of Pre-Spliners has not been tested and is merely a prototype!" << endl;
		cout << "It is advised that you change your arguments to select one and run again." << endl;
		//return !SUCCESS;
	}

	/*
	//This is commented out because I am calling the current calls
	cout << endl;
	cout << "WARNING!!!!!" << endl;
	cout << "Monte Carlo and it's subsequent functions are the original functions from MergeBathy v3.6." << endl;
	cout << "Therefore, all changes, improvements, bug fixes, etc. including Kalman, hU and vU, scalecInterp, and scalecInterpTile will not be applied. << end;
	cout << "Take care when using Monte Carlo until it is tested for correctness!!!" << endl;
	cout << "Press any key to continue... " << endl;
	*/
	//cin.get(); //listens for key press

	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int returnValue = SUCCESS;
	int i, stdLoc;
	OUTPUT_DATA xyzOut;

	double UTMEasting;
	double UTMNorthing;

	RNG randomNumber;
	stringstream ss;
	string numStr;
	string fileName;
	string extInterpFileName;

	vector<double> xMC = vector<double>((*inputDataX).size());
	vector<double> yMC = vector<double>((*inputDataY).size());
	vector<double> zMC = vector<double>((*inputDataX).size());
	vector<double> eMC = vector<double>((*inputDataY).size());
	//SJZ
	vector<double> hMC = vector<double>((*inputDataY).size());
	vector<double> vMC = vector<double>((*inputDataY).size());

	dgrid xMeshGridMC;
	dgrid yMeshGridMC;
	vector<double> xtMC;
	vector<double> ytMC;
	vector<double> xMeshVectorMC;
	vector<double> yMeshVectorMC;

	dvector ziToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
	dvector eiToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
	dvector neiToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
	dvector reiToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);

	//************************************************************************************
	// I. Iterate over the number of Monte Carlo runs to be done
	//************************************************************************************
	stdLoc = 0;
	for (int mcRunNum = 0; mcRunNum < numMCRuns; mcRunNum++)
	{
		printf("\n\nEntering Monte Carlo Run Number: %d\n",(mcRunNum + 1));

		//A. Preinit the file numbers
		ss << mcRunNum;
		ss >> numStr;
		ss.clear();

		fileName = outputFileName;

		//B. Depending on the number of data runs pre pad the output file names with the run number
		if (mcRunNum < 10)
			fileName.append("_00");
		else if(mcRunNum < 100)
			fileName.append("_0");
		else
			fileName.append("_");
		fileName.append(numStr);
		//Done with preInit

		//C. Set disposable variables for each loop
		for (i = 0; i < (const int)xMC.size(); i++){
			xMC[i] = (*inputDataX)[i] + (*inputDataHErr)[i]*randomNumber.normal();
			yMC[i] = (*inputDataY)[i] + (*inputDataHErr)[i]*randomNumber.normal();
			zMC[i] = (*inputDataZ)[i];
			eMC[i] = (*inputDataE)[i];

			//SJZ
			hMC[i] = (*inputDataHErr)[i];
			vMC[i] = (*inputDataVErr)[i];
		}
		
		#pragma region -- Gridding
		if(bathyGrid->GriddingFlag)
		{
			//We intend to ensemble instead of re-splining already splined results (iteratively splining).
			vector<double> copyInputX;
			vector<double> copyInputY;
			vector<double> copyInputZ;
			vector<double> copyInputE;
			vector<double> copyInputH;
			vector<double> copyInputV;
			if(sumG > 1)
			{
				copyInputX.assign(xMC.begin(),xMC.end());
				copyInputY.assign(yMC.begin(),yMC.end());
				copyInputZ.assign(zMC.begin(),zMC.end());
				copyInputE.assign(eMC.begin(),eMC.end());
				copyInputH.assign(hMC.begin(),hMC.end());
				copyInputV.assign(vMC.begin(),vMC.end());
			}
			//************************************************************************************
			// II. Check and do MB_ZGrid if necessary
			//************************************************************************************
			if (additionalOptions.find("-ZGrid")->second == 1)
			{
				#pragma region --MBZ
				cout << "********************************************************" << endl;
				cout << "* Computing MB_ZGrid" << endl;
				cout << "********************************************************" << endl;

				//A. Set iterative filenames
				extInterpFileName = (*MB_ZGridInput).z_OutputFileName;
				if (mcRunNum < 10)
					extInterpFileName.append("_00");
				else if(mcRunNum < 100)
					extInterpFileName.append("_0");
				else
					extInterpFileName.append("_");
				extInterpFileName.append(numStr);

				extInterpFileName.append((".txt"));

				//B. Now call mb_zgrid
				externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);
				//bool mbzReturn = extInterp.run_MB_ZGrid_ORIGINAL(&xMC, &yMC, &zMC, &eMC, inputDataHErr, inputDataVErr, x0, y0, x1, y1, (*MB_ZGridInput).spacingX, (*MB_ZGridInput).spacingY, (*MB_ZGridInput).tension, extInterpFileName, (*MB_ZGridInput).usage);
				bool mbzReturn = extInterp.run_MB_ZGrid(&xMC, &yMC, &zMC, &eMC, &hMC, &vMC, x0, y0, x1, y1, additionalOptions, (*MB_ZGridInput).spacingX, (*MB_ZGridInput).spacingY, (*MB_ZGridInput).tension, extInterpFileName, (*MB_ZGridInput).usage, bathyGrid);

				if (!mbzReturn)
				{
					cerr << "MB_ZGrid FAILED! ABORTING!!!" << endl;
					return MB_ZGRID_ERROR;
				}
				cout << "Done Computing MB_ZGrid" << endl << endl;
				if(abs((*MB_ZGridInput).usage) == 2)
					return returnValue;
				if(sumG > 1)
				{
					xMC.assign(copyInputX.begin(),copyInputX.end());
					yMC.assign(copyInputY.begin(),copyInputY.end());
					zMC.assign(copyInputZ.begin(),copyInputZ.end());
					eMC.assign(copyInputE.begin(),copyInputE.end());
					hMC.assign(copyInputH.begin(),copyInputH.end());
					vMC.assign(copyInputV.begin(),copyInputV.end());
				}
				#pragma endregion
			}
			//Done with MB_ZGrid

			//************************************************************************************
			// III. Check and do GMT Surface if necessary
			//************************************************************************************
			if (additionalOptions.find("-GMTSurface")->second == 1)
			{
				#pragma region --GMTSurface
				cout << "********************************************************" << endl;
				cout << "* Computing GMT Surface" << endl;
				cout << "********************************************************" << endl;

				//A. Set iterative filenames
				extInterpFileName = (*GMTSurfaceInput).z_OutputFileName;
				if (mcRunNum < 10)
					extInterpFileName.append("_00");
				else if(mcRunNum < 100)
					extInterpFileName.append("_0");
				else
					extInterpFileName.append("_");
				extInterpFileName.append(numStr);

				extInterpFileName.append((".txt"));

				//B. Now call GMT Surface
				externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);
				//bool mbzReturn = extInterp.run_Surface_ORIGINAL(&xMC, &yMC, &zMC, &eMC, inputDataHErr, inputDataVErr, (*GMTSurfaceInput).spacingX, (*GMTSurfaceInput).spacingY, (*GMTSurfaceInput).tension, extInterpFileName, (*GMTSurfaceInput).scaleFactor, (*GMTSurfaceInput).alpha, (*GMTSurfaceInput).usage);
				bool gmtReturn = extInterp.run_Surface(&xMC, &yMC, &zMC, &eMC, &hMC, &vMC, additionalOptions, (*GMTSurfaceInput).spacingX, (*GMTSurfaceInput).spacingY, (*GMTSurfaceInput).tension, extInterpFileName, (*GMTSurfaceInput).scaleFactor, (*GMTSurfaceInput).alpha, (*GMTSurfaceInput).usage, bathyGrid);

				if (!gmtReturn)
				{
					cerr << "GMT Surface FAILED! ABORTING!!!" << endl;
					return GMT_SURF_ERROR;
				}
				cout << "Done Computing GMT Surface" << endl << endl;
				if(abs((*GMTSurfaceInput).usage) == 2)
					return returnValue;
				if(sumG > 1)
				{
					xMC.assign(copyInputX.begin(),copyInputX.end());
					yMC.assign(copyInputY.begin(),copyInputY.end());
					zMC.assign(copyInputZ.begin(),copyInputZ.end());
					eMC.assign(copyInputE.begin(),copyInputE.end());
					hMC.assign(copyInputH.begin(),copyInputH.end());
					vMC.assign(copyInputV.begin(),copyInputV.end());
				}
				#pragma endregion
			}
			//Done with GMT Surface
		
			//************************************************************************************
			// III. Check and do ALGSpline Surface if necessary
			//************************************************************************************
			if (additionalOptions.find("-ALGSpline")->second == 1)
			{
				#pragma region -- ALGSpline
				cout << "********************************************************" << endl;
				cout << "* Computing ALG Spline" << endl;
				cout << "********************************************************" << endl;
			
				//A. Set iterative filenames
				extInterpFileName = (*GMTSurfaceInput).z_OutputFileName;
				if (mcRunNum < 10)
					extInterpFileName.append("_00");
				else if(mcRunNum < 100)
					extInterpFileName.append("_0");
				else
					extInterpFileName.append("_");
				extInterpFileName.append(numStr);

				extInterpFileName.append((".txt"));

				double ySize = (*xMeshGrid).rows();
				double xSize = (*xMeshGrid).cols();
				externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);

				//Must accept a regular grid.  Problem sam.
				bool algReturn = extInterp.run_ALGSpline(inputDataX, inputDataY, inputDataZ, x0,  x1, y0, y1, xSize, ySize, &eMC, &hMC, &vMC, additionalOptions, (*ALGSplineInput).spacingX, (*ALGSplineInput).spacingY, (*ALGSplineInput).tension, (*ALGSplineInput).z_OutputFileName, (*ALGSplineInput).scaleFactor, (*ALGSplineInput).alpha, (*ALGSplineInput).usage, bathyGrid);

				if (!algReturn)
				{
					cerr << "ALG Spline FAILED! ABORTING!!!" << endl;
					return ALG_SPLINE_ERROR;
				}
				cout << "Done Computing ALG Spline." << endl << endl;
				if(abs((*ALGSplineInput).usage) == 1)
					return returnValue;
				if(sumG > 1)
				{
					xMC.assign(copyInputX.begin(),copyInputX.end());
					yMC.assign(copyInputY.begin(),copyInputY.end());
					zMC.assign(copyInputZ.begin(),copyInputZ.end());
					eMC.assign(copyInputE.begin(),copyInputE.end());
					hMC.assign(copyInputH.begin(),copyInputH.end());
					vMC.assign(copyInputV.begin(),copyInputV.end());
				}
				#pragma endregion
			}

			if(bathyGrid->getGrids().size() > 1){
				bathyGrid->ensemble();
				vector<double> ensembleResultsX = bathyGrid->getEnsemble_X();
				vector<double> ensembleResultsY = bathyGrid->getEnsemble_Y();
				vector<double> ensembleResultsZ = bathyGrid->getEnsemble_Z();
				vector<double> ensembleResults  = bathyGrid->getEnsemble_U();
				vector<double> ensembleResultsH = bathyGrid->getEnsemble_H();
				vector<double> ensembleResultsV = bathyGrid->getEnsemble_V();

				xMC.assign(ensembleResultsX.begin(),ensembleResultsX.end());
				yMC.assign(ensembleResultsY.begin(),ensembleResultsY.end());
				zMC.assign(ensembleResultsZ.begin(),ensembleResultsZ.end());
				eMC.assign(ensembleResults.begin(),ensembleResults.end());
				//SJZ
				hMC.assign(ensembleResultsH.begin(),ensembleResultsH.end());
				vMC.assign(ensembleResultsV.begin(),ensembleResultsV.end());
			}
		}
		bathyGrid->clear();
		#pragma endregion


		xMeshVectorMC = vector<double>((*xInterpVector));
		yMeshVectorMC = vector<double>((*yInterpVector));

		//************************************************************************************
		// IV. If we are computing the data points then run scalecInterpTile.  Otherwise go the longer route of scalecInterp.
		//************************************************************************************
		if (additionalOptions.find("-preInterpolatedLocations")->second == 0)
		{
			//A. Reassign for manipulation
			xMeshGridMC = dgrid((*xMeshGrid));
			yMeshGridMC = dgrid((*yMeshGrid));
			xtMC = vector<double>((*xSingleVector));
			ytMC = vector<double>((*ySingleVector));

			//returnValue = bathyTool_ORIGINAL(&xMC, &yMC, &zMC, &eMC, &xMeshGridMC, &yMeshGridMC, &xtMC, &ytMC, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, 1.0, true, true, NEITOL, &xyzOut);
			returnValue = bathyTool(&xMC, &yMC, &zMC, &eMC, &hMC, &vMC, &xMeshGridMC, &yMeshGridMC, &xtMC, &ytMC, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, 1.0, true, true, NEITOL, &xyzOut);

			//B. Clear so it can be used again
			xMeshGridMC.clear();
			yMeshGridMC.clear();
			xtMC.clear();
			ytMC.clear();
		}
		else
		{
			//returnValue = bathyToolPreDefined_ORIGINAL(&xMC, &yMC, inputDataZ, inputDataE, &xMeshVectorMC, &yMeshVectorMC, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, NEITOL, &xyzOut);
			returnValue = bathyToolPreDefined(&xMC, &yMC, &zMC, &eMC, &hMC, &vMC, &xMeshVectorMC, &yMeshVectorMC, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, NEITOL, &xyzOut);
		}

		//C. Clear up the variables if something went wrong
		if(returnValue != 0)
		{
			xyzOut.depth.clear();
			xyzOut.error.clear();
			xyzOut.nEi.clear();
			xyzOut.rEi.clear();
			
			if(abs(additionalOptions.find("-propUncert")->second) == 1)
			{
				xyzOut.depth0.clear();
				xyzOut.error0.clear();
			}
			if(abs(additionalOptions.find("-kalman")->second) == 1)
			{
				xyzOut.depthK.clear();
				xyzOut.errorK.clear();
			}
			return returnValue;
		}

		//************************************************************************************
		// V. Convert from UTM to Lat/Lon so we can output the data.
		//************************************************************************************
		double newLat, newLon;
		if ((additionalOptions.find("-inputInMeters")->second == 0) || (additionalOptions.find("-outputRasterFile")->second == 0))
		{
			for(i = 0; i < (const int)xMeshVectorMC.size(); i++){
				UTMEasting = (xMeshVectorMC[i])*cos(deg2rad*(-rotAngle)) - (yMeshVectorMC[i])*sin(deg2rad*(-rotAngle));
				UTMNorthing = (xMeshVectorMC[i])*sin(deg2rad*(-rotAngle)) + (yMeshVectorMC[i])*cos(deg2rad*(-rotAngle));

				UTMNorthing += UTMNorthingRef;
				UTMEasting += UTMEastingRef;

				UTMtoLL(RefEllip, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
				yMeshVectorMC[i] = newLat;
				xMeshVectorMC[i] = newLon;
			}
		}else if (additionalOptions.find("-outputRasterFile")->second == 1)
		{
			//1. Remove the central point offset from the data for raster format grids
			for(i = 0; i < (const int)xMeshVectorMC.size(); i++){
				UTMEasting = (xMeshVectorMC[i])*cos(deg2rad*(-rotAngle)) - (yMeshVectorMC[i])*sin(deg2rad*(-rotAngle));
				UTMNorthing = (xMeshVectorMC[i])*sin(deg2rad*(-rotAngle)) + (yMeshVectorMC[i])*cos(deg2rad*(-rotAngle));

				UTMNorthing += UTMNorthingRef;
				UTMEasting += UTMEastingRef;

				xMeshVectorMC[i] = (UTMEasting)*cos(deg2rad*rotAngle) - (UTMNorthing)*sin(deg2rad*rotAngle);
				yMeshVectorMC[i] = (UTMEasting)*sin(deg2rad*rotAngle) + (UTMNorthing)*cos(deg2rad*rotAngle);
			}
		}

		//************************************************************************************
		// VI. Write the data to an output file and clean up.
		//************************************************************************************
		if (additionalOptions.find("-outputRasterFile")->second == 0)
		{
			//A. Write the data to a flat file
			writeFile(fileName, &xMeshVectorMC, &yMeshVectorMC, &xyzOut.depth, &xyzOut.error, &xyzOut.nEi, &xyzOut.rEi, &additionalOptions);
		}else
		{
			//B. Write the data to a raster file
			double xSize = (*xMeshGrid).rows();
			double ySize = (*xMeshGrid).cols();
			writeRasterFile(fileName, &xMeshVectorMC, &yMeshVectorMC, &xyzOut.depth, &xyzOut.error, &xyzOut.nEi, &xyzOut.rEi, gridSpacingX, gridSpacingY, xSize, ySize, &additionalOptions, UTMZoneRef);
		}

		//C. Do the standard deviation stuff
		for (i = 0; i < (int) xyzOut.depth.size(); i++)
		{
			ziToStandardDeviation[stdLoc] = xyzOut.depth[i];
			eiToStandardDeviation[stdLoc] = xyzOut.error[i];
			neiToStandardDeviation[stdLoc] = xyzOut.nEi[i];
			reiToStandardDeviation[stdLoc] = xyzOut.rEi[i];
			stdLoc++;
		}

		//D. Clean up
		cout << "Done Creating Output File" << endl;
		xyzOut.depth.clear();
		xyzOut.error.clear();
		xyzOut.nEi.clear();
		xyzOut.rEi.clear();
		xMeshVectorMC.clear();
		yMeshVectorMC.clear();
	}

	//Take the standard deviation
	double ziSTD = standardDeviation(&ziToStandardDeviation);
	double eiSTD = standardDeviation(&eiToStandardDeviation);
	double neiSTD = standardDeviation(&neiToStandardDeviation);
	double reiSTD = standardDeviation(&reiToStandardDeviation);

	printf("Computed Standard Deviation:\n");
	printf("\tDepth: %.8f\n", ziSTD);
	printf("\tError: %.8f\n", eiSTD);
	printf("\tNormalized Error: %.8f\n", neiSTD);
	printf("\tResidual Error: %.8f\n", reiSTD);

	xMC.clear();
	yMC.clear();
	zMC.clear();
	eMC.clear();
	hMC.clear();
	vMC.clear();

	return returnValue;
}








int run(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, dgrid *xMeshGrid, dgrid *yMeshGrid, vector<double> *xSingleVector, vector<double> *ySingleVector, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double smoothingScaleX, double smoothingScaleY, double &x0, double &y0, double &x1, double &y1, double meanXSingle, double meanYSingle, string &kernelName,  map<string, int> additionalOptions, string outputFileName, double subDataMulitplier, double UTMNorthingRef, double UTMEastingRef, double rotAngle, int RefEllip, char UTMZoneRef[4], int numMCRuns, MB_ZGRID_DATA *MB_ZGridInput, GMT_SURFACE_DATA *GMTSurfaceInput, ALG_SPLINE_DATA *ALGSplineInput, Bathy_Grid* bathyGrid, int USE_UTM)
{
	//MergeBathy used to be able to only run one gridding algorithm at a time.  
	//These gridding algorithms spline sparse data in order to obtain more 
	//values before continuing.
	//If multiple gridding algorithms have 2/-2 USAGE values
	//set, the new up-sampled pre-splined values replace the input
	//x,y,z,e.  These values would then be input in the next gridding algorithm
	//and replaced with those new up-sampled pre-splined values, etc. before continuing on
	//to mergeBathy.
	//This is no longer the case.  Now MergeBathy will send the original input
	//to each spline routine for separate calculations and then combine them with
	//into an ensemble average afterwards.
	int sumG = additionalOptions.find("-ZGrid")->second + additionalOptions.find("-GMTSurface")->second
		+ additionalOptions.find("-ALGSpline")->second;// + additionalOptions.find("-Sibson_Natural_Neighbor_Interp")->second;
	
	if(sumG > 1)
	{
		cout << "Ensembling of Pre-Spliners has not been tested and is merely a prototype!" << endl;
		cout << "It is advised that you exit, change your arguments to select one and run again." << endl;
		//return !SUCCESS;
	}

	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int returnValue = SUCCESS;
	int i, stdLoc;
	OUTPUT_DATA xyzOut;

	double UTMEasting;
	double UTMNorthing;
	double newLon;
	double newLat;

	RNG randomNumber;
	stringstream ss;
	string numStr;
	string fileName, outputFileNameT;
	string extInterpFileName;

	vector<double> xMC = vector<double>((*inputDataX).size());
	vector<double> yMC = vector<double>((*inputDataY).size());
	vector<double> zMC = vector<double>((*inputDataX).size());
	vector<double> eMC = vector<double>((*inputDataY).size());
	//SJZ
	vector<double> hMC = vector<double>((*inputDataY).size());
	vector<double> vMC = vector<double>((*inputDataY).size());

	dgrid xMeshGridMC;
	dgrid yMeshGridMC;
	vector<double> xtMC;
	vector<double> ytMC;
	vector<double> xMeshVectorMC;
	vector<double> yMeshVectorMC;
		
	string ensembleFName = outputFileName;
	ensembleFName.append("_Ensemble.txt");
	dvector ziToStandardDeviation;
	dvector eiToStandardDeviation;
	dvector neiToStandardDeviation;
	dvector reiToStandardDeviation;
	dvector zi0ToStandardDeviation;
	dvector ei0ToStandardDeviation;
	dvector ziKToStandardDeviation;
	dvector eiKToStandardDeviation;
	bool MCFlag = true;
	if(numMCRuns == -1)
	{
		MCFlag = false;
		numMCRuns = 1;
	}
	else{
		ziToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
		eiToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
		neiToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
		reiToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
		zi0ToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
		ei0ToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
		ziKToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
		eiKToStandardDeviation = dvector((*xInterpVector).size()*numMCRuns);
	}
	//************************************************************************************
	// I. Iterate over the number of Monte Carlo runs to be done
	//************************************************************************************
	stdLoc = 0;
	for (int mcRunNum = 0; mcRunNum < numMCRuns; mcRunNum++)
	{
		if(MCFlag)
		{
			printf("\n\nEntering Monte Carlo Run Number: %d\n",(mcRunNum + 1));
		
			//A. Preinit the file numbers
			ss << mcRunNum;
			ss >> numStr;
			ss.clear();

			outputFileNameT = outputFileName;

			//B. Depending on the number of data runs pre pad the output file names with the run number
			if (mcRunNum < 10)
				outputFileNameT.append("_00");
			else if(mcRunNum < 100)
				outputFileNameT.append("_0");
			else
				outputFileNameT.append("_");
			outputFileNameT.append(numStr);
			//Done with preInit
		}else 
			outputFileNameT = outputFileName;

		//C. Set disposable variables for each loop
		for (i = 0; i < (const int)xMC.size(); i++){
			if(MCFlag){
				xMC[i] = (*inputDataX)[i] + (*inputDataHErr)[i]*randomNumber.normal();
				yMC[i] = (*inputDataY)[i] + (*inputDataHErr)[i]*randomNumber.normal();
			}else{
				xMC[i] = (*inputDataX)[i];
				yMC[i] = (*inputDataY)[i];
			}
			zMC[i] = (*inputDataZ)[i];
			eMC[i] = (*inputDataE)[i];
			hMC[i] = (*inputDataHErr)[i];
			vMC[i] = (*inputDataVErr)[i];
		}
		
		if(MCFlag)
		{
			/*The ensembler does not work with MC!!! The grids are unable to align!!!
			This would crash because perturbing the input x,y changes the mins and maxes upon which MBZ and GMT build 
			their grids. This cause grids that cannot be aligned to generate. 
			To fix this, it needs to be determined whether MC is happening at the proper place.  If x,y were perturbed closer towards initial ingestion, there would be no problem as all the manipulations to x,y and variables computed from x,y would then take place and the subsequent grids and variables generated would follow the same conventions of a single run. It may be decided though that perturbing x,y here was done deliberately, and if this is the case, either recalculate the mins, maxes or modify how GMT and MBZ build their grids so that they may align.  For example, MBZ will build its grid from the original mins, maxes but this may be determined to be incorrect implementation.  Conversely, this maybe the correct implementation and GMT should be passed the values and them to build its grid.  --I was not able to investigate into why GMT and MBZ built their grids differently.
			We will assume for now that this implementation is correct as is and will quit the ensembler.
			*/
			//cout << "Ensembling with Monte Carlo does not work!" << endl;
			//Find new mins and maxes, gets the floors and ceilings similar to what is done if not perturbing.
			/*x0 = floor(*max((*inputDataX).begin(),(*inputDataX).end()));
			x1 = ceil(*min((*inputDataX).begin(),(*inputDataX).end()));
			y0 = floor(*max((*inputDataY).begin(),(*inputDataY).end()));
			y1 = ceil(*min((*inputDataY).begin(),(*inputDataY).end()));*/

		}


		#pragma region -- Gridding
		if(bathyGrid->GriddingFlag)
		{
			//We intend to ensemble instead of re-splining already splined results (iteratively splining).
			vector<double> copyInputX;
			vector<double> copyInputY;
			vector<double> copyInputZ;
			vector<double> copyInputE;
			vector<double> copyInputH;
			vector<double> copyInputV;
			if(sumG > 1)
			{
				copyInputX.assign(xMC.begin(),xMC.end());
				copyInputY.assign(yMC.begin(),yMC.end());
				copyInputZ.assign(zMC.begin(),zMC.end());
				copyInputE.assign(eMC.begin(),eMC.end());
				copyInputH.assign(hMC.begin(),hMC.end());
				copyInputV.assign(vMC.begin(),vMC.end());
			}

			if(MCFlag)
			{
				//Rebuild if we perturbed our x,y inputs
				bathyGrid->clear();
				bathyGrid->Construct_Tin(&xMC, &yMC, &zMC, &hMC, &vMC);
			}

			//************************************************************************************
			// II. Check and do MB_ZGrid if necessary
			//************************************************************************************
			if (additionalOptions.find("-ZGrid")->second == 1)
			{
				#pragma region --MBZ
				cout << "********************************************************" << endl;
				cout << "* Computing MB_ZGrid" << endl;
				cout << "********************************************************" << endl;

				//A. Set iterative filenames
				extInterpFileName = (*MB_ZGridInput).z_OutputFileName;
				if(MCFlag){
					if (mcRunNum < 10)
						extInterpFileName.append("_00");
					else if(mcRunNum < 100)
						extInterpFileName.append("_0");
					else
						extInterpFileName.append("_");
					extInterpFileName.append(numStr);
				}

				//B. Now call mb_zgrid
				externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);
				//bool mbzReturn = extInterp.run_MB_ZGrid_ORIGINAL(&xMC, &yMC, &zMC, &eMC, inputDataHErr, inputDataVErr, x0, y0, x1, y1, (*MB_ZGridInput).spacingX, (*MB_ZGridInput).spacingY, (*MB_ZGridInput).tension, extInterpFileName, (*MB_ZGridInput).usage);
				bool mbzReturn = extInterp.run_MB_ZGrid(&xMC, &yMC, &zMC, &eMC, &hMC, &vMC, x0, y0, x1, y1, additionalOptions, (*MB_ZGridInput).spacingX, (*MB_ZGridInput).spacingY, (*MB_ZGridInput).tension, extInterpFileName, (*MB_ZGridInput).usage, bathyGrid);

				if (!mbzReturn)
				{
					cerr << "MB_ZGrid FAILED! ABORTING!!!" << endl;
					return MB_ZGRID_ERROR;
				}
				cout << "Done Computing MB_ZGrid" << endl << endl;
				/*if(abs((*MB_ZGridInput).usage) == 1)
					return returnValue;*/
				if(sumG > 1)
				{
					xMC.assign(copyInputX.begin(),copyInputX.end());
					yMC.assign(copyInputY.begin(),copyInputY.end());
					zMC.assign(copyInputZ.begin(),copyInputZ.end());
					eMC.assign(copyInputE.begin(),copyInputE.end());
					hMC.assign(copyInputH.begin(),copyInputH.end());
					vMC.assign(copyInputV.begin(),copyInputV.end());
				}
				#pragma endregion
			}
			//Done with MB_ZGrid

			//************************************************************************************
			// III. Check and do GMT Surface if necessary
			//************************************************************************************
			if (additionalOptions.find("-GMTSurface")->second == 1)
			{
				#pragma region --GMTSurface
				cout << "********************************************************" << endl;
				cout << "* Computing GMT Surface" << endl;
				cout << "********************************************************" << endl;

				//A. Set iterative filenames
				extInterpFileName = (*GMTSurfaceInput).z_OutputFileName;
				if(MCFlag){
					if (mcRunNum < 10)
						extInterpFileName.append("_00");
					else if(mcRunNum < 100)
						extInterpFileName.append("_0");
					else
						extInterpFileName.append("_");
					extInterpFileName.append(numStr);
				}

				//B. Now call GMT Surface
				externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);
				//bool mbzReturn = extInterp.run_Surface_ORIGINAL(&xMC, &yMC, &zMC, &eMC, inputDataHErr, inputDataVErr, (*GMTSurfaceInput).spacingX, (*GMTSurfaceInput).spacingY, (*GMTSurfaceInput).tension, extInterpFileName, (*GMTSurfaceInput).scaleFactor, (*GMTSurfaceInput).alpha, (*GMTSurfaceInput).usage);
				bool gmtReturn = extInterp.run_Surface(&xMC, &yMC, &zMC, &eMC, &hMC, &vMC, additionalOptions, (*GMTSurfaceInput).spacingX, (*GMTSurfaceInput).spacingY, (*GMTSurfaceInput).tension, extInterpFileName, (*GMTSurfaceInput).scaleFactor, (*GMTSurfaceInput).alpha, (*GMTSurfaceInput).usage, bathyGrid);

				if (!gmtReturn)
				{
					cerr << "GMT Surface FAILED! ABORTING!!!" << endl;
					return GMT_SURF_ERROR;
				}
				cout << "Done Computing GMT Surface" << endl << endl;
				/*if(abs((*GMTSurfaceInput).usage) == 1)
					return returnValue;*/
				if(sumG > 1)
				{
					xMC.assign(copyInputX.begin(),copyInputX.end());
					yMC.assign(copyInputY.begin(),copyInputY.end());
					zMC.assign(copyInputZ.begin(),copyInputZ.end());
					eMC.assign(copyInputE.begin(),copyInputE.end());
					hMC.assign(copyInputH.begin(),copyInputH.end());
					vMC.assign(copyInputV.begin(),copyInputV.end());
				}
				#pragma endregion
			}
			//Done with GMT Surface
		
			//************************************************************************************
			// III. Check and do ALGSpline Surface if necessary
			//************************************************************************************
			if (additionalOptions.find("-ALGSpline")->second == 1)
			{
				#pragma region -- ALGSpline
				cout << "********************************************************" << endl;
				cout << "* Computing ALG Spline" << endl;
				cout << "********************************************************" << endl;
				//A. Set iterative filenames
				extInterpFileName = (*GMTSurfaceInput).z_OutputFileName;
				if(MCFlag){
					if (mcRunNum < 10)
						extInterpFileName.append("_00");
					else if(mcRunNum < 100)
						extInterpFileName.append("_0");
					else
						extInterpFileName.append("_");
					extInterpFileName.append(numStr);
				}

				double ySize = (*xMeshGrid).rows();
				double xSize = (*xMeshGrid).cols();
				externalInterpolators extInterp(UTMNorthingRef, UTMEastingRef, rotAngle, RefEllip, UTMZoneRef);

				//Must accept a regular grid.  Problem sam.
				bool algReturn = extInterp.run_ALGSpline(inputDataX, inputDataY, inputDataZ, x0,  x1, y0, y1, xSize, ySize, &eMC, &hMC, &vMC, additionalOptions, (*ALGSplineInput).spacingX, (*ALGSplineInput).spacingY, (*ALGSplineInput).tension, (*ALGSplineInput).z_OutputFileName, (*ALGSplineInput).scaleFactor, (*ALGSplineInput).alpha, (*ALGSplineInput).usage, bathyGrid);

				if (!algReturn)
				{
					cerr << "ALG Spline FAILED! ABORTING!!!" << endl;
					return ALG_SPLINE_ERROR;
				}
				cout << "Done Computing ALG Spline." << endl << endl;
				//if(abs((*ALGSplineInput).usage) == 1)
				//	return returnValue;
				if(sumG > 1)
				{
					xMC.assign(copyInputX.begin(),copyInputX.end());
					yMC.assign(copyInputY.begin(),copyInputY.end());
					zMC.assign(copyInputZ.begin(),copyInputZ.end());
					eMC.assign(copyInputE.begin(),copyInputE.end());
					hMC.assign(copyInputH.begin(),copyInputH.end());
					vMC.assign(copyInputV.begin(),copyInputV.end());
				}
				#pragma endregion
			}

			if(bathyGrid->getGrids().size() > 1){
				//An interpGrid is not created when the spline-gridder output is not used as input (-1/1).
				//Therefore, only the ones with (-2/2) set will be ensembled.
				cout << "********************************************************" << endl;
				cout << "* Ensembling pre-splined grids to use as input" << endl;
				cout << "********************************************************" << endl;
				bool ensembleReturn = bathyGrid->ensemble();
				if (!ensembleReturn)
				{
					if(MCFlag)
						cout << "Ensembling with Monte Carlo does not work!" << endl;
					cerr << "Ensemble FAILED! ABORTING!!!" << endl;
					return !SUCCESS;
				}
				bathyGrid->printensemble(ensembleFName);
				vector<double> ensembleResultsX = bathyGrid->getEnsemble_X();
				vector<double> ensembleResultsY = bathyGrid->getEnsemble_Y();
				vector<double> ensembleResultsZ = bathyGrid->getEnsemble_Z();
				vector<double> ensembleResults  = bathyGrid->getEnsemble_U();
				vector<double> ensembleResultsH = bathyGrid->getEnsemble_H();
				vector<double> ensembleResultsV = bathyGrid->getEnsemble_V();

				xMC.assign(ensembleResultsX.begin(),ensembleResultsX.end());
				yMC.assign(ensembleResultsY.begin(),ensembleResultsY.end());
				zMC.assign(ensembleResultsZ.begin(),ensembleResultsZ.end());
				eMC.assign(ensembleResults.begin(),ensembleResults.end());
				//SJZ
				hMC.assign(ensembleResultsH.begin(),ensembleResultsH.end());
				vMC.assign(ensembleResultsV.begin(),ensembleResultsV.end());
			}
		}
		bathyGrid->clear();
		#pragma endregion

		xMeshVectorMC = vector<double>((*xInterpVector));
		yMeshVectorMC = vector<double>((*yInterpVector));

		//************************************************************************************
		// IV. If we are computing the data points then run scalecInterpTile.  Otherwise go the longer route of scalecInterp.
		//************************************************************************************
		if (additionalOptions.find("-preInterpolatedLocations")->second == 0)
		{
			//A. Reassign for manipulation
			xMeshGridMC = dgrid((*xMeshGrid));
			yMeshGridMC = dgrid((*yMeshGrid));
			xtMC = vector<double>((*xSingleVector));
			ytMC = vector<double>((*ySingleVector));

			//returnValue = bathyTool_ORIGINAL(&xMC, &yMC, &zMC, &eMC, &xMeshGridMC, &yMeshGridMC, &xtMC, &ytMC, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, 1.0, true, true, NEITOL, &xyzOut);
			returnValue = bathyTool(&xMC, &yMC, &zMC, &eMC, &hMC, &vMC, &xMeshGridMC, &yMeshGridMC, &xtMC, &ytMC, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, 1.0, true, true, NEITOL, &xyzOut);

			//B. Clear so it can be used again
			xMeshGridMC.clear();
			yMeshGridMC.clear();
			xtMC.clear();
			ytMC.clear();
		}
		else
		{
			//returnValue = bathyToolPreDefined_ORIGINAL(&xMC, &yMC, inputDataZ, inputDataE, &xMeshVectorMC, &yMeshVectorMC, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, NEITOL, &xyzOut);
			returnValue = bathyToolPreDefined(&xMC, &yMC, &zMC, &eMC, &hMC, &vMC, &xMeshVectorMC, &yMeshVectorMC, smoothingScaleX, smoothingScaleY, x0, y0, meanXSingle, meanYSingle, kernelName, additionalOptions, NEITOL, &xyzOut);
		}

		//C. Clear up the variables if something went wrong
		if(returnValue != 0)
		{
			xyzOut.depth.clear();
			xyzOut.error.clear();
			xyzOut.nEi.clear();
			xyzOut.rEi.clear();
			
			if(abs(additionalOptions.find("-propUncert")->second) == 1)
			{
				xyzOut.depth0.clear();
				xyzOut.error0.clear();
			}
			if(abs(additionalOptions.find("-kalman")->second) == 1)
			{
				xyzOut.depthK.clear();
				xyzOut.errorK.clear();
			}
			return returnValue;
		}

		for(int k = 0; k < (const int)xyzOut.nEi.size(); k++)
		{
			xyzOut.error[k] = sqrt(xyzOut.error[k]);
			xyzOut.nEi[k] = sqrt(xyzOut.nEi[k]);
			xyzOut.rEi[k] = sqrt(xyzOut.rEi[k]);
			xyzOut.standardDev[k] = sqrt(xyzOut.standardDev[k]);
			if(abs(additionalOptions.find("-propUncert")->second) == 1)
			{
				xyzOut.error0[k] = sqrt(xyzOut.error0[k]);
			}
		}

		vector<double> xInterpTemp;
		vector<double> yInterpTemp;
		vector<double> nEiTemp;
		vector<double> rEiTemp;
		vector<double> depthTemp;
		vector<double> errorTemp;
		vector<double> depthKTemp;
		vector<double> errorKTemp;
		vector<double> depth0Temp;
		vector<double> error0Temp;

		//same as above but creates new vars so both raster and txt can be printed for debugging
		/*vector<double> alldepthTemp(xyzOut.depth);
		vector<double> alldepth0Temp(xyzOut.depth0);
		vector<double> alldepthKTemp(xyzOut.depthK);
		vector<double> allerrorTemp(xyzOut.error);
		vector<double> allerror0Temp(xyzOut.error0);
		vector<double> allerrorKTemp(xyzOut.errorK);
		vector<double> allnEiTemp(xyzOut.nEi);
		vector<double> allrEiTemp(xyzOut.rEi);*/

		// Post-processing
		double nescale;
		vector<double>::iterator temp;
		vector<double> bar;
		vector<double>::iterator it;
		if(additionalOptions.find("-outputRasterFile")->second == 0 && additionalOptions.find("-outputBagFile")->second == 0)
		{
			for(it = xyzOut.nEi.begin(); it < xyzOut.nEi.end(); it++)
			{
				if(*it != NaN)
					bar.push_back(*it);
			}

			 nescale = median(bar) + 0.1;

			if (nescale > 0.95)
				nescale = 0.95;

			// removes edge extrapolation effects
			nescale = max(nescale,0.75);
			temp = std::min_element(xyzOut.nEi.begin(),xyzOut.nEi.end());
			if (*temp > nescale)
				nescale = *temp + 0.1;

			//remove edge extrapolation effects
			if(additionalOptions.find("-modelflag")->second == 0)
			{
				if(additionalOptions.find("-nonegdepth")->second==0)
				{
					//remove only edge extrapolation effects 
					for(int cnt = 0; cnt < (const int)xyzOut.nEi.size(); cnt++)
					{
						if(xyzOut.nEi[cnt] > nescale )
						{ //really poor samples
							xyzOut.depth[cnt] = NaN;
							xyzOut.error[cnt] = NaN;
							xyzOut.nEi[cnt] = NaN;
							xyzOut.rEi[cnt] = NaN;
							xyzOut.depth0[cnt] = NaN;
							xyzOut.error0[cnt] = NaN;
							xyzOut.depthK[cnt] = NaN;
							xyzOut.errorK[cnt] = NaN;
						}
					}
				}
				else
				{
					//remove edge extrapolation effects and negative depths
					for(int cnt = 0; cnt < (const int)xyzOut.nEi.size(); cnt++)
					{
						if(xyzOut.nEi[cnt] > nescale || xyzOut.depth[cnt] < 0.1)//SJZ changed it from && 1/21/16
						{//This is weird... I think it should be || instead of &&
							xyzOut.depth[cnt] = NaN;
							xyzOut.error[cnt] = NaN;
							xyzOut.nEi[cnt] = NaN;
							xyzOut.rEi[cnt] = NaN;
							xyzOut.depth0[cnt] = NaN;
							xyzOut.error0[cnt] = NaN;
							xyzOut.depthK[cnt] = NaN;
							xyzOut.errorK[cnt] = NaN;
						}
					}
				}
				cout << "Points with normalized errors greater than "<< nescale << " removed." << endl;
			}
			else
				//Model Mode; Keep all points
				cout << "Nothing is bad!" << endl;

			//Keep only data points for regular output; use smaller data set
			for(int i = 0; i < (const int)(xMeshVectorMC).size(); i++)
			{
				if(xyzOut.depth[i] != NaN)
					xInterpTemp.push_back((xMeshVectorMC)[i]);
				if(xyzOut.depth[i] != NaN)
					yInterpTemp.push_back((yMeshVectorMC)[i]);
				if(xyzOut.depth[i] != NaN)
					depthTemp.push_back(xyzOut.depth[i]);
				if(xyzOut.error[i] != NaN)
					errorTemp.push_back(xyzOut.error[i]);
				if(xyzOut.nEi[i] != NaN)
					nEiTemp.push_back(xyzOut.nEi[i]);
				if(xyzOut.rEi[i] != NaN)
					rEiTemp.push_back(xyzOut.rEi[i]);
				if(xyzOut.depth0[i] != NaN)
					depth0Temp.push_back(xyzOut.depth0[i]);
				if(xyzOut.error0[i] != NaN)
					error0Temp.push_back(xyzOut.error0[i]);
				if(xyzOut.depthK[i] != NaN)
					depthKTemp.push_back(xyzOut.depthK[i]);
				if(xyzOut.errorK[i] != NaN)
					errorKTemp.push_back(xyzOut.errorK[i]);
			}
		}
		else 
		{	//Raster output
			xInterpTemp		= xMeshVectorMC;
			yInterpTemp		= yMeshVectorMC;
			depthTemp		= xyzOut.depth;
			depth0Temp		= xyzOut.depth0;
			depthKTemp		= xyzOut.depthK;
			errorTemp		= xyzOut.error;
			error0Temp		= xyzOut.error0;
			errorKTemp		= xyzOut.errorK;
			nEiTemp			= xyzOut.nEi;
			rEiTemp			= xyzOut.rEi;
		}

		vector<double> xMeshVector;
		vector<double> yMeshVector;
		vector<double> depthTemp2;
		vector<double> errorTemp2;
		vector<double> neiTemp2;
		vector<double> reiTemp2;
		vector<double> depth0Temp2;
		vector<double> error0Temp2;
		vector<double> depthKTemp2;
		vector<double> errorKTemp2;

	//	vector<double> UTMNorthings((*xInterpVector).size());
	//	vector<double> UTMEastings((*xInterpVector).size());
		vector<double> UTMNorthings(xInterpTemp.size());
		vector<double> UTMEastings(xInterpTemp.size());
		//vector<double> xInterpTemp2(xInterpTemp.size());
		//vector<double> yInterpTemp2(xInterpTemp.size());

		//SJZ clean up to try to run FOUO on 32bit**********************************
		//vector<double> xInterpTemp3;
		//vector<double> yInterpTemp3;
		//vector<double> depthTemp3;
		//vector<double> errorTemp3;
		vector<double> neiTemp3;
		//vector<double> reiTemp3;
		//vector<double> depth0Temp3;
		//vector<double> error0Temp3;
		//vector<double> depthKTemp3;
		//vector<double> errorKTemp3;
		//	
		//vector<double> xInterpTemp4;
		//vector<double> yInterpTemp4;
		//vector<double> depthTemp4;
		//vector<double> errorTemp4;
		//vector<double> neiTemp4;
		//vector<double> reiTemp4;
		//vector<double> depth0Temp4;
		//vector<double> error0Temp4;
		//vector<double> depthKTemp4;
		//vector<double> errorKTemp4;
		//end clean up
		int xSize = 0, ySize = 0;
		vector<double> xt;
		vector<double> yt;
		vector<double> x;
		vector<double> y;
		dgrid newxx;
		dgrid newyy;
		vector<double> newxxVec;
		vector<double> newyyVec;
		//************************************************************************************
		// IV. Convert from UTM to Lat/Lon so we can output the data.
		//************************************************************************************
		if ((additionalOptions.find("-inputInMeters")->second == 0) && (additionalOptions.find("-outputRasterFile")->second == 0 && additionalOptions.find("-outputBagFile")->second == 0))
		{
			for(int i = 0; i < (const int)xInterpTemp.size(); i++)
			{
				if(!USE_UTM) // SJZ
				{
					UTMEasting = (xInterpTemp[i])*cos(deg2rad*(-rotAngle)) - (yInterpTemp[i])*sin(deg2rad*(-rotAngle));
					UTMNorthing = (xInterpTemp[i])*sin(deg2rad*(-rotAngle)) + (yInterpTemp[i])*cos(deg2rad*(-rotAngle));

					UTMNorthing += UTMNorthingRef;
					UTMEasting += UTMEastingRef;
				}
				else
				{ // SJZ
					UTMEasting	= xInterpTemp[i];
					UTMNorthing = yInterpTemp[i];
				}
				UTMtoLL(RefEllip, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
				if(newLon >= 180) // SJZ
					newLon -= 360;
				yInterpTemp[i] = newLat;
				xInterpTemp[i] = newLon;
			}
			if (additionalOptions.find("-printMatlabMatch")->second == 1)
			{
				if (((rotAngle > 45) & (rotAngle <= 135)) | ((rotAngle > 225) & (rotAngle <= 315)))
					cerr << "Rotation Angle must be (<= 45 & >315) OR (> 135 & <= 225) to match Matlab output otherwise the data is transposed!" << endl;
			}
		}
		else if (additionalOptions.find("-outputRasterFile")->second == 1 || additionalOptions.find("-outputBagFile")->second == 1)
		{ 
			dgrid xx;
			dgrid yy;
			vector<double> zz;//, zz0, zzK, ee, ee0, eeK, neei, reei; // SJZ
			double minNorthing	= (double)MAX_INT;
			double minEasting	= (double)MAX_INT;
			double maxNorthing	= (double)MIN_INT;
			double maxEasting	= (double)MIN_INT;
			/*double minXInterp	= (double)MAX_INT;
			double minYInterp	= (double)MAX_INT;
			double maxXInterp	= (double)MIN_INT;
			double maxYInterp	= (double)MIN_INT;*/
			int INIT_TEMPS = 1;

			//1. Remove the central point offset from the data for raster format grids
			//A new raster grid will be created by interpolating from current data to raster grid locations
			//Convert only points with data to UTM.
			//for(int i = 0; i < (const int)xInterpVector->size(); i++)
			for(int i = 0; i < (const int)xInterpTemp.size(); i++)
			{
				if(!USE_UTM) // SJZ
				{ 
					UTMEasting = (xInterpTemp)[i]*cos(deg2rad*(-rotAngle)) - (yInterpTemp)[i]*sin(deg2rad*(-rotAngle));
					UTMNorthing = (xInterpTemp)[i]*sin(deg2rad*(-rotAngle)) + (yInterpTemp)[i]*cos(deg2rad*(-rotAngle));

					//UTMEasting = (*xInterpVector)[i]*cos(deg2rad*(-rotAngle)) - (*yInterpVector)[i]*sin(deg2rad*(-rotAngle));
					//UTMNorthing = (*xInterpVector)[i]*sin(deg2rad*(-rotAngle)) + (*yInterpVector)[i]*cos(deg2rad*(-rotAngle));
			
					UTMNorthing += UTMNorthingRef;
					UTMEasting  += UTMEastingRef;
				}
				else
				{ // SJZ
					UTMEasting	= xInterpTemp[i];
					UTMNorthing = yInterpTemp[i];
				}
				if(INIT_TEMPS)
				{
					minNorthing	= UTMNorthing;
					minEasting	= UTMEasting;
					maxNorthing	= UTMNorthing;
					maxEasting	= UTMEasting;
				/*	minXInterp	= MAX_INT;
					minYInterp	= MAX_INT;
					maxXInterp	= MIN_INT;
					maxYInterp	= MIN_INT;*/
					INIT_TEMPS = 0;
				}
				UTMNorthings[i] = UTMNorthing;	//y
				UTMEastings [i] = UTMEasting;	//x
				//yInterpTemp2[i] = UTMNorthing;
				//xInterpTemp2[i] = UTMEasting;
				//(*yInterpVector)[i] = UTMNorthing;
				//(*xInterpVector)[i] = UTMEasting; 
				UTMtoLL(RefEllip, UTMNorthing, UTMEasting, UTMZoneRef, newLat, newLon);
				if(newLon >= 180) // SJZ
					newLon -= 360;
				yInterpTemp[i] = newLat;
				xInterpTemp[i] = newLon;

				minEasting  = min(minEasting, UTMEasting);
				maxEasting  = max(maxEasting, UTMEasting);
				minNorthing = min(minNorthing, UTMNorthing);
				maxNorthing = max(maxNorthing, UTMNorthing);
			}
		
			int Ni = (const int)UTMNorthings.size();
			string interpMethod;
			Bathy_Grid new_bathyGrid=Bathy_Grid();
				
		//	new_bathyGrid.Construct_TinRaster(&UTMEastings, &UTMNorthings, &alldepthTemp, &allerrorTemp, &allnEiTemp, &allrEiTemp, &alldepth0Temp, &allerror0Temp, &alldepthKTemp, &allerrorKTemp);
			new_bathyGrid.Construct_TinRaster(&UTMEastings, &UTMNorthings, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &depth0Temp, &error0Temp, &depthKTemp, &errorKTemp);//sjz

			//A. Obtain a Grid
			// Calculate x and y for our new grid where we go from
			createMeshXYDims(minEasting, minNorthing, maxEasting, maxNorthing, gridSpacingX, gridSpacingY, &xt, &yt);

			// get our new x,y dimensions
			xSize = (const int)xt.size();
			ySize = (const int)yt.size();
			
			// Resize vectors to match our new dimensions
			xMeshVector.resize(xSize*ySize, 0.00);
			yMeshVector.resize(xSize*ySize, 0.00);
			xx.resize(ySize, xSize, 0.00);
			yy.resize(ySize, xSize, 0.00);
		
			// Create grid with those dimensions
			createMeshGrid(&xt, &yt, &xMeshVector, &yMeshVector, &xx, &yy);
			
			// Reshape grid
			newxx.resize(xx.rows(), xx.cols(),0.00);
			newyy.resize(yy.rows(), yy.cols(),0.00);
			reshapeGrid(&xx, &newxx);
			reshapeGrid(&yy, &newyy);

			//B.Interpolate to get z, e, nmsei, rei, kz, kvar values at raster grid locations.
			InterpGrid new_gmt = InterpGrid(GMT);
			zz.resize((newxx.vec()).size(), 0.00);
			(new_gmt).estimate(&(newxx.vec()), &(newyy.vec()), &zz, 1.96, 2.00, 0, new_bathyGrid.getTin(), "bilinear", "bilinear", "none");

			depthTemp2 = *(new_gmt).getZ();
			errorTemp2 = (new_gmt).getE();
			depth0Temp2 = (new_gmt).getZ0();
			error0Temp2 = (new_gmt).getE0();
			depthKTemp2 = (new_gmt).getZK();
			errorKTemp2 = (new_gmt).getEK();
			neiTemp2 = (new_gmt).getNEI();
			reiTemp2 = (new_gmt).getREI();

			new_gmt.clear();
			new_bathyGrid.clear();
			xx.clear();
			yy.clear();
			xt.clear();
			yt.clear();
		
			newxxVec = newxx.vec();
			newyyVec = newyy.vec();
		
			//Keep only points in our range; set bad points to NaN
			for(int i = 0; i < (const int)(newxxVec).size(); i++)
			{
				if((newxxVec[i] < minEasting || newxxVec[i] > maxEasting) || (newyyVec[i] < minNorthing || newyyVec[i] > maxNorthing))
				{
					depthTemp2[i] = NaN;
					errorTemp2[i] = NaN;
					neiTemp2[i] = NaN;
					reiTemp2[i] = NaN;
					if(abs(additionalOptions.find("-propUncert")->second) == 1)
					{
						depth0Temp2[i] = NaN;
						error0Temp2[i] = NaN;
					}
					if(abs(additionalOptions.find("-kalman")->second) == 1)
					{
						depthKTemp2[i] = NaN;
						errorKTemp2[i] = NaN;
					}
				}
				else
				{
					/*xInterpTemp3.push_back(newxxVec[i]); 
					yInterpTemp3.push_back(newyyVec[i]); 
					depthTemp3.push_back(depthTemp2[i]); 
					errorTemp3.push_back(errorTemp2[i]);*/
					neiTemp3.push_back(neiTemp2[i]);
					/*reiTemp3.push_back(reiTemp2[i]); 
					if(abs(additionalOptions.find("-propUncert")->second) == 1)
					{
						depth0Temp3.push_back(depth0Temp2[i]);
						error0Temp3.push_back(error0Temp2[i]);
					}
					if(abs(additionalOptions.find("-kalman")->second) == 1)
					{
						depthKTemp3.push_back(depthKTemp2[i]);
						errorKTemp3.push_back(errorKTemp2[i]);
					}*/
				}
			}

			bar.clear();
			// Post-processing
			for(it = neiTemp3.begin(); it < neiTemp3.end(); it++)
			{
				if(*it != NaN)
					bar.push_back(*it);
			}

			nescale = median(bar) + 0.1;

			if (nescale > 0.95)
				nescale = 0.95;

			// removes edge extrapolation effects
			nescale = max(nescale,0.75);
			temp = std::min_element(neiTemp3.begin(),neiTemp3.end());
			if (*temp > nescale)
				nescale = *temp + 0.1;
	
			//Remove edge extrapolation effects
			if(additionalOptions.find("-modelflag")->second == 0)
			{
				if(additionalOptions.find("-nonegdepth")->second==0)
				{
					//remove only edge extrapolation effects
					for(int cnt = 0; cnt < (const int)neiTemp2.size(); cnt++)
					{
						if(neiTemp2[cnt] < 0 )
							int stop;
						if(neiTemp2[cnt] > nescale || neiTemp2[cnt] < 0 )
						{ //really poor samples
							depthTemp2[cnt] = NaN;
							errorTemp2[cnt] = NaN;
							//Keep all nei
							//neiTemp2[cnt] = NaN;
							reiTemp2[cnt] = NaN;
							if(abs(additionalOptions.find("-propUncert")->second) == 1)
							{
								depth0Temp2[cnt] = NaN;
								error0Temp2[cnt] = NaN;
							}
							if(abs(additionalOptions.find("-kalman")->second) == 1)
							{
								depthKTemp2[cnt] = NaN;
								errorKTemp2[cnt] = NaN;
							}
						}
					}
				}
				else
				{
					//remove edge extrapolation effects and negative depths
					for(int cnt = 0; cnt < (const int)neiTemp2.size(); cnt++)
					{ 
						if((neiTemp2[cnt] > nescale || neiTemp2[cnt] < 0) || depthTemp2[cnt] < 0.1)//SJZ changed it from && 1/21/16
						{
							depthTemp2[cnt] = NaN;
							errorTemp2[cnt] = NaN;
							//Keep all nei
							//neiTemp2[cnt] = NaN;
							reiTemp2[cnt] = NaN;
							if(abs(additionalOptions.find("-propUncert")->second) == 1)
							{
								depth0Temp2[cnt] = NaN;
								error0Temp2[cnt] = NaN;
							}
							if(abs(additionalOptions.find("-kalman")->second) == 1)
							{
								depthKTemp2[cnt] = NaN;
								errorKTemp2[cnt] = NaN;
							}
						}
					}
				}
				cout << "Points with normalized errors greater than "<< nescale << " removed." << endl;
			}
			else //Model Mode; Keep all points!
				cout << "Nothing is bad!" << endl;

			//for(int i = 0; i < (const int)(newxxVec).size(); i++) 
			//{
			//	if(depthTemp2[i] != NaN)
			//		xInterpTemp4.push_back((newxxVec)[i]);
			//	if(depthTemp2[i] != NaN)
			//		yInterpTemp4.push_back((newyyVec)[i]);
			//	if(depthTemp2[i] != NaN)
			//		depthTemp4.push_back(depthTemp2[i]);
			//	if(errorTemp2[i] != NaN)
			//		errorTemp4.push_back(errorTemp2[i]);
			//	//Keep all nei
			///*	if(neiTemp2[i] != NaN)
			//		neiTemp4.push_back(neiTemp2[i]);*/
			//	if(reiTemp2[i] != NaN)
			//		reiTemp4.push_back(reiTemp2[i]);
			//	if(abs(additionalOptions.find("-propUncert")->second) == 1)
			//	{
			//		if(depth0Temp2[i] != NaN)
			//			depth0Temp4.push_back(depth0Temp2[i]);
			//		if(error0Temp2[i] != NaN)
			//			error0Temp4.push_back(error0Temp2[i]);
			//	}
			//	if(abs(additionalOptions.find("-kalman")->second) == 1)
			//	{
			//		if(depthKTemp2[i] != NaN)
			//			depthKTemp4.push_back(depthKTemp2[i]);
			//		if(errorKTemp2[i] != NaN)
			//			errorKTemp4.push_back(errorKTemp2[i]);
			//	}
			//}
		}

		//************************************************************************************
		// VI. Write the data to an output file and clean up.
		//************************************************************************************
		string fileName = outputFileNameT;
		string fname;
		string f;

		if(additionalOptions.find("-appendFilename")->second == 1)
		{
			//size_t found = fileName.back();
			size_t found = *fileName.rbegin(); // UNIX
		
			//Print results to file and attach the window name used to filename.
			if(kernelName.compare("hanning") == 0 || kernelName.compare("hann") == 0)
				fname = fileName.substr(0,found+1).append("hann");
			else if (kernelName.compare("boxcar") == 0)
				fname = fileName.substr(0,found+1).append("boxcar");
			else if (kernelName.compare("quadloess") == 0)
				fname = fileName.substr(0,found+1).append("quadloess");
			else if (kernelName.compare("loess") == 0)
				fname = fileName.substr(0,found+1).append("loess");

			if (additionalOptions.find("-nnInterp")->second == 1)
				fname = fname.append("_NN");
		}
		else 
			fname = fileName;

		if (additionalOptions.find("-outputRasterFile")->second == 0 && additionalOptions.find("-outputBagFile")->second == 0)
		{
			//A. Write the data to a flat file
			if((additionalOptions.find("-mse")->second == 1 || additionalOptions.find("-propUncert")->second == 1) || additionalOptions.find("-kalman")->second == 1)
			{
				f = fname;
				writeFiles(f, &xInterpTemp, &yInterpTemp, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions, &depth0Temp, &error0Temp, &depthKTemp, &errorKTemp);
			}
			//i. Print MSE with Kalman results tagged on the end.
			if(additionalOptions.find("-printMSEwK")->second == 1)
			{
				f = fname;
				f = f.append("MSEwK");
				writeFile(f, &xInterpTemp, &yInterpTemp, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions, &depthKTemp, &errorKTemp);
			}
			//ii. Print matlab matched output;
			if (additionalOptions.find("-printMatlabMatch")->second == 1)
			{
				f = fname;
				f = f.append("Matlab");
				writeFileMatch(f, &xInterpTemp, &yInterpTemp, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions,&depthKTemp, &errorKTemp);
			}
			//iii. Each estimator prints to its own output file.
		/*	if(additionalOptions.find("-mse")->second == 1)
			{
				f = fname;
				writeFile(f, &xInterpTemp, &yInterpTemp, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions);
			}
			if(additionalOptions.find("-propUncert")->second == 1)
			{
				f = fname;
				f = f.append("P");
				writeFile(f, &xInterpTemp, &yInterpTemp, &depth0Temp, &error0Temp, &nEiTemp, &rEiTemp, &additionalOptions);
			}
			if(additionalOptions.find("-kalman")->second == 1)
			{
				f = fname;
				f = f.append("K");
				writeFile(f, &xInterpTemp, &yInterpTemp, &depthKTemp, &errorKTemp, &nEiTemp, &rEiTemp, &additionalOptions);
			}*/
		}
		else if (additionalOptions.find("-outputRasterFile")->second == 0)
		{
			//B. Write the data to a bag file
			xSize = newxx.rows();
			ySize = newxx.cols();

			f = fname;
			writeBagFile(f, &newxxVec, &newyyVec, &depthTemp2, &errorTemp2, &neiTemp2, &reiTemp2,(int)gridSpacingX, (int)gridSpacingY, xSize, ySize, &additionalOptions, UTMZoneRef, &depth0Temp2, &error0Temp2, &depthKTemp2, &errorKTemp2);

			//Debug by printing bag to raster
		//	f = fname;
		//	readBagFile(f, &newxxVec, &newyyVec, &depthTemp2, &errorTemp2, &neiTemp2, &reiTemp2, (int)gridSpacingX, (int)gridSpacingY, xSize, ySize, &additionalOptions);
		}
		else
		{
		
			//B. Write the data to a raster file
			xSize = newxx.rows();
			ySize = newxx.cols();

			////Print MSE raster files 
			////write mse file for debugging raster--delete when done
			//f = fname;
			//f = f.append("mse");
			////original mse values no nans
			////writeFile(f, &xInterpTemp2, &yInterpTemp2, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions);
			//writeFile(f, &UTMEastings, &UTMNorthings, &depthTemp, &errorTemp, &nEiTemp, &rEiTemp, &additionalOptions);
			//
			////f = fname;
			////f = f.append("mseALL");
			////writeFile(f, &UTMEastings, &UTMNorthings, &alldepthTemp, &allerrorTemp, &allnEiTemp, &allrEiTemp, &additionalOptions);
			//	
			////new raster values with before its been naned
			//f = fname;
			//f = f.append("mse3");
			//writeFile(f, &xInterpTemp3, &yInterpTemp3,  &depthTemp3, &errorTemp3, &neiTemp3, &reiTemp3, &additionalOptions);
			
			f = fname;
			writeRasterFile(f, &newxxVec, &newyyVec, &depthTemp2, &errorTemp2, &neiTemp2, &reiTemp2, gridSpacingX, gridSpacingY, xSize, ySize, &additionalOptions, UTMZoneRef, &depth0Temp2, &error0Temp2, &depthKTemp2, &errorKTemp2);

		}

		if(MCFlag)
		{
			//C. Do the standard deviation stuff
			for (i = 0; i < (int) xyzOut.depth.size(); i++)
			{
				ziToStandardDeviation[stdLoc] = xyzOut.depth[i];
				eiToStandardDeviation[stdLoc] = xyzOut.error[i];
				neiToStandardDeviation[stdLoc] = xyzOut.nEi[i];
				reiToStandardDeviation[stdLoc] = xyzOut.rEi[i];

				if(abs(additionalOptions.find("-propUncert")->second) == 1)
				{
					zi0ToStandardDeviation[stdLoc] = xyzOut.depth0[i];
					ei0ToStandardDeviation[stdLoc] = xyzOut.error0[i];
				}
				if(abs(additionalOptions.find("-kalman")->second) == 1 )
				{
					ziKToStandardDeviation[stdLoc] = xyzOut.depthK[i];
					eiKToStandardDeviation[stdLoc] = xyzOut.errorK[i];
				}
				stdLoc++;
			}
		}
		//D. Clean up
		cout << "Done Creating Output File" << endl;
		xyzOut.depth.clear();
		xyzOut.error.clear();
		xyzOut.nEi.clear();
		xyzOut.rEi.clear();
		if(abs(additionalOptions.find("-propUncert")->second) == 1)
		{
			xyzOut.depth0.clear();
			xyzOut.error0.clear();
		}
		if(abs(additionalOptions.find("-kalman")->second) == 1 )
		{
			xyzOut.depthK.clear();
			xyzOut.errorK.clear();
		}

		xMeshVectorMC.clear();
		yMeshVectorMC.clear();
	}

	if(MCFlag)
	{
		//Take the standard deviation
		double ziSTD = standardDeviation(&ziToStandardDeviation);
		double eiSTD = standardDeviation(&eiToStandardDeviation);
		double neiSTD = standardDeviation(&neiToStandardDeviation);
		double reiSTD = standardDeviation(&reiToStandardDeviation);
		
		printf("Computed Standard Deviation:\n");
		printf("\tDepth: %.8f\n", ziSTD);
		printf("\tError: %.8f\n", eiSTD);
		printf("\tNormalized Error: %.8f\n", neiSTD);
		printf("\tResidual Error: %.8f\n", reiSTD);

		if(abs(additionalOptions.find("-propUncert")->second) == 1)
		{
			double zi0STD = standardDeviation(&zi0ToStandardDeviation);
			double ei0STD = standardDeviation(&ei0ToStandardDeviation);
			printf("\tDepth (PU): %.8f\n", zi0STD);
			printf("\tError (PU): %.8f\n", ei0STD);
		}
		if(abs(additionalOptions.find("-kalman")->second) == 1 )
		{
			double ziKSTD = standardDeviation(&ziKToStandardDeviation);
			double eiKSTD = standardDeviation(&eiKToStandardDeviation);
			printf("\tDepth (K): %.8f\n", ziKSTD);
			printf("\tError (K): %.8f\n", eiKSTD);
		}

	}
	xMC.clear();
	yMC.clear();
	zMC.clear();
	eMC.clear();
	hMC.clear();
	vMC.clear();

	return returnValue;
}

#ifdef _MSC_VER
//Restore warning state -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( pop )
#endif 
#endif
