#include "fileReader.h"
#include "supportedFileTypes.h"
#include <string>
#include <cstring>	//UNIX
#include <cmath>
#include <string.h>
#include "standardOperations.h"
#include "constants.h" //UNIX for INT_MIN

#ifdef _MSC_VER
//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	//#include "../WarningStates.h"	//Disable all Warnings!!!
	#pragma warning ( disable : 4996 )	//Deprecated call for strcpy
#endif
#else 
//UNIX
#define strtok_s strtok_r
#define sscanf_s sscanf
#define _countof(a) (sizeof(a)/sizeof(*(a)))
#endif

//************************************************************************************
// SUBROUTINE I: Function call for reading input files. The starting point for file reader.
//************************************************************************************
int readFile(string &fileName, vector<TRUE_DATA> *inputData, BOUNDING_BOX bbox, int noerr, int nonegdepth, int unScaledAvgInputs)
{
	int retVal = SUCCESS;

	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	ifstream inListFile;
	string inputFile;
	char inputFileChar[512];
	vector<string> inputFiles;

	//************************************************************************************
	// I. Try to open the input file
	//************************************************************************************
	inListFile.open(fileName.c_str(), ifstream::in);
	if (!inListFile.is_open())
	{
		cerr << "Unable to open input file list!" << endl;
		cerr << fileName.c_str() << endl;
		return LIST_FILE_ERROR;
	}
	while(!inListFile.eof())
	{
		//Read in each input file to be used
		//Discard if it's a blank at the end of the file
		inListFile.getline(inputFileChar, 512);
		if (inListFile.gcount() > 2)
		{
			inputFile = inputFileChar;
			size_t found = inputFile.find_first_of("#");
			if(found!=string::npos)
				inputFile.erase(found);
			inputFile.erase(inputFile.find_last_not_of(" \r\t")+1);
			if(!inputFile.empty())
				inputFiles.push_back(inputFile);
		}
	}
	inListFile.clear();
	inListFile.close();
	cout << "Number of Input Files to Read: " << inputFiles.size() << endl;
	(*inputData) = vector<TRUE_DATA>(inputFiles.size());

	//************************************************************************************
	// II. Check the supported file types header file and process accordingly
	//************************************************************************************
	//If the same value is read-in for all e, h, or v, error we assume this value is an average
	//and scale it as a function of each sounding's depth for a new value.  The may be disabled
	//for testing purposes with the -useUnScaledAvgInputs flag.
	//If no value is given for e, we use 1% of the depth unless the depth is 0, then we use 0.01.
	//The sonar resolution (SR) is a function of the speed of sound and the sampling speed.
	//Specifically, SR = 1500 m/s/30,000 hz = 0.05 m.
	//Here we are using 0.01 as a default which corresponds to a sampling frequency of 150,000 hz
	//We need to check the average frequency for sonar.
	int retval =!SUCCESS;
	for (int i = 0; i < (const int)inputFiles.size(); i++)
	{
		//A. No Error
		if ((inputFiles[i].find(no_error_DAT) != string::npos) || (inputFiles[i].find(no_error_TXT) != string::npos))
		{
			retVal = readXYZ(inputFiles[i], bbox, inputData, i, noerr, nonegdepth);
			if (retVal != SUCCESS)
				break;
		}//B. Error
		else if ((inputFiles[i].find(error_DAT) != string::npos) || (inputFiles[i].find(error_TXT) != string::npos))
		{
			retVal = readXYZE(inputFiles[i], bbox, inputData, i, noerr, nonegdepth, unScaledAvgInputs);
			if (retVal != SUCCESS)
				break;
		}//C. HV Error
		else if ((inputFiles[i].find(hv_error_DAT) != string::npos) || (inputFiles[i].find(hv_error_TXT) != string::npos))
		{
			retVal = readXYZHV(inputFiles[i], bbox, inputData, i, noerr, nonegdepth, unScaledAvgInputs);
			if (retVal != SUCCESS)
				break;
		}//D. GSF
		else if ((inputFiles[i].find(GSF) != string::npos) || (inputFiles[i].find(GSF2) != string::npos))
		{
			retVal = readGSF(inputFiles[i], bbox, inputData, i, noerr, nonegdepth, unScaledAvgInputs);
			if (retVal != SUCCESS)
				break;
		}//E. ARC_ASCII_RASTER
		else if ((inputFiles[i].find(raster) != string::npos) || (inputFiles[i].find(raster2) != string::npos))
		{
			retVal = readARC_ASCII_RASTER(inputFiles[i], bbox, inputData, i, noerr, nonegdepth);
			if (retVal != SUCCESS)
				break;
		}
		else 
		{
			retval = LIST_FILE_EXT_ERROR;
			cerr << "Unable to open input file list!" << endl;
			cerr << "File type not supported (no identifiable extension found)!" << endl;
			cerr << fileName.c_str() << endl;
			break; //error opening files
		}
	}
	
	printf("Done Reading Input Files\n\n");
	inputFiles.clear();
	return retVal;
}

//************************************************************************************
// SUBROUTINE II: Function call for reading pre-interpolated file locations.
//************************************************************************************
int readLocationsFile(string &fileName, FORCED_LOCATIONS *inputData, int usagePreInterpLocsLatLon)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int retVal = SUCCESS;
	char inputFileChar[128];
	char delim[] = " \t";
	double lon, lat;
	char header[128];
	ifstream inFile;
	(*inputData).longitudeSum = 0.00;
	(*inputData).latitudeSum  = 0.00;
	char *next_token;

	//************************************************************************************
	// I. Try to open the input file
	//************************************************************************************
	inFile.open(fileName.c_str(), ios::in);
	if (!inFile.is_open())
	{
		cerr << "Unable to open input file: " << fileName << endl;
		return IN_FILE_ERROR;
	}

	printf("Reading Interpolation Locations Input File: %s ......", fileName.c_str());

	//Check for header
	if ((inFile.gcount() > 3) && !inFile.eof())
	{
		next_token = NULL;
		if(! (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%s", &header, _countof(header) ) <= 0))
			inFile.getline(inputFileChar, 128);
	}

	//************************************************************************************
	// II. Read the file locations
	//************************************************************************************
	inFile.getline(inputFileChar, 128);
	while ((inFile.gcount() > 3))
	{
		next_token = NULL;
		
		if(usagePreInterpLocsLatLon < 0) //1
		{
			if (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%lf", &lat) <= 0)
			{
				retVal = IN_FILE_ERROR;
				break;
			}

			if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &lon) <= 0)
			{
				retVal = IN_FILE_ERROR;
				break;
			}
		}
		else if(usagePreInterpLocsLatLon > 0) //-1
		{
			if (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%lf", &lon) <= 0)
			{
				retVal = IN_FILE_ERROR;
				break;
			}

			if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &lat) <= 0)
			{
				retVal = IN_FILE_ERROR;
				break;
			}
		}
		else if(abs(usagePreInterpLocsLatLon) != 2)
		{
			cout<< "Wrong -PreInterpolatedLocations usage value." << endl;
			return !SUCCESS;
		}
		(*inputData).forcedLonCoord.push_back(lon);
		(*inputData).forcedLatCoord.push_back(lat);
		if(abs(usagePreInterpLocsLatLon) != 2)
		{
			(*inputData).longitudeSum += (lon+180.00) - int((lon+180.00)/360.00)*360.00-180.00;
			(*inputData).latitudeSum  += (lat+180.00) - int((lat+180.00)/360.00)*360.00-180.00;
		}
		if(inFile.eof())
			break;
		inFile.getline(inputFileChar, 128);
		
	}
	printf(" %d Positions Read Successfully.\n", (int)(*inputData).forcedLonCoord.size());

	inFile.clear();
	inFile.close();

	printf("Done Reading Interpolation Location Input File\n\n");

	return retVal;
}

//************************************************************************************
// SUBROUTINE III: Function call for reading GSF Files
//************************************************************************************
int readGSF(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth, int unScaledAvgInputs)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int retVal = SUCCESS;

	int gsfHandle = -1;
	char *fName = (char*)"";

	double lon, lat, depth, hError, vError, error, wMultiplier;
	double tideCorrection = 0.00;
	double maximumDataOffset = INT_MAX;
	int i = 0;
	int beamon = -2;
	int stat;
	double z_avg;
	double z_max = MIN_INT, z_min = MAX_INT;

	gsfRecords gsfRec;
	gsfDataID id;

	string fileNameTemp;
	char *next_token;
	(*inputData)[pos].longitudeSum = 0.00;
	(*inputData)[pos].latitudeSum  = 0.00;

	//A. Parse out the user weighting
	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "&\t\n", &next_token);
	char* wM = strtok_s(NULL, "& \t\n", &next_token);
	wMultiplier = ((wM==NULL)||(sscanf_s(wM, "%lf", &wMultiplier) <= 0))? 0.00:(1/wMultiplier);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "^\t\n", &next_token);
	wM = strtok_s(NULL, "^ \t\n", &next_token);
	tideCorrection = ((wM==NULL)||(sscanf_s(wM, "%lf", &tideCorrection) <= 0))? 0.00:(tideCorrection);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "#\t\n", &next_token);
	wM = strtok_s(NULL, "# \t\n", &next_token);
	maximumDataOffset = ((wM==NULL)||(sscanf_s(wM, "%lf", &maximumDataOffset) <= 0))? INT_MAX:(maximumDataOffset);

	next_token = NULL;
	fileName = (char*)strtok_s((char*)fileName.c_str(), "&^#\t\n", &next_token);

	//************************************************************************************
	// I. Try to open the input file
	//************************************************************************************
	strcpy(fName,fileName.c_str());

	gsfHandle = open_gsffile(fName);
	if(gsfHandle == -1)
	{
		cerr << "Unable to open input file: " << fileName << endl;
		return IN_FILE_ERROR;
	}

	//************************************************************************************
	// II. Read the file data
	//************************************************************************************
	printf("Reading Input File: %s ......", fileName.c_str());
	while (!gsfReader(&lon, &lat, &depth, &hError, &vError,	0, gsfHandle, &beamon, &stat, &gsfRec, &id))
	{
		if ((lon == NONSENSE) && (lat == NONSENSE) && (depth == NONSENSE))
			continue;
		//Cut the bounding box
		if (bbox.doBoundingBox)
		{
			if ( (lat > bbox.bboxTop) || (lat < bbox.bboxBottom) || (lon > bbox.bboxRight) || (lon < bbox.bboxLeft))
				continue;
		}
		depth += tideCorrection;
		//Keep only water
		if (!nonegdepth || depth >= 0)
		{
			(*inputData)[pos].lon.push_back(lon/60.00);
			(*inputData)[pos].lat.push_back(lat/60.00);
			(*inputData)[pos].depth.push_back(depth);
			(*inputData)[pos].h_Error.push_back(abs(hError));
			(*inputData)[pos].v_Error.push_back(abs(vError));
			if(z_max < depth)
				z_max = depth;
			else if(z_min > depth)
				z_min = depth;

			//A. Compute the longitude sum for use later when converting to UTM
			(*inputData)[pos].longitudeSum += (lon+180.00) - int((lon+180.00)/360.00)*360.00-180.00;
			(*inputData)[pos].latitudeSum  += (lat+180.00) - int((lat+180.00)/360.00)*360.00-180.00;
		}
	}
	hError = standardDeviation(&(*inputData)[pos].h_Error, False);
	vError = standardDeviation(&(*inputData)[pos].v_Error, False);
	z_avg = (z_max-z_min)/2;
	
	if((hError == 0 && vError == 0) && ((*inputData)[pos].h_Error[0] ==  0 && (*inputData)[pos].v_Error[0] ==  0))
	{
		//No e provided; all 0's.
		//Calculate error as 1% of depth and find h and v.
		//Default 0.01 if depth is 0.
		double hvError;
		for(int cnt=0; cnt<(*inputData)[pos].h_Error.size();cnt++)
		{
			if((*inputData)[pos].depth[cnt] == 0)
				error = 0.01;
			else
				error = abs((*inputData)[pos].depth[cnt]*0.01);
			//hvError = sqrt(pow(error,2)/2);
			(*inputData)[pos].error.push_back(error);
			(*inputData)[pos].h_Error[cnt] = 0.0;//hvError;
			(*inputData)[pos].v_Error[cnt] = error;//hvError;
		}
	}
	else
	{
		//if (hError == 0 && (*inputData)[pos].h_Error[0] == 0)
		//{
		//	//Double v if no h (std dev of h is 0)
		//	(*inputData)[pos].h_Error = vector<double>((*inputData)[pos].v_Error);
		//}
		//else if(vError == 0 && (*inputData)[pos].v_Error[0] == 0)
		//{
		//	//Double h if no v (std dev of v is 0)
		//	(*inputData)[pos].v_Error = vector<double>((*inputData)[pos].h_Error);
		//}

		//Replace missing h and v data with their std devs, and Calculate e 
		for(int cnt = 0; cnt < (*inputData)[pos].v_Error.size(); cnt++)
		{
			if(!unScaledAvgInputs)
			{
				if(hError == 0) 
				{
					if(z_avg != 0)
						(*inputData)[pos].h_Error[cnt] = (*inputData)[pos].h_Error[cnt] * (1 + ((z_avg - (*inputData)[pos].depth[cnt])/z_avg));
				}
				if(vError == 0)
				{
					if(z_avg != 0)
						(*inputData)[pos].v_Error[cnt] = (*inputData)[pos].v_Error[cnt] * (1 + ((z_avg - (*inputData)[pos].depth[cnt])/z_avg));
				}
			}
			if((*inputData)[pos].h_Error[cnt] == 0)
				(*inputData)[pos].h_Error[cnt] = hError;
			if((*inputData)[pos].v_Error[cnt] == 0)
				(*inputData)[pos].v_Error[cnt] = vError;
			
			//error = sqrt(pow(abs((*inputData)[pos].h_Error[cnt]),2) + pow(abs((*inputData)[pos].v_Error[cnt]),2));
			(*inputData)[pos].error.push_back(abs((*inputData)[pos].h_Error[cnt]));//error);
		}	
	}
	//B. Zero out errors so we don't crash later since none were provided
	(*inputData)[pos].x = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].y = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].maximumDataOffset = maximumDataOffset;

	close_gsffile(&gsfHandle);
	printf(" %d Records Read Successfully.\n", (int)(*inputData)[pos].lon.size());
	if (wMultiplier != 0)
		cout << "\tUsing Multiplier: " << wMultiplier << endl;
	if (tideCorrection != 0)
		cout << "\tUsing Tide Correction: " << tideCorrection << endl;
	if (maximumDataOffset != INT_MAX)
		cout << "\tUsing Maximum Data Offset: " << maximumDataOffset << endl;

	return retVal;
}

//************************************************************************************
// SUBROUTINE IV: Function call for reading XYZ Files
//************************************************************************************
int readXYZ(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int retVal = SUCCESS;
	double lon, lat, depth, wMultiplier, error, hvError;
	double tideCorrection = 0.00;
	double maximumDataOffset = INT_MAX;
	char inputFileChar[128];
	char delim[] = " \t";
	string fileNameTemp;
	ifstream inFile;
	char *next_token;
	char header[128];

	//A. Parse out the user weighting
	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "&\t\n", &next_token);
	char* wM = strtok_s(NULL, "& \t\n", &next_token);
	wMultiplier = ((wM==NULL)||(sscanf_s(wM, "%lf", &wMultiplier) <= 0))? 0.00:(1/wMultiplier);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "^\t\n", &next_token);
	wM = strtok_s(NULL, "^ \t\n", &next_token);
	tideCorrection = ((wM==NULL)||(sscanf_s(wM, "%lf", &tideCorrection) <= 0))? 0.00:(tideCorrection);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "#\t\n", &next_token);
	wM = strtok_s(NULL, "# \t\n", &next_token);
	maximumDataOffset = ((wM==NULL)||(sscanf_s(wM, "%lf", &maximumDataOffset) <= 0))? INT_MAX:(maximumDataOffset);

	//Check if format is at the end of the filename and remove it
	string typeTemp = (fileNameTemp.find(no_error_DAT) != string::npos) ? no_error_DAT: no_error_TXT;
	size_t found = fileNameTemp.find(typeTemp);
	size_t foundFmt = fileNameTemp.substr(found).find_first_not_of(typeTemp);
	if(foundFmt!=string::npos)
	{
		fileNameTemp.erase(found + typeTemp.length(), fileNameTemp.length());
		fileName = fileNameTemp;
	}
	next_token = NULL;
	fileName = (char*)strtok_s((char*)fileName.c_str(), "&^#\t\n", &next_token);

	//************************************************************************************
	// I. Try to open the input file
	//************************************************************************************
	inFile.open(fileName.c_str(), ios::in);
	if (!inFile.is_open())
	{
		cerr << "Unable to open input file: " << fileName << endl;
		return IN_FILE_ERROR;
	}

	(*inputData)[pos].longitudeSum = 0.00;
	(*inputData)[pos].latitudeSum  = 0.00;

	printf("Reading Input File: %s ......", fileName.c_str());
	
	//Check for header
	if ((inFile.gcount() > 3) && !inFile.eof())
	{
		next_token = NULL;
		if(! (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%s", &header, _countof(header) ) <= 0))
			inFile.getline(inputFileChar, 128);
	}

	//************************************************************************************
	// II. Read the file data
	//************************************************************************************
	inFile.getline(inputFileChar, 128);
	while (inFile.gcount() > 3)
	{
		next_token = NULL;
		if (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%lf", &lon) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &lat) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &depth) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}
		//Cut the bounding box
		if (bbox.doBoundingBox)
		{
			if ( (lat > bbox.bboxTop) || (lat < bbox.bboxBottom) || (lon > bbox.bboxRight) || (lon < bbox.bboxLeft))
			{
				if(inFile.eof())
					break;
				inFile.getline(inputFileChar, 128);
				continue;
			}
		}
		depth += tideCorrection;
		//Keep only water
		if (!nonegdepth || depth >= 0)
		{
			(*inputData)[pos].lon.push_back(lon);
			(*inputData)[pos].lat.push_back(lat);
			(*inputData)[pos].depth.push_back(depth);
			//B. Compute e, h, and v since none were provided.
			//Calculate error as 1% of depth and find h and v.
			//Default 0.01 if depth is 0
			if(depth == 0)
				error = 0.01;
			else
				error = abs(depth*.01);
			//hvError = sqrt(pow(error,2)/2);
			(*inputData)[pos].error.push_back(error);
			(*inputData)[pos].h_Error.push_back(0.0);//hvError);
			(*inputData)[pos].v_Error.push_back(error);//hvError);
			//C. Compute the longitude sum for use later when converting to UTM
			//Make sure the longitude is between -180.00 .. 179.9 and then add it to the sum
			(*inputData)[pos].longitudeSum += (lon+180.00) - int((lon+180.00)/360.00)*360.00-180.00;
			(*inputData)[pos].latitudeSum  += (lat+180.00) - int((lat+180.00)/360.00)*360.00-180.00;
		}
		if(inFile.eof())
			break;
		inFile.getline(inputFileChar, 128);
	}
	
	//D. Zero out errors so we don't crash later since none were provided
	(*inputData)[pos].x = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].y = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].maximumDataOffset = maximumDataOffset;
	printf(" %d Records Read Successfully.\n", (int)(*inputData)[pos].lon.size());
	if (wMultiplier != 0)
		cout << "\tUsing Multiplier: " << wMultiplier << endl;
	if (tideCorrection != 0)
		cout << "\tUsing Tide Correction: " << tideCorrection << endl;
	if (maximumDataOffset != INT_MAX)
		cout << "\tUsing Maximum Data Offset: " << maximumDataOffset << endl;

	inFile.clear();
	inFile.close();
	//vector<double>::iterator minIt= min_element((*inputData)[pos].lon.begin(),(*inputData)[pos].lon.end());
	//vector<double>::iterator maxIt= max_element((*inputData)[pos].lon.begin(),(*inputData)[pos].lon.end());
	//if(*minIt<0 & *maxIt>0)
	//{
	//	if(0-*minIt > -180-*maxIt)
	//	{
	//		for(int i =0;i<(*inputData)[pos].lon.size();i++)
	//		{
	//			if((*inputData)[pos].lon[i] < 0)
	//				(*inputData)[pos].lon[i]+=360;
	//		}
	//	}
	//}


	return retVal;
}

//************************************************************************************
// SUBROUTINE V: Function call for reading XYZE Files
//************************************************************************************
int readXYZE(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth, int unScaledAvgInputs)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int retVal = SUCCESS;
	char inputFileChar[128];
	char delim[] = " \t";
	double lon, lat, depth, error, hvError, wMultiplier; double z_avg = 0;
	double tideCorrection = 0.00;
	double maximumDataOffset = INT_MAX;
	string fileNameTemp;
	ifstream inFile;
	char *next_token;
	char header[128];
	double z_max = MIN_INT, z_min = MAX_INT;

	//A. Parse out the user weighting
	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "&\t\n", &next_token);
	char* wM = strtok_s(NULL, "& \t\n", &next_token);
	wMultiplier = ((wM==NULL)||(sscanf_s(wM, "%lf", &wMultiplier) <= 0))? 1.00:(1/wMultiplier);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "^\t\n", &next_token);
	wM = strtok_s(NULL, "^ \t\n", &next_token);
	tideCorrection = ((wM==NULL)||(sscanf_s(wM, "%lf", &tideCorrection) <= 0))? 0.00:(tideCorrection);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "#\t\n", &next_token);
	wM = strtok_s(NULL, "# \t\n", &next_token);
	maximumDataOffset = ((wM==NULL)||(sscanf_s(wM, "%lf", &maximumDataOffset) <= 0))? INT_MAX:(maximumDataOffset);

	next_token = NULL;
	fileName = (char*)strtok_s((char*)fileName.c_str(), "&^#\t\n", &next_token);

	//************************************************************************************
	// I. Try to open the input file
	//************************************************************************************
	inFile.open(fileName.c_str(), ios::in);
	if (!inFile.is_open())
	{
		cerr << "Unable to open input file: " << fileName << endl;
		return IN_FILE_ERROR;
	}

	(*inputData)[pos].longitudeSum = 0.00;
	(*inputData)[pos].latitudeSum  = 0.00;

	printf("Reading Input File: %s ......", fileName.c_str());
	
	//Check for header
	if ((inFile.gcount() > 3) && !inFile.eof())
	{
		next_token = NULL;
		if(! (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%s", &header, _countof(header) ) <= 0))
			inFile.getline(inputFileChar, 128);
	}

	//************************************************************************************
	// II. Read the file data
	//************************************************************************************
	inFile.getline(inputFileChar, 128);
	while ((inFile.gcount() > 3))
	{
		next_token = NULL;
		if (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%lf", &lon) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &lat) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &depth) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &error) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}
		//Cut the bounding box
		if (bbox.doBoundingBox)
		{
			if ( (lat > bbox.bboxTop) || (lat < bbox.bboxBottom) || (lon > bbox.bboxRight) || (lon < bbox.bboxLeft))
			{
				if(inFile.eof())
					break;
				inFile.getline(inputFileChar, 128);
				continue;
			}
		}
		depth += tideCorrection;
		//Keep only water
		if (!nonegdepth || depth >= 0 )
		{
			(*inputData)[pos].lon.push_back(lon);
			(*inputData)[pos].lat.push_back(lat);
			(*inputData)[pos].depth.push_back(depth);
			if(z_max < depth)
				z_max = depth;
			else if(z_min > depth)
				z_min = depth;

			(*inputData)[pos].error.push_back(abs(error));
			//B. Compute the longitude sum for use later when converting to UTM
			//Make sure the longitude is between -180.00 .. 179.9 and then add it to the sum
			(*inputData)[pos].longitudeSum += (lon+180.00) - int((lon+180.00)/360.00)*360.00-180.00;
			(*inputData)[pos].latitudeSum  += (lat+180.00) - int((lat+180.00)/360.00)*360.00-180.00;
		}
		if(inFile.eof())
			break;
		inFile.getline(inputFileChar, 128);
	}
	error = standardDeviation(&(*inputData)[pos].error, False);
	z_avg = (z_max - z_min)/2;
	
	for(int cnt = 0; cnt < (*inputData)[pos].error.size(); cnt++)
	{
		if(error == 0)
		{
			if((*inputData)[pos].error[0] != 0 && !unScaledAvgInputs)
			{
				//All e the same value.
				//We assume this is the average and scale by the point's depth.
				if(z_avg != 0)
					(*inputData)[pos].error[cnt] = (*inputData)[pos].error[cnt] * (1 + ((z_avg-(*inputData)[pos].depth[cnt])/z_avg));
			}
			else 
			{
				//All 0s.
				//Calculate e as 1% of depth and find h and v.
				//Default 0.01 if depth is 0.
				if((*inputData)[pos].depth[cnt] == 0)
					(*inputData)[pos].error[cnt] = 0.01;
				else
					(*inputData)[pos].error[cnt] = abs((*inputData)[pos].depth[cnt]*0.01);
			}
		}
		else
		{
			if((*inputData)[pos].error[cnt] == 0)
				(*inputData)[pos].error[cnt] = error;
		}
		//hvError = sqrt(pow((*inputData)[pos].error[cnt],2)/2);
		(*inputData)[pos].h_Error.push_back(0.0);//hvError);
		(*inputData)[pos].v_Error.push_back(abs((*inputData)[pos].error[cnt]));//hvError);
	}
	
	//C. Zero out errors so we don't crash later since none were provided
	(*inputData)[pos].x = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].y = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].maximumDataOffset = maximumDataOffset;
	printf(" %d Records Read Successfully.\n", (int)(*inputData)[pos].lon.size());
	if (wMultiplier != 0)
		cout << "\tUsing Multiplier: " << wMultiplier << endl;
	if (tideCorrection != 0)
		cout << "\tUsing Tide Correction: " << tideCorrection << endl;
	if (maximumDataOffset != INT_MAX)
		cout << "\tUsing Maximum Data Offset: " << maximumDataOffset << endl;

	inFile.clear();
	inFile.close();

	return retVal;
}

//************************************************************************************
// SUBROUTINE VI: Function call for reading XYZHV Files
//**********************************************************************************
int readXYZHV(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth, int unScaledAvgInputs)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int retVal = SUCCESS;
	char inputFileChar[128];
	char delim[] = " \t";
	double lon, lat, depth, hError, vError, error, wMultiplier;
	double tideCorrection = 0.00;
	double maximumDataOffset = INT_MAX;
	string fileNameTemp;
	ifstream inFile;
	char *next_token;
	char header[128];
	double z_avg;
	double z_max = MIN_INT, z_min = MAX_INT;

	//A. Parse out the user weighting
	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "&\t\n", &next_token);
	char* wM = strtok_s(NULL, "& \t\n", &next_token);
	wMultiplier = ((wM==NULL)||(sscanf_s(wM, "%lf", &wMultiplier) <= 0))? 0.00:(1/wMultiplier);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "^\t\n", &next_token);
	wM = strtok_s(NULL, "^ \t\n", &next_token);
	tideCorrection = ((wM==NULL)||(sscanf_s(wM, "%lf", &tideCorrection) <= 0))? 0.00:(tideCorrection);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "#\t\n", &next_token);
	wM = strtok_s(NULL, "# \t\n", &next_token);
	maximumDataOffset = ((wM==NULL)||(sscanf_s(wM, "%lf", &maximumDataOffset) <= 0))? INT_MAX:(maximumDataOffset);

	next_token = NULL;
	fileName = (char*)strtok_s((char*)fileName.c_str(), "&^#\t\n", &next_token);

	//************************************************************************************
	// I. Try to open the input file
	//************************************************************************************
	inFile.open(fileName.c_str(), ios::in);
	if (!inFile.is_open())
	{
		cerr << "Unable to open input file: " << fileName << endl;
		return IN_FILE_ERROR;
	}

	(*inputData)[pos].longitudeSum = 0.00;
	(*inputData)[pos].latitudeSum  = 0.00;

	printf("Reading Input File: %s ......", fileName.c_str());

	//Check for header
	if ((inFile.gcount() > 3) && !inFile.eof())
	{
		next_token = NULL;
		if(! (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%s", &header, _countof(header) ) <= 0))
			inFile.getline(inputFileChar, 128);
	}

	//************************************************************************************
	// II. Read the file data
	//************************************************************************************
	inFile.getline(inputFileChar, 128);
	while ((inFile.gcount() > 3))
	{
		next_token = NULL;
		if (sscanf_s(strtok_s(inputFileChar, delim, &next_token), "%lf", &lon) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &lat) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &depth) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &hError) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}

		if (sscanf_s(strtok_s(NULL, delim, &next_token), "%lf", &vError) <= 0)
		{
			retVal = IN_FILE_ERROR;
			break;
		}
		//Cut the bounding box
		if (bbox.doBoundingBox)
		{
			if ( (lat > bbox.bboxTop) || (lat < bbox.bboxBottom) || (lon > bbox.bboxRight) || (lon < bbox.bboxLeft))
			{
				if(inFile.eof())
					break;
				inFile.getline(inputFileChar, 128);
				continue;
			}
		}
		depth += tideCorrection;
		//Keep only water
		if (!nonegdepth || depth >= 0)
		{
			(*inputData)[pos].lon.push_back(lon);
			(*inputData)[pos].lat.push_back(lat);
			(*inputData)[pos].depth.push_back(depth);
			(*inputData)[pos].h_Error.push_back(abs(hError));
			(*inputData)[pos].v_Error.push_back(abs(vError));
			if(z_max < depth)
				z_max = depth;
			else if(z_min > depth)
				z_min = depth;

			//A. Compute the longitude sum for use later when converting to UTM
			//Make sure the longitude is between -180.00 .. 179.9 and then add it to the sum
			(*inputData)[pos].longitudeSum += (lon+180.00) - int((lon+180.00)/360.00)*360.00-180.00;
			(*inputData)[pos].latitudeSum  += (lat+180.00) - int((lat+180.00)/360.00)*360.00-180.00;
		}
		if(inFile.eof()) // SJZ
			break;
		inFile.getline(inputFileChar, 128);
	}
	hError = standardDeviation(&(*inputData)[pos].h_Error, False);
	vError = standardDeviation(&(*inputData)[pos].v_Error, False);
	z_avg = (z_max-z_min)/2;
	
	if((hError == 0 && vError == 0) && ((*inputData)[pos].h_Error[0] ==  0 && (*inputData)[pos].v_Error[0] ==  0))
	{
		//No e provided; all 0s.
		//Calculate error as 1% of depth and find h and v.
		//Default 0.01 if depth is 0.
		double hvError;
		for(int cnt=0; cnt<(*inputData)[pos].h_Error.size();cnt++)
		{
			if((*inputData)[pos].depth[cnt] == 0)
				error = 0.01;
			else
				error = abs((*inputData)[pos].depth[cnt]*0.01);
			//hvError = sqrt(pow(error,2)/2);
			(*inputData)[pos].error.push_back(error);
			(*inputData)[pos].h_Error[cnt] = 0;//hvError;
			(*inputData)[pos].v_Error[cnt] = error;//hvError;
		}
	}
	else
	{
		//if (hError == 0 && (*inputData)[pos].h_Error[0] == 0)
		//{
		//	//Double v if no h (std dev of h is 0)
		//	(*inputData)[pos].h_Error = vector<double>((*inputData)[pos].v_Error);
		//}
		//else if(vError == 0 && (*inputData)[pos].v_Error[0] == 0)
		//{
		//	//Double h if no v (std dev of v is 0)
		//	(*inputData)[pos].v_Error = vector<double>((*inputData)[pos].h_Error);
		//}

		//Replace missing h and v data with their std devs, and Calculate e 
		for(int cnt = 0; cnt < (*inputData)[pos].v_Error.size(); cnt++)
		{
			if(!unScaledAvgInputs)
			{
				if(hError == 0) 
				{
					if(z_avg != 0)
						(*inputData)[pos].h_Error[cnt] = (*inputData)[pos].h_Error[cnt] *( 1 + ((z_avg - (*inputData)[pos].depth[cnt])/z_avg));
				}
				if(vError == 0)
				{
					if(z_avg != 0)
						(*inputData)[pos].v_Error[cnt] = (*inputData)[pos].v_Error[cnt] * (1 + ((z_avg - (*inputData)[pos].depth[cnt])/z_avg));
				}
			}
			if((*inputData)[pos].h_Error[cnt] == 0)
				(*inputData)[pos].h_Error[cnt] = hError;
			if((*inputData)[pos].v_Error[cnt] == 0)
				(*inputData)[pos].v_Error[cnt] = vError;

			//error = sqrt(pow(abs((*inputData)[pos].h_Error[cnt]),2) + pow(abs((*inputData)[pos].v_Error[cnt]),2));
			(*inputData)[pos].error.push_back(abs((*inputData)[pos].v_Error[cnt]));//error);
		}	
	}
	
	//B. Zero out errors so we don't crash later since none were provided
	(*inputData)[pos].x = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].y = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].maximumDataOffset = maximumDataOffset;
	printf(" %d Records Read Successfully.\n", (int)(*inputData)[pos].lon.size());
	if (wMultiplier != 0)
		cout << "\tUsing Multiplier: " << wMultiplier << endl;
	if (tideCorrection != 0)
		cout << "\tUsing Tide Correction: " << tideCorrection << endl;
	if (maximumDataOffset != INT_MAX)
		cout << "\tUsing Maximum Data Offset: " << maximumDataOffset << endl;

	inFile.clear();
	inFile.close();

	return retVal;
}

//************************************************************************************
// SUBROUTINE VII: Function call for reading ARC ASCII Raster Files
//************************************************************************************
int readARC_ASCII_RASTER(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int retVal = SUCCESS;
	string temp;
	double noData = -9999; //Is this suppose to be 999999? SJZ
	double depth, wMultiplier;
	double tideCorrection = 0.00;
	double maximumDataOffset = INT_MAX;
	int nx, ny, lowerLeftX, lowerLeftY, lx;
	string fileNameTemp;
	ifstream inFile;
	char *next_token;
	double error,hvError;

	//A. Parse out the user weighting
	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "&\t\n", &next_token);
	char* wM = strtok_s(NULL, "& \t\n", &next_token);
	wMultiplier = ((wM==NULL)||(sscanf_s(wM, "%lf", &wMultiplier) <= 0))? 0.00:(1/wMultiplier);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "^\t\n", &next_token);
	wM = strtok_s(NULL, "^ \t\n", &next_token);
	tideCorrection = ((wM==NULL)||(sscanf_s(wM, "%lf", &tideCorrection) <= 0))? 0.00:(tideCorrection);

	next_token = NULL;
	fileNameTemp = fileName;
	fileNameTemp = (char*)strtok_s((char*)fileNameTemp.c_str(), "#\t\n", &next_token);
	wM = strtok_s(NULL, "# \t\n", &next_token);
	maximumDataOffset = ((wM==NULL)||(sscanf_s(wM, "%lf", &maximumDataOffset) <= 0))? INT_MAX:(maximumDataOffset);

	next_token = NULL;
	fileName = (char*)strtok_s((char*)fileName.c_str(), "&^#\t\n", &next_token);

	//************************************************************************************
	// I. Try to open the input file
	//************************************************************************************
	inFile.open(fileName.c_str(), ios::in);
	if (!inFile.is_open())
	{
		cerr << "Unable to open input file: " << fileName << endl;
		return IN_FILE_ERROR;
	}

	(*inputData)[pos].longitudeSum = 0.00;
	(*inputData)[pos].latitudeSum  = 0.00;

	printf("Reading Input File: %s ......", fileName.c_str());
	//************************************************************************************
	// II. Read the file data
	//************************************************************************************

	//A. Read the header data
	inFile >> temp >> nx;
	inFile >> temp >> ny;
	inFile >> temp >> lowerLeftX;
	inFile >> temp >> lowerLeftY;
	inFile >> temp >> lx;
	inFile >> temp >> noData;

	double xLeft = lowerLeftX + (0.5*lx);
	double yTop = lowerLeftY + (ny*lx) - (0.5*lx);
	double xPos = xLeft;
	double yPos = yTop;

	//B. Read the grid data
	for (int i = 0; i < ny; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			inFile >> depth;
			if (!(depth == noData || depth == 0))
			{
				//Cut the bounding box
				if (bbox.doBoundingBox)
				{
					if ( (yPos > bbox.bboxTop) || (yPos < bbox.bboxBottom) || (xPos > bbox.bboxRight) || (xPos < bbox.bboxLeft))
					{
						xPos = xPos + lx;
						continue;
					}
				}
				depth += tideCorrection;
				//Keep only water
				if (!nonegdepth || depth >= 0)
				{
					(*inputData)[pos].lon.push_back(xPos);
					(*inputData)[pos].lat.push_back(yPos);
					(*inputData)[pos].depth.push_back(depth);
					//C. Compute e, h, and v since none were provided.
					//Calculate error as 1% of depth and find h and v
					//Default 0.01 if depth is 0.
					if(depth == 0)
						error = 0.01;
					else
						error = abs(depth * 0.01);
					//hvError = sqrt(pow(error,2)/2);
					(*inputData)[pos].error.push_back(error);
					(*inputData)[pos].h_Error.push_back(0.0);//hvError);
					(*inputData)[pos].v_Error.push_back(error);//hvError);
					//D. Compute the longitude sum for use later when converting to UTM
					//Make sure the longitude is between -180.00 .. 179.9 and then add it to the sum
					(*inputData)[pos].longitudeSum += (xPos+180.00) - int((xPos+180.00)/360.00)*360.00-180.00;
					(*inputData)[pos].latitudeSum  += (yPos+180.00) - int((yPos+180.00)/360.00)*360.00-180.00;
				}
			}
			xPos = xPos + lx;
		}
		xPos = xLeft;
		yPos = yPos - lx;
	}
	(*inputData)[pos].x = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].y = vector<double>((*inputData)[pos].lon.size(), 0.00);
	(*inputData)[pos].maximumDataOffset = maximumDataOffset;
	printf(" %d Records Read Successfully.\n", (int)(*inputData)[pos].lon.size());
	if (wMultiplier != 0)
		cout << "\tUsing Multiplier: " << wMultiplier << endl;
	if (tideCorrection != 0)
		cout << "\tUsing Tide Correction: " << tideCorrection << endl;
	if (maximumDataOffset != INT_MAX)
		cout << "\tUsing Maximum Data Offset: " << maximumDataOffset << endl;

	inFile.clear();
	inFile.close();
	return retVal;
}


#ifdef _MSC_VER
//Restore warning state -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( pop )
#endif 
#endif