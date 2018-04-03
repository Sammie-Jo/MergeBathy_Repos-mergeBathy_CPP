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
#include "mergeBathy.h"
#include <time.h>
#include <string>
#include <string.h> //UNIX
#ifdef _WIN32
#include <direct.h>
#define getcwd _getcwd // stupid MSFT "deprecation" warning
#else
#include <unistd.h>
#define _MAX_PATH FILENAME_MAX
#define _getcwd getcwd
//PATH_MAX
#endif

using namespace std;

int main(char argc, char *argv[])
{
	//************************************************************************************
	//0. Declare variables and data structures:
	//************************************************************************************
	int returnValue = SUCCESS;

	//1. Input Arguments Variables
	string outputFileName;
	string inputFileList;
	string kernelName;
	string interpolationLocationsFileName;
	int usagePreInterpLocsLatLon;
	//int metersPreInterpLocsLatLon;
	double gridSpacingX_Lon;
	double gridSpacingY_Lat;
	double smoothingScaleX = -1;
	double smoothingScaleY = -1;
	double refLon;
	double refLat;
	double rotationAngle;
	int argLocation = 0;
	int numMCRuns = -1;
	double start, stop;

	map<string, int> additionalOptions;
	vector<TRUE_DATA> inputData;

	//2. Data structures for external interpolates and irregular grid points
	MB_ZGRID_DATA mbzData;
	GMT_SURFACE_DATA GMTSurfaceData;
	ALG_SPLINE_DATA ALGsplineData;
	FORCED_LOCATIONS forcedLocPositions;
	BOUNDING_BOX bbox;
	bbox.doBoundingBox = false;

	//4. Any additional options to be used must be set to their default values here before they can be used.
	additionalOptions["-noerr"] = 0;
	additionalOptions["-nmsei"] = 0;
	additionalOptions["-msri"] = 0;
	additionalOptions["-inputInMeters"] = 0;
	additionalOptions["-kriging"] = 0;
	additionalOptions["-ZGrid"] = 0;
	additionalOptions["-GMTSurface"] = 0;
	additionalOptions["-ALGSpline"] = 0;
	additionalOptions["-preInterpolatedLocations"] = 0;
	additionalOptions["-boundingBox"] = 0;
	additionalOptions["-computeOffset"] = 0;
	additionalOptions["-outputRasterFile"] = 0;
	additionalOptions["-outputBagFile"] = 0;
	additionalOptions["-multiThread"] = 0;
	additionalOptions["-numMCRuns"] = 0;
	additionalOptions["-modelflag"] = 0;
	additionalOptions["-nonegdepth"] = 0;
	additionalOptions["-printMatlabMatch"] = 0;
	additionalOptions["-nnInterp"] = 0;
	additionalOptions["-mse"] = 1;//0; 0 will not compute mse; we always compute by default.
	additionalOptions["-propUncert"] = 0;
	additionalOptions["-kalman"] = 0;
	additionalOptions["-printMSEwK"] = 0;
	additionalOptions["-appendFilename"] = 0;
	additionalOptions["-useUnscaledAvgInputs"] = 0;
	vector<string> unrecognizedParams;
	//************************************************************************************
	//I. Check input arguments and store data. If not then display the usage instructions.
	//************************************************************************************
	if (argc < 9)
	{
		cerr << "MergeBathy Build:" << OBF_VERSION_NUMBER << endl;
		cerr << "Usage: mergeBathy.exe <output_depth_file> <grid_spacing> [grid_spacing_Y]" << endl;
		cerr << "                   <kernel_name> <input_file_list> " << endl;
		cerr << "                   <ref_lon> <ref_lat> <rotation_angle> <num_MC_runs(-1 if not MC)>" << endl;
		cerr << "                   [-noerr] [-nmsei] [-msri] [-modelflag] [-nonegdepth] [-inputInMeters] [-kriging] [-msmooth <smoothing_scale_x> <smoothing_scale_y>]" << endl;
		cerr << "                   [-llsmooth <smoothing_scale_longitude (X)> <smoothing_scale_latitude (Y)>] [-llgrid]" << endl;
		cerr << "					[-computeOffset] [-outputRasterFile] [-outputBagFile] [-multiThread <num_threads>]" << endl;
		cerr << "                   [-ZGrid <grid_spacing_X> <grid_spacing_Y> <Z_Grid_Output_File_Name> <Tension_Factor (Typically 1e10)> <Usage: (1: Do not use as input. 2: Use as input. Negate the value to include error in the computation)> ]" << endl;
		cerr << "                   [-GMTSurface <grid_spacing_X> <grid_spacing_Y> <GMT_Surface_Output_File_Name> <Tension_Factor (Between 0 and 1)> <scale_factor> <alpha> <Usage: (1: Do not use as input. 2: Use as input. Negate the value to include error in the computation)> ]" << endl;
		cerr << "[-ALGSpline <grid_spacing_X> <grid_spacing_Y> <ALG_Surface_Output_File_Name> ]" << endl;
		cerr << "					[-preInterpolatedLocations <interpolation_location_file_name> <Usage: (1: Read in Lat,Lon. Negate to read in Lon,Lat.>]" << endl;
		cerr << "					[-boundingBox <upper_bound> <lower_bound> <right_bound> <left_bound>]" << endl;
		cerr << "					[-printMatlabMatch] [-nnInterp] [-useUnScaledAvgInputs]" << endl; 
		cerr << "					[-mse <Print: (1: Print output file. Negate to disable)>]" << endl;
		cerr << "					[-propUncert <Print: (1: Print output file. Negate to disable)>]" << endl;
		cerr << "					[-kalman <Print: (1: Print output file. Negate to disable)>]" << endl; 
		cerr << "					[-printMSEwK ] [-appendFilename]" << endl;
		return ARGS_ERROR;
	}else
	{
		outputFileName = argv[++argLocation];
		gridSpacingX_Lon = atof(argv[++argLocation]);

		//1. Check if Y_Lat spacing is provided
		if (isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0])
		{
			gridSpacingY_Lat = atof(argv[++argLocation]);
		}else
		{
			gridSpacingY_Lat = gridSpacingX_Lon;
		}
		//Added 4/17/15 SJZ
		if(gridSpacingX_Lon < 0 || gridSpacingY_Lat < 0)
		{
			cout << "Grid spacing cannot be negative." << endl;
			return ARGS_ERROR;
		}
		if(gridSpacingX_Lon == 0 && gridSpacingY_Lat != 0)
			gridSpacingX_Lon = gridSpacingY_Lat;
		else if(gridSpacingY_Lat == 0 && gridSpacingX_Lon != 0)
			gridSpacingY_Lat = gridSpacingX_Lon;

		kernelName = argv[++argLocation];
		if (!(kernelName == "hann" || kernelName == "hanning" || kernelName == "boxcar" || kernelName == "loess" || kernelName == "quadloess")){
			cout << "Improper argument for kernel passed. Exiting!" << endl;
			return ARGS_ERROR;
		}

		inputFileList = argv[++argLocation];

		refLon = atof(argv[++argLocation]);
		refLat = atof(argv[++argLocation]);
		rotationAngle = atof(argv[++argLocation]);
		numMCRuns = atoi(argv[++argLocation]);
		//SJZ forces single run interpolation
		if(numMCRuns == 0)
			numMCRuns = -1;
		additionalOptions["-numMCRuns"] = numMCRuns;

		//2. Check additional arguments
		for (argLocation=argLocation+1; argLocation < argc; ++argLocation)
		{
			//additionalOptions[argv[argLocation]] = 1;
			//a. ZGrid will be used
			if (strcmp(argv[argLocation], "-ZGrid") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -ZGrid. Exiting!" << endl;
					return ARGS_ERROR;
				}
				mbzData.spacingX = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -ZGrid. Exiting!" << endl;
					return ARGS_ERROR;
				}
				mbzData.spacingY = atof(argv[++argLocation]);

				mbzData.z_OutputFileName = argv[++argLocation];

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -ZGrid. Exiting!" << endl;
					return ARGS_ERROR;
				}
				mbzData.tension = atof(argv[++argLocation]);

				mbzData.usage = atoi(argv[++argLocation]);
				if (!(mbzData.usage >= -2 && mbzData.usage <= 2)){
					cout << "Improper argument passed to -ZGrid. Exiting!" << endl;
					return ARGS_ERROR;
				}
			}

			//b. GMT Surface will be used
			else if (strcmp(argv[argLocation], "-GMTSurface") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -GMTSurface. Exiting!" << endl;
					return ARGS_ERROR;
				}
				GMTSurfaceData.spacingX = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -GMTSurface. Exiting!" << endl;
					return ARGS_ERROR;
				}
				GMTSurfaceData.spacingY = atof(argv[++argLocation]);

				GMTSurfaceData.z_OutputFileName = argv[++argLocation];

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -GMTSurface. Exiting!" << endl;
					return ARGS_ERROR;
				}
				GMTSurfaceData.tension = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -GMTSurface. Exiting!" << endl;
					return ARGS_ERROR;
				}
				GMTSurfaceData.scaleFactor = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -GMTSurface. Exiting!" << endl;
					return ARGS_ERROR;
				}
				GMTSurfaceData.alpha = atof(argv[++argLocation]);

				GMTSurfaceData.usage = atoi(argv[++argLocation]);
				if (!(GMTSurfaceData.usage >= -2 && GMTSurfaceData.usage <= 2)){
					cout << "Improper argument passed to -GMTSurface. Exiting!" << endl;
					return ARGS_ERROR;
				}
			}

			//b. ALGSpline Surface will be used
			else if (strcmp(argv[argLocation], "-ALGSpline") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -ALGSpline. Exiting!" << endl;
					return ARGS_ERROR;
				}
				ALGsplineData.spacingX = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -ALGSpline. Exiting!" << endl;
					return ARGS_ERROR;
				}
				ALGsplineData.spacingY = atof(argv[++argLocation]);

				ALGsplineData.z_OutputFileName = argv[++argLocation];

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -ALGSpline. Exiting!" << endl;
					return ARGS_ERROR;
				}
				ALGsplineData.tension = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -ALGSpline. Exiting!" << endl;
					return ARGS_ERROR;
				}
				ALGsplineData.scaleFactor = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -ALGSpline. Exiting!" << endl;
					return ARGS_ERROR;
				}
				ALGsplineData.alpha = atof(argv[++argLocation]);

				ALGsplineData.usage = atoi(argv[++argLocation]);
				if (!(ALGsplineData.usage >= -2 && ALGsplineData.usage <= 2)){
					cout << "Improper argument passed to -ALGSpline. Exiting!" << endl;
					return ARGS_ERROR;
				}
			}

			//c. Irregular grid spacing will be used
			else if (strcmp(argv[argLocation], "-preInterpolatedLocations") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				interpolationLocationsFileName = argv[++argLocation];
				usagePreInterpLocsLatLon = atoi(argv[++argLocation]);
				if (abs(usagePreInterpLocsLatLon) != 1 && abs(usagePreInterpLocsLatLon) != 2){
					cout << "Improper argument passed to -preInterpolatedLocations. Exiting!" << endl;
					return ARGS_ERROR;
				}
				/*metersPreInterpLocsLatLon = atoi(argv[++argLocation]);
				if (abs(usagePreInterpLocsLatLon) != 1){
					cout << "Improper argument passed to -preInterpolatedLocations. Exiting!" << endl;
					return ARGS_ERROR;
				}*/
			}

			//d. Multi-Threading
			else if (strcmp(argv[argLocation], "-multiThread") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				if (!isdigit(argv[argLocation+1][0])){
					cout << "Improper argument passed to -multiThread. Exiting!" << endl;
					return ARGS_ERROR;
				}
				additionalOptions["-multiThread"] = atoi(argv[++argLocation]);
			}

			//e. Meter Smoothing
			else if (strcmp(argv[argLocation], "-msmooth") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -msmooth. Exiting!" << endl;
					return ARGS_ERROR;
				}
				smoothingScaleX = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -msmooth. Exiting!" << endl;
					return ARGS_ERROR;
				}
				smoothingScaleY = atof(argv[++argLocation]);
			}

			//f. Lat, Lon Smoothing
			else if (strcmp(argv[argLocation], "-llsmooth") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -llsmooth. Exiting!" << endl;
					return ARGS_ERROR;
				}
				smoothingScaleX = atof(argv[++argLocation]);

				if (!isdigit(argv[argLocation+1][0]) && '.' != argv[argLocation+1][0]){
					cout << "Improper argument passed to -llsmooth. Exiting!" << endl;
					return ARGS_ERROR;
				}
				smoothingScaleY = atof(argv[++argLocation]);

				smoothingScaleX = (deg2km*smoothingScaleX)*1000.0*cos(PI*refLat/180.0);
				smoothingScaleY = (deg2km*smoothingScaleY)*1000.0;
			}

			//g. Using Lat/Long Grid Spacing
			else if (strcmp(argv[argLocation], "-llgrid") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				gridSpacingX_Lon = (deg2km*gridSpacingX_Lon)*1000.0*cos(PI*refLat/180.0);
				gridSpacingY_Lat = (deg2km*gridSpacingY_Lat)*1000.0;
			}

			//h. Bounding Box
			else if (strcmp(argv[argLocation], "-boundingBox") == 0)
			{
				additionalOptions[argv[argLocation]] = 1;
				bbox.doBoundingBox = true;

				bbox.bboxTop = atof(argv[++argLocation]);
				bbox.bboxBottom = atof(argv[++argLocation]);
				bbox.bboxRight = atof(argv[++argLocation]);
				bbox.bboxLeft = atof(argv[++argLocation]);

				if (bbox.bboxTop <= bbox.bboxBottom)
				{
					cout << "Improper argument passed to -boundingBox. Exiting!" << endl;
					return ARGS_ERROR;
				}

				if (bbox.bboxRight <= bbox.bboxLeft)
				{
					cout << "Improper argument passed to -boundingBox. Exiting!" << endl;
					return ARGS_ERROR;
				}
			}
			
			////h. Don't compute offset; Default will compute offset.
			//else if (strcmp(argv[argLocation], "-noOffset") == 0)
			//{
			//	additionalOptions["-computeOffset"] = 0;
			//}

			//i. Compute offset; Default will not compute offset.
			else if (strcmp(argv[argLocation], "-computeOffset") == 0)
				additionalOptions["-computeOffset"] = 1;

			//j. Don't use data errors; Default will use errors.
			else if (strcmp(argv[argLocation], "-noerr") == 0)
				additionalOptions["-noerr"] = 1;

			//k. nmsei
			else if (strcmp(argv[argLocation], "-nmsei") == 0)
				additionalOptions["-nmsei"] = 1;
			
			//l. msri
			else if (strcmp(argv[argLocation], "-msri") == 0)
				additionalOptions["-msri"] = 1;

			//m. nmsei
			else if (strcmp(argv[argLocation], "-inputInMeters") == 0)
				additionalOptions["-inputInMeters"] = 1;

			//n. kriging
			else if (strcmp(argv[argLocation], "-kriging") == 0)
				additionalOptions["-kriging"] = 1;

			//o. Use all points; Default will remove bad points.
			else if (strcmp(argv[argLocation], "-modelflag") == 0)
				additionalOptions["-modelflag"] = 1;

			//p. Don't use negative depths, water only; Default uses negative depths
			else if (strcmp(argv[argLocation], "-nonegdepth") == 0)
				additionalOptions["-nonegdepth"] = 1;
			
			//q. Print output file with MSE followed by K
			else if (strcmp(argv[argLocation], "-useUnscaledAvgInputs") == 0)
				additionalOptions["-useUnscaledAvgInputs"] = 1;
	
			//r. Use Nearest Neighbor interpolation; Default uses bilinear interpolation
			else if (strcmp(argv[argLocation], "-nnInterp") == 0)
				additionalOptions["-nnInterp"] = 1;
			
			//s. Use Mean Square Error (linear) Estimator
			else if (strcmp(argv[argLocation], "-mse") == 0)
			{
				additionalOptions["-mse"] = atoi(argv[++argLocation]);
				if (!(abs(additionalOptions.find("-mse")->second) == 1))
				//if (!(additionalOptions.find("-mse")->second >= -1 && additionalOptions.find("-mse")->second <= 1))
				{
					cout << "Improper argument passed to -mse. Exiting!" << endl;
					return ARGS_ERROR;
				}
			}

			//t. Use Propagated Uncertainty Estimator
			else if (strcmp(argv[argLocation], "-propUncert") == 0)
			{
				additionalOptions["-propUncert"] = atoi(argv[++argLocation]);
				if (!(additionalOptions.find("-propUncert")->second >= -1 && additionalOptions.find("-propUncert")->second <= 1))
				{
					cout << "Improper argument passed to -propUncert. Exiting!" << endl;
					return ARGS_ERROR;
				}
			}
			
			//u. Use Kalman Estimator
			else if (strcmp(argv[argLocation], "-kalman") == 0)
			{
				additionalOptions["-kalman"] = atoi(argv[++argLocation]);
				if (!(additionalOptions.find("-kalman")->second >= -1 && additionalOptions.find("-kalman")->second <= 1))
				{
					cout << "Improper argument passed to -kalman. Exiting!" << endl;
					return ARGS_ERROR;
				}
			}
				
			//v. Print output file with MSE followed by K
			else if (strcmp(argv[argLocation], "-printMSEwK") == 0)
			{
				additionalOptions["-printMSEwK"] = 1;
				if (additionalOptions.find("-mse")->second == 0)
				{
					additionalOptions["-mse"] = -1;
				}
				if (additionalOptions.find("-kalman")->second == 0)
				{
					additionalOptions["-kalman"] = -1;
				}
			}
			
			//w. Format output Files to match MATLAB; Default uses original output format
			else if (strcmp(argv[argLocation], "-printMatlabMatch") == 0)
				additionalOptions["-printMatlabMatch"] = 1;

			//x. Print output file with MSE followed by K
			else if (strcmp(argv[argLocation], "-printMatlabMatch") == 0)
				additionalOptions["-printMatlabMatch"] = 1;

			//y. Print output file with MSE followed by K
			else if (strcmp(argv[argLocation], "-outputRasterFile") == 0)
				additionalOptions["-outputRasterFile"] = 1;

			//z. Print output file with MSE followed by K
			else if (strcmp(argv[argLocation], "-outputBagFile") == 0)
				additionalOptions["-outputBagFile"] = 1;
			
			//aa. Append filenames
			else if (strcmp(argv[argLocation], "-appendFilename") == 0)
				additionalOptions["-appendFilename"] = 1;
			
			//bb.Unrecognized parameter
			else 
				unrecognizedParams.push_back(argv[argLocation]);
		}
		if(!unrecognizedParams.empty())
		{
			for(vector<string>::iterator it = unrecognizedParams.begin(); it != unrecognizedParams.end();  it++)
				cout << "ERROR: Unrecognized parameter: " << *it << endl;
			return ARGS_ERROR;
		}
	}
	//************************************************************************************
	//II. Print Banner and options.
	//************************************************************************************
	cout << "\nInitializing .... " << OBF_VERSION_NUMBER << endl << endl;;

	cout << "Using the " << kernelName << " Weighting Window" << endl;
	cout << "Using Reference Longitude: " << refLon << endl;
	cout << "Using Reference Latitude: " << refLat << endl;
	cout << "Using X Grid Spacing (meters): " << gridSpacingX_Lon << endl;
	cout << "Using Y Grid Spacing (meters): " << gridSpacingY_Lat << endl;

	if (smoothingScaleX == -1)
	{
		smoothingScaleX = gridSpacingX_Lon;
		smoothingScaleY = gridSpacingY_Lat;
	}
	else
	{
		//Added 4/17/15 SJZ
		if(smoothingScaleX < gridSpacingX_Lon || smoothingScaleY < gridSpacingY_Lat)
		{
			cout << "Smoothing scales cannot be less than grid spacing." << endl;
			return ARGS_ERROR;
		}
	}

	cout << "Using X Smoothing Scale (meters): " << smoothingScaleX << endl;
	cout << "Using Y Smoothing Scale (meters): " << smoothingScaleY << endl;

	if (additionalOptions.find("-kriging")->second == 1)
	{
		cout << "Computing using Kriging" << endl;
	}
	if (additionalOptions.find("-ZGrid")->second == 1)
	{
		cout << "Using ZGrid in Interpolation" << endl;
	}
	if (additionalOptions.find("-GMTSurface")->second == 1)
	{
		cout << "Using GMT Surface in Interpolation" << endl;
	}
	if (additionalOptions.find("-ALGSpline")->second == 1)
	{
		cout << "Using ALG Spline Surface in Interpolation" << endl;
	}
	if (additionalOptions.find("-preInterpolatedLocations")->second == 1)
	{
		cout << "Computing to Pre-Defined Spacing" << endl;
	}
	if (additionalOptions.find("-multiThread")->second != 0)
	{
		cout << "Using " << additionalOptions.find("-multiThread")->second << " Cores in Interpolation" << endl;
	}
	if (additionalOptions.find("-inputInMeters")->second == 1)
	{
		cout << "Input data in (x,y) meters instead of (lon, lat); no UTM conversions will be computed." << endl;
	}
	if (additionalOptions.find("-computeOffset")->second == 1)
	{
		cout << "Compute offset." << endl;
	}
	if (additionalOptions.find("-noerr")->second == 1)
	{
		cout << "Don't use data errors in computation; Default uses data errors." << endl;
	}
	if (additionalOptions.find("-modelflag")->second == 1)
	{
		cout << "Do not NaN depths in the end. No bad points!" << endl;
	}
	if (additionalOptions.find("-nonegdepth")->second == 1)
	{
		cout << "Do not allow negative depths." << endl;
	}
	if (additionalOptions.find("-nnInterp")->second == 1)
	{
		cout << "Using Nearest Neighbor Interpolation for pre-splining depths." << endl;
	}
	if (abs(additionalOptions.find("-mse")->second) == 1)
	{
		cout << "Using MSE (linear) Estimator." << endl;
		if( additionalOptions.find("-mse")->second == -1)
			cout << "...MSE output file disabled." << endl;
	}
	if (abs(additionalOptions.find("-propUncert")->second) == 1)
	{
		cout << "Use Propagated Uncertainty Estimator." << endl;
		if( additionalOptions.find("-propUncert")->second == -1)
			cout << "...Propagated Uncertainty output file disabled." << endl;
	}
	if (abs(additionalOptions.find("-kalman")->second) == 1)
	{
		cout << "Using Kalman Estimator." << endl;
		if( additionalOptions.find("-kalman")->second == -1)
			cout << "...Kalman output file disabled." << endl;
	}
	//Additional printouts 
	/*if (additionalOptions.find("-printMatlabMatch")->second == 1)
	{
		cout << "Match MATLAB output." << endl;
	}
	if (additionalOptions.find("-printMSEwK")->second == 1)
	{
		cout << "Print MSE output file with K results tagged on the end ." << endl;
	}*/

	//1. Let the user know how many data runs are to be used
	if (numMCRuns == -1)
	{
		cout << "Performing Single Interpolation" << endl << endl;
	}else
	{
		cout << "Performing " << numMCRuns << " Monte Carlo Interpolations" << endl << endl;
	}

	//Alert user of warnings in using mergeBathy 4.0
	if(additionalOptions.find("-kriging")->second == 1 || numMCRuns != -1){
		string ans;
		cout << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "! WARNING!!!!" << endl;
		cout << "! You are entering UNSAFE code!" << endl;
		cout << "! The code path taken for Kriging and Monte Carlo have not been tested!" << endl; 
		cout << "! Enter at you own risk! "<< endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

		/*cout << "Do you want to continue? [(n)/y]" << endl;
		cin >> ans;
		if(ans != "y"){
			"Exiting...";
			return !SUCCESS;
		}*/
	}

	//Alert user of warnings in using mergeBathy 4.0
	if(additionalOptions.find("-outputBagFile")->second == 1){
		string ans;
		cout << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
		cout << "! WARNING!!!!" << endl;
		cout << "! BAG format not available on x64 Debug!" << endl;
		cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
	}


	//Set default estimator to MSE if not specified
	if((additionalOptions.find("-mse")->second == 0 && additionalOptions.find("-propUncert")->second == 0) && additionalOptions.find("-kalman")->second == 0)
		additionalOptions["-mse"] = 1;

	/*if(abs(usagePreInterpLocsLatLon)==2 && additionalOptions.find("-inputInMeters")->second == 0)
	{
		cout << "Improper argument passed to -preInterpolatedLocations. Input must be in meters too. Exiting!" << endl;
		return ARGS_ERROR;	
	}*/
	//Set default to performing Ensemble Tests for dissertation.  This was a quick fix so that prespline algorithms will interpolate to the same grid as specified in the preinterpLocs file
	char buffer[_MAX_PATH]; string wrkdir;
	/* Get the current working directory: */
	if (_getcwd(buffer, _MAX_PATH) == NULL)
		perror("_getcwd error");
	else {
		//printf("%s\n", buffer);
		cout << "Working dir: " << buffer << endl << endl;
		wrkdir = buffer;
	}
	std::size_t found = inputFileList.find("Ensemble");
	std::size_t found2 = wrkdir.find("Ensemble");
	if (found != std::string::npos || found2 != std::string::npos)
		additionalOptions["-ensembleTests"] = 1;
	else
		additionalOptions["-ensembleTests"] = 0;


	//************************************************************************************
	//II. Read the files
	//************************************************************************************
	//1. Read the input files
	start = clock();
	returnValue = readFile(inputFileList, &inputData, bbox, additionalOptions.find("-noerr")->second, additionalOptions.find("-nonegdepth")->second,additionalOptions.find("-useUnscaledAvgInputs")->second);

	if (returnValue != SUCCESS)
	{
		cerr << "FAILED TO READ FILES! ABORTING!!!" << endl;
		return returnValue;
	}

	//2. Read in the interpolated locations
	if (additionalOptions.find("-preInterpolatedLocations")->second == 1)
	{
		returnValue = readLocationsFile(interpolationLocationsFileName, &forcedLocPositions, usagePreInterpLocsLatLon);
		if (returnValue != SUCCESS)
		{
			cerr << "FAILED TO READ INTERPOLATION LOCATION FILE! ABORTING!!!" << endl;
			return returnValue;
		}
	}

	//************************************************************************************
	//III. Start timer and main subroutine for data interpolation.
	//************************************************************************************
	start = clock();

	returnValue = mergeBathy_PreCompute(&inputData, refLon, refLat, rotationAngle, gridSpacingX_Lon, gridSpacingY_Lat, smoothingScaleX, smoothingScaleY, kernelName, outputFileName, additionalOptions, numMCRuns, &mbzData, &GMTSurfaceData, &ALGsplineData, &forcedLocPositions,usagePreInterpLocsLatLon);

	if (returnValue != SUCCESS)
	{
		cerr << "INTERPOLATION FAILED! ABORTING!!!" << endl;
	}

	//************************************************************************************
	//IV. Clean up input structures.
	//************************************************************************************
	for (int i = 0; i < (const int)inputData.size(); i++)
	{
		inputData[i].lon.clear();
		inputData[i].lat.clear();
		inputData[i].depth.clear();
		inputData[i].error.clear();
		inputData[i].time.clear();
		inputData[i].h_Error.clear();
		inputData[i].v_Error.clear();
		inputData[i].x.clear();
		inputData[i].y.clear();
	}
	stop = clock();
	cout << "\nTotal Computation Time: " << (double(stop)-double(start))/CLOCKS_PER_SEC << " seconds" << endl;

	return returnValue;
}
