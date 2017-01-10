/** 
* @file			inFileStructs.h
* @brief		Define input data structures going from mergeBathy to mergeBathyOld.
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/

#ifndef INFILESTRUCTS_H_
#define INFILESTRUCTS_H_

#include <vector>
#include <string>

using namespace std;

/**
* This structure is used to define the input data points read in from the files to be interpolated.
*/
typedef struct
{
	/**
	* Found utm zone.
	*/
	char utmzone[4]; // SJZ 
//	char *utmzone;

	/**
	* Known longitude value read from a file.
	*/
	vector<double> lon;

	/**
	* Known latitude value read from a file.
	*/
	vector<double> lat;

	/**
	* Known depth value read from a file.
	*/
	vector<double> depth;

	/**
	* Known combined error value read from a file.
	*/
	vector<double> error;

	/**
	* Known time value read from a file.
	*/
	vector<double> time;

	/**
	* Known horizontal error value read from a file.
	*/
	vector<double> h_Error;

	/**
	* Known vertical error value read from a file.
	*/
	vector<double> v_Error;

	/**
	* Longitude coordinate converted to meters.
	*/
	vector<double> x;

	/**
	* Latitude coordinate converted to meters.
	*/
	vector<double> y;

	/**
	* The sum of all the longitude points.
	*/
	double longitudeSum;
	
	/**
	* The sum of all the latitude points.
	*/
	double latitudeSum;// SJZ

	/**
	* The maximum offset that may be applied to a data set.  User defined value that may or may not be passed.
	*/
	double maximumDataOffset;

} TRUE_DATA;

/**
* This structure is used to define all values necessary to run MB_ZGrid.
*/
typedef struct
{
	/**
	* User defined output file name for the MB_ZGrid data.
	*/
	string z_OutputFileName;
	
	/**
	* MB_ZGrid specific spacing in the X direction.
	*/
	double spacingX;

	/**
	* MB_ZGrid specific spacing in the Y direction.
	*/
	double spacingY;
	
	/**
	*  Tension factor for MB_ZGrid computation (1e10 is default).
	*/
	double tension;
	
	/**
	* Just do the computation (0) or use as input to mergeBathy (1).
	*/
	int usage;

} MB_ZGRID_DATA;

/**
* This structure is used to define all values necessary to run GMT Surface.
*/
typedef struct
{
	/**
	* User defined output file name for the GMT_Surface data.
	*/
	string z_OutputFileName;

	/**
	* GMT_Surface specific spacing in the X direction.
	*/
	double spacingX;

	/**
	* GMT_Surface specific spacing in the Y direction.
	*/
	double spacingY;

	/**
	* Tension factor for GMT_Surface computation (Between 0 and 1).
	*/
	double tension;

	/**
	* The scale factor value for computing error (1.96 is default).
	*/
	double scaleFactor;

	/**
	* The alpha value for computing error (2.0 is default).
	*/
	double alpha;

	/**
	* Just do the computation (0) or use as input to mergeBathy (1).
	*/
	int usage;

} GMT_SURFACE_DATA;
/**
* This structure is used to define all values necessary to run ALGSpline.
*/
typedef struct
{
	/**
	* User defined output file name for the GMT_Surface data.
	*/
	string z_OutputFileName;

	/**
	* GMT_Surface specific spacing in the X direction.
	*/
	double spacingX;

	/**
	* GMT_Surface specific spacing in the Y direction.
	*/
	double spacingY;

	/**
	* Tension factor for GMT_Surface computation (Between 0 and 1).
	*/
	double tension;

	/**
	* The scale factor value for computing error (1.96 is default).
	*/
	double scaleFactor;

	/**
	* The alpha value for computing error (2.0 is default).
	*/
	double alpha;

	/**
	* Just do the computation (0) or use as input to mergeBathy (1).
	*/
	int usage;

} ALG_SPLINE_DATA;
/**
* This structure is used to define all values necessary to compute an irregular grid.
* This structure is used to define a specific set of data points to be interpolated
* These points will be used instead of computed points
* Especially useful for re-gridding to known points to a known data file
*/
typedef struct
{
	/**
	* Found utm zone.
	*/
	char utmzone[4]; // SJZ 
//	char *utmzone; 

	/**
	* Vector of longitude coordinates.
	*/
	vector<double> forcedLonCoord;

	/**
	* Vector of latitude coordinates.
	*/
	vector<double> forcedLatCoord; 

	/**
	* The sum of all the longitude points.
	*/
	double longitudeSum;
	
	/**
	* The sum of all the latitude points.
	*/
	double latitudeSum; // SJZ

} FORCED_LOCATIONS;

/**
* This structure is used to define a bounding box used to cut out data from an input grid.
*/
typedef struct
{
	/**
	* Boolean to determine if we do the boundingBox.
	*/
	bool doBoundingBox;

	/**
	* Top of the bounding box.
	*/
	double bboxTop;

	/**
	* Bottom of the bounding box.
	*/
	double bboxBottom;

	/**
	* Right of the bounding box.
	*/
	double bboxRight;

	/**
	* Left of the bounding box.
	*/
	double bboxLeft;
	
} BOUNDING_BOX;

#endif
