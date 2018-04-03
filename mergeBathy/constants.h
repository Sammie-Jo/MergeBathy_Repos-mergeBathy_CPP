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
/**
* @file			constants.h
* @brief		Defines all standardized constants and return codes to be used throughout mergeBathy.
* @author		Kevin Duvieilh
* @date			02 August 2011
*
*/

#pragma once
#include <fstream>

//#ifdef WIN32
//const double NaN = std::numeric_limits<double>::quiet_NaN();
//#else
//unsigned long nan2[2]={0xffffffff, 0x7fffffff};
//const double NaN = *( double* )nan2;
////const int NaN = MAX_INT; //ALP for UNIX
//#endif


#ifndef INT_MIN //limits.h microsoft specific define for linux SJZ 11/21/14
#define INT_MIN     (-2147483647 - 1)  /* minimum (signed) int value */
#define INT_MAX       2147483647    /* maximum (signed) int value */
#endif

#ifndef MAX_INT
#define MAX_INT 9999999	//Denson says this is platform specific be careful
#define MIN_INT -9999999
#endif

//Interpolation method for depth and uncertainty
const int DEBUG_DISABLE_PTLOFFSET = 1;

const int NaN = MAX_INT; //ALP for UNIX
const double RCOND_TOL = 1e-6;  //Reciprocal condition test threshold

/**
* Version Number.
*/
const static char *OBF_VERSION_NUMBER = "BUILD 5.0.2: July 23, 2015";

/**
* PI.
*/
const double PI = 3.1415926535897932384626433832795028841971693993751058209749;

/**
* Value used to convert degrees to radians.
*/
const double deg2rad = PI / 180.00000000000000000000000;

/**
* Value used to convert radians to degrees.
*/
const double rad2deg = 180.00000000000000000000000 / PI;

/**
* Value used to convert degrees to kilometers.
*/
const double deg2km = 111.1949266445587;

/**
* divide lambda by this for scale.
*/
const double Dscale = 4.0000; 

/**
* Tolerance of minimum normalize error acceptable to not expand the smoothing scale. Otherwise, expand it.
*/
const double NEITOL = 0.5000;

/**
* Tolerance of minimum normalize error acceptable to not expand the smoothing scale. Otherwise, expand it.
*/
const double NEITOL_COMPUTE_OFFSET = 1.0000;

const double RMAX = 1.0000;   
const double DR = 0.0100;
const double eps = 2.220446049250313e-016;

/**
* Sets the maximum amount of data points to Krig before having to further subtile.
*/
const int KRIGED_SIZE_THRESHOLD = 64;

/**
* Input argument error. -1.
*/
const int ARGS_ERROR = -1;

/**
* Successful return value. 0.
*/
const int SUCCESS = 0;

/**
* There was a problem opening or reading the input list file. 1.
*/
const int LIST_FILE_ERROR = 1;

/**
* There was a problem opening or reading the input list file. 11.
*/
const int LIST_FILE_EXT_ERROR = 11; //Missing identifiable extension

/**
* There was a problem opening or reading an individual data input file. 2.
*/
const int IN_FILE_ERROR = 2;

/**
* There was a problem opening or writing the output file. 3.
*/
const int OUT_FILE_ERROR = 3;

/**
* There was a problem computing the offset between input data sets. 4.
*/
const int COMPUTE_OFFSET_ERROR = 4;

/**
* There was a problem using MB ZGrid to pre-spline the data. 5.
*/
const int MB_ZGRID_ERROR = 5;

/**
* There was a problem using GMT Surface or the GMT Surface Error Estimator to pre-spline the data. 6.
*/
const int GMT_SURF_ERROR = 6;

/**
* Error return code for ALG external interpolator; currently bilinear and bicubic spline. 7.
*/
const int ALG_SPLINE_ERROR = 7;

/**
* Default E value to use if E = 0 for divide by 0 exception handling in scalecInterpTile, consistentWeights, and subSampleData.
*/
const double DEFAULT_E = 1e-6;//0.000001;

