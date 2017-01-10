/** 
* @file			externalInterpolators.h
* @brief		Define class used for calling external interpolators such as MB_ZGrid and GMT Surface.
* @author		Kevin Duvieilh
* @date			16 June 2011
*
*/

#pragma once
#include <string>
#include <vector>
#include <map>
#include "Error_Estimator/Bathy_Grid.h"
using namespace std;

struct COLUMN_RNG_STRUCT
{
	int begin;
	int end;
};

class externalInterpolators
{

public:
	/**
	* A constructor for the external interpolator.
	* @param NorthingRef - The UTM Northing reference point.
	* @param EastingRef - The UTM Easting reference point.
	* @param rotAngle - An angle [in degrees] of arbitrary rotation.  This allows the computation X and Y space to undergo a rotation, such that the along-shore or cross-shore direction can be oriented along a given axis.  This allows the user to more easily select unequal arbitrary smoothing length scales.
	* @param RefEllip - The reference ellipsoid used for converting from Longitude and Latitude to UTM.
	* @param ZoneRef - The UTM Zone referece character.
	*/
	externalInterpolators(double NorthingRef, double EastingRef, double rotAngle, int RefEllip, char ZoneRef[4]);

	/*SIBSON non working!
	double vector_max(std::vector<double>);
	double vector_min(std::vector<double>);
	void grid_segmentor(double num_columns,double num_threads,std::vector<double>& col_strip_start,std::vector<double>& col_strip_end);
	int build_input_grid_locs(const std::vector<double>&,const std::vector<double>&,const std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,double, double);
	int build_output_grid_locs(const std::vector<double>&,const std::vector<double>&,double,double);
	DWORD WINAPI compute_sibson_ntrlnbr_cont_inter_col_thread(LPVOID);
	*/

	/**
	* Run the MB_ZGrid routine.
	* @param x - Vector of input X data points.
	* @param y - Vector of input Y data points.
	* @param z - Vector of input Depth data points.
	* @param e - Vector of input Error data points.
	* @param x0 - Minimum value contained in the x vector.
	* @param y0 - Minimum value contained in the y vector.
	* @param x1 - Maximum value contained in the x vector.
	* @param y1 - Maximum value contained in the y vector.
	* @param additionalOptions - std::map container of additional options flags set in the main function. 
	* @param spacingX - Spacing of the MB_ZGrid computational area in the X direction.
	* @param spacingY - Spacing of the MB_ZGrid computational area in the Y direction.
	* @param tension - Sets the tension of the interpolation.  A value of 0.0 yields a pure Laplace (minimum curvature) solution and a value of infinity yields a pure thin plate spline solution. A value of 1e10 value has commonly been used to yield spline solutions.  
	* @param z_OutputFileName - MB_ZGrid output file.  Each line consists of three columns, corresponding to Longitude, Latitude, and Depth respectively.
	* @param usage - Sets how MB_ZGrid output will be used. A value of 0 will perform the MB_ZGrid interpolation and write the results to the specified output file.  A value of 1 will perform the MB_ZGrid interpolation, write the results to the specified output file, and use the computed X,Y, and Z values as input for mergeBathy.  If used as input then the data will take the place of the data read from the input files.  This allows for a pre-smoothing effect before mergeBathy is run.
	* @return Success or failure value.
	*/
	bool run_MB_ZGrid(vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *hError, vector<double> *v, double x0, double y0, double x1, double y1, map<string, int> additionalOptions, double spacingX, double spacingY, double tension, string z_OutputFileName, int usage, Bathy_Grid* bathyGrid);

	/**
	* Run the MB_ZGrid routine.
	* @param x - Vector of input X data points.
	* @param y - Vector of input Y data points.
	* @param z - Vector of input Depth data points.
	* @param e - Vector of input Error data points.
	* @param hError - Vector of input Horizontal Error data points.
	* @param v - Vector of input Vertical Error data points.
	* @param x0 - Minimum value contained in the x vector.
	* @param y0 - Minimum value contained in the y vector.
	* @param x1 - Maximum value contained in the x vector.
	* @param y1 - Maximum value contained in the y vector.
	* @param additionalOptions - std::map container of additional options flags set in the main function. 
	* @param spacingX - Spacing of the MB_ZGrid computational area in the X direction.
	* @param spacingY - Spacing of the MB_ZGrid computational area in the Y direction.
	* @param tension - Sets the tension of the interpolation.  A value of 0.0 yields a pure Laplace (minimum curvature) solution and a value of infinity yields a pure thin plate spline solution. A value of 1e10 value has commonly been used to yield spline solutions.  
	* @param z_OutputFileName - MB_ZGrid output file.  Each line consists of three columns, corresponding to Longitude, Latitude, and Depth respectively.
	* @param scaleFactor - The multiplier value for a Confidence Interval to be used in error calculation. A value of 1.96 is typically used for a 95% Confidence Interval.
	* @param alpha - This is the alpha value for error computation. Typically 2.0.
	* @param usage - Sets how MB_ZGrid output will be used. A value of 0 will perform the MB_ZGrid interpolation and write the results to the specified output file.  A value of 1 will perform the MB_ZGrid interpolation, write the results to the specified output file, and use the computed X,Y, and Z values as input for mergeBathy.  If used as input then the data will take the place of the data read from the input files.  This allows for a pre-smoothing effect before mergeBathy is run.
	* @return Success or failure value.
	*/
	bool run_Surface(vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *hError, vector<double> *v, map<string, int> additionalOptions, double spacingX, double spacingY, double tension, string z_OutputFileName, double scaleFactor, double alpha, int usage, Bathy_Grid* bathyGrid);

	/**
	* Run the ALG Spline routine.
	* 
		vector<double> *pxin_vect	- pointer to input vector x(longitude)
		vector<double> *pyin_vect	- pointer to input vector y(latitude)
		vector<double> *pzin_vect	- pointer to input vector z(depth)
		double spacingX             - input x grid spacing
		double spacingY             - input y grid spacing
		double xsize                - number of input columns
		double ysize                - number of input rows
		string z_OutputFileName     - output filename 
		int usage                   - usage parameter
	*/
	/*bool run_ALGSpline(vector<double>	*pxin_gvect, vector<double>	*pyin_gvect, vector<double>	*pzin_gvect, double utm_xmin,  double utm_xmax, double utm_ymin,  double utm_ymax, double xin_size, double yin_size, double spacingX, double spacingY, string z_OutputFileName, int usage);*/

	bool run_ALGSpline(vector<double>	*pxin_gvect, vector<double>	*pyin_gvect, vector<double>	*pzin_gvect, double utm_xmin,  double utm_xmax, double utm_ymin,  double utm_ymax, double xin_size, double yin_size, vector<double> *e, vector<double> *hError, vector<double> *v, map<string, int> additionalOptions, double spacingX, double spacingY, double tension, string z_OutputFileName, double scaleFactor, double alpha, int usage, Bathy_Grid* bathyGrid);

	bool run_Sibson_Natural_Neighbor_Interp(vector<double> *pxin_gvect, vector<double>	*pyin_gvect, vector<double>	*pzin_gvect, double utm_xmin,  double utm_xmax, double utm_ymin,  double utm_ymax, double xin_size, double yin_size, vector<double> *e, vector<double> *hError, vector<double> *v, map<string, int> additionalOptions, double spacingX, double spacingY, double tension, string z_OutputFileName, double scaleFactor, double alpha, int usage, Bathy_Grid* bathyGrid);

	int write_output_file(string chAbsFilename, const vector<double>& grid_xv, const vector<double>& grid_yv, const vector<double>& grid_zv);



	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// ORIGINAL FUNCTIONS FOR MONTE CARLO RUNS AND KRIGING!!!
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
	/**
	* Run the MB_ZGrid routine.
	* @param x - Vector of input X data points.
	* @param y - Vector of input Y data points.
	* @param z - Vector of input Depth data points.
	* @param e - Vector of input Error data points.
	* @param x0 - Minimum value contained in the x vector.
	* @param y0 - Minimum value contained in the y vector.
	* @param x1 - Maximum value contained in the x vector.
	* @param y1 - Maximum value contained in the y vector.
	* @param spacingX - Spacing of the MB_ZGrid computational area in the X direction.
	* @param spacingY - Spacing of the MB_ZGrid computational area in the Y direction.
	* @param tension - Sets the tension of the interpolation.  A value of 0.0 yields a pure Laplace (minimum curvature) solution and a value of infinity yields a pure thin plate spline solution. A value of 1e10 value has commonly been used to yield spline solutions.  
	* @param z_OutputFileName - MB_ZGrid output file.  Each line consists of three columns, corresponding to Longitude, Latitude, and Depth respectively.
	* @param usage - Sets how MB_ZGrid output will be used. A value of 0 will perform the MB_ZGrid interpolation and write the results to the specified output file.  A value of 1 will perform the MB_ZGrid interpolation, write the results to the specified output file, and use the computed X,Y, and Z values as input for mergeBathy.  If used as input then the data will take the place of the data read from the input files.  This allows for a pre-smoothing effect before mergeBathy is run.
	* @return Success or failure value.
	*/
	//bool run_MB_ZGrid_ORIGINAL(vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *hError, vector<double> *v, double x0, double y0, double x1, double y1, double spacingX, double spacingY, double tension, string z_OutputFileName, int usage);

	/**
	* Run the MB_ZGrid routine.
	* @param x - Vector of input X data points.
	* @param y - Vector of input Y data points.
	* @param z - Vector of input Depth data points.
	* @param e - Vector of input Error data points.
	* @param hError - Vector of input Horizontal Error data points.
	* @param v - Vector of input Vertical Error data points.
	* @param x0 - Minimum value contained in the x vector.
	* @param y0 - Minimum value contained in the y vector.
	* @param x1 - Maximum value contained in the x vector.
	* @param y1 - Maximum value contained in the y vector.
	* @param spacingX - Spacing of the MB_ZGrid computational area in the X direction.
	* @param spacingY - Spacing of the MB_ZGrid computational area in the Y direction.
	* @param tension - Sets the tension of the interpolation.  A value of 0.0 yields a pure Laplace (minimum curvature) solution and a value of infinity yields a pure thin plate spline solution. A value of 1e10 value has commonly been used to yield spline solutions.  
	* @param z_OutputFileName - MB_ZGrid output file.  Each line consists of three columns, corresponding to Longitude, Latitude, and Depth respectively.
	* @param scaleFactor - The multiplier value for a Confidence Interval to be used in error calculation. A value of 1.96 is typically used for a 95% Confidence Interval.
	* @param alpha - This is the alpha value for error computation. Typically 2.0.
	* @param usage - Sets how MB_ZGrid output will be used. A value of 0 will perform the MB_ZGrid interpolation and write the results to the specified output file.  A value of 1 will perform the MB_ZGrid interpolation, write the results to the specified output file, and use the computed X,Y, and Z values as input for mergeBathy.  If used as input then the data will take the place of the data read from the input files.  This allows for a pre-smoothing effect before mergeBathy is run.
	* @return Success or failure value.
	*/
	//bool run_Surface_ORIGINAL(vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *hError, vector<double> *v, double spacingX, double spacingY, double tension, string z_OutputFileName, double scaleFactor, double alpha, int usage);
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




private:
	/**
	* The UTM Northing reference point.
	*/
	double UTMNorthingRef;

	/**
	* The UTM Easting reference point.
	*/
	double UTMEastingRef;

	/**
	* An angle [in degrees] of arbitrary rotation.  This allows the computation X and Y space to undergo a rotation, such that the along-shore or cross-shore direction can be oriented along a given axis.  This allows the user to more easily select unequal arbitrary smoothing length scales.
	*/
	double rotationAngle;

	/**
	* The reference ellipsoid used for converting from Longitude and Latitude to UTM.
	*/
	int refEllipsoid;

	/**
	* The UTM Zone reference character.
	*/
	char UTMZoneRef[4];

};

