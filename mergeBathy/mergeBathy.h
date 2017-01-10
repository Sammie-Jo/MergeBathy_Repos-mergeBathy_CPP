/** 
* @file			mergeBathy.h
* @brief		Header file for the main function of mergeBathy.
* @author		Kevin Duvieilh
* @date			04 August 2011
*
*/
#pragma once
#include <iostream>
#include <string>
#include <map>
#include "constants.h"
#include "fileReader.h"
#include "mergeBathyOld.h"

/**
* @fn main
* mergeBathy main function.  The main function takes 8 mandatory command line arguments and is capable of having additional arguments for separate processing.
* Arguments encapsulated by <> are mandatory arguments. 
* Arguments encapsulated by [] are optional arguments.
* Some optional arguments have mandatory arguments that must be passed after the optional argument is invoked.  These mandatory arguments are indented below the optional argument that requries them.
* 
* <output_depth_file> - A text file consisting of one line for each output point at which a valid depth was computed.  Each line consists of four columns, corresponding to Longitude, Latitude, Depth, and Uncertainty (the mean square interpolation error estimate) respectively.  Additional columns can be present depending on other additional command line options provided at runtime.
* <grid_spacing> - Computational grid spacing in the X direction defined in meters.
* [grid_spacing_Y] - Computational grid spacing in the Y direction defined in meters.  This argument is optional and if it is not provided the grid spacing in the Y direction will be assumed to be the same as the spacing in the X direction.
* <kernel_name> - Interpolation window (i.e. filter window) to use, the choices are 'quadloess', 'linloess', 'hanning' or 'hann', and 'boxcar'. Hanning is the standard smoothing window.
* <input_file_list> - Name of a textfile containing names of individual data files containing data to be merged.  Each line of the listfile should have a single filename.
* <ref_lon> - Center longitude of the given computational area.
* <ref_lat> - Center latitude of the given computational area.
* <rotation_angle> - An angle [in degrees] of arbitrary rotation.  This allows the computation X and Y space to undergo a rotation, such that the along-shore or cross-shore direction can be oriented along a given axis.  This allows the user to more easily select unequal arbitrary smoothing length scales.
* <num_MC_runs(-1 if not MC)> - Integer number of Monte Carlo iterations to perform. Use -1 if no Monte Carlo simulations are to be run.
* [-noerr] - Do not output the uncertainty estimate (the mean square interpolation error estimate) in the output file.
* [-nmsei] - Write out the normalized mean square error estimate to the output file. This appears as an additional column in the ASCII file.
* [-msri] - Write out the mean square of the residuals to the output file. This appears as an additional column in the ASCII file. It is given after the msei-column if both are given.
* [-inputInMeters] - Specifies that all input files are using an X/Y coordinate system in meters instead of the typical Longitude/Latitude in degrees.
* [-kriging] - Indicates that kriging will be used on to correct oversmoothing of the interpolation surface at the input data points. This step restores finer details lost from smoothing alone, but will add to the computation time required.
* [-ZGrid] - Forces the included Z_Grid software package to run.  This function was NOT created by the design team of newMergeBathy but is included under the GNU Public License Agreement.
*		<grid_spacing_X> - Computational grid spacing in the X direction defined in meters that will be used in MB_ZGrid computation.
*		<grid_spacing_Y> - Computational grid spacing in the Y direction defined in meters that will be used in MB_ZGrid computation.
*		<Z_Grid_Output_File_Name> - MB_ZGrid output file.  Each line consists of three columns, corresponding to Longitude, Latitude, and Depth respectively.
*		<Tension_Factor> - Sets the tension of the interpolation.  A value of 0.0 yields a pure Laplace (minimum curvature) solution and a value of infinity yields a pure thin plate spline solution. A value of 1e10 value has commonly been used to yield spline solutions.  
*		<Usage> - A value of 1 will perform the MB_ZGrid interpolation and write the results to the specified output file.  A value of 2 will perform the MB_ZGrid interpolation, write the results to the specified output file, and use the computed X,Y, and Z values as input for mergeBathy.  If used as input then the data will take the place of the data read from the input files.  This allows for a pre-smoothing effect before mergeBathy is run.  Negate the usage values (-1 or -2) to compute the error associated with the pre-smoothing.
* [-GMTSurface] - Forces the included GMT Surface software package to run.  This function was NOT created by the design team of newMergeBathy.
*		<grid_spacing_X> - Computational grid spacing in the X direction defined in meters that will be used in GMT Surface computation.
*		<grid_spacing_Y> - Computational grid spacing in the Y direction defined in meters that will be used in GMT Surface computation.
*		<GMT_Surface_Output_File_Name> - MB_ZGrid output file.  Two versions of this file will be created.  The first creates a file where each line consists of three columns, corresponding to Longitude, Latitude, and Depth respectively.  The second creates a file where each line consists of four columns, corresponding to Longitude, Latitude, Depth, and Error respectively.  The second file has "_includeError.txt" appended to the end of the specified file name.
*		<Tension_Factor> - The Laplacian tension operator between 0 and 1.  Typcially 0.1
*		<scale_factor> - The multiplier value for a Confidence Interval to be used in error calculation. A value of 1.96 is typically used for a 95% Confidence Interval.
*		<alpha> - This is the alpha value for error computation. Typically 2.0.
*		<Usage> - A value of 1 will perform the GMT Surface interpolation and write the results to the specified output file.  A value of 2 will perform the GMT Surface interpolation, write the results to the specified output file, and use the computed X,Y,Z, and E values as input for mergeBathy.  If used as input then the data will take the place of the data read from the input files.  This allows for a pre-smoothing effect before mergeBathy is run.  Negate the usage values (-1 or -2) to compute the error associated with the pre-smoothing.
* [-preInterpolatedLocations] - Allows a file to be provided to determine the exact Longitude and Latitude locations for the interpolation.  This allows irregularly gridded output at set points instead of points determined from a user defined grid spacing.  Using this option slows down computation speed significantly and therefore it is only recommended for small data sets; however, this option does function properly on large data sets.
*		<interpolation_location_file_name> - The file name of a file containing Longitude and Latitude columns of data points that specify output interpolation grid locations.
* [-computeOffset] - Computes the offset between multiple data sets and normalizes depth.
* [-mse] - Perform Mean Square Error Estimator.
*		<Print> - A value of 1 to write results to an output file.  A value of -1 to disable output file printout.
* [-propUncert] - Perform Propagated Uncertainty Estimator.
*		<Print> - A value of 1 to write results to an output file.  A value of -1 to disable output file printout.
* [-kalman] - Perform Kalman Estimator.
*		<Print> - A value of 1 to write results to an output file.  A value of -1 to disable output file printout.
* [-nnInterp] - Perform Nearest Neighbor interpolation when pre-splining.  Default performs bilinear interpolation.
* [-printMatlabMatch] - print output file with results formatted to match Matlab's output file.
*/

