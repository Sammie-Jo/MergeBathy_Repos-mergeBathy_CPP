/** 
* @file			supportedFileTypes.h
* @brief		Defines all filetypes that are supported by fileReader
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/

#pragma once
#include <string>

/**
* Define flat text tile extensions. Allow _xyze.dat files.
* Contents of the files are:
*	@n error_DAT and error_TXT:
*		Longitude	Latitude	Depth	Error
*/
const string error_DAT = "_xyze.dat";

/**
* Define flat text tile extensions. Allow _xyzhv_mc.dat files.
* Contents of the files are:
*	@n hv_error_DAT and hv_error_TXT:
*		Longitude	Latitude	Depth	Horizontal_Error	Vertical_Error
*/
const string hv_error_DAT = "_xyzhv_mc.dat";

/**
* Define flat text tile extensions. Allow _xyz.dat files.
* Contents of the files are:
*	@n no_error_DAT and no_error_TXT:
*		Longitude	Latitude	Depth
*/
const string no_error_DAT = "_xyz.dat";

/**
* Define flat text tile extensions. Allow _xyde.txt files.
* Contents of the files are:
*	@n error_DAT and error_TXT:
*		Longitude	Latitude	Depth	Error
*/
const string error_TXT = "_xyde.txt";

/**
* Define flat text tile extensions. Allow _xydhv_mc.txt files.
* Contents of the files are:
*	@n hv_error_DAT and hv_error_TXT:
*		Longitude	Latitude	Depth	Horizontal_Error	Vertical_Error
*/
const string hv_error_TXT = "_xydhv_mc.txt";

/**
* Define flat text tile extensions.  Allow _xyd.txt files.
* Contents of the files are:
*	@n no_error_DAT and no_error_TXT:
*		Longitude	Latitude	Depth
*/
const string no_error_TXT = "_xyd.txt";

/**
* Define GSF extensions. Allow .d files.
* Contents of the files follow GSF file standards.
*/
const string GSF = ".d";

/**
* Define GSF extensions. Allow .gsf files.
* Contents of the files follow GSF file standards.
*/
const string GSF2 = ".gsf";

/**
* Define ARC ASCII Raster extensions. Allow .grd files.
* Contents of the file follow ARC ASCII Raster standards.
* @n The header of the file should be established as follows.
*	@n ncols = number of columns
*	@n nrows = number of rows
*	@n xllcorner = lower left X coordinate
*	@n yllcorner = lower left Y coordinate
*	@n cellsize = size of the data cell 
*	@n nodata = flag used to signify no data
*/
const string raster = ".grd";

/**
* Define ARC ASCII Raster extensions. Allow .asc files.
* Contents of the file follow ARC ASCII Raster standards.
* @n The header of the file should be established as follows.
*	@n ncols = number of columns
*	@n nrows = number of rows
*	@n xllcorner = lower left X coordinate
*	@n yllcorner = lower left Y coordinate
*	@n cellsize = size of the data cell 
*	@n nodata = flag used to signify no data
*/
const string raster2 = ".asc";

