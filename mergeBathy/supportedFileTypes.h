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

