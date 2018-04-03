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
* @file			fileReader.h
* @brief		Header file for reading xyz, xyze, xyz_hv, gsf, and ARC ASCII Raster file formats for mergeBathy.
* @author		Kevin Duvieilh
* @date			04 August 2011
*
*/

#ifndef FILEREADER_H_
#define FILEREADER_H_

#include <iostream>
#include <fstream>
#include <string>
#include "constants.h"
#include "inFileStructs.h"

extern "C" {
	#include "GSF/gsfReader.h"
}

using namespace std;

/**
* Read a file list and pass the individual files to be read to the appropriate functions.  Each line of the input file must contain the path to a single input file.  Additional input files must be placed on their own separate row in this file.  Any file read from the file list must match one of the supported filetypes in "supportedFileTypes.h".
* File Format:
*	${PATH_TO_FILES}/file_name_1.txt
*	${PATH_TO_FILES}/file_name_2.txt
*	${PATH_TO_FILES}/file_name_3.txt
*	...
* @param fileName - File name of the list file that contains each individual data set.
* @param inputData - Vector TRUE_DATA that contains the values read from each file.  The length of the vector corresponds to the number of input files.  Each index of the vector corresponds the TRUE_DATA read from a specific input file.  (Returned).
* @param bbox - The bounding extents of the input data if it is to be cut.
* @return Success or failure boolean.
*/
int readFile(string &fileName, vector<TRUE_DATA> *inputData, BOUNDING_BOX bbox, int noerr, int nonegdepth, int unScaledAvgInputs);

/**
* Read a file containing Longitude and Latitude data columns.  Used to grid data to an irregular grid at the Longitude and Latitude coordinates specified in this file.
* File Format:
*	Longitude	Latitude
*	Longitude	Latitude
*	Longitude	Latitude
*	...
* @param fileName - zFile name of the file that contains the Longitude and Latitude coordinates.
* @param inputData - FORCED_LOCATIONS that contains the Longitude and Latitude values read from the input file.  (Returned).
* @return Success or failure boolean.
*/
int readLocationsFile(string &fileName, FORCED_LOCATIONS *inputData,int usagePreInterpLocsLatLon);

/**
* Read a file containing GSF data.
* File Format:
*	Standard GSF file format
* @param fileName - File name of the GSF file to be read.
* @param bbox - The bounding extents of the input data if it is to be cut.
* @param inputData - Vector TRUE_DATA that contains the values read from each file.  Only the index specified by pos is modified.  The TRUE_DATA at inputData[pos] will contain all of data read from the inputFile.  (Returned).
* @param pos - Integer representing the location of inputData where the data from the input file should be stored.
* @return Success or failure boolean.
*/
int readGSF(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth, int unScaledAvgInputs);

/**
* Read a file containing Longitude, Latitude, and Depth data columns. 
* File Format:
*	Longitude	Latitude	Depth
*	Longitude	Latitude	Depth
*	Longitude	Latitude	Depth
*	...
* @param fileName - File name of the XYZ file to be read.
* @param bbox - The bounding extents of the input data if it is to be cut.
* @param inputData - Vector TRUE_DATA that contains the values read from each file.  Only the index specified by pos is modified.  The TRUE_DATA at inputData[pos] will contain all of data read from the inputFile.  (Returned).
* @param pos - Integer representing the location of inputData where the data from the input file should be stored.
* @return Success or failure boolean.
*/
int readXYZ(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth);

/**
* Read a file containing Longitude, Latitude, Depth, and Error data columns. 
* File Format:
*	Longitude	Latitude	Depth	Error
*	Longitude	Latitude	Depth	Error
*	Longitude	Latitude	Depth	Error
*	...
* @param fileName - File name of the XYZE file to be read.
* @param bbox - The bounding extents of the input data if it is to be cut.
* @param inputData - Vector TRUE_DATA that contains the values read from each file.  Only the index specified by pos is modified.  The TRUE_DATA at inputData[pos] will contain all of data read from the inputFile.  (Returned).
* @param pos - Integer representing the location of inputData where the data from the input file should be stored.
* @return Success or failure boolean.
*/
int readXYZE(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth, int unScaledAvgInputs);

/**
* Read a file containing Longitude, Latitude, Depth, Horizontal Error, and Vertical Error data columns. 
* File Format:
*	Longitude	Latitude	Depth	Horizontal_Error	Vertical_Error
*	Longitude	Latitude	Depth	Horizontal_Error	Vertical_Error
*	Longitude	Latitude	Depth	Horizontal_Error	Vertical_Error
*	...
* @param fileName - File name of the XYZ_HV file to be read.
* @param bbox - The bounding extents of the input data if it is to be cut.
* @param inputData - Vector TRUE_DATA that contains the values read from each file.  Only the index specified by pos is modified.  The TRUE_DATA at inputData[pos] will contain all of data read from the inputFile.  (Returned).
* @param pos - Integer representing the location of inputData where the data from the input file should be stored.
* @return Success or failure boolean.
*/
int readXYZHV(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth, int unScaledAvgInputs);

/**
* Read a containing ARC ASCII Raster coordinate system. 
* File Format:
*	ncols = number of columns
*	nrows = number of rows
*	xllcorner = lower left X coordinate
*	yllcorner = lower left Y coordinate
*	cellsize = size of the data cell 
*	nodata = flag used to signify no data
*	dataPoint[1][1]			dataPoint[1][2]			dataPoint[1][3]		...		dataPoint[1][ncols]
*	dataPoint[2][1]			dataPoint[2][2]			dataPoint[2][3]		...		dataPoint[2][ncols]
*	dataPoint[3][1]			dataPoint[3][2]			dataPoint[3][3]		...		dataPoint[3][ncols]
*	...						...						...					...		...
*	dataPoint[nrows][1]		dataPoint[nrows][2]		dataPoint[nrows][3] ...		dataPoint[nrows][ncols]
*
* @param fileName - File name of the XYZ_HV file to be read.
* @param bbox - The bounding extents of the input data if it is to be cut.
* @param inputData - Vector TRUE_DATA that contains the values read from each file.  Only the index specified by pos is modified.  The TRUE_DATA at inputData[pos] will contain all of data read from the inputFile.  (Returned).
* @param pos - Integer representing the location of inputData where the data from the input file should be stored.
* @return Success or failure boolean.
*/
int readARC_ASCII_RASTER(string &fileName, BOUNDING_BOX bbox, vector<TRUE_DATA> *inputData, int &pos, int noerr, int nonegdepth);

#endif
