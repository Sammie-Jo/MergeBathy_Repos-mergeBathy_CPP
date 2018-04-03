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
* @file			fileWriter.h
* @brief		Header file for writing text and ARC ASCII Raster file formats for mergeBathy.
* @author		Kevin Duvieilh
* @date			04 August 2011
*
*/

#pragma once
//#include <iostream>
//#include <fstream>
#include <string>
#include <vector>
#include <map>
#include "constants.h"
//#include "grid.h"

using namespace std;

/**
* Write a raster output file.
* @param fileName - File name of the list file that contains each individual data set.
* @param x - The x coordinates.
* @param y - The y coordinates.
* @param z - The depth values.
* @param spacingX - Grid spacing in the X dimension.
* @param spacingY - Grid spacing in the X dimension.
* @param xSize - Grid size in the X dimension.
* @param ySize - Grid size in the Y dimension.
* @param additionalOptions - A map of the additional output file options.
* @param ZoneRef - The UTM reference zone.
* @return Success or failure boolean.
*/


int writeBagFile(string &fName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, int gridSpacingX, int gridSpacingY, int xSize, int ySize, map<string, int> *additionalOptions, char ZoneRef[4], vector<double> *z0, vector<double> *e0, vector<double> *zK, vector<double> *eK);

//int baginit(string XMLFile, string fName, vector<double> *z, vector<double> *e, vector<double> *x, vector<double> *y);   

int readBagFile(string &fName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, int gridSpacingX, int gridSpacingY, int xSize, int ySize, map<string, int> *additionalOptions);

//int ProcessCommandInput( int argc, char **argv, char *gisFile, char *xmlFile, char *bagFile, int *summaryOnly );
