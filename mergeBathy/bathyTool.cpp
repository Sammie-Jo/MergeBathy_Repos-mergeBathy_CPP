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
#include "bathyTool.h"
#include "subSampleData.h"
#include <fstream>

//************************************************************************************
// SUBROUTINE I: Function call for running standard BathyTool
//************************************************************************************
int bathyTool(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, dgrid *xMeshGrid, dgrid *yMeshGrid, vector<double> *xSingleVector, vector<double> *ySingleVector, double gridSpacingX, double gridSpacingY, double &x0, double &y0, double meanXSingle, double meanYSingle, string &kernelName, map<string, int> additionalOptions, double subDataMulitplier, bool dispIntermResults, bool useDscale, const double neitol, OUTPUT_DATA *xyzOut)
{
	//************************************************************************************
	//0. Declare and initialize local variables and objects for use with bathyTool.
	//************************************************************************************
	int returnValue = 0;
	//Apply smoothing multiplier.  Necessary for compute offset
	double subSpacingX = gridSpacingX * subDataMulitplier;
	double subSpacingY = gridSpacingY * subDataMulitplier;
	if (useDscale)
	{
		subSpacingX = subSpacingX / Dscale;
		subSpacingY = subSpacingY / Dscale;
	}

	vector< vector <double> > subsampledData = vector< vector<double> >(7);

	//************************************************************************************
	//I. Subsample the data for use in interpolating .
	//************************************************************************************
	returnValue = subsampleData(inputDataX, inputDataY, inputDataZ, inputDataE, inputDataHErr, inputDataVErr, subSpacingX, subSpacingY, x0, y0, meanXSingle, meanYSingle, dispIntermResults, &subsampledData);
	if (returnValue != 0)
	{
		return returnValue;
	}

	//************************************************************************************
	//II. Use subsampled data in regular grid interpolation subroutine.
	//************************************************************************************
	returnValue = scalecInterpTile(&subsampledData, xMeshGrid, yMeshGrid, xSingleVector, ySingleVector, gridSpacingX, gridSpacingY, meanXSingle, meanYSingle, kernelName, additionalOptions, dispIntermResults, neitol, xyzOut);
	if (returnValue != 0)
	{
		return returnValue;
	}

	//************************************************************************************
	//III. Clean up data structures and return.
	//************************************************************************************
	subsampledData[0].clear();
	subsampledData[1].clear();
	subsampledData[2].clear();
	subsampledData[3].clear();
	subsampledData[4].clear();
	subsampledData[5].clear();
	subsampledData[6].clear();

	return returnValue;
}

//************************************************************************************
// SUBROUTINE II: Function call for running Pre-Interpolated BathyTool
//************************************************************************************
int bathyToolPreDefined(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, double &x0, double &y0, double meanXSingle, double meanYSingle, string &kernelName,  map<string, int> additionalOptions, const double neitol, OUTPUT_DATA *xyzOut)
{
	//************************************************************************************
	//0. Declare and initialize local variables and objects for use with bathyTool.
	//************************************************************************************
	int returnValue = 0;
	vector< vector <double> > subsampledData = vector< vector<double> >(7);

	//************************************************************************************
	//I. Subsample the data for use in interpolating .
	//************************************************************************************
	returnValue = subsampleData(inputDataX, inputDataY, inputDataZ, inputDataE, inputDataHErr, inputDataVErr, gridSpacingX/Dscale, gridSpacingY/Dscale, x0, y0, meanXSingle, meanYSingle, true, &subsampledData);
	if (returnValue != 0)
	{
		return returnValue;
	}

	//************************************************************************************
	//II. Use subsampled data in regular grid interpolation subroutine.
	//************************************************************************************
	returnValue = scalecInterp(&subsampledData, xInterpVector, yInterpVector, gridSpacingX, gridSpacingY, meanXSingle, meanYSingle, kernelName, additionalOptions, neitol, xyzOut);
	if (returnValue != 0)
	{
		return returnValue;
	}

	//************************************************************************************
	//III. Clean up data structures and return.
	//************************************************************************************
	subsampledData[0].clear();
	subsampledData[1].clear();
	subsampledData[2].clear();
	subsampledData[3].clear();
	subsampledData[4].clear();
	subsampledData[5].clear();
	subsampledData[6].clear();

	return returnValue;
}

