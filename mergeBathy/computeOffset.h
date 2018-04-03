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
* @file			computeOffset.h
* @brief		Header file for computing the offset of various input data sets in mergeBathy.
* @author		Kevin Duvieilh
* @date			08 July 2011
*
*/

#pragma once
//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	#include "WarningStates.h"		//Disable all Warnings!!!
#endif

//Third-Party Includes
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning( pop )			//Restore warning state
#endif

//Project Includes
#include <string>
#include <vector>
#include "grid.h"
#include <map>
#include "inFileStructs.h"
#include "constants.h"
#include <cmath>
#include <algorithm>

/**
* Computes the offset between input data sets.  Each data set is smoothed using mergeBathy and then meshed back together.
* @param inputData - Vector TRUE_DATA that contains the values read from each file.  The length of the vector corresponds to the number of input files.  Each index of the vector corresponds the TRUE_DATA read from a specific input file.
* @param x0 - Minimum value contained in the x vector.
* @param y0 - Minimum value contained in the y vector.
* @param x1 - Maximum value contained in the x vector.
* @param y1 - Maximum value contained in the y vector.
* @param xInterpVector - Vector of the X points to be interpolated.
* @param yInterpVector - Vector of the Y points to be interpolated.
* @param gridSpacingX - Computational grid spacing in the X direction defined in meters.
* @param gridSpacingY - Computational grid spacing in the Y direction defined in meters.
* @param zDataIn - The vector of of original Depth values. The computed offset is applied to this vector. (Returned).
* @param eDataIn - The vector of of original Error values. The computed offset is applied to this vector. (Returned).
* @return Success or failure value.
*/
void computeOffset(vector<TRUE_DATA> *inputData, double x0, double y0, double x1, double y1, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, vector<double> *zDataIn, vector<double> *eDataIn, vector<double> *hEDataIn, vector<double> *vEDataIn, string &kernel,  map<string, int> addOpts);

template<class T>
bool InvertMatrix (const boost::numeric::ublas::compressed_matrix<T>& input, boost::numeric::ublas::compressed_matrix<T>& inverse);
template <typename size_type, typename A> // SJZ
int determinant(const boost::numeric::ublas::permutation_matrix<size_type,A>& pm);

