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
* @file			subSampleData.h
* @brief		Produce a subsampled data grid using a boxcar window
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/
#pragma once
#include <iostream>
#include <vector>
#include "grid.h"

/**
* This subsamples and spaces the data points by interpolating data into regular sample bins using boxcar window.  Will remove duplicate data points and sort data for faster interpolation.
* @param inputDataX - Vector of input X data points.
* @param inputDataY - Vector of input Y data points.
* @param inputDataZ - Vector of input Depth data points.
* @param inputDataE - Vector of input Error data points.
* @param inputDataH - Vector of input Horizontal Error data points.
* @param inputDataV - Vector of input Vertical Error data points.
* @param gridSpacingX - Spacing of the computational area in the X direction.
* @param gridSpacingY - Spacing of the computational area in the Y direction.
* @param inX0 - Minimum value contained in the inputDataX.
* @param inY0 - Minimum value contained in the inputDataY.
* @param meanX - Mean value of inputDataX.
* @param meanX - Mean value of inputDataY.
* @param dispIntermResults - Determine if intermediate output should be displayed to the command line.
* @param subsampledData - A n by 5 vector containing the subsampled data.  Index 0 contains the X data normalized to the mean of X. Index 1 contains the Y data normalized to the mean of Y. Index 2 contains the Depth data. Index 3 contains the Error data.  Index 4 contains the Error data squared. (Returned).
* @return Success or failure value.
*/
int subsampleData(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, double gridSpacingX, double gridSpacingY, double &inX0, double &inY0, double meanX, double meanY, bool dispIntermResults, vector< vector<double> > *subsampledData);

/**
* Round a double value up or down.
* >= .5 rounds up.
* < .5 rounds down.
*
* @param x - Double value to be rounded.
* @return Rounded value.
*/
double round(double x);

/**
* Partition to be used on a dgrid, called from QuickSort.  Column 1 of the grid contains the data values, column 2 contains the original index location of the data point in column 1 in the same row.
* Used to partition the grid for another quick sort.
*
* @param d - dgrid containing points to be sorted in column 1 and original indexes of those coordinate in column 2.  When sorting the entire row is sorted d[1][1] and d[1][2] are moved together to maintain knowledge of where the data point was in the original grid. (Returned).
* @param left - Left index to quick sort around.
* @param right - Right index to quick sort around.
*/
int Partition( dgrid *d, int left, int right);

/**
* Selection sort to be used on a dgrid, called from QuickSort.  Column 1 of the grid contains the data values, column 2 contains the original index location of the data point in column 1 in the same row.
* Only used when the number of elements is too small to continue partitioning.
*
* @param d - dgrid containing points to be sorted in column 1 and original indexes of those coordinate in column 2.  When sorting the entire row is sorted d[1][1] and d[1][2] are moved together to maintain knowledge of where the data point was in the original grid. (Returned).
* @param left - Left index to quick sort around.
* @param right - Right index to quick sort around.
*/
void SelectionSort(dgrid *data, int left, int right);

/**
* Quicksort to be used on a dgrid.  Column 1 of the grid contains the data values, column 2 contains the original index location of the data point in column 1 in the same row.
*
* @param d - dgrid containing points to be sorted in column 1 and original indexes of those coordinate in column 2.  When sorting the entire row is sorted d[1][1] and d[1][2] are moved together to maintain knowledge of where the data point was in the original grid. (Returned).
* @param left - Left index to quick sort around.
* @param right - Right index to quick sort around.
*/
void Quicksort( dgrid *d, int left, int right);

/**
* merge_helper to be used on a dgrid.  Column 1 of the grid contains the data values, column 2 contains the original index location of the data point in column 1 in the same row.
*
* @param input - dgrid containing points to be sorted in column 1 and original indexes of those coordinate in column 2.  When sorting the entire row is sorted d[1][1] and d[1][2] are moved together to maintain knowledge of where the data point was in the original grid. (Returned).
* @param left - Left index to merge sort around.
* @param right - Right index to merge sort around.
* @param scratch - dgrid for sorting.
*/
void merge_helper(dgrid *input, int left, int right, dgrid *scratch);

/**
* MergeSort to be used on a dgrid.  Column 1 of the grid contains the data values, column 2 contains the original index location of the data point in column 1 in the same row.
*
* @param input - dgrid containing points to be sorted in column 1 and original indexes of those coordinate in column 2.  When sorting the entire row is sorted d[1][1] and d[1][2] are moved together to maintain knowledge of where the data point was in the original grid. (Returned).
* @param size - Left index to quick sort around.
*/
int MergeSort(dgrid *input, int size);

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ORIGINAL FUNCTIONS FOR MONTE CARLO RUNS!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

/**
* This subsamples and spaces the data points by interpolating data into regular sample bins using boxcar window.  Will remove duplicate data points and sort data for faster interpolation.
* @param inputDataX - Vector of input X data points.
* @param inputDataY - Vector of input Y data points.
* @param inputDataZ - Vector of input Depth data points.
* @param inputDataE - Vector of input Error data points.
* @param gridSpacingX - Spacing of the computational area in the X direction.
* @param gridSpacingY - Spacing of the computational area in the Y direction.
* @param inX0 - Minimum value contained in the inputDataX.
* @param inY0 - Minimum value contained in the inputDataY.
* @param meanX - Mean value of inputDataX.
* @param meanX - Mean value of inputDataY.
* @param dispIntermResults - Determine if intermediate output should be displayed to the command line.
* @param subsampledData - A n by 5 vector containing the subsampled data.  Index 0 contains the X data normalized to the mean of X. Index 1 contains the Y data normalized to the mean of Y. Index 2 contains the Depth data. Index 3 contains the Error data.  Index 4 contains the Error data squared. (Returned).
* @return Success or failure value.
*/
//int subsampleData_ORIGINAL(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, double gridSpacingX, double gridSpacingY, double &inX0, double &inY0, double meanX, double meanY, bool dispIntermResults, vector< vector<double> > *subsampledData);
