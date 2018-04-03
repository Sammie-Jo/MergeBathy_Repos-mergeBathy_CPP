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
#include "standardOperations.h"

double findMin(double &a, double &b)
{
	return (((a) < (b)) ? (a) : (b));
}

double findMax(double &a, double &b)
{
	return (((a) > (b)) ? (a) : (b));
}

int findMin(int &a, int &b)
{
	return (((a) < (b)) ? (a) : (b));
}

int findMax(int &a, int &b)
{
	return (((a) > (b)) ? (a) : (b));
}

double roundDouble (double a)
{
	return (((ceil(a)-a) < .5) ? (ceil(a)) : (floor(a)));
}

double min (const double& x, const double& y) { if(x < y) { return x; } return y; }

double max (const double& x, const double& y) { if(x > y) { return x; } return y; }
/**
* Compute the mean from a vector.
* @param toMean - Vector of data points.
* @return The mean of the vector toMean.
*/
double mean(dvector *toMean, bool useZeros)
{
	//Initialize local variables
	double meanOut = 0.00;
	size_t size = (*toMean).size();
	int i;

	//Compute the mean
	for (i = 0; i < (const int)size; i++)
	{
		if(!useZeros && (*toMean)[i] == 0)
			size--;
		meanOut += (*toMean)[i];
	}
	meanOut = meanOut / (double)size;

	return meanOut;
}

/**
* Compute the standard deviation from a vector.
* @param toStd - Vector of data points.
* @return The standard deviation of the vector toStd.
*/
double standardDeviation(dvector *toStd, bool useZeros)
{
	//Initialize local variables
	double meanStd = 0.00;
	double stdOut = 0.00;
	size_t size = (*toStd).size();
	int i;

	//Compute the mean
	meanStd = mean(toStd, useZeros);

	//Compute standard deviation
	for (i = 0; i < (const int)size; i++)
	{
		if(!useZeros && (*toStd)[i] == 0)
			size--;
		else
			stdOut += pow( ( (*toStd)[i] - meanStd ), 2);
	}
	stdOut = stdOut / (double)(size - 1.00);
	stdOut = sqrt(stdOut);

	return stdOut;
}

/**
* Compute the median from a vector.
* Comparable to matlab, if 2 medians return the average
* @param median - Vector of data points.
* @return The median of the vector bar.
*/
double median(vector<double> bar)
{
	vector<double>::iterator temp;

	sort(bar.begin(),bar.end());
	size_t n = bar.size();

	double m1, m2;
	if(n%2==1)
		m1 = bar[(n+1)/2-1];
	else{
		m1 = bar[n/2-1];
		m2 = bar[n/2+1-1];
		m1 += abs(m1-m2)/2;
	}
	return m1;
}

double fix(double x)//trunc
{
    return (x>0) ? floor(x) : ceil(x);
}

void createMeshXYDims(double x0, double y0, double x1, double y1, double gridSpacingX, double gridSpacingY, vector<double> *xt, vector<double> *yt)
{
	double meanXt = 0.00;
	double meanYt = 0.00;
	double locationValue = x0;
	int xtSize = 0;
	while (locationValue <= (const double)(x1 )) {
//	while (locationValue <= (const double)(x1 + gridSpacingX)) {
		(*xt).push_back(locationValue);
		meanXt += locationValue;
		locationValue += gridSpacingX;
		xtSize += 1;
	}
	//Get the mean for later use
	meanXt = meanXt / (double)xt->size();

	//B. Calculate Y
	locationValue = y0;
	int ytSize = 0;

	while (locationValue <= (const double)(y1 )) {
//	while (locationValue <= (const double)(y1 + gridSpacingY)) {
		yt->push_back(locationValue);
		meanYt += locationValue;
		locationValue += gridSpacingY;
		ytSize += 1;
	}
	//Get the mean for later use
	meanYt = meanYt / (double)yt->size();
}

void createMeshGrid(vector<double> *xt, vector<double> *yt, vector<double> *xMeshVector, vector<double> *yMeshVector, dgrid *xMeshGrid, dgrid *yMeshGrid)
{
	/*int temp = xtSize*ytSize;
	xMeshVector->resize(temp);
	yMeshVector->resize(temp);
	xMeshGrid->resize(ytSize, xtSize);
	yMeshGrid->resize(ytSize, xtSize);*/

	//if(ytSize*xtSize != xMeshVector->size())
	//	cerr<<"Error: Dimensions of parameters do not match calculated dimensions!"<<endl;

	size_t xtSize = xt->size();
	size_t ytSize = yt->size();
	//C. Reform the data to a grid that we can use for calculations
	int clx = 0;
	int cly = 0;
	int currentLoc = 0;
	for (int i = 0; i < (const int)ytSize; i++)
	{
		currentLoc = i;
		for (int j = 0; j < (const int)xtSize; j++)
		{
			(*xMeshGrid)(i,j) = (*xt)[clx];
			(*yMeshGrid)(i,j) = (*yt)[cly];
			(*xMeshVector)[currentLoc] = (*xt)[clx];
			(*yMeshVector)[currentLoc] = (*yt)[cly];
			currentLoc = currentLoc + (int) ytSize;
			clx += 1;
		}
		clx = 0;
		cly += 1;
	}
}

void reshapeGrid(dgrid *g, dgrid *newg)
{
	int l = 0;
	int k = 0;
	int g_rows = g->rows();
	int g_cols = g->cols();
	int newg_rows = newg->rows();
	int newg_cols = newg->cols();

	int s1 = g_rows * g_cols;
	int s2 = newg_rows * newg_cols;
	if(s1!=s2)
		std::cerr << "Number of elements must not change!" << std::endl;
	else{
		for (int i = 0; i < g_rows; i++)
		{
			for (int j = 0; j < g_cols; j++)
			{
				(*newg)(k,l) = (*g)(i,j);
				l++;
				if(l == newg_cols){ l = 0; k++;}
			}
			//k++;
			//if(k > newg_rows){ k = 0; l++;}
		}
	}
}


