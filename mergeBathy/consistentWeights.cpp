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
#include "consistentWeights.h"
#include "subSampleData.h"
/*
--Notes taken from consistenWeight.m--
Over-Relaxation (David Young) - Determines how much of a frequency to take out. This is determined by 0.

Compute weights consistent with observations and expected error variance for each observation.
The weights need to be consistent (consistently increase) in order to be smooth and continuous.  As the errors approach 0, so should the weightings. (Paul said this was Bayesian approach)

Inputs:
z, the observations
e2, the error VARIANCE (SQUARED DAMMIT!) at each observation
wtol, (optional) convergence criteria (suggest 0.1 to 0.01)

Output:
w, the consistent weights
s2, the estimated true variance

Assume (wj^2) = s2/(s2+e2j), where s2 is true variance and e2j is error variance at observation j.
Solve problem by guessing wj, then compute s2, then update wj.
Assumes that true variance, s2, is constant over the data.
--Notes takes from consistentWeight.m--
*/

void consistentWeights(const vector<double> *z, const vector<double> *e, double *wtol, vector<double> *w, double *s3Value)
{
//************************************************************************************
// 0. Declare local variables and objects
//************************************************************************************
	int i;
	int cnt = 0;
	double tempMax, tempMax2;
	double mu, s2, s3;
	double sumWinit, sumZ;
	int zSize = (const int)(*z).size();
	vector<double> winit (zSize,1);

//************************************************************************************
//I. Initialize variables, with weights all set to 1 as an initial guess.
//************************************************************************************
	tempMax = abs((*w)[0]-winit[0]);
	mu = 0;

//************************************************************************************
// II. Iterate to find weights. Exit criteria: max(abs(w-winit)) <= wtol
// or 10 loops, whichever comes first.
//************************************************************************************
	while((cnt < 10) && (tempMax > (*wtol))){
		cnt += 1;

		//A. Reassign "winit" after first loop s.t. large weights do not
		//dominate too soon
		if (cnt > 1){
			for (i = 0; i < zSize; i++){
				winit[i] = ((*w)[i] + winit[i]) / 2.00;
			}
		}

		//B. Compute terms of the iteration
		//1. Weighted mean
		sumWinit = 0;
		sumZ = 0;
		for (i = 0; i < zSize; i++){
			sumWinit = sumWinit + winit[i];
			sumZ += winit[i] * (*z)[i];
		}

		mu = sumZ/sumWinit;

		//2a. Weighted mean squared residual-> true variance if weights are accurate
		sumZ = 0;
		sumWinit = 0;
		for (i = 0; i < zSize; i++){
			sumZ = (*z)[i] - mu;
			sumWinit += pow(winit[i] * sumZ,2);
		}
		s2 = sumWinit / double(zSize);

		//2b. Use this instead when nu is small, the a priori error
		sumZ = 0;
		sumWinit = 0;
		for (i = 0; i < zSize; i++){
			sumZ += ((zSize-1.00)*s2 + (*e)[i])/double(zSize);
		}
		s3 = sumZ / double(zSize);
		(*s3Value) = s3;

		//3. Weights...
		sumZ = 0;
		sumWinit = 0;

		//Handle division by 0.
		//The weighting should not be far off from the error.
		//As the error approaches 0, so should the weighting
		//in order to have a smooth continuity at 0.
		double eTemp;
		if((*e)[0] == 0)
			eTemp = DEFAULT_E;
		else eTemp = (*e)[0];

		tempMax = std::abs(sqrt( s3 / (s3 + eTemp))-winit[0]);
		tempMax2 = 0;
		for (i = 0; i < zSize; i++){
			if((*e)[i] == 0)
				eTemp = DEFAULT_E;
			else eTemp = (*e)[i];

			(*w)[i] = sqrt(s3 / (s3 + eTemp));
			tempMax2 = std::abs((*w)[i]-winit[i]);
			if (tempMax2 > tempMax){
				tempMax = tempMax2;
			}
		}
	}
	winit.clear();
}

void consistentWeights3(const vector<double> *z, const vector<double> *e, double *wtol, vector<double> *w, double *weightsMax)
{
//************************************************************************************
// 0. Declare local variables and objects
//************************************************************************************
	int i;
	int cnt = 0;
	double tempMax, tempMax2;
	double mu, s2, s3;
	double sumWinit, sumZ;
//	double returnValue;
	int zSize = (const int)(*z).size();
	vector<double> winit (zSize,1);
	(*weightsMax) = -9999999.0;

//************************************************************************************
//I. Initialize variables, with weights all set to 1 as an initial guess.
//************************************************************************************
	tempMax = abs((*w)[0]-winit[0]);
	mu = 0;

//************************************************************************************
// II. Iterate to find weights. Exit criteria: max(abs(w-winit)) <= wtol
// or 10 loops, whichever comes first.
//************************************************************************************
	while((cnt < 10) && (tempMax > (*wtol))){
		cnt += 1;

		//A. Reassign "winit" after first loop s.t. large weights do not
		//dominate too soon
		if (cnt > 1){
			for (i = 0; i < zSize; i++){
				winit[i] = ((*w)[i] + winit[i]) / 2.00;
			}
		}

		//B. Compute terms of the iteration
		//1. Weighted mean
		sumWinit = 0;
		sumZ = 0;
		for (i = 0; i < zSize; i++){
			sumWinit = sumWinit + winit[i];
			sumZ += (winit[i] * (*z)[i]);
		}

		mu = sumZ/sumWinit;

		//2a. Weighted mean squared residual-> true variance if weights are accurate
		sumZ = 0;
		sumWinit = 0;
		for (i = 0; i < zSize; i++){
			sumZ = ((*z)[i] - mu);
			sumWinit += pow(winit[i] * sumZ,2);
		}
		s2 = sumWinit / double(zSize);

		//2b. Use this instead when nu is small, the a priori error
		sumZ = 0;
		sumWinit = 0;
		for (i = 0; i < zSize; i++){
			sumZ += ((zSize-1.00)*s2 + (*e)[i])/double(zSize);
		}
		s3 = sumZ / double(zSize);

		//3. Weights...
		sumZ = 0;
		sumWinit = 0;
		tempMax = abs(sqrt(s3 / (s3 + (*e)[0]))-winit[0]);
		tempMax2 = 0;
		for (i = 0; i < zSize; i++){
			(*w)[i] = sqrt(s3 / (s3 + (*e)[i]));
			tempMax2 = abs((*w)[i]-winit[i]);
			if ((*w)[i] > (*weightsMax))
			{
				(*weightsMax) = (*w)[i];
			}
			if (tempMax2 > tempMax){
				tempMax = tempMax2;
			}
		}
	}
	winit.clear();
}