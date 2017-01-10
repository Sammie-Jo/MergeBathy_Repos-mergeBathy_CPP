#include "scalecInterp.h"
#include "regr_xzw.h"
#include "standardOperations.h"
#include <fstream>
#include <time.h>
#include <math.h>

//************************************************************************************
// SUBROUTINE I: Function call for doing 2D linear interpolation.
//************************************************************************************
double interp1(double x, const vector<double> *xi, const vector<double> *yi){
	double y;
	unsigned int j;

	//A. if x is outside the xi[] interval take a boundary value (left or right)
	if (x <= (*xi)[0])
		return y = (*yi)[0];
	if (x >= (*xi)[(*xi).size()-1])
		return y = (*yi)[(*yi).size()-1];

	//B. loop to find j so that x[j-1] < x < x[j]
	j = 0;
	while (j <= ((*xi).size() - 1))
	{
		if ((*xi)[j] >= x)
			break;
		j = j + 1;
	}
	y = (*yi)[j-1] + ((*yi)[j] - (*yi)[j-1])*(x - (*xi)[j-1])/((*xi)[j]-(*xi)[j-1]);
	return y;
}

//************************************************************************************
// SUBROUTINE II: Function call for doing loess weighting.
// This function was just transplanted from the old version of
// mergeBathy It needs to be streamlined and vectorized to
// match the rest of the Optimized code.
//************************************************************************************
void loessKernel(int dtemp, int p, vector<double> *ri, vector<double> *ai){
/*
 Input
   d, the dimension of the kernel
   p, the power of the loess (1=linear, 2=quadratic)

 Output
   r, the radial distance
   ar, the weights, sum to one, are scaled by number of inputs
*/

//************************************************************************************
//0. Declare and initialize variables
//************************************************************************************
		double DX = 0.05;
		double rItem = -(1 + DX);

		vector<double> x = vector<double>((uint)(((1+DX)+(1+DX))/DX) + 1, 0);
		dgrid X = dgrid((uint)(x.size()*x.size()), 2, 0);
		dgrid W;
		dgrid w;
		dgrid r;
		dgrid tmpGrid1, tmpGrid2;
		dgrid X2tmp;
		dgrid r2tmp;
		unsigned int i, j, k;
		ofstream outfile;
		uint d = (uint) dtemp;

//************************************************************************************
//I. Initialize computational grid for building the loess windows.
//************************************************************************************
		int cnt = 0;
		for (i = 0; i < x.size(); i++){
			x[i] = rItem;
			rItem = rItem + DX;
		}
		for (i = 0; i < x.size(); i++){
			for (j = 0; j < x.size(); j++){
				X(cnt, 0) =  x[j];
				X(cnt, 1) =  x[i];
				cnt += 1;
			}
		}

//************************************************************************************
//II. Get weights of the tri-cube function (Cleveland, 1979).
//************************************************************************************
		//A. Compute Euclidean tri-cube function
		tmpGrid1 = dgrid(X.rows(), 1, 1);
		for (i = 0; i < X.rows(); i++){
			tmpGrid1(i, 0) = pow(X(i, 0), 2) + pow(X(i, 1), 2);
		}

		W = dgrid( tmpGrid1.rows(), 1, 0);
		for (i = 0; i < tmpGrid1.rows(); i++){
			if (tmpGrid1.get(i,0) < 1){
				W(i, 0) = pow((1- pow(tmpGrid1(i,0), 3.00)), 3);
			}else{
				W(i, 0) = 0;
			}
		}

		r = dgrid((uint)x.size(), d, 0);
		for (i = 0; i < x.size(); i++){
			r(i,0) = x[i];
		}

		//B. Compute radial tri-cube function
		tmpGrid2 = dgrid(r.rows(), 1, 0);
		for (i = 0; i < r.rows(); i++){
			tmpGrid2(i, 0) = pow(r(i, 0), 2) + pow(r(i, 1), 2);
		}
		w = dgrid( tmpGrid2.rows(), 1, 0);

		for (i = 0; i < tmpGrid2.rows(); i++){
			if (tmpGrid2.get(i,0) < 1){
				w(i, 0) = pow((1- pow(tmpGrid2(i,0), 3)), 3);
			}else{
				w(i, 0) = 0;
			}
		}

//************************************************************************************
//III. Get weights of the tri-cube function (Cleveland, 1979).
//************************************************************************************
		uint m = d+1;
		double n = pow((double)x.size(), (double)d);
		uint q = (uint)(0.5 * m * (m+1));

		//A. Build GRID Objects for holding weights
		//1. Linear loess window
		if (p == 1){
			X2tmp = dgrid(X.rows(), 3, 0);
			r2tmp = dgrid(r.rows(), 3, 0);
			for (i = 0; i < X.rows(); i++){
				X2tmp(i,0) = 1;
				X2tmp(i,1) = X(i,0);
				X2tmp(i,2) = X(i,1);
			}
			for (i = 0; i < r.rows(); i++){
				r2tmp(i,0) = 1;
				r2tmp(i,1) = r(i,0);
				r2tmp(i,2) = r(i,1);
			}

		//2. Quadratic loess window
		}else if (p == 2){
			X2tmp = dgrid(X.rows(),1 + X.cols() + (q-m), 0);
			r2tmp = dgrid(r.rows(), 1 + r.cols() + (q-m), 0);

			for (i = 0; i < X.rows(); i++){
				X2tmp(i,0) = 1;
				X2tmp(i,1) = X(i,0);
				X2tmp(i,2) = X(i,1);
			}
			for (i = 0; i < r.rows(); i++){
				r2tmp(i,0) = 1;
				r2tmp(i,1) = r(i,0);
				r2tmp(i,2) = r(i,1);
			}

			for (i = 1; i < (d+1); i++){
				for (j = 1; j < i+1; j++){
					for (k = 0; k < X2tmp.rows(); k++){
						X2tmp(k, m) = X2tmp(k,i) * X2tmp(k,j);
					}
					for (k = 0; k < r2tmp.rows(); k++){
						r2tmp(k, m) = r2tmp(k,i) * r2tmp(k,j);
					}
					m  = m+1;
				}
			}
		}else{
			printf("FATAL ERROR:\nENDING DATA RUN\n");
			exit(1);
		}

		//B. Apply tri-00 weights
		//1. Euclidean case
		dgrid XW(W.rows(), X2tmp.cols(), 0);
		for (i = 0; i < W.rows(); i++){
			for (j = 0; j < X2tmp.cols(); j++){
				XW(i,j) =  X2tmp(i,j)*W(i, 0);
			}
		}
		//2. Polar case
		dgrid rw(w.rows(),r2tmp.cols(), 0);
		for (i = 0; i < w.rows(); i++){
			for (j = 0; j < r2tmp.cols(); j++){
				rw(i,j) = r2tmp(i,j)*w(i,0);
			}
		}

//************************************************************************************
//IV. Calculate the Loess windows via linear regression.
//************************************************************************************
		//1. Euclidean case
		dgrid XX( matmult(trans(XW),XW));
		XX = XX / ( pow((double)x.size(), (double)d));

		dgrid XX_inv = dgrid(inv(XX));
		dgrid XX_inv2 = dgrid( XX_inv.rows(), 1, 0);
		for (i = 0; i < XX_inv.rows(); i++){
			XX_inv2(i,0) = XX_inv.get(i,0);
		}
		w = w / (double)x.size();

		tmpGrid1.clear();

		//2. Polar case
		tmpGrid1 = dgrid(matmult(rw,XX_inv2));

		dgrid ar = dgrid(tmpGrid1.rows(), 1, 0);
		for (i = 0; i < tmpGrid1.rows(); i++){
			ar(i,0) = tmpGrid1(i,0)*w(i,0);
		}

//************************************************************************************
//V. Assign loess windows weights to output variables.
//************************************************************************************

		for (i = 0; i < x.size(); i++){
			(*ai)[i] = ar(i,0);
			(*ri)[i] = x[i];
		}
}

//************************************************************************************
// SUBROUTINE III: Preliminary Function call used for weighting data and applying window weights.
//************************************************************************************
void scalecInterpPerturbations_PreCompute(const vector<double> *subDataZ, const vector<double> *subDataE, vector<double> *perturbWeights, PERTURBS *perturb)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int i;
	double wtol = 0.1;
	double weightsMax;

	//************************************************************************************
	// I. Compute weighting scales across the entire interpolation window
	//************************************************************************************
	consistentWeights3(subDataZ, subDataE, &wtol, perturbWeights, &weightsMax);

	for (i = 0; i < (const int)(*subDataZ).size(); i++){
		(*perturbWeights)[i] = ((*perturbWeights)[i] + eps)/(eps + weightsMax);
	}
	//************************************************************************************
	// III. Generate ri (radial distances) and ai (weights) based
	// on the specific smoothing scale defined by the user
	//************************************************************************************
	if ((*perturb).kernelName== "hanning" || (*perturb).kernelName == "hann"){
		double rItem = 0;
		(*perturb).riVector = vector<double>((int)((RMAX + DR) / DR) + 1, 0);
		(*perturb).aiVector = vector<double>((int)((RMAX + DR) / DR) + 1, 0);

		for (i = 0; i < (const int)((RMAX + DR) / DR) + 1; i++)
		{
			(*perturb).riVector[i] = rItem;
			rItem = rItem + DR;

			//A. This is the hanning weight routine from MATLAB w = hanning_wt(x)
			// Input: x are nxm inputs, weights will be centered on x=0
			// Output: w are weights 0<=w<=1

			//Convert to radial distance
			(*perturb).aiVector[i] = abs((*perturb).riVector[i]);

			//hanning window
			if ((*perturb).aiVector[i] <= 1){
				(*perturb).aiVector[i] = (1 - pow( cos(PI * (0.5 + (0.5 * (*perturb).aiVector[i]))), 2));
			}else{
				(*perturb).aiVector[i] = 0;
			}
		}
	//B. Boxcar window
	}else if ((*perturb).kernelName == "boxcar"){
		double rItem = 0;
		(*perturb).riVector = vector<double>((int)((RMAX + DR) / DR) + 1, 0);
		(*perturb).aiVector = vector<double>((int)((RMAX + DR) / DR) + 1, 1);

		for (i = 0; i < (const int)((RMAX + DR) / DR) + 1; i++){
			(*perturb).riVector[i] = rItem;
			rItem = rItem + DR;
		}
	//C. Quadratic loess window
	}else if ((*perturb).kernelName == "quadloess"){
		(*perturb).riVector = vector<double>((int)((RMAX + DR) / DR) + 1, 1);
		(*perturb).aiVector = vector<double>((int)((RMAX + DR) / DR) + 1, 1);

		loessKernel(2, 2, &(*perturb).riVector, &(*perturb).aiVector);

	//D. Linear loess window
	}else if ((*perturb).kernelName == "loess"){//linloess
		(*perturb).riVector = vector<double>((int)((RMAX + DR) / DR) + 1, 1);
		(*perturb).aiVector = vector<double>((int)((RMAX + DR) / DR) + 1, 1);

		loessKernel(2, 1, &(*perturb).riVector, &(*perturb).aiVector);

	//E. Error message
	}else{
		printf("Error: Unrecognized smoothing window.\n");
		printf("Aborting data run.\n");
		exit(1);
	}
}

//************************************************************************************
// SUBROUTINE IV: Primary computation routine.
//************************************************************************************
void scalecInterpPerturbations_Compute(const vector<double> *subDataX, const vector<double> *subDataY, const vector<double> *subDataZ, const vector<double> *subDataE, const vector<double> *subDataH, const vector<double> *subDataV, double *xGridValue, double *yGridValue, const vector<double> *weights, const double neitol, double dmin, double slope, const vector<double> *subDataX0, const vector<double> *subDataY0, dgrid *Xiii, PERTURBS *perturb, bool MSE, bool PROP_UNCERT, bool KALMAN, bool KRIGING)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	int count;
	const int idySize = (const int)(*subDataX).size();
	vector<int> aid; 
	vector<double> nWeights;	//a 
	vector<double> nWeights0;	//a0 
	vector<double> r = vector<double>(idySize, 0.00);
	vector<double> rx0;
	vector<double> newrx0;

	double p;
	double q;
	int na;
	int naTmp = 0;
	double sumCount = 0;
	double sumNormWeights_init;	//suma
	double sumNormWeights;
	double valueCalculated1, valueCalculated2;
	double rCompute, rCompute0;
	double matMultValue, matMultValue2;
	vector<double> grow;

	//Safety Check vU and hU
	if((*subDataV).empty())
		cerr << "scalecInterpPerturbations Error: vU is empty.";
	else if((*subDataV).size() != (*subDataZ).size())
		cerr << "scalecInterpPerturbations Error: vU is wrong size.";

	if((*subDataH).empty())
		cerr << "scalecInterpPerturbations Error: hU is empty.";
	else if((*subDataH).size() != (*subDataZ).size())
		cerr << "scalecInterpPerturbations Error: hU is wrong size.";

	//************************************************************************************
	//I. Begin interpolation loop over the Ni interpolation points. Establish the values to be used in r
	//************************************************************************************
	for (int j = 0; j < idySize; j++)
	{
		sumCount = 0;
		valueCalculated1 = ((*subDataX)[j]) - ((*xGridValue));
		valueCalculated2 = ((*subDataY)[j]) - ((*yGridValue));
		sumCount += pow(valueCalculated1,2);
		sumCount += pow(valueCalculated2,2);
		r[j] = sqrt(sumCount);
	}

	count = 0;
	naTmp = 0;
	na = 0;
	aid = vector<int>();
	nWeights = vector<double>();
	nWeights0 = vector<double>();

	//************************************************************************************
	//II. THIS PART IS IMPORTANT! We must expand smoothing scales if error
	//	tolerances are not met until either we met error tolerances or have
	//	expanded the smoothing scale 10 times.
	//	If the maximum Normalized Mean Squared Error of Interpolation
	//	tolerance (nmseitol) is not passed in or is set to [],
	//	we had set nmseitol = inf so that tolerance is never invoked.
	//  ri (radial distances) and ai (weights) based
	//  on the specific smoothing scale defined by the user
	//************************************************************************************
	while((((*perturb).perturbationNEi > neitol) && (count < 10)) || (count == 0))
	{
		aid.clear();
		nWeights.clear();
		nWeights0.clear();
		if(KRIGING)
			count = count + 1;
		p = pow(2.00,count);
		if(!KRIGING)
			count = count + 1;
		sumNormWeights = 0.0;
		sumNormWeights_init = 0;
		matMultValue = 0.0;

		//A. Find the weights inside the linear regression
		for (int j = 0; j < idySize; j++){
			if (r[j] < p){
				rCompute0 = interp1((r[j] / p), &(*perturb).riVector, &(*perturb).aiVector);
				rCompute = rCompute0 * (*weights)[j]; 
				aid.push_back(j); 
				nWeights.push_back(rCompute); 
				nWeights0.push_back(rCompute0); 
				sumNormWeights_init = sumNormWeights_init + rCompute;
			}
		}
		na = (const int)aid.size(); 

		//cout << nWeights << endl;

		if (na == 0)
			continue;

		//B. Normalize if we have some weights
		if (abs(sumNormWeights_init) > 0) 
		{
			for (int j = 0; j < (const int)na; j++){
				nWeights[j] /= sumNormWeights_init;
				sumNormWeights += nWeights[j];
				matMultValue = matMultValue + nWeights[j] * nWeights[j]; //nmsei
			}
		}else
		{
			//weights not normalized
			for (int j = 0; j < (const int)na; j++){
				sumNormWeights += nWeights[j];
				matMultValue = matMultValue + nWeights[j] * nWeights[j];
			}
		}

		//C. Check weights for error calculations. NP note: actual error is
		//	the noised passed (nmsei*s) + signal reduced (1-nmsei)*s. In the
		//	first case, we want to distinguish points with at least one nonzero
		//	weight

		if(sumNormWeights > 0){
			double item = matMultValue;
			item = item * (1.00 - eps);
			(*perturb).perturbationNEi = item; 
		}else{
			(*perturb).perturbationNEi = 1; 
		}
	}

	double S_H = 1.0;

	//Define a alternate dmin to take care of cases of irregularly
	//space output points. We will use this dmin instead of
	//from a grid. dmin is defined to be as follows. PAE - 07 May 2014
	for (int j = 0; j < na; j++)
	{
		sumCount = 0;
		valueCalculated1 = (*subDataX0)[aid[j]] - (*Xiii)(0,0); // *xGridValue;
		valueCalculated2 = (*subDataY0)[aid[j]] - (*Xiii)(0,1); // *yGridValue;
		sumCount += pow(valueCalculated1,2);
		sumCount += pow(valueCalculated2,2);
		rx0.push_back(sqrt(sumCount));
		if(rx0[j] != 0){
			newrx0.push_back(rx0[j]);
		}
	}

	if(dmin == NaN)
		if(newrx0.empty())
			cerr<<"newrx0 is empty!!"<<endl;
		else dmin = median(newrx0);

	grow = vector<double>((*subDataE).size(),1);
	double tempVal = 0.0;
	double matMultVal3 = 0.0;

	if(PROP_UNCERT || MSE)
	{
		//************************************************************************************
		//III. Finally!!! Convolve window against the data to get interpolated
		//	depths.  This is a matrix multiplication operation, so we need C-code
		//for doing that.
		//************************************************************************************
		matMultValue = 0;
		for (int j = 0; j < (const int)na; j++)
		{
			matMultValue += ((*subDataZ)[aid[j]]  *  nWeights[j]);
		}
		(*perturb).perturbationZ = matMultValue; 

		//Linear estimators hann, loess, quadloess, boxcar
		if(PROP_UNCERT)
		{
			for(int j = 0; j < na; j++)
			{
				matMultVal3 += nWeights0[j] * (*subDataZ)[ aid[j] ]; 
				// Variance growth
				grow[j] = 1 + pow(((rx0[j] + S_H * (*subDataH)[aid[j] ])/dmin),2);
				tempVal += pow(nWeights0[j], 2) * ( grow[j] * pow((*subDataV)[ aid[j] ],2 ) + pow( (*subDataH)[aid[j]] * tan((PI/180) * slope), 2)); //SJZ 9/20/16
			}
			(*perturb).perturbationZ0 = matMultVal3; 
			(*perturb).perturbationE0 = tempVal;//sqrt(tempVal); SJZ 9/20/16

			if(0>((*perturb).perturbationE0))
			{
				cerr << "Nan";
			}
		}

		if(MSE)
		{
			q = pow(p, -2); 
			(*perturb).perturbationNEi = 1.00 - q * (1.00- (*perturb).perturbationNEi);

			//printf("PD: %.6f\t%.6f\n", (*perturbationZ), (*perturbationNEi));

			//A. Calculate residuals (difference between smoothed surface and observations) if possible
			sumCount = 0;
			if (  ( (*perturb).perturbationNEi < 1 ) && ( na > 0 )  )
			{
				//B. Calculate weighted residuals
				matMultValue = 0;
				matMultValue2 = 0;

				for (int j = 0; j < (const int)na; j++)
				{
					sumCount = ((*subDataZ)[aid[j]] -  (*perturb).perturbationZ)  *  nWeights[j]; 
					matMultValue = matMultValue + sumCount*sumCount; 
					matMultValue2 = matMultValue2 + (nWeights[j] * (*subDataE)[ aid[j] ] );
				}

				//C. Calculate weighted mean square residual, msri
				(*perturb).perturbationREi = matMultValue / (*perturb).perturbationNEi;

				//D. Recalculate if small degrees of freedom, add estimate that converges to observation error for (na-1)->0
				(*perturb).perturbationREi = (((na-1) * ((*perturb).perturbationREi)) + matMultValue2 ) / ((double)na); 

				//E. Calculate weighted mean square error of the interpolation, msei
				(*perturb).perturbationE = (*perturb).perturbationREi  *   ((*perturb).perturbationNEi) / double(1.00 - (*perturb).perturbationNEi);

				// Insert correct calculations here
				(*perturb).standardDev2 = standardDeviation((dvector*)(subDataZ));

				//printf("PD: %.6f\t%.6f\n", (*perturbationREi), (*perturbationE));
			}
			else
			{
				//F. Set nmsei to one if
				//(((*perturbationData).nmsei.get(i,0) < 1) && (na > 0)) is false
				(*perturb).perturbationNEi = 1; //nmsei
			}
		}
	}//end propUncert/mse

	int tKALMAN = 1;
	if(aid.empty())
	{
		tKALMAN = 0;
		(*perturb).perturbationZK = NaN;
		(*perturb).perturbationEK = NaN;
	}

	//Kalman recursive estimator
	//Doesn't matter which window used because they will all converge
	//to the same value since weights are not used
	//and the same points are used since they fall in the same radius.
	if(tKALMAN)
	{
		// Start by using the first values to initialize
		double d_meas;
		double K;
		double var_est_new;
		double var_meas;
		double z_est_new;
		double innov;

		double grows = 1 + pow(( ( rx0[0] + S_H * (*subDataH)[ aid[0] ]) /dmin), 2);
		double z_est_past = (*subDataZ)[ aid[0] ];
		double var_est_past = pow((*subDataV)[aid[0]], 2) * grows + pow( (*subDataH)[aid[0]] * tan((PI/180) * slope), 2);

		if(na == 1)
		{
			(*perturb).perturbationZK = z_est_past;
			(*perturb).perturbationEK = var_est_past;
		}
		else
		{
			// Then we loop through the remaining data
			for(int j = 1; j < na; j++)
			{
				grows = 1+ pow(( ( rx0[j] + S_H * (*subDataH)[ aid[j] ]) /dmin),2);
				// get new measurement values
				d_meas = (*subDataZ)[ aid[j] ];

				// Increase variance with distance
				var_meas = pow((*subDataV)[aid[j]], 2) * grows + pow( (*subDataH)[aid[j]]*tan((PI/180) * slope), 2);

				// Compute the Kalman gain
				if(var_est_past + var_meas == 0)
					K = 0;
				else
					K = var_est_past/(var_est_past + var_meas);

				// Compute depth innovations
				innov = d_meas - z_est_past;
				// Compute updated depth estimate
				z_est_new = z_est_past + K*innov;
				// Compute updated variance estimate
				var_est_new = K*var_meas;
				// Substitute for next loop
				z_est_past = z_est_new;
				var_est_past = var_est_new;
			} 

			// Save the results
			(*perturb).perturbationZK = z_est_new;
			(*perturb).perturbationEK = var_est_new;
		}
	}//end KALMAN
	aid.clear();
	nWeights.clear();
	nWeights0.clear();
	grow.clear();
	rx0.clear();
	newrx0.clear();

	r.clear();
}


//From _Serial and is the same as mergeBathy_v3.7.1_Paul which is modified to pass standardDev in the Perturbations call
void scalecInterpPerturbations_Compute_ForKriging(const vector<double> *subDataX, const vector<double> *subDataY, const vector<double> *subDataZ, const vector<double> *subDataE, double *xGridValue, double *yGridValue, const vector<double> *weights, const vector<double> *ri, const vector<double> *ai, const double neitol, double *perturbationZ, double *perturbationE, double *perturbationNEi, double *perturbationREi, double *standardDev)
{
	//************************************************************************************
	// 0. Declare local variables and objects
	//************************************************************************************
	double p;
	int j, count;
	const size_t idySize = (*subDataX).size();

	vector<double> r = vector<double>(idySize, 0.00);

	vector<int> aid;
	vector<double> nWeights;

	int naTmp = 0;
	int na = 0;
	double q;
	double sumCount = 0;
	double sumNormWeights_init, sumNormWeights;
	double valueCalculated1, valueCalculated2;
	double rCompute;
	double matMultValue, matMultValue2;

	//************************************************************************************
	//I. Begin interpolation loop over the Ni interpolation points. Establish the values to be use in r
	//************************************************************************************

	for (j = 0; j < (const int)idySize; j++)
	{
		sumCount = 0;
		valueCalculated1 = ((*subDataX)[j]) - ((*xGridValue));
		valueCalculated2 = ((*subDataY)[j]) - ((*yGridValue));
		sumCount += pow(valueCalculated1,2);
		sumCount += pow(valueCalculated2,2);
		r[j] = sqrt(sumCount);
	}

	count = 0;
	naTmp = 0;
	na = 0;

	//************************************************************************************
	//II. THIS PART IS IMPORTANT! We must expand smoothing scales if error
	//	tolerances are not met until either we met error tolerances or have
	//	expanded the smoothing scale 10 times.
	//	If the maximum Normalized Mean Squared Error of Interpolation
	//	tolerance (nmseitol) is not passed in or is set to [],
	//	we had set nmseitol = inf so that tolerance is never invoked.
	//************************************************************************************
	while((((*perturbationNEi) > neitol) && (count < 10)) || (count == 0))
	{
		aid.clear();
		nWeights.clear();
		//Increment count before pow() since its for kriging
		count = count + 1;
		p = pow(2.00,count);
		sumNormWeights = 0;
		sumNormWeights_init = 0;
		matMultValue = 0;

		//A. Find the weights inside the linear regression
		for (int j = 0; j < (const int)idySize; j++){
			if (r[j] < p){
				rCompute = interp1((r[j] / p), ri, ai) * (*weights)[j];
				aid.push_back(j);
				nWeights.push_back(rCompute);
				sumNormWeights_init = sumNormWeights_init + rCompute;
			}
		}
		na = (const int)aid.size();

		if (na == 0)
		{
			continue;
		}

		//B. Normalize if we have some weights
		if (abs(sumNormWeights_init) > 0)
		{
			for (int j = 0; j < (const int)na; j++){
				nWeights[j] /= sumNormWeights_init;
				sumNormWeights += nWeights[j];
				matMultValue = matMultValue + nWeights[j]*nWeights[j];
			}
		}else
		{
			//weights not normalized
			for (int j = 0; j < (const int)na; j++){
				sumNormWeights += nWeights[j];
				matMultValue = matMultValue + nWeights[j]*nWeights[j];
			}
		}

		//C. Check weights for error calculations. NP note: actual error is
		//	the noised passed (nmsei*s) + signal reduced (1-nmsei)*s. In the
		//	first case, we want to distinguish points with at least one nonzero
		//	weight

		if(sumNormWeights > 0){
			double item = matMultValue;
			item = item * (1.00 - eps);
			(*perturbationNEi) = item;
		}else{
			(*perturbationNEi) = 1;
		}
	}

	q = pow(p, -2);
	(*perturbationNEi) = 1.00 - q * (1.00-(*perturbationNEi));

	//************************************************************************************
	//III. Finally!!! Convolve window against the data to get interpolated
	//	depths.
	//************************************************************************************
	matMultValue = 0;
	for (j = 0; j < (const int)na; j++)
	{
		matMultValue += ((*subDataZ)[aid[j]] * nWeights[j]);
	}
	(*perturbationZ) = matMultValue;

	//A. Calculate residuals (difference between smoothed surface and observations) if possible
	sumCount = 0;
	if (((*perturbationNEi) < 1) && (na > 0))
	{
		//B. Calculate weighted residuals
		matMultValue = 0;
		matMultValue2 = 0;

		for (j = 0; j < (const int)na; j++){
			sumCount = ((*subDataZ)[aid[j]] - (*perturbationZ)) * nWeights[j];
			matMultValue = matMultValue + sumCount*sumCount;
			matMultValue2 = matMultValue2 + ((nWeights[j]) * ((*subDataE)[aid[j]]));
		}

		//C. Calculate weighted mean square residual, msri
		(*perturbationREi) = (matMultValue / (*perturbationNEi));


		//D. Recalculate if small degrees of freedom, add estimate that converges to observation error for (na-1)->0
		(*perturbationREi) = ( ( (na - 1) * (*perturbationREi)) + matMultValue2 ) / (double)na;

		//E. Calculate weighted mean square error of the interpolation, msei
		(*perturbationE) = ((*perturbationREi) * (*perturbationNEi)) / double(1.00 - (*perturbationNEi));

		// Insert correct calculations here
		(*standardDev) = standardDeviation((dvector*)(subDataZ));

	}else
	{
		//F. Set nmsei to one if
		//(((*perturbationData).nmsei.get(i,0) < 1) && (na > 0)) is false
		(*perturbationNEi) = 1;
	}

	r.clear();
	aid.clear();
	nWeights.clear();
}

/*This function is the original scalecInterpPerturbations_Compute and only differs from Compute_ForKriging in 1 line where count is incremented before pow().  This function does not contain any of the new functionality added or bug fixes.  It is only here in case the developer wishes to implement the original Monte Carlo processing path.*/
//void scalecInterpPerturbations_Compute(const vector<double> *subDataX, const vector<double> *subDataY, const vector<double> *subDataZ, const vector<double> *subDataE, double *xGridValue, double *yGridValue, const vector<double> *weights, const vector<double> *ri, const vector<double> *ai, const double neitol, double *perturbationZ, double *perturbationE, double *perturbationNEi, double *perturbationREi, double *standardDev)
//{
//	//************************************************************************************
//	// 0. Declare local variables and objects 
//	//************************************************************************************
//	double p;
//	int j, count;
//	const int idySize = (*subDataX).size();
//
//	vector<double> r = vector<double>(idySize, 0.00);
//
//	vector<int> aid;
//	vector<double> nWeights;
//	
//	int naTmp = 0;
//	int na = 0;
//	double q;
//	double sumCount = 0;
//	double sumNormWeights_init, sumNormWeights;
//	double valueCalculated1, valueCalculated2;
//	double rCompute;
//	double matMultValue, matMultValue2;
//	
//	//************************************************************************************
//	//I. Begin interpolation loop over the Ni interpolation points. Establish the values to be use in r
//	//************************************************************************************
//
//	for (j = 0; j < idySize; j++)
//	{
//		sumCount = 0;
//		valueCalculated1 = ((*subDataX)[j]) - ((*xGridValue));
//		valueCalculated2 = ((*subDataY)[j]) - ((*yGridValue));
//		sumCount += pow(valueCalculated1,2);
//		sumCount += pow(valueCalculated2,2);
//		r[j] = sqrt(sumCount);
//	}
//
//	//cout << r << endl;
//
//	count = 0;
//	naTmp = 0;
//	na = 0;
//
//	//************************************************************************************
//	//II. THIS PART IS IMPORTANT! We must expand smoothing scales if error 
//	//	tolerances are not met until either we met error tolerances or have
//	//	expanded the smoothing scale 10 times. 
//	//	If the maximum Normalized Mean Squared Error of Interpolation 
//	//	tolerance (nmseitol) is not passed in or is set to [], 
//	//	we had set nmseitol = inf so that tolerance is never invoked.
//	//************************************************************************************
//	while((((*perturbationNEi) > neitol) && (count < 10)) || (count == 0))
//	{
//		aid.clear();
//		nWeights.clear();
//		p = pow(2.00,count);
//		count = count + 1;
//		sumNormWeights = 0;
//		sumNormWeights_init = 0;
//		matMultValue = 0;
//
//		//A. Find the weights inside the linear regression
//		for (j = 0; j < idySize; j++){
//			if (r[j] < p){
//				rCompute = interp1((r[j] / p), ri, ai) * (*weights)[j];
//				aid.push_back(j);
//				nWeights.push_back(rCompute);
//				sumNormWeights_init = sumNormWeights_init + rCompute;
//			}
//		}
//		na = aid.size();
//
//		//cout << nWeights << endl;
//		
//		if (na == 0)
//		{
//			continue;
//		}
//
//		//B. Normalize if we have some weights
//		if (abs(sumNormWeights_init) > 0)
//		{
//			for (j = 0; j < (const int)na; j++){
//				nWeights[j] /= sumNormWeights_init;
//				sumNormWeights += nWeights[j];
//				matMultValue = matMultValue + nWeights[j]*nWeights[j];
//			}
//		}else
//		{   
//			//weights not normalized
//			for (j = 0; j < (const int)na; j++){
//				sumNormWeights += nWeights[j];
//				matMultValue = matMultValue + nWeights[j]*nWeights[j];
//			}
//		}
//
//		//C. Check weights for error calculations. NP note: actual error is
//		//	the noised passed (nmsei*s) + signal reduced (1-nmsei)*s. In the
//		//	first case, we want to distinguish points with at least one nonzero
//		//	weight
//
//		if(sumNormWeights > 0){
//			double item = matMultValue;
//			item = item * (1.00 - eps);
//			(*perturbationNEi) = item; 
//		}else{ 
//			(*perturbationNEi) = 1;
//		}
//	}
//
//	q = pow(p, -2);
//	(*perturbationNEi) = 1.00 - q * (1.00-(*perturbationNEi));
//
//	//************************************************************************************
//	//III. Finally!!! Convolve window against the data to get interpolated 
//	//	depths.
//	//************************************************************************************
//	matMultValue = 0;
//	for (j = 0; j < (const int)na; j++)
//	{
//		matMultValue += ((*subDataZ)[aid[j]] * nWeights[j]);
//	}
//	(*perturbationZ) = matMultValue;
//
//	//printf("PD: %.6f\t%.6f\n", (*perturbationZ), (*perturbationNEi));
//	
//	//A. Calculate residuals (difference between smoothed surface and observations) if possible
//	sumCount = 0;
//	if (((*perturbationNEi) < 1) && (na > 0))
//	{
//		//B. Calculate weighted residuals
//		matMultValue = 0;
//		matMultValue2 = 0;
//
//		for (j = 0; j < (const int)na; j++){
//			sumCount = ((*subDataZ)[aid[j]] - (*perturbationZ)) * nWeights[j];
//			matMultValue = matMultValue + sumCount*sumCount;
//			matMultValue2 = matMultValue2 + ((nWeights[j]) * ((*subDataE)[aid[j]]));
//		}
//
//		//C. Calculate weighted mean square residual, msri
//		(*perturbationREi) = (matMultValue / (*perturbationNEi));
//
//		//D. Recalculate if small degrees of freedom, add estimate that converges to observation error for (na-1)->0
//		(*perturbationREi) = ( ( (na - 1) * (*perturbationREi)) + matMultValue2 ) / (double)na;
//
//		//E. Calculate weighted mean square error of the interpolation, msei
//		(*perturbationE) = ((*perturbationREi) * (*perturbationNEi)) / double(1.00 - (*perturbationNEi));
//
//		// Insert correct calculations here
//		(*standardDev) = standardDeviation((dvector*)(subDataZ));
//
//		//printf("PD: %.6f\t%.6f\n", (*perturbationREi), (*perturbationE));
//
//	}else
//	{
//		//F. Set nmsei to one if 
//		//(((*perturbationData).nmsei.get(i,0) < 1) && (na > 0)) is false
//		(*perturbationNEi) = 1;
//	}
//
//	r.clear();
//	aid.clear();
//	nWeights.clear();
//}
//
//
