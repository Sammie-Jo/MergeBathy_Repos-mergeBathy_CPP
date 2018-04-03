/**********************************************************************
* CC0 License
**********************************************************************
* MergeBathy - Tool to combine one or more bathymetric data files onto a single input grid.
* written in 2015 by Samantha J.Zambo(samantha.zambo@gmail.com) while employed by the U.S.Naval Research Laboratory.
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
#include "regr_xzw.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>
#include "constants.h"
#include "ALG/linalg.h"
#include "ALG/ap.h"

//extern "C" {
//	#include "f2c.h"
//    #include "blaswrap.h"
//    #include "clapack.h"
//}
//extern int dgecon_(const char *norm, const int *n, double *a, const int *lda, const double *anorm, double *rcond, double *work, int *iwork, int *info, int len);
//extern int dgetrf_(const int *m, const int *n, double *a, const int *lda, int *lpiv, int *info);
//extern double dlange_(const char *norm, const int *m, const int *n, const double *a, const int *lda, double *work, const int norm_len);

void regr_xzw(std::vector<double> *regrX, dvector *subDataX, dvector *subDataY, dvector *subDataZ, dvector *subDataE, vector<double> *weights, dvector *btrend, dvector *bi)
{
//************************************************************************************
//0. Declare and initialize local variables and objects.
//************************************************************************************
	double  n;
	int m = 3;
	int i;

//************************************************************************************
//I. Compute the trend surface.
//************************************************************************************
	//A. Compute the degrees of freedom
	double sumW = 0;
	double maxW = (*weights)[0];
	int subDataXSize = (const int)(*subDataX).size();
	int weightsSize = (const int)(*weights).size();

	dgrid xNew(subDataXSize, 3, 0);

	dvector zNew(subDataXSize,0);
	dvector x1(subDataXSize,0);
	dvector x2(subDataXSize,0);
	dvector x3(subDataXSize,0);

	//B. Loop and apply the weights
	for (i = 0; i < weightsSize; i++){
		sumW = sumW + (*weights)[i];
		if ((*weights)[i] > maxW){
			maxW = (*weights)[i];
		}
		zNew[i] = (*subDataZ)[i] * (*weights)[i];
		xNew(i,0) = (*regrX)[i] * (*weights)[i];
		xNew(i,1) = (*subDataX)[i] * (*weights)[i];
		xNew(i,2) = (*subDataY)[i] * (*weights)[i];

		x1[i] = (*regrX)[i] * (*weights)[i];
		x2[i] = (*subDataX)[i] * (*weights)[i];
		x3[i] = (*subDataY)[i] * (*weights)[i];
	}
	// number of dof
	n = sumW / maxW;

	//C. Compute covariances
	dgrid XX((matmult(trans(xNew), xNew)/n));
	dgrid XZ((matmult(trans(zNew), xNew)/n)); // XZ(((trans(zNew) * xNew)/n));

	//D. Solve the weighted least squares problem
	//dgrid XXinv;//(inv(XX, 0, 'U', true, 'N'));

	//Compute an estimate for the reciprocal condition of XX in 1-norm.
	//If XX is well conditioned, the reciprocal condition number is near 1.0.
	//If XX is badly conditioned and nearly singular, the reciprocal condition number is near 0.
	//and the system is sensitive to perturbations
	alglib::ae_int_t N = XX.rows();
	alglib::real_2d_array XX_Temp = alglib::real_2d_array();
	double* d = &(XX.vec())[0];
	XX_Temp.setcontent(XX.rows(), XX.cols(), d);

	//Debug Printout XX
	//for(int ii = 0; ii < XX.rows(); ii++) {
	//	for(int j = 0; j < XX.cols(); j++) cout << XX_Temp[ii][j] << "  ";
	//	cout << endl;
	//}

	double rcondValue = alglib::rmatrixrcond1(XX_Temp, N);
	if(rcondValue <= RCOND_TOL)
	{
	/*	cerr << endl;
		cerr << "Reciprocal Condition Number = " << rcondValue << " < " << RCOND_TOL << endl;*/
		return;
	}

	dgrid XXinv(inv(XX, 0, 'U', true, 'N'));

	//E. compute parameters
	dgrid b(matmult(XXinv, trans(XZ)));//(XXinv * trans(XZ)));
	for(i = 0; i < (const int)b.rows(); i++){
		(*btrend)[i] = b(i,0);
	}

	//F. model residuals
	dgrid mszG( (matmult(trans(zNew),zNew) / n));
	double msz = mszG(0,0);

	dgrid msrG( mszG - (matmult(matmult(trans(b),XX),b)));
	double msr = msrG(0,0);

	double sk( 1 - (msr/msz));
	XX.clear();
	XZ.clear();
	b.clear();
	mszG.clear();
	msrG.clear();
	xNew.clear();
	zNew.clear();

	//G. Compute variance estimates
	for (i = 0; i < (const int) XXinv.rows(); i++){
		(*bi)[i] = sqrt((XXinv(i,i) * msr) / (n - m));
	}
	XXinv.clear();
}

/*
Notes from mergeBathy developers:
 To get the right error estimate, use this estimator.
 Read priestly, page 368: mse/msr ~ chi-sq(m)/(N-m)
 Here is argument:

 var_observed = var_model + var_residual = var_true + var_noise

 and

 var_model = var_true + var_artificial

    var_true is variance of perfect model

    var_noise is additive white noise

    var_artificial is due to random correlations, i.e., aliased

 a linear model will pick up this much of noise

 var_a = (m/N) var_n, our extra bit of information

 Thus,

 var_n = var_r + var_a = var_r + (m/N) var_n = var_r /(1-m/N)
        = var_r*N/(N-m), expected value

		and

 var_a = var_r (m/N)/(1-m/N) = (m/(N-m)) var_r

 thus: mse_of_model = ( msr*m/(sumw-m) );

 This is plausible error in model, assumed uniform over range of data
 BUT, we are interested in error of estimate of b(1), the value at x=0

 See priestly page on regression models, which states that
 b are normally dist. around b_true, with

 var(b) = var_noise*diag(XX_inv)

 msn = msr*sumw/(sumw-m);

 mse_of_b(1) = XX_inv(1)*msn/N;

 -- Look to the F-dist to explain ratio of variances
 after testing synthetic examples, conclude that this is robust
 even leaving large mean values in!
*/