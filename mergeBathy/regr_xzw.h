/** 
* @file			regr_xzw.h
* @purpose		Gneral linear regression call
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/

#ifndef REGRXZW_H_
#define REGRXZW_H_

#include <vector>
#include "grid.h"

using namespace std;

/**
* General linear regression call.
* @param regrX - Vector of dependant variables.
* @param subDataX - Vector of subsampled X data points.
* @param subDataY - Vector of subsampled Y data points.
* @param subDataZ - Vector of subsampled Z data points.
* @param subDataE - Vector of subsampled E data points.
* @param weights - Vector of weight values for each depth data point.
* @param btrend - dVector of estimated parameters: z^ = X*b. (Returned).
* @param bi - dVector stimated variances (root mean square error) of parameter(s) (confidence intervals assume gaussian white-noise dist, with bmse estimated variance). (Returned).
*/
void regr_xzw(vector<double> *regrX, dvector *subDataX, dvector *subDataY, dvector *subDataZ, dvector *subDataE, vector<double> *weights, dvector *btrend, dvector *bi);

#endif
