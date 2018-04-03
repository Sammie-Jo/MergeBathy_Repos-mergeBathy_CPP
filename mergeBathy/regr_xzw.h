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
