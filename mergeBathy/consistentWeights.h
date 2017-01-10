/**
* @file			consistentWeights.h
* @brief		Compute weights consistent with measurements and expected errors variance for each measurement.
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/

#pragma once
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include "constants.h"
#include "grid.h"

/**
* Compute weights consistent with measurements and expected errors variance for each measurement.
* Assume (wj^2) = s2/(s2 + e2j), where s2 is true variance and e2j is error
* variance at observation j. Solve problem by guessing wj, then compute s2,
* then update wj. Assumes that true variance, s2, is constant over the data
* Algorithm quits after 10 loops or when max(abs(w-winit)) < wtol.
* Putting in e = 0 and wtol = 0.1 will result in w = 1.
* @param z - Vector of depth values.
* @param e - Vector of the square of the error values.
* @param wtol - Convergence criteria (NP suggests 0.1 to 0.01).
* @param w - Vector of weights to be applied to each depth.  Index locations of the weights correspond to the index locations of the depths. (Returned).
* @param s3Value - The estimated true variance of the observations. (Returned).
*/
void consistentWeights(const vector<double> *z, const vector<double> *e, double *wtol, vector<double> *w, double *s3Value);

/**
* Compute weights consistent with measurements and expected errors variance for each measurement.
* Assume (wj^2) = s2/(s2 + e2j), where s2 is true variance and e2j is error
* variance at observation j. Solve problem by guessing wj, then compute s2,
* then update wj. Assumes that true variance, s2, is constant over the data
* Algorithm quits after 10 loops or when max(abs(w-winit)) < wtol.
* Putting in e = 0 and wtol = 0.1 will result in w = 1.
* @param z - Vector of depth values.
* @param e - Vector of the square of the error values.
* @param wtol - Convergence criteria (NP suggests 0.1 to 0.01).
* @param w - Vector of weights to be applied to each depth.  Index locations of the weights correspond to the index locations of the depths. (Returned).
* @param weightsMax - The maximum value located in the weights (w) vector. (Returned).
*/
void consistentWeights3(const vector<double> *z, const vector<double> *e, double *wtol, vector<double> *w, double *weightsMax);

