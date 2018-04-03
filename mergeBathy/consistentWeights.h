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

