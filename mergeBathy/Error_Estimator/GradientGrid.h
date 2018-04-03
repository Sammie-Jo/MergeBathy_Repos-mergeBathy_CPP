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
#pragma once
#include "../grid.h"

class GradientGrid
{
	private:
		
		dgrid dzdxOut;
		dgrid dzdyOut;
		dgrid magOut;
		dgrid slopeOut;
		dgrid aspectOut;
		bool existFlag;

	public:
		GradientGrid(){ existFlag = false; };
		~GradientGrid(){ clear();}
		void calc_GradientGrid(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn);
		void clear();
		bool exist(){ return existFlag; }
		vector<double> getDzdxOut_Vector(){ return dzdxOut.vec(); }
		vector<double> getDzdyOut_Vector(){ return dzdyOut.vec(); }
		vector<double> getMagOut_Vector(){ return magOut.vec(); }
		vector<double>* getSlopeOut_Vector(){ return &(slopeOut.vec()); }
		vector<double> getAspectOut_Vector(){ return aspectOut.vec(); }
		dgrid getDzdxOut_Grid(){ return dzdxOut; }
		dgrid getDzdyOut_Grid(){ return dzdyOut; }
		dgrid getMagOut_Grid(){ return magOut; }
		dgrid* getSlopeOut(){ return &slopeOut; }
		dgrid getAspectOut(){ return aspectOut; }
};

