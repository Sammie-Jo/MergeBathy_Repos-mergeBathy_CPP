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

