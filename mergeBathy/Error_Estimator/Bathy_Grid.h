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
#include <list>
#include <vector>
#include "geom.h"
#include "pointList.h"
#include "sHullDelaunay.h"
#include "GradientGrid.h"

enum GridType{ MBZ, GMT, ALGSpline, Sibson_Natural_Neighbor_Interp };

class InterpGrid
{
	private:
		// uncertainty estimation vectors
		std::vector<double> e; 
		std::vector<double> e2; 
		std::vector<double> e3; 
		std::vector<double> e4;
		std::vector<double> e5;
	
		// depth vector
		std::vector<double> z;
		// slope vector
		std::vector<double> s;
		//raster
		std::vector<double> z0;
		std::vector<double> e0;
		std::vector<double> zK;
		std::vector<double> eK;
		std::vector<double> nei; 
		std::vector<double> rei; 
	

		// query list
		PointList* ptlSurface;
		// ptlSurface's gradients
		GradientGrid* grads;
		// ptlSurface's point locations in Tin
		std::vector<Triangle > triangles;
		// grid type: MBZ, GMT
		GridType gtype;
		// constants for computing uncertainties
		double scaleFactor;
		double alpha;
		double deltaMin;
		
	public:	
		// getters
		std::vector<double> getE()const { return e; }
		std::vector<double> getE2()const { return e2; }
		std::vector<double> getE3()const { return e3; }
		std::vector<double> getE4()const { return e4; }
		std::vector<double> getE5()const { return e5; }

		std::vector<double>* getZ() { return &z; }
		//raster/bag
		std::vector<double> getZ0() { return z0; }
		std::vector<double> getE0() { return e0; }
		std::vector<double> getZK() { return zK; }
		std::vector<double> getEK() { return eK; }
		std::vector<double> getNEI() { return nei; }
		std::vector<double> getREI() { return rei; }



		////Setters for Gridding Ensembler
		//void setPtl (vector<double> *xSurf, vector<double> *ySurf, vector<double> *zSurf, double xOffset, double yOffset) { ptlSurface = new PointList(); ptlSurface->setFromVectors(*xSurf, *ySurf, *zSurf, xOffset, yOffset); }
		//void setE (vector<double> eTemp) { e.assign(eTemp.begin(), eTemp.end()); }
		//void setE2(vector<double> eTemp) { e2.assign(eTemp.begin(), eTemp.end()); }
		//void setE3(vector<double> eTemp) { e3.assign(eTemp.begin(), eTemp.end()); }
		//void setE4(vector<double> eTemp) { e4.assign(eTemp.begin(), eTemp.end()); }
		//void setE5(vector<double> eTemp) { e5.assign(eTemp.begin(), eTemp.end()); }

		GradientGrid* getGrads() const { return grads; }

		// constructors
		InterpGrid(GridType g) { gtype = g; ptlSurface = NULL; grads = NULL; }
		// destructors
		~InterpGrid() { clear(); }
		void clear();

		// Calculate uncertainties function
		void estimate(std::vector<double> *xSurf, std::vector<double> *ySurf, std::vector<double> *zSurf, const double& sH, const double& alpha, const double& deltaMin, SHullDelaunay* tin, string depthInterpMethod, string errorInterpMethod, string extrapMethod);

		//void InterpGrid::estimate(vector<double> *xSurf, vector<double> *ySurf, vector<double> *zSurf,vector<double> *eSurf,vector<double> *hSurf,vector<double> *vSurf, vector<double> *nmseiSurf,vector<double> *reiSurf,const double& sH, const double& alpha, const double& deltaMin, SHullDelaunay* tin, string interpMethod);
		friend class Bathy_Grid;
};




class Bathy_Grid
{
	private:
		// input data
		PointList* ptl;
		// Delaunay triangulation created from ptl
		SHullDelaunay* tin;
		// interpGrids queried against tin
		std::vector<InterpGrid*> interpGrids;
		// ensemble estimation vectors
		std::vector<double> ensemble_X;
		std::vector<double> ensemble_Y;
		std::vector<double> ensemble_Z;
		std::vector<double> ensemble_U;
		std::vector<double> ensemble_U2;
		std::vector<double> ensemble_U3;
		std::vector<double> ensemble_U4;
		std::vector<double> ensemble_U5;
		std::vector<double> ensemble_H;
		std::vector<double> ensemble_V;

	public: 
		// flag to perform gridding computations. 
		int GriddingFlag;
		// getters
		PointList* getPointList() const{ return ptl; }
		SHullDelaunay* getTin() const { return tin; }
		InterpGrid* getGrid(GridType g);
		std::vector<InterpGrid*> getGrids();
		// get ensemble estimation vectors
		std::vector<double> getEnsemble_X () const{ return ensemble_X; }
		std::vector<double> getEnsemble_Y () const{ return ensemble_Y; }
		std::vector<double> getEnsemble_Z () const{ return ensemble_Z; }
		std::vector<double> getEnsemble_U () const{ return ensemble_U; }
		std::vector<double> getEnsemble_U2() const{ return ensemble_U2; }
		std::vector<double> getEnsemble_U3() const{ return ensemble_U3; }
		std::vector<double> getEnsemble_U4() const{ return ensemble_U4; }
		std::vector<double> getEnsemble_U5() const{ return ensemble_U5; }
		std::vector<double> getEnsemble_H() const{ return ensemble_H; }
		std::vector<double> getEnsemble_V() const{ return ensemble_V; }
		
		void Construct_TinRaster(vector<double> *xInit, vector<double> *yInit, vector<double> *zInit, vector<double> *eInit);
		void Construct_TinRaster(vector<double> *xInit, vector<double> *yInit, vector<double> *zInit, vector<double> *eInit, vector<double> *neiInit, vector<double> *reiInit, vector<double> *z0Init, vector<double> *e0Init, vector<double> *zKInit, vector<double> *eKInit);

		// construct Delaunay tin
		void Construct_Tin(std::vector<double> *xInit, std::vector<double> *yInit, std::vector<double> *zInit, std::vector<double> *hInit, std::vector<double> *vInit);
		
		// perform ensemble of interpGrids
		bool ensemble();
		
		// add interpGrid to list of interpGrids
		void addToList( InterpGrid* g);
		
		// compare two interpGrids
		void operator=(const Bathy_Grid& g);
		
		// print output to file
		void printensemble(std::string z_OutputFileName);

		// constructor
		Bathy_Grid() {ptl=NULL; tin=NULL;}
		
		// destructors
		~Bathy_Grid() { clear(); }
		void clear();

};






