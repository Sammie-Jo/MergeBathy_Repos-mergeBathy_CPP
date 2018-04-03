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
#include "GradientGrid.h"
#include <list>
#include <vector>
#include "geom.h"
#include "mesh.h"
#include "pointList.h"
#include <string>
#include <cstring> //UNIX

class SHullDelaunay
{
	private:
		Edge *startingEdge;
		Edge *hullEdge;
		//std::list < Edge * > edges;
		vector < Edge * > edges; // SJZ

		vector <Point*> ptPtrs;//keep up with allocated new pts for deletion later
		
		PointList *pL;
		
		void insert( const Point& p);
		int checkEdge( Edge* e );
		bool pointOnSurface(const Point& p);
		
		Triangle locateNearestNeighbor(const Point& p, string extrapMethod);
		Triangle locateNearestNeighborTri(const Point& p, string extrapMethod);
		Edge *locate(const Point& p);
		vector< Edge* > locate3NN(const Point& p);
		vector< Edge* > findNearest3(const Point& p, Edge* e, Edge* e3, double distance, Edge* oNextLoop, Edge* u);
		Edge* nearestDest (Edge* e);
		Edge* nearestDestnTri2P (Edge* e, const Point& p, Point& pn);
		Point* computeNearestNeighborWeights(const Point& p, const Triangle& t, const double& sH, const double& alpha, const double& deltaMin, const double& slope);
		Point* computeDepth(const Point& p, const Triangle& t, const double& sH, const double& alpha, const double& deltaMin, string interpMethod);

		//Interpolation Methods
		double bilinearInterp(double x1, double y1, double w1, double x2, double y2, double w2, double x3, double y3, double w3, double x, double y);		
		
	public:

		SHullDelaunay() { pL=NULL; }
		~SHullDelaunay() { clear(); }

		void insert(PointList& pl);

		double sample(const Point& p);
		Mesh* getMesh();

		vector<Point> determinePingLocationInTriangle(PointList& pl, const double& sH, const double& alpha, const double& deltaMin, GradientGrid& grads, vector<Triangle>& triangles, string depthInterpMethod, string uncertInterpMethod, string extrapMethod);
		//void gradientGrid(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn, vector<Gradient*>& gradients);
	
		double getOffsetX() const { return pL->getOffsetX(); }
		double getOffsetY() const { return pL->getOffsetY(); }

		void clear();
};

