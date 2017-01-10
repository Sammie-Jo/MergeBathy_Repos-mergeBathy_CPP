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

