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
//Copyright (c) 2010 Dr David Sinclair http://s-hull.org/

#include <iomanip>
#include <exception>
#include<iostream>
#include<list>
#include<math.h>
#include <algorithm>
#include "geom.h"
#include "mesh.h"
#include "pointList.h"
#include "../regr_xzw.h"
#include "../consistentWeights.h"
//#include "s_hull.h"
#include "../grid.h"
#include "sHullDelaunay.h"
#include "GradientGrid.h"

using namespace std;

void SHullDelaunay::clear()
{
	//list< Edge * >::iterator it;
	vector< Edge * >::iterator it;
	// Clear list of Edges
	for(it=edges.begin(); it != edges.end(); it++)
		delete *it;
	//for(int i=0; i < (int)edges.size(); i++)
	//	delete edges[i];
	edges.clear();

	//Clear list of Points on Edges
	for(int i = 0; i < (int)ptPtrs.size(); i++)
		delete ptPtrs[i];
	ptPtrs.clear();
	
	if(pL != NULL) // SJZ
	{
		delete pL;
		pL = NULL;
	}
}

void SHullDelaunay::insert(PointList& pl)
{
	Point a1, b2, c3;
//	pL = &pl;
	//SHullDelaunay gets its own copy. This is necessary for multi-threading.
	//If not, every thread will delete the same referenced bathyGrid pointlist, an error.
	pL = new PointList(pl);

	//Select Seed - a1
	//This can be any point, I just pick one in the middle of the list
	//This doesn't work because it assumes a1 is index 0
	//Seeding implementation needs to be revisited.
	//For now, a1 is index 0.
	int index = pL->size() / 2;
	//a1 = (*pL)[index];
	a1 = (*pL)[0];

	//Sort based on distance from a1
	pL->sort(a1);
	pL->unique();

	//Closest to Seed - b2
	// 0 is going to be a1, so index 1 is the next closest point
	b2 = (*pL)[1];

	//Smallest Circle - c3
	double radius = 0.0;
	double currentMin = 99999999999.0;
	int i;
	for( i = 2; i < pL->size(); i++)
	{
		radius = radiusOfCircumCircle(a1, b2, (*pL)[i]);
		if(radius < currentMin)
		{
			currentMin = radius;
			c3 = (*pL)[i];
		}
	}

	//init Triangle
	Point *da, *db, *dc;
	da = new Point(a1), db = new Point(b2), dc = new Point(c3);
	ptPtrs.push_back(da);
	ptPtrs.push_back(db);
	ptPtrs.push_back(dc);

	Edge* ea = makeEdge(); // Make first edge
	ea->endPoints(da, db); // Set endpoints of first edge to a and b

	Edge* eb = makeEdge(); // Make second Edge
	splice(ea->sym(), eb); // Connect first to second edge
	eb->endPoints(db, dc); // Set endpoints of second edge to b and c

	Edge* ec = makeEdge(); // Make third Edge
	splice(eb->sym(), ec); // Connect second and third edge
	ec->endPoints(dc, da);
	splice(ec->sym(), ea);

	startingEdge = ea;

	// Ensure that my starting hulledge is facing counter-clockwise
	if(rightOf(*dc, ea))
		hullEdge = ea->sym();
	else
		hullEdge = ea;

	// Add starting edges to edge list
	edges.push_back(ea);
	edges.push_back(eb);
	edges.push_back(ec);

	//Sort based on distance from center of circumscribed circle formed by a1, b2, and c3
	Point cp = centerOfCircumCircle(a1, b2, c3);
	pL->sort(cp);

	//pL->assignID(); //TIN debug
	// Insert all points building radially from Hull
	for(int i = 3; i < pL->size(); i++)
		insert((*pL)[i]);

	// Insert 3 outer edge points, this can be skipped in many cases, but helps
	// in the case of sampling whenever some samples can lie outside the hull of the
	// structure. This will  prevent infinite loops in searching, or bad data samples
	// along the hull of the structure. Samples within triangles containing these points
	// should be ignored.
	Point oe1(0.0, 9999999.0, 0.0);
	Point oe2(9999999.0, -9999999.0, 0.0);
	Point oe3(-9999999.0, -9999999.0, 0.0);
	oe1.hU = 0.0;
	oe1.vU = 0.0;
	oe2.hU = 0.0;
	oe2.vU = 0.0;
	oe3.hU = 0.0;
	oe3.vU = 0.0;
	//oe1.e  = 0.0;
	//oe2.e  = 0.0;
	//oe3.e  = 0.0;
	//oe1.nei  = 0.0;
	//oe2.nei  = 0.0;
	//oe3.nei  = 0.0;
	//oe1.rei  = 0.0;
	//oe2.rei  = 0.0;
	//oe3.rei  = 0.0;
	//oe1.z0  = 0.0;
	//oe2.z0  = 0.0;
	//oe3.z0  = 0.0;
	//oe1.e0  = 0.0;
	//oe2.e0  = 0.0;
	//oe3.e0  = 0.0;
	//oe1.zK  = 0.0;
	//oe2.zK  = 0.0;
	//oe3.zK  = 0.0;
	//oe1.eK  = 0.0;
	//oe2.eK  = 0.0;
	//oe3.eK  = 0.0;

	insert(oe1);
	insert(oe2);
	insert(oe3);

	// Perform edge flipping on all edges until no flips occur.
	//list< Edge * >::iterator it; // SJZ
	vector< Edge * >::iterator it; // SJZ
	int numFlips = 0;
	int oldNumFlips = 0, veryOldNumFlips = 0, count1 = 0, count2 = 0;
	vector<double> pattern, oldpattern;
	int plimit = 20;		//pattern size limit
	int mcnt = 0;			//counts pattern match attempts
	int mlimit = 20;		//pattern match attempt limit, has to be larger than verifylimit; if exceeded increase plimit to attempt to capture pattern
	int verifycnt = 0;		//counts consecutive pattern matches
	int verifylimit = 10;	//consecutive pattern matches needed to verify
	do
	{
		if(numFlips==oldNumFlips && (numFlips!=0 && count2==0))
			count1++;
		else if(numFlips==veryOldNumFlips && (numFlips!=0 && count1==0))
			count2++;
		else {count1=0;count2=0;}
		if(count1==20 || count2==20)
			break;
		if((int)pattern.size() > plimit)
		{
			mcnt++;
			if(!oldpattern.empty() && (std::equal(pattern.begin(),pattern.end(),oldpattern.begin())))
			{
				verifycnt++;
				if(verifycnt == verifylimit)
					break;
			}else verifycnt = 0;
			oldpattern.assign(pattern.begin(),pattern.end());

			if(mcnt > mlimit)
			{
				mcnt = 0;
				plimit++;
				oldpattern.clear();
			}
			pattern.clear();
		}
		pattern.push_back(numFlips);

		veryOldNumFlips = oldNumFlips;
		oldNumFlips = numFlips;
		numFlips = 0;
		for(it = edges.begin(); it != edges.end(); it++)
			numFlips += checkEdge(*it);
	}while(numFlips!=0);
}

void SHullDelaunay::insert(const Point& p)
{
	Edge* temp = hullEdge;
	// Find first edge in counter-clockwise motion that sees the point.
	if(rightOf(p, hullEdge))
		while(rightOf(p, hullEdge->oPrev()->sym())) // Go Backward to first sight
			hullEdge = hullEdge->oPrev()->sym();
	else
		while(!rightOf(p, hullEdge)) // Go Forward to first sight
		{
			hullEdge = hullEdge->sym()->oNext();
			if(temp == hullEdge)
				return;
		}

	// Make first new edge
	temp = hullEdge->oPrev();
	Edge* base = makeEdge();
	edges.push_back(base);

	Point* pPtr = new Point(p);
	ptPtrs.push_back(pPtr);
	
	base->endPoints(hullEdge->org(), pPtr);//new Point(p)); //SJZ 12/12/14
	splice(temp, base);

	// Make additional edges for all vertices of the hull edges
	// that are visible to point p.
	while(rightOf(p, hullEdge)) //startingEdge can see Point
	{
		base = connect(base, hullEdge->sym());
		edges.push_back(base);
		base = base->sym();
		hullEdge = hullEdge->sym()->oNext()->oNext();
	}
}

vector<Point> SHullDelaunay::determinePingLocationInTriangle(PointList& pl, const double& sH, const double& alpha, const double& deltaMin, GradientGrid& grads, vector<Triangle>& triangles, string depthInterpMethod, string uncertInterpMethod, string extrapMethod)
{
	vector<Point> pos = vector<Point>(pl.size());

	Point* recomputedPoint;
	Triangle tri;
	triangles.reserve(pl.size());

	//A. Find our depths and gradients if not already done.
	//If we are pre-splining and need to estimate uncertainty for a pre-splined zGrid or GMT grid
	//then we have already computed our depths and gradients during those routines,
	//so skip down and compute our uncertainties and triangles.
	//If our gradients don't exist, then we don't know our depths either
	//because we aren't pre-splining, they weren't provided in a pre-interpolated file or we are interpolating to a raster.
	if(!grads.exist())
	{
		double ztemp, etemp;
		double avgZ = 0.0;
		double minZ = (double)MAX_INT;
		double maxZ = (double)MIN_INT;
		vector<double> x;
		vector<double> y;
		vector<double> z;

		for (int i = 0; i < pl.size(); i++)
		{
			//Get the triangle from the list of known pings
			if(depthInterpMethod == "NN"){
				//Nearest Neighbor Interpolation
				//Find 3 nearest neighbors
				tri = locateNearestNeighbor(pl[i], extrapMethod);
			}
			else{
				//BILINEAR Interpolation
				//Find encompassing tri and replace
				//invalid points with nearest neighbors
				tri = locateNearestNeighborTri(pl[i], extrapMethod);
			}

			triangles.push_back(tri);
			// Compute Depth via interpMethod chosen
			recomputedPoint = computeDepth(pl[i], tri, sH, alpha, deltaMin, depthInterpMethod);
			pos[i] = *recomputedPoint;
			x.push_back(pos[i].x);
			y.push_back(pos[i].y);
			z.push_back(pos[i].z);
			ztemp = pos[i].z;
			minZ = min(ztemp, minZ);
			maxZ = max(ztemp, maxZ);
			avgZ += ztemp;
			pl.setZ(i,ztemp);
			//Check if interpolating for raster output
			if(uncertInterpMethod == "bilinear")
			{
				//************************************************************
				//etemp = pos[i].e;
				etemp = pos[i].u;
				pl.setE(i,etemp);
				
				etemp = pos[i].nei;
				pl.setNEI(i,etemp);
				
				etemp = pos[i].rei;
				pl.setREI(i,etemp);
				
				ztemp = pos[i].z0;
				pl.setZ0(i,ztemp);
				
				etemp = pos[i].e0;
				pl.setE0(i,etemp);
				
				ztemp = pos[i].zK;
				pl.setZK(i,ztemp);
				
				etemp = pos[i].eK;
				pl.setEK(i,etemp);
				
			}
			delete recomputedPoint;
			recomputedPoint = NULL;
		}

		avgZ /= (double)pl.size();
		pl.setAvgZ(avgZ);
		pl.setMaxZ(maxZ);
		pl.setMinZ(minZ);
		if(uncertInterpMethod == "bilinear")
			return pos;
		grads.calc_GradientGrid(x, y, z);
		x.clear();
		y.clear();
		z.clear();
	}

	vector<double>* slopes = grads.getSlopeOut_Vector();
	
	//B. Find our uncertainty. 
	//Check to see if we have already located our triangles.
	//If not, do that now.  If so, then we just need to
	//compute our weighted uncertainties.
	int FINDTRIS_FLAG = 0;
	if(triangles.empty())
		FINDTRIS_FLAG = 1;
	for (int i = 0; i < pl.size(); i++)
	{
		if(FINDTRIS_FLAG)
		{
			//Get the triangle from the list of known pings
			if(depthInterpMethod == "NN"){
				//Nearest Neighbor Interpolation
				//Find 3 nearest neighbors of all points.
				tri = locateNearestNeighbor(pl[i], extrapMethod);
			}
			else{
				//3 DT NN with Bilinear for Depth and with IDW for Uncertainty
				//Find encompassing triangle and replace
				//invalid points with nearest neighbors
				tri = locateNearestNeighborTri(pl[i], extrapMethod);
			}
			triangles.push_back(tri);
		}
		//Compute weights for interpolation points
		recomputedPoint = computeNearestNeighborWeights(pl[i], triangles[i], sH, alpha, deltaMin, (*slopes)[i]);
		pos[i] = *recomputedPoint;
		delete recomputedPoint;
		recomputedPoint = NULL;
	}
	slopes = NULL;
	triangles.clear();
	return pos;
}

//Computes depth using vertices of triangle t.
//Bilinear interpolation is the intended INTERP_METHOD.
//INTERP_METHOD may be NN for debugging and MATLAB comparisons.
//This is because an exact replication of how MATLAB handles
//points that fall outside the convex hull was not possibly.
//Therefore, only when p is inside the convex hull will bilinear
//interpolation results match.
Point* SHullDelaunay::computeDepth(const Point& p, const Triangle& t, const double& sH, const double& alpha, const double& deltaMin, string interpMethod)
{
	Point *computedPing = new Point();

	double z, e = 0;
	double nei = 0;
	double rei = 0;
	double z0 = 0;
	double e0 = 0;
	double zK = 0;
	double eK = 0;
	int nn = 0;

	//If we have bad points then we 
	//have decided not to extrapolate or
	//locateNearestNeighborTri could not find
	//a valid points to replace the bad points which
	//happens when interpolating from a regular grid
	//to another grid for raster output, but for 
	//this case extrapolation is not allowed anyway.

	// Invalid Outside Boundary Points
	Point oe1(0.0, 9999999.0, 0.0);
	Point oe2(9999999.0, -9999999.0, 0.0);
	Point oe3(-9999999.0, -9999999.0, 0.0);

	//Check for bad points
	//so that it is not used in extrapolation.
	//Attempts to find a valid point failed.
	//This will happen for Raster output if we
	//create a tin from a grid. This will give
	//collinear nearest neighbors.
	int p1Flag=0;
	int p2Flag=0;
	int p3Flag=0;
	if((t.v1==oe1 || t.v1==oe2) || t.v1==oe3)
		p1Flag=1;
	if((t.v2==oe1 || t.v2==oe2) || t.v2==oe3)
		p2Flag=1;
	if((t.v3==oe1 || t.v3==oe2) || t.v3==oe3)
		p3Flag=1;

	//See if p is at an vertex.
	//If we are then its an exact interpolation.
	if(p.x == t.v1.x && p.y == t.v1.y)
	{
		z = t.v1.z;
		//Raster output uses only bilinear so these are only here for that
		nei = t.v1.nei;
		rei = t.v1.rei;
		z0	= t.v1.z0;
		e0	= t.v1.e0;
		zK	= t.v1.zK;
		eK	= t.v1.eK;
	}
	else if(p.x == t.v2.x && p.y == t.v2.y) 
	{
		z = t.v2.z;
		//Raster output uses only bilinear so these are only here for that
		nei = t.v2.nei;
		rei = t.v2.rei;
		z0	= t.v2.z0;
		e0	= t.v2.e0;
		zK	= t.v2.zK;
		eK	= t.v2.eK;
	}
	else if(p.x == t.v3.x && p.y == t.v3.y)
	{
		z = t.v3.z;
		//Raster output uses only bilinear so these are only here for that
		nei = t.v3.nei;
		rei = t.v3.rei;
		z0	= t.v3.z0;
		e0	= t.v3.e0;
		zK	= t.v3.zK;
		eK	= t.v3.eK;
	}
	//Nearest Neighbor interpolation for depth z
	//Find the Nearest Neighbor of the 3 NN
	else if(interpMethod == "NN")
	{
		if((p1Flag || p2Flag) || p3Flag)
			z	= (double)NaN;
		else
		{
			double distv1 = sqrt(t.v1.distance2d(p));
			double distv2 = sqrt(t.v2.distance2d(p));
			double distv3 = sqrt(t.v3.distance2d(p));

			if(distv1 <= distv2 && distv1 <= distv3)
			{
				nn = 1;
				z = t.v1.z;
			}
			else if(distv2 <= distv3)
			{
				nn = 2;
				z = t.v2.z;
			}
			else
			{
				nn = 3;
				z = t.v3.z;
			}
		}
	}
	//Bilinear Interpolation to find depth z;
	else
	{
		if((p1Flag || p2Flag) || p3Flag)
		{ //set pts to Nan; no extrapolation
			z	= (double)NaN;
			e	= (double)NaN;

			//Raster output uses only bilinear so these are only here for that
			nei = (double)NaN;
			rei = (double)NaN;
			z0	= (double)NaN;
			e0	= (double)NaN;
			zK	= (double)NaN;
			eK	= (double)NaN;
		}
		else
		{ //bilinear interpolation
			double x1 = t.v1.x;
			double x2 = t.v2.x;
			double x3 = t.v3.x;
		
			double y1 = t.v1.y;
			double y2 = t.v2.y;
			double y3 = t.v3.y;
		
			//depths
			double w1 = t.v1.z;
			double w2 = t.v2.z;
			double w3 = t.v3.z;
		
			//error
			double w1_e = t.v1.u;
			double w2_e = t.v2.u;
			double w3_e = t.v3.u;
		
			//double w1_e = t.v1.e;
			//double w2_e = t.v2.e;
			//double w3_e = t.v3.e;

			double w1_nei = t.v1.nei;
			double w2_nei = t.v2.nei;
			double w3_nei = t.v3.nei;
			double w1_rei = t.v1.rei;
			double w2_rei = t.v2.rei;
			double w3_rei = t.v3.rei;
		
			double w1_z0 = t.v1.z0;
			double w2_z0 = t.v2.z0;
			double w3_z0 = t.v3.z0;
			double w1_e0 = t.v1.e0;
			double w2_e0 = t.v2.e0;
			double w3_e0 = t.v3.e0;

			double w1_zK = t.v1.zK;
			double w2_zK = t.v2.zK;
			double w3_zK = t.v3.zK;
			double w1_eK = t.v1.eK;
			double w2_eK = t.v2.eK;
			double w3_eK = t.v3.eK;

			//Query point
			double x = p.x;
			double y = p.y;
		
			z	= bilinearInterp(x1, y1, w1, x2, y2, w2, x3, y3, w3, x, y);
			e	= bilinearInterp(x1, y1, w1_e, x2, y2, w2_e, x3,y3,w3_e, x, y);
			nei = bilinearInterp(x1, y1, w1_nei, x2, y2, w2_nei, x3, y3, w3_nei, x, y);
			rei = bilinearInterp(x1, y1, w1_rei, x2, y2, w2_rei, x3, y3, w3_rei, x, y);
			z0	= bilinearInterp(x1, y1, w1_z0, x2, y2, w2_z0, x3, y3, w3_z0, x, y);
			e0	= bilinearInterp(x1, y1, w1_e0, x2, y2, w2_e0, x3, y3, w3_e0, x, y);
			zK	= bilinearInterp(x1, y1, w1_zK, x2, y2, w2_zK, x3, y3, w3_zK, x, y);
			eK	= bilinearInterp(x1, y1, w1_eK, x2, y2, w2_eK, x3, y3, w3_eK, x, y);
		}
	}
	computedPing->x = p.x;
	computedPing->y = p.y;
	computedPing->z = z;
	computedPing->hU = p.hU;
	computedPing->vU = p.vU;
	//************************************************************
	//computedPing->e = e;//Error calculated via bilinear interp
	computedPing->u = e;//Error calculated via bilinear interp

	computedPing->nei	= nei;//Error calculated via bilinear interp
	computedPing->rei	= rei;//Error calculated via bilinear interp
	computedPing->z0	= z0;//Error calculated via bilinear interp
	computedPing->e0	= e0;//Error calculated via bilinear interp
	computedPing->zK	= zK;//Error calculated via bilinear interp
	computedPing->eK	= eK;//Error calculated via bilinear interp

	return computedPing;
}

double SHullDelaunay::bilinearInterp(double x1, double y1, double w1, double x2, double y2, double w2, double x3, double y3, double w3, double x, double y)
{
	//Bilinear Interpolation
	double w;
	double DET	= x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3;
	if(DET == 0 || DET == NaN)
		w = (double)NaN;
	else{
		double A 	= ((y2-y3)*w1+(y3-y1)*w2+(y1-y2)*w3) / DET;
		double B 	= ((x3-x2)*w1+(x1-x3)*w2+(x2-x1)*w3) / DET;
		double C 	= ((x2*y3-x3*y2)*w1+(x3*y1-x1*y3)*w2+(x1*y2-x2*y1)*w3) / DET;
		w 	= A*x+B*y+C;
	}
	return w;
}

//Computes the weighted uncertainty from the 3 vertices of the containing tri or nearest neighbor
Point* SHullDelaunay::computeNearestNeighborWeights(const Point& p, const Triangle& t, const double& sH, const double& alpha, const double& deltaMin, const double& slope)
{
	Point *computedPing = new Point();
	double sigma=0,sigma2=0,sigma3=0,sigma4=0,sigma5=0;

	//if(p == t.v1)
	//	sigma = t.v1.e;	
	//else if(p == t.v2)
	//	sigma = t.v2.e;	
	//else if(p == t.v3)
	//	sigma = t.v3.e;	

	double minNodeSpacing = min(t.v1.distance2d(t.v2),t.v2.distance2d(t.v3));
	minNodeSpacing = min(minNodeSpacing, t.v3.distance2d(t.v1));

	double distv1 = sqrt(t.v1.distance2d(p));
	double distv2 = sqrt(t.v2.distance2d(p));
	double distv3 = sqrt(t.v3.distance2d(p));
	double slopeP = tan((PI/180.0)*slope);

	//These numbers match equations in a MATLAB function Paul' developed.
	//1.Regular inverse weighted distance
	double weightN1_C1 = pow(t.v1.vU,2) * (1.00 + pow( (distv1 + sH*t.v1.hU) / deltaMin, alpha) );
	double weightN2_C1 = pow(t.v2.vU,2) * (1.00 + pow( (distv2 + sH*t.v2.hU) / deltaMin, alpha) );
	double weightN3_C1 = pow(t.v3.vU,2) * (1.00 + pow( (distv3 + sH*t.v3.hU) / deltaMin, alpha) );

	//3.Regular inverse weighted distance
	double weightN1_C3 = pow(t.v1.vU,2) + pow(slopeP*sH*t.v1.hU,alpha);
	double weightN2_C3 = pow(t.v2.vU,2) + pow(slopeP*sH*t.v2.hU,alpha);
	double weightN3_C3 = pow(t.v3.vU,2) + pow(slopeP*sH*t.v3.hU,alpha);

	//2.Regular inverse weighted distance
 	double weightN1_C2 = pow(t.v1.vU,2) * (1.00 + pow( (distv1 + slopeP*sH*t.v1.hU) / deltaMin, alpha) );
	double weightN2_C2 = pow(t.v2.vU,2) * (1.00 + pow( (distv2 + slopeP*sH*t.v2.hU) / deltaMin, alpha) );
	double weightN3_C2 = pow(t.v3.vU,2) * (1.00 + pow( (distv3 + slopeP*sH*t.v3.hU) / deltaMin, alpha) );

	//4.Regular inverse weighted distance
	double weightN1_C4 = weightN1_C1 + pow(slopeP*sH*t.v1.hU, alpha);
	double weightN2_C4 = weightN2_C1 + pow(slopeP*sH*t.v2.hU, alpha);
	double weightN3_C4 = weightN3_C1 + pow(slopeP*sH*t.v3.hU, alpha);

	//5.Regular inverse weighted distance if iDW then needs to move down and sqr
	double weightN1_C5 = weightN1_C2 + pow(slopeP*sH*t.v1.hU, alpha);
	double weightN2_C5 = weightN2_C2 + pow(slopeP*sH*t.v2.hU, alpha);
	double weightN3_C5 = weightN3_C2 + pow(slopeP*sH*t.v3.hU, alpha);

	int nn = 0;
	////Find NN
	////Comment out for IDW. 
	////To test the function comment in distv1 ifelses and comment out the other.
	////if(distv1<=distv2 && distv1<=distv3) 
	//if(p == t.v1)
	//{
	//	nn = 1;
	//	sigma= sqrt(weightN1_C1);
	//	sigma2= sqrt(weightN1_C2);
	//	sigma3= sqrt(weightN1_C3);
	//	sigma4= sqrt(weightN1_C4);
	//	sigma5= sqrt(weightN1_C5);
	//}
	////else if(distv2<=distv3)
	//else if(p == t.v2)
	//{
	//	nn = 2;
	//	sigma= sqrt(weightN2_C1);
	//	sigma2= sqrt(weightN2_C2);
	//	sigma3= sqrt(weightN2_C3);
	//	sigma4= sqrt(weightN2_C4);
	//	sigma5= sqrt(weightN2_C5);
	//}
	////else
	//else if(p == t.v3)
	//{
	//	nn = 3;
	//	sigma= sqrt(weightN3_C1);
	//	sigma2= sqrt(weightN3_C2);
	//	sigma3= sqrt(weightN3_C3);
	//	sigma4= sqrt(weightN3_C4);
	//	sigma5= sqrt(weightN3_C5);
	//}//Comment out for IDW
	////if((p != t.v1 && p != t.v2) && p !=  t.v3)
	//else
	{
		//Added to handle division by 0.  
		//Don't know why this wasn't caught before therefore it needs 
		//to be double checked to make sure scenarios where division 
		//by 0 occur are allowed (ie. if this is because we are at 
		//an grid location, make sure that we can use grid location 
		//points in computations). --SJZ 9/7/16
		
		//Compute the inverse weighted distance from all three neighbors
		if(abs(distv1) < eps)
		{
			sigma = sqrt(weightN1_C1);
			sigma2 = sqrt(weightN1_C2);
			sigma3 = sqrt(weightN1_C3);
			sigma4 = sqrt(weightN1_C4);
			sigma5 = sqrt(weightN1_C5);
		}
		else if(distv2 == 0)
		{
			sigma = sqrt(weightN2_C1);
			sigma2 = sqrt(weightN2_C2);
			sigma3 = sqrt(weightN2_C3);
			sigma4 = sqrt(weightN2_C4);
			sigma5 = sqrt(weightN2_C5);
		}
		else if(distv3 == 0)
		{
			sigma = sqrt(weightN3_C1);
			sigma2 = sqrt(weightN3_C2);
			sigma3 = sqrt(weightN3_C3);
			sigma4 = sqrt(weightN3_C4);
			sigma5 = sqrt(weightN3_C5);
		}
		else
		{
			double d1_C1 = weightN1_C1/distv1;
			double d2_C1 = weightN2_C1/distv2;
			double d3_C1 = weightN3_C1/distv3;

			double d1_C2 = weightN1_C2/distv1;
			double d2_C2 = weightN2_C2/distv2;
			double d3_C2 = weightN3_C2/distv3;

			double d1_C3 = weightN1_C3/distv1;
			double d2_C3 = weightN2_C3/distv2;
			double d3_C3 = weightN3_C3/distv3;

			double d1_C4 = weightN1_C4/distv1;
			double d2_C4 = weightN2_C4/distv2;
			double d3_C4 = weightN3_C4/distv3;

			double d1_C5 = weightN1_C5/distv1;
			double d2_C5 = weightN2_C5/distv2;
			double d3_C5 = weightN3_C5/distv3;
			double distSum = (1.00/distv1) + (1.00/distv2) + (1.00/distv3);
			sigma = sqrt((d1_C1 + d2_C1 + d3_C1)/distSum); //Comment In for IDW
			sigma2 = sqrt((d1_C2 + d2_C2 + d3_C2)/distSum); //Comment In for IDW
			sigma3 = sqrt((d1_C3 + d2_C3 + d3_C3)/distSum); //Comment In for IDW
			sigma4 = sqrt((d1_C4 + d2_C4 + d3_C4)/distSum); //Comment In for IDW
			sigma5 = sqrt((d1_C5 + d2_C5 + d3_C5)/distSum); //Comment In for IDW
		//	sigma = sqrt((weightN1 + weightN2 + weightN3)/3);
		}
	}
	computedPing->x = p.x;
	computedPing->y = p.y;
	computedPing->z = p.z;
	computedPing->hU = p.hU;
	computedPing->vU = p.vU;
	computedPing->u = sigma4;//sigma; swap cube's with curve's
	computedPing->u2 = sigma2;
	computedPing->u3 = sigma3;
	computedPing->u4 = sigma;//sigma4;
	computedPing->u5 = sigma5;
	computedPing->s = slopeP;

	return computedPing;
}

// Finds the 3 points closest to my query point p regardless
// if they are a vertex of the encompassing triangle.
Triangle SHullDelaunay::locateNearestNeighbor(const Point& p, string extrapMethod)
{
	if(extrapMethod == "none")//Raster
	{
		//Find which triangle p is in first to see if it lies 
		//outside the convex hull.  If it does we return 
		//since we don't want to extrapolate.
		Triangle t = locateNearestNeighborTri(p, extrapMethod);
				
		// Invalid Outside Boundary Points
		Point oe1(0.0, 9999999.0, 0.0);
		Point oe2(9999999.0, -9999999.0, 0.0);
		Point oe3(-9999999.0, -9999999.0, 0.0);

		//Check for bad points again! Won't necessarily return bad pts.
		//We will have bad points when p is outside the convex hull.
		//Did not attempt to replace bad points.
		int p1Flag=0;
		int p2Flag=0;
		int p3Flag=0;
		if((t.v1==oe1 || t.v1==oe2) || t.v1==oe3)
			p1Flag=1;
		if((t.v2==oe1 || t.v2==oe2) || t.v2==oe3)
			p2Flag=1;
		if((t.v3==oe1 || t.v3==oe2) || t.v3==oe3)
			p3Flag=1;

		//Return if we are outside the hull
		if(((p1Flag || p2Flag) || p3Flag))
			return t;
	}

	//Find 3 nearest neighbors of all points
	vector<Edge*> e = locate3NN(p);
	Triangle t;

	if(e[0]->id == -1)
	{
		t = Triangle();
		return t;
	}

	Point p1, p2, p3;
	p1 = e[0]->org2d();
	p2 = e[0]->dest2d();
	p3 = e[1]->dest2d();

	t = Triangle(p1, p2, p3);
	return t;
}

// Finds the encompassing triangle of query point p.
// If p is outside the natural convex hull,
// a triangle with invalid boundary points will be selected.
// These n-invalid points are replaced with query point
// p's n-nearest neighbors of all points which are valid and not already selected.
// Invalid boundary points are:
// Point oe1(0.0, 9999999.0, 0.0);
// Point oe2(9999999.0, -9999999.0, 0.0);
// Point oe3(-9999999.0, -9999999.0, 0.0);
Triangle SHullDelaunay::locateNearestNeighborTri(const Point& p, string extrapMethod)
{
	int PRINT_WARNINGS = 0; //degenerate cases?
	int printDebugFlag = 0;
	//Find triangle containing p
	Edge* e = locate(p);
	Triangle t;

	if(e->id == -1)
	{
		t = Triangle();
		return t;
	}

	Point p1, p2, p3;
	p1 = e->org2d();
	p2 = e->dest2d();
	p3 = e->oNext()->dest2d();

	//We don't need to replace invalid points
	//because we don't want to extrapolate.
	//Since the point lies outside the convex hull
	if(extrapMethod == "none") //Raster
	{
		t = Triangle(p1, p2, p3);
		return t;
	}

	// Invalid Outside Boundary Points
	Point oe1(0.0, 9999999.0, 0.0);
	Point oe2(9999999.0, -9999999.0, 0.0);
	Point oe3(-9999999.0, -9999999.0, 0.0);

	//Set flag when point falls outside convex hull
	//so that it is not used in extrapolation.
	//A valid point needs to be found to replace the boundary point.
	int p1Flag=0;
	int p2Flag=0;
	int p3Flag=0;
	if((p1==oe1 || p1==oe2) || p1==oe3)
		p1Flag=1;
	if((p2==oe1 || p2==oe2) || p2==oe3)
		p2Flag=1;
	if((p3==oe1 || p3==oe2) || p3==oe3)
		p3Flag=1;

	#pragma region --p1 & p2 Invalid
	if((p1Flag && p2Flag) && !p3Flag)
	{
		//*************************************************************************
		//p1 and p2 are invalid points.
		//We only have 1 valid point p3.
		//Find the 2 non-collinear nearest neighboring points to p3.
		//*************************************************************************
		Edge* e2 = e->sym()->oPrev();		//get edge containing valid pt p2 and p3

		//I. Replace invalid point with query point's nearest neighbor if it is not already one of the valid points
		vector<Edge*> en = locate3NN(p);
		Point p1n = en[0]->org2d();
		Point p2n = en[0]->dest2d();
		Point p3n = en[1]->dest2d();
		vector<Point> nntemp;
		nntemp.push_back(p1n);
		nntemp.push_back(p2n);
		nntemp.push_back(p3n);

		//Sort by distance from p
		MinDistance md;
		md.c = p;
		std::sort(nntemp.begin(), nntemp.end(), md);

		//Find if already a valid pt
		int p1nSameFlag = 0;
		int p2nSameFlag = 0;
		int p3nSameFlag = 0;
		if(p1n == p3)
			p1nSameFlag = 1;
		if(p2n == p3)
			p2nSameFlag = 1;
		if(p3n == p3)
			p3nSameFlag = 1;
		if((p1nSameFlag && p2nSameFlag) && p3nSameFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are the same as what we started with!"<<endl;

		//Find if invalid pt
		int p1nBadFlag = 0;
		int p2nBadFlag = 0;
		int p3nBadFlag = 0;
		if((p1n==oe1 || p1n==oe2) || p1n==oe3)
			p1nBadFlag = 1;
		if((p2n==oe1 || p2n==oe2) || p2n==oe3)
			p2nBadFlag = 1;
		if((p3n==oe1 || p3n==oe2) || p3n==oe3)
			p3nBadFlag = 1;
		if((p1nBadFlag && p2nBadFlag) && p3nBadFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are invalid boundary points!  Where's the hull?!"<<endl;

		//Find if collinear pt
		int p1nCoFlag = 0;
		int p2nCoFlag = 0;
		int p3nCoFlag = 0;
		if(collinear(p1n,e2))
			p1nCoFlag = 1;
		if(collinear(p2n, e2))
			p2nCoFlag = 1;
		if(collinear(p3n, e2))
			p3nCoFlag = 1;
		if((p1nCoFlag && p2nCoFlag) && p3nCoFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are collinear with initial valid points!"<<endl;

		//Select the nearest neighbor to replace p1 that is valid, new, and non-collinear
		int p1nUsedFlag = 0;
		int p2nUsedFlag = 0;
		int p3nUsedFlag = 0;
		int cnt=2;
		Point tempP = p1;
		while(cnt){
			if(((!p1nSameFlag && !p1nBadFlag) && !p1nCoFlag) && !p1nUsedFlag){
				tempP = p1n;
				p1nUsedFlag = 1;
			}else if(((!p2nSameFlag && !p2nBadFlag) && !p2nCoFlag) && !p2nUsedFlag){
				tempP = p2n;
				p2nUsedFlag = 1;
			}else if(((!p3nSameFlag && !p3nBadFlag) && !p3nCoFlag) && !p3nUsedFlag){
				tempP = p3n;
				p3nUsedFlag = 1;
			}else if((!p1nUsedFlag && !p2nUsedFlag) && !p3nUsedFlag){
				if(PRINT_WARNINGS)
					cerr<<"WARNING: Point is far from hull.  All nearest neighbor points found are the same as what we started with and the nearest neighbor is invalid."<<endl;
			}else if(p1nUsedFlag + p2nUsedFlag + p3nUsedFlag == 1){
				if(PRINT_WARNINGS)
					cerr<<"WARNING: Only replaced 1 invalid point out of 2!"<<endl;
			}

			cnt--;
			if(cnt)
			{
				p1 = tempP;
				tempP = p2;
			}
			else p2 = tempP;
		}

		if(printDebugFlag){
			cout<<"**P1 & P2 Replaced by p's NNs out of all points: "<<endl;
			p1.print();		//new p1
			p2.print();		//new p2
			p3.print();		//same p3
		}
	}
	#pragma endregion p1 & p2 Invalid

	#pragma region --p1 & p3 Invalid
	else if((p1Flag && p3Flag) && !p2Flag)
	{
		//*************************************************************************
		//p1 and p3 are bad points
		//Find the 2 non-collinear nearest neighboring points to p2.
		//*************************************************************************
		Edge* e2 = e->sym()->oPrev();		//get edge containing valid pt p2 and p3

		//I. Replace invalid point with query point's nearest neighbor if it is not already one of the valid points
		vector<Edge*> en = locate3NN(p);
		Point p1n = en[0]->org2d();
		Point p2n = en[0]->dest2d();
		Point p3n = en[1]->dest2d();
		vector<Point> nntemp;
		nntemp.push_back(p1n);
		nntemp.push_back(p2n);
		nntemp.push_back(p3n);

		//Sort by distance from p
		MinDistance md;
		md.c = p;
		std::sort(nntemp.begin(), nntemp.end(), md);

		//Find if already a valid pt
		int p1nSameFlag = 0;
		int p2nSameFlag = 0;
		int p3nSameFlag = 0;
		if(p1n == p2)
			p1nSameFlag = 1;
		if(p2n == p2)
			p2nSameFlag = 1;
		if(p3n == p2)
			p3nSameFlag = 1;
		if((p1nSameFlag && p2nSameFlag) && p3nSameFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are the same as what we started with!"<<endl;

		//Find if invalid pt
		int p1nBadFlag = 0;
		int p2nBadFlag = 0;
		int p3nBadFlag = 0;
		if((p1n==oe1 || p1n==oe2) || p1n==oe3)
			p1nBadFlag = 1;
		if((p2n==oe1 || p2n==oe2) || p2n==oe3)
			p2nBadFlag = 1;
		if((p3n==oe1 || p3n==oe2) || p3n==oe3)
			p3nBadFlag = 1;
		if((p1nBadFlag && p2nBadFlag) && p3nBadFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are invalid boundary points!  Where's the hull?!"<<endl;

		//Find if collinear pt
		int p1nCoFlag = 0;
		int p2nCoFlag = 0;
		int p3nCoFlag = 0;
		if(collinear(p1n,e2))
			p1nCoFlag = 1;
		if(collinear(p2n, e2))
			p2nCoFlag = 1;
		if(collinear(p3n, e2))
			p3nCoFlag = 1;
		if((p1nCoFlag && p2nCoFlag) && p3nCoFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are collinear with initial valid points!"<<endl;

		//Select the nearest neighbor to replace p1 that is valid, new, and non-collinear
		int p1nUsedFlag = 0;
		int p2nUsedFlag = 0;
		int p3nUsedFlag = 0;
		int cnt=2;
		Point tempP = p1;
		while(cnt){
			if(((!p1nSameFlag && !p1nBadFlag) && !p1nCoFlag) && !p1nUsedFlag){
				tempP = p1n;
				p1nUsedFlag = 1;
			}else if(((!p2nSameFlag && !p2nBadFlag) && !p2nCoFlag) && !p2nUsedFlag){
				tempP = p2n;
				p2nUsedFlag = 1;
			}else if(((!p3nSameFlag && !p3nBadFlag) && !p3nCoFlag) && !p3nUsedFlag){
				tempP = p3n;
				p3nUsedFlag = 1;
			}else if((!p1nUsedFlag && !p2nUsedFlag) && !p3nUsedFlag){	//==0
				if(PRINT_WARNINGS)
					cerr<<"WARNING: Point is far from hull.  All nearest neighbor points found are the same as what we started with and the nearest neighbor is invalid."<<endl;
			}else if(p1nUsedFlag + p2nUsedFlag + p3nUsedFlag == 1){
				if(PRINT_WARNINGS)
					cerr<<"WARNING: Only replaced 1 invalid point out of 2!"<<endl;
			}

			cnt--;
			if(cnt)
			{
				p1 = tempP;	
				tempP = p3;
			}
			else p3 = tempP;
		}

		if(printDebugFlag){
			cout<<"***P1 & P3 Replaced by p's NNs out of all points: "<<endl;
			p1.print();		//new p1
			p2.print();		//same p2
			p3.print();		//new p3
		}
	}
	#pragma endregion p1 & p3 Invalid

	#pragma region --p2 & p3 Invalid
	else if((p2Flag && p3Flag) && !p1Flag)
	{
		//*************************************************************************
		//p2 and p3 are bad points
		//Find the 2 non-collinear nearest neighboring points to p1.
		//*************************************************************************
		Edge* e2 = e->sym()->oPrev();		//get edge containing valid pt p2 and p3

		//I. Replace invalid point with query point's nearest neighbor if it is not already one of the valid points
		vector<Edge*> en = locate3NN(p);
		Point p1n = en[0]->org2d();
		Point p2n = en[0]->dest2d();
		Point p3n = en[1]->dest2d();
		vector<Point> nntemp;
		nntemp.push_back(p1n);
		nntemp.push_back(p2n);
		nntemp.push_back(p3n);

		//Sort by distance from p
		MinDistance md;
		md.c = p;
		std::sort(nntemp.begin(), nntemp.end(), md);

		//Find if already a valid pt
		int p1nSameFlag = 0;
		int p2nSameFlag = 0;
		int p3nSameFlag = 0;
		if(p1n == p1)
			p1nSameFlag = 1;
		if(p2n == p1)
			p2nSameFlag = 1;
		if(p3n == p1)
			p3nSameFlag = 1;
		if((p1nSameFlag && p2nSameFlag) && p3nSameFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are the same as what we started with!"<<endl;

		//Find if invalid pt
		int p1nBadFlag = 0;
		int p2nBadFlag = 0;
		int p3nBadFlag = 0;
		if((p1n==oe1 || p1n==oe2) || p1n==oe3)
			p1nBadFlag = 1;
		if((p2n==oe1 || p2n==oe2) || p2n==oe3)
			p2nBadFlag = 1;
		if((p3n==oe1 || p3n==oe2) || p3n==oe3)
			p3nBadFlag = 1;
		if((p1nBadFlag && p2nBadFlag) && p3nBadFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are invalid boundary points!  Where's the hull?!"<<endl;

		//Find if collinear pt
		int p1nCoFlag = 0;
		int p2nCoFlag = 0;
		int p3nCoFlag = 0;
		if(collinear(p1n,e2))
			p1nCoFlag = 1;
		if(collinear(p2n, e2))
			p2nCoFlag = 1;
		if(collinear(p3n, e2))
			p3nCoFlag = 1;
		if((p1nCoFlag && p2nCoFlag) && p3nCoFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are collinear with initial valid points!"<<endl;

		//Select the nearest neighbor to replace p1 that is valid, new, and non-collinear
		int p1nUsedFlag = 0;
		int p2nUsedFlag = 0;
		int p3nUsedFlag = 0;
		int cnt=2;
		Point tempP = p2;
		while(cnt){
			if(((!p1nSameFlag && !p1nBadFlag) && !p1nCoFlag) && !p1nUsedFlag){
				tempP = p1n;
				p1nUsedFlag = 1;
			}else if(((!p2nSameFlag && !p2nBadFlag) && !p2nCoFlag) && !p2nUsedFlag){
				tempP = p2n;
				p2nUsedFlag = 1;
			}else if(((!p3nSameFlag && !p3nBadFlag) && !p3nCoFlag) && !p3nUsedFlag){
				tempP = p3n;
				p3nUsedFlag = 1;
			}else if((!p1nUsedFlag && !p2nUsedFlag) && !p3nUsedFlag){
				if(PRINT_WARNINGS)
					cerr<<"WARNING: Point is far from hull.  All nearest neighbor points found are the same as what we started with and the nearest neighbor is invalid."<<endl;
			}else if(p1nUsedFlag + p2nUsedFlag + p3nUsedFlag == 1){
				if(PRINT_WARNINGS)
					cerr<<"WARNING: Only replaced 1 invalid point out of 2!"<<endl;
			}

			cnt--;
			if(cnt)
			{
				p2 = tempP;
				tempP = p3;
			}
			else p3 = tempP;
		}

		if(printDebugFlag){
			cout<<"***P2 & P3 Replaced by p's NNs out of all points: "<<endl;
			p1.print();		//same p1
			p2.print();		//new p2
			p3.print();		//new p3
		}
	}
	#pragma endregion p2 & p3 Invalid

	#pragma region --p1 Invalid
	else if(p1Flag && (!p2Flag && !p3Flag))
	{
		//*************************************************************************
		//p1 is invalid only
		//Find 1 non-collinear nearest neighbor for p1.
		//***************************************************************************
		Edge* e2 = e->sym()->oPrev();		//get edge containing valid pt p2 and p3

		//I. Replace invalid point with query point's nearest neighbor if it is not already one of the valid points
		vector<Edge*> en = locate3NN(p);
		Point p1n = en[0]->org2d();
		Point p2n = en[0]->dest2d();
		Point p3n = en[1]->dest2d();
		vector<Point> nntemp;
		nntemp.push_back(p1n);
		nntemp.push_back(p2n);
		nntemp.push_back(p3n);

		//Sort by distance from p
		MinDistance md;
		md.c = p;
		std::sort(nntemp.begin(), nntemp.end(), md);

		//Find if already a valid pt
		int p1nSameFlag = 0;
		int p2nSameFlag = 0;
		int p3nSameFlag = 0;
		if(p1n == p2 || p1n == p3)
			p1nSameFlag = 1;
		if(p2n == p2 || p2n == p3)
			p2nSameFlag = 1;
		if(p3n == p2 || p3n ==p3)
			p3nSameFlag = 1;
		if((p1nSameFlag && p2nSameFlag) && p3nSameFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are the same as what we started with!"<<endl;

		//Find if invalid pt
		int p1nBadFlag = 0;
		int p2nBadFlag = 0;
		int p3nBadFlag = 0;
		if((p1n==oe1 || p1n==oe2) || p1n==oe3)
			p1nBadFlag = 1;
		if((p2n==oe1 || p2n==oe2) || p2n==oe3)
			p2nBadFlag = 1;
		if((p3n==oe1 || p3n==oe2) || p3n==oe3)
			p3nBadFlag = 1;
		if((p1nBadFlag && p2nBadFlag) && p3nBadFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are invalid boundary points!  Where's the hull?!"<<endl;

		//Find if collinear pt
		int p1nCoFlag = 0;
		int p2nCoFlag = 0;
		int p3nCoFlag = 0;
		if(collinear(p1n,e2))
			p1nCoFlag = 1;
		if(collinear(p2n, e2))
			p2nCoFlag = 1;
		if(collinear(p3n, e2))
			p3nCoFlag = 1;
		if((p1nCoFlag && p2nCoFlag) && p3nCoFlag)//Raster extrapolation
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are collinear with initial valid points!"<<endl;

		//Select the nearest neighbor to replace p1 that is valid, new, and non-collinear
		if((!p1nSameFlag && !p1nBadFlag) && !p1nCoFlag)
			p1 = p1n;
		else if((!p2nSameFlag && !p2nBadFlag) && !p2nCoFlag)
			p1 = p2n;
		else if((!p3nSameFlag && !p3nBadFlag) && !p3nCoFlag)
			p1 = p3n;
		else if(PRINT_WARNINGS)
			cerr<<"WARNING: Point is far from hull.  All nearest neighbor points found are the same as what we started with and the nearest neighbor is invalid."<<endl;

		if(printDebugFlag){
			cout<<"***P1 Replaced by p's NN out of all points: "<<endl;
			p1.print();		//new p1 to replace bad org point
			p2.print();		//same p2
			p3.print();		//same p3
		}
	}
	#pragma endregion p1 Invalid

	#pragma region --p2 Invalid
	else if(p2Flag && (!p1Flag && !p3Flag))
	{
		//*************************************************************************
		//p2 is only bad point
		//Find 1 non-collinear nearest neighbor for p2.
		//*************************************************************************
		Edge* e2 = e->sym()->oPrev();		//get edge containing valid pt p2 and p3

		//I. Replace invalid point with query point's nearest neighbor if it is not already one of the valid points
		vector<Edge*> en = locate3NN(p);
		Point p1n = en[0]->org2d();
		Point p2n = en[0]->dest2d();
		Point p3n = en[1]->dest2d();
		vector<Point> nntemp;
		nntemp.push_back(p1n);
		nntemp.push_back(p2n);
		nntemp.push_back(p3n);

		//Sort by distance from p
		MinDistance md;
		md.c = p;
		std::sort(nntemp.begin(), nntemp.end(), md);

		//Find if already a valid pt
		int p1nSameFlag = 0;
		int p2nSameFlag = 0;
		int p3nSameFlag = 0;
		if(p1n == p1 || p1n == p3)
			p1nSameFlag = 1;
		if(p2n == p1 || p2n == p3)
			p2nSameFlag = 1;
		if(p3n == p1 || p3n ==p3)
			p3nSameFlag = 1;
		if((p1nSameFlag && p2nSameFlag) && p3nSameFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are the same as what we started with!"<<endl;

		//Find if invalid pt
		int p1nBadFlag = 0;
		int p2nBadFlag = 0;
		int p3nBadFlag = 0;
		if((p1n==oe1 || p1n==oe2) || p1n==oe3)
			p1nBadFlag = 1;
		if((p2n==oe1 || p2n==oe2) || p2n==oe3)
			p2nBadFlag = 1;
		if((p3n==oe1 || p3n==oe2) || p3n==oe3)
			p3nBadFlag = 1;
		if((p1nBadFlag && p2nBadFlag) && p3nBadFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are invalid boundary points!  Where's the hull?!"<<endl;

		//Find if collinear pt
		int p1nCoFlag = 0;
		int p2nCoFlag = 0;
		int p3nCoFlag = 0;
		if(collinear(p1n,e2))
			p1nCoFlag = 1;
		if(collinear(p2n, e2))
			p2nCoFlag = 1;
		if(collinear(p3n, e2))
			p3nCoFlag = 1;
		if((p1nCoFlag && p2nCoFlag) && p3nCoFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are collinear with initial valid points!"<<endl;

		//Select the nearest neighbor to replace p1 that is valid, new, and non-collinear
		if((!p1nSameFlag && !p1nBadFlag) && !p1nCoFlag)
			p2 = p1n;
		else if((!p2nSameFlag && !p2nBadFlag) && !p2nCoFlag)
			p2 = p2n;
		else if((!p3nSameFlag && !p3nBadFlag) && !p3nCoFlag)
			p2 = p3n;
		else if(PRINT_WARNINGS)
			cerr<<"WARNING: Point is far from hull.  All nearest neighbor points found are the same as what we started with and the nearest neighbor is invalid."<<endl;

		if(printDebugFlag){
			cout<<"***P2 Replaced by p's NN out of all points: "<<endl;
			p1.print();		//same p1
			p2.print();		//new p2
			p3.print();		//same p3
		}
	}
	#pragma endregion p2 Invalid

	#pragma region --p3 Invalid
	else if(p3Flag && (!p1Flag && !p2Flag))
	{
		//*************************************************************************
		//p3 is an extremum only
		//Find 1 non-collinear nearest neighbor for p3.
		//*************************************************************************
		Edge* e2 = e->sym()->oPrev();		//get edge containing valid pt p2 and p3

		//I. Replace invalid point with query point's nearest neighbor if it is not already one of the valid points
		vector<Edge*> en = locate3NN(p);
		Point p1n = en[0]->org2d();
		Point p2n = en[0]->dest2d();
		Point p3n = en[1]->dest2d();
		vector<Point> nntemp;
		nntemp.push_back(p1n);
		nntemp.push_back(p2n);
		nntemp.push_back(p3n);

		//Sort by distance from p
		MinDistance md;
		md.c = p;
		std::sort(nntemp.begin(), nntemp.end(), md);

		//Find if already a valid pt
		int p1nSameFlag = 0;
		int p2nSameFlag = 0;
		int p3nSameFlag = 0;
		if(p1n == p2 || p1n == p1)
			p1nSameFlag = 1;
		if(p2n == p2 || p2n == p1)
			p2nSameFlag = 1;
		if(p3n == p2 || p3n ==p1)
			p3nSameFlag = 1;
		if((p1nSameFlag && p2nSameFlag) && p3nSameFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are the same as what we started with!"<<endl;

		//Find if invalid pt
		int p1nBadFlag = 0;
		int p2nBadFlag = 0;
		int p3nBadFlag = 0;
		if((p1n==oe1 || p1n==oe2) || p1n==oe3)
			p1nBadFlag = 1;
		if((p2n==oe1 || p2n==oe2) || p2n==oe3)
			p2nBadFlag = 1;
		if((p3n==oe1 || p3n==oe2) || p3n==oe3)
			p3nBadFlag = 1;
		if((p1nBadFlag && p2nBadFlag) && p3nBadFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are invalid boundary points!  Where's the hull?!"<<endl;

		//Find if collinear pt
		int p1nCoFlag = 0;
		int p2nCoFlag = 0;
		int p3nCoFlag = 0;
		if(collinear(p1n,e2))
			p1nCoFlag = 1;
		if(collinear(p2n, e2))
			p2nCoFlag = 1;
		if(collinear(p3n, e2))
			p3nCoFlag = 1;
		if((p1nCoFlag && p2nCoFlag) && p3nCoFlag)
			if(PRINT_WARNINGS)
				cerr<<"WARNING: All nearest neighbor points found are collinear with initial valid points!"<<endl;

		//Select the nearest neighbor to replace p1 that is valid, new, and non-collinear
		if((!p1nSameFlag && !p1nBadFlag) && !p1nCoFlag)
			p3 = p1n;
		else if((!p2nSameFlag && !p2nBadFlag) && !p2nCoFlag)
			p3 = p2n;
		else if((!p3nSameFlag && !p3nBadFlag) && !p3nCoFlag)
			p3 = p3n;
		else if(PRINT_WARNINGS)
			cerr<<"WARNING: Point is far from hull.  All nearest neighbor points found are the same as what we started with and the nearest neighbor is invalid."<<endl;

		if(printDebugFlag){
			cout<<"***P3 Replaced by p's NN out of all points: "<<endl;
			p1.print();		//same p1
			p2.print();		//same p2
			p3.print();		//new p3
		}
	}
	#pragma endregion p3 Invalid

	#pragma region --All Invalid
	else if((p1Flag && p2Flag) && p3Flag)
	{
		//All points are invalid.  This should never happen.
		cerr<<"This should never happen but let's add a safety net and check just in case."<<endl;
	 	cerr<<"NO Points in Delaunay Tri! If using preInterpolatedLocations, make sure lat/lon."<<endl;
	//preInterpolatedLocations hits this everytime, need to investigate and fix. //SJZ 9/19/16
		exit(0);
	}
	#pragma endregion All Invalid

	#pragma region --None Invalid
	else
	{
		if(printDebugFlag)
			cout << "All points valid!" << endl;
	}
	#pragma endregion None Invalid

	t = Triangle(p1, p2, p3);
	return t;
}

//Finds the nearest destination (AKA shortest edge) of a given edge by pivoting
//its origin and also finds the nearest vertex (either the left or the right
//of the nearest dest) that will give the fattest tri.
//Used to replace invalid boundary points for extrapolation when a point falls outside the convex hull
//Should pass the edge with the nearest pt at the origin if replacing a bad point.
Edge* SHullDelaunay::nearestDestnTri2P (Edge* e, const Point& p, Point& pn)
{
	double distance2, distance3;
	double distance;
	double distanceOld;
	Edge* uOld;
	Edge* u = e;
	distanceOld = p.distance2d(u->dest2d());

	while(TRUE)
	{
		//get next ccw edge around origin
		u = u->oNext();
		distance = p.distance2d(u->dest2d());
		if(e == u)
		{
			//Find the triangle.  Pick the edge that will make it the fattest.
			distance2 = p.distance2d(uOld->dPrev()->org2d());	//get prev cw edge around dest
			distance3 = p.distance2d(uOld->dNext()->org2d());	//get next ccw edge around dest

			if(distance2 < distance3)
				pn = uOld->dPrev()->org2d();
			else
				pn = uOld->dNext()->org2d();
			return uOld;		//we've checked all our edges on the pivot
		}
		if(distance < distanceOld)
		{
			distanceOld = distance;
			uOld = u;
		}
	}
}

//Finds the nearest destination (AKA shortest edge) of a given edge by pivoting its origin
//Used to replace invalid boundary points for extrapolation when a point falls outside the convex hull
Edge* SHullDelaunay::nearestDest (Edge* e)
{
	double distance;
	double distanceOld;
	Edge* uOld;
	Edge* u = e;
	distance = u->org2d().distance2d(u->dest2d());	//get length of current edge e
	distanceOld = distance;

	while(TRUE)
	{
		//get next ccw edge around origin
		u = u->oNext();
		distance = u->org2d().distance2d(u->dest2d());
		if(e == u)
			return uOld;			//we've checked all our edges on the pivot
		if(distance < distanceOld)
		{
			distanceOld = distance;
			uOld = u;
		}
	}
}

int SHullDelaunay::checkEdge( Edge* e)
{
	Edge* t = e->oPrev();

	// Examine edges to ensure that the Delaunay condition is satisfied
	if(rightOf(t->dest2d(), e) && inCircle(e->org2d(), t->dest2d(), e->dest2d(), e->oNext()->dest2d())==1)
	{
		swap(e);
		e = e->oPrev();
		return 1;
	}
	return 0;
}

bool SHullDelaunay::pointOnSurface(const Point& p)
{
	Edge* temp = hullEdge;
	while(!rightOf(p, hullEdge))
	{
		hullEdge = hullEdge->sym()->oNext();
		if(temp == hullEdge) return true;
	}
	return false;
}

//Original Locate to find containing triangle
Edge* SHullDelaunay::locate(const Point& p)
{
	Edge* e = startingEdge;
	int i = 0;
	while (TRUE)
	{
		if (e->org2d().eq2d(p) || e->dest2d().eq2d(p))
			return e;
		else if (rightOf(p, e))
			e = e->sym();
		else if (!rightOf(p, e->oNext()))
			e = e->oNext();
		else if (!rightOf(p, e->dPrev()))
			e = e->dPrev();
		else
			return e;

		// This is bad, but I need to solve some floating point errors before I can get rid of it.
		if(i++ > 2000000)
			return e;
	}
}

//This function is used by locateNearestNeighbor
//Start by finding the containing triangle
//because it the quickest path  to get to the nearest neighbors
vector<Edge*> SHullDelaunay::locate3NN(const Point& p)
{
	Edge* e3=NULL;
	vector< Edge * >::iterator ite;
	Edge* e = startingEdge;
	Edge* u=NULL;
	Edge* oNextLoop=startingEdge;
	double distance=0;
	double i = 0;

	while (TRUE)
	{
		if (e->org2d().eq2d(p) || e->dest2d().eq2d(p))
		{
			//If p is a vertex on our starting edge
			e = e->sym();
			vector<Edge*> edges = findNearest3(p, e, e3, distance, oNextLoop, u);
			return edges;
		}
		else if (rightOf(p, e))
		{
			e = e->sym();
			oNextLoop=e;
		}

		else if (!rightOf(p, e->oNext()))
		{
			e = e->oNext();
			//Find nn when infinitely cyclic and cannot find a containing tri
			//This happens when outside boundary pts are not added
			//and p is outside the convex hull.
			if (e == oNextLoop)
			{
				vector<Edge*> edges = findNearest3(p, e, e3, distance, oNextLoop, u);
				return edges;
			}
		}
		else if (!rightOf(p, e->dPrev()))
		{
			e = e->dPrev();
			oNextLoop=e;
		}
		else
		{
			vector<Edge*> edges = findNearest3(p, e, e3, distance, oNextLoop, u);
			return edges;
		}

		if(i++ > 99999)
		{
			vector<Edge*> edges = findNearest3(p, e, e3, distance, oNextLoop, u);
			return edges;
		}
	}
}

//Added to find the 3 Nearest Neighbors by cycling edges of a point and following the path of the closest edge
//This function is used by locate3NN()
//Beginning with edge e, pivot the org to cycle through all edges until edge e3 reached
vector<Edge*> SHullDelaunay::findNearest3(const Point& p, Edge* e, Edge* e3, double distance, Edge* oNextLoop, Edge* u)
{
	double distance2;
	double distance3;
	double uDistance;
	Edge* e2;

	//Loop through all edges
	while(TRUE)
	{
		//Determine which end of the current edge is closer to the point.
		//Initialize values, then compare which is closer.
		u = e;
		e3 = u->oNext();
		distance3 = p.distance2d(u->oNext()->dest2d());	//Distance from point to next edge about origin destination
		distance  = p.distance2d(u->dest2d()); //Distance from point to current edge destination
		distance2 = p.distance2d(u->org2d()); //Distance from point to current edge origin
		oNextLoop = u; //Next edge about current edge origin to check
		e2 = u->sym(); //Current edge flipped

		//If the origin of the current edge is closer than the destination, then swap values.
		//We want the closest to be the destination of the edge.  Now we can pivot the rest of the edges
		//around the origin to compare with our closest to make sure it is the closest out of all edges around the origin.
		//If we find a closer destination than our current we will make this our new edge and repeat the process.
		//In doing this, our new origin will be the previous closest point found.
		if(distance2 < distance)
		{
			double distancet = distance2;
			distance2 = distance;
			distance  = distancet;
			e2 = e;
			e  = e->sym();
			u  = e;
			oNextLoop = u;
		}
		while(TRUE)
		{
			//We will then pivot edges extending from this origin
			//looking for any destinations that are closer to our point.
			while(u->oNext() != oNextLoop) //Cycle until we have pivoted all edges and are back to our beginning edge
			{
				//Get our next edge pivoted about origin
				u = u->oNext();
				uDistance = p.distance2d(u->dest2d());

				//If this edge's destination is closer, replace as our closest.
				if(uDistance < distance)
				{ //is closer than our closest
					e3 = e2;
					e2 = e;
					e  = u;
					distance3 = distance2;
					distance2 = distance;
					distance  = uDistance;
				}
				else if(uDistance < distance2 )
				{ //is closer than our 2nd closest
					e3 = e2;
					distance3 = distance2;
					e2 = u;
					distance2 = uDistance;
				}
				else if((uDistance == distance2) & (e2->org2d().eq2d(e->org2d())))
				{	//Isosceles triangle formed where the 2nd closest is the destination
					//Therefore, haven't found the 2nd closest in our current cycle
					//Therefore, the 2nd closest is saved on a edge from the previous origin
					//and needs to be updated to the edge from the current.
					e2 = u;
					distance2 = uDistance;
				}
				else if(uDistance <= distance3 )
				{ //is closer or as close as our 3rd closest
					e3 = u;
					distance3 = uDistance;
				}
			}
			if( oNextLoop->sym() == e)
			{	//We have pivoted all edges about the origin and didn't find any destinations closer than our closest.
				//Therefore, we have the closest.
				//Now we pivot our 2nd closest to check to see that we found the 3rd closest
				u = e2->sym();
				oNextLoop = e2->sym();
				while(u->oNext() != oNextLoop)
				{
					u = u->oNext();
					uDistance = p.distance2d(u->dest2d());
					if((uDistance < distance3) & (u != e))//*
					{ //is closer than 3rd closest and not the closest
						e3 = u;
						distance3 = uDistance;
					}
				}
				if(e->org2d().eq2d(e2->dest2d()))
				{//if the closest and 2nd closest are the same edge
					vector<Edge*> edges;
					edges.push_back(e);
					edges.push_back(e3);
					return edges; // e;
				}
				vector<Edge*> edges;
				edges.push_back(e2);
				edges.push_back(e3);
				return edges; // e2;
			}
			//Get next origin to pivot
			u=e->sym();
			oNextLoop=e->sym();
		}
	}
}

Mesh* SHullDelaunay::getMesh()
{
	Mesh *m = new Mesh();//where is this getting deleted?
	m->insertPoints(*pL);

	Triangle t;
	Point p1, p2, p3;

	list<Triangle> tri;
	list<Triangle>::iterator it;
	//list< Edge * >::iterator ite;// SJZ
	vector< Edge * >::iterator ite;// SJZ
	for(ite = edges.begin(); ite != edges.end(); ite++)
	{
		if(!rightOf((*ite)->oNext()->dest2d(), (*ite)))
		{
			p1 = (*ite)->org2d();
			p2 = (*ite)->dest2d();
			p3 = (*ite)->oNext()->dest2d();
			t.set(p1, p2, p3);
			//t.set(p1, p3, p2);
			tri.push_back(t);
		}

		if(rightOf((*ite)->oPrev()->dest2d(), (*ite)))
		{
			p1 = (*ite)->org2d();
			p2 = (*ite)->oPrev()->dest2d();
			p3 = (*ite)->dest2d();
			t.set(p1, p2, p3);
			//t.set(p1, p3, p2);
			tri.push_back(t);
		}
	}

	tri.sort();
	tri.unique();

	for(it = tri.begin(); it != tri.end(); it++)
		m->insertIndices(*it);

	return m;
}

double SHullDelaunay::sample(const Point& p)
{
	Edge* e;
	if(pointOnSurface(p))
		e = locate(p);
	else
		return 0;

	if ((e->org2d()).eq2d(p))
		return (e->org2d()).z;
	if ((e->dest2d().eq2d(p))) // point is already in
		return (e->dest2d()).z;

	Point p1, p2, p3;
	p1 = e->org2d();
	p2 = e->dest2d();
	p3 = e->oNext()->dest2d();

	return samplePlane(p, p1, p2, p3);
}

