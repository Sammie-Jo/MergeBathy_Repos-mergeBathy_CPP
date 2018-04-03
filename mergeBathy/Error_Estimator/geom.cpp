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
* Written in 2015 by Eric Mixon while employed by the U.S.Naval Research Laboratory.
* To the extent possible under law, the author(s) and the U.S.Naval Research Laboratory have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide.This software is distributed without any warranty.
* You should have received a copy of the CC0 Public Domain Dedication along with this software.If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
**********************************************************************/
//	Author: Eric Mixon		Date:7/14/10
//	Description:
//	- Library of geometric objects and functions
//	- Edge, Ping, Vector2d, Line, Triangle, QuadEdge
//	- Also some other math oriented functions
//    for various algorithms such as splice for
//	  the incremental Delaunay algorithm

// Prevent redefinition of Bathy geom classes and functions
//#pragma once

// Some c++ libraries
#include<math.h>
#include"geom.h"
#include<list>
#include<iostream>
#include<iomanip>
#include<assert.h>
#include<vector>
#include "../standardOperations.h"
//#include"pointList.h"

// Standard namespace
using namespace std;

#pragma region --Vector2d
/*
* Class Vector2d
*
* Member Variables:
* x, y - positions x and y
*
* Member Functions:
* 	Constructors:
*	 Vector2d() - sets x and y to 0
*	 Vector2d(double a, double b) sets x and y to a and b
*
*	norm() returns sqrt(x^2 + y^2)
*	normalize() - normalizes the vector
*	operator+(Vector2d) - performs vector addition between two vectors
*		returning the resulting vector
*	operator-(Vector2d) - performs vector subtraction between two vectors
*		returning the resulting vector
*	operator*(double, Vector2d) - performs multiplication of a double and a vector
*		returning the resulting vector
*	dot(vector2d, Vector2d) - returns the dot product of two 2d vectors which
*		is a double
*/
Vector2d::Vector2d() { x = 0; y = 0; }

Vector2d::Vector2d(const double& a, const double& b) { x = a; y = b; }

double Vector2d::norm() const { return sqrt(x * x + y * y); }

void Vector2d::normalize()
{
	double len;
	if (abs((len = sqrt(x*x + y*y)) <=.0000001))//APL
//	if ((len = sqrt(x * x + y * y)) == 0.0)
		cerr << "Vector2d::normalize: Division by 0\n";
	else
	{
		x /= len;
		y /= len;
	}
}

Vector2d Vector2d::operator+(const Vector2d& v) const
{
	return Vector2d(x + v.x, y + v.y);
}

Vector2d Vector2d::operator-(const Vector2d& v) const
{
	return Vector2d(x - v.x, y - v.y);
}

inline Vector2d operator*(const double& c, const Vector2d& v)
{
	return Vector2d(c * v.x, c * v.y);
}

inline double dot(const Vector2d& u, const Vector2d& v)
{
	return u.x * v.x + u.y * v.y;
}
#pragma endregion

#pragma region -- Ping Class

/*
* Class Ping
*/

Ping::Ping(double lon, double lat, double dep) : longitude(lon), latitude(lat), depth(dep) {}

Ping::Ping( const Ping& p )
{
	longitude = p.getLong();
	latitude = p.getLat();
	depth = p.getDepth();
	uncertainty = p.getUncert();
	horizontalUncertainty = p.getHorizontalUncert();
	verticalUncertainty = p.getVerticalUncert();
}

// setLong takes in and sets the longitude
void Ping::setLong(const double& lon) { longitude = lon; }

// setLat takes in and sets the latitude
void Ping::setLat(const double& lat) { latitude = lat; }

// setDepth takes in and sets the depth
void Ping::setDepth(const double& dep) { depth = dep; }

// set lat long and depth at once
void Ping::set(const double& lon, const double& lat, const double& dep)
{
	longitude = lon;
	latitude = lat;
	depth = dep;
}

// Set the horizontal uncertainty
void Ping::setHorizontalUncert( double h ){ horizontalUncertainty = h; }

// Get the horizontal uncertainty
double Ping::getHorizontalUncert() const{ return horizontalUncertainty; }

// Set the vertical uncertainty
void Ping::setVerticalUncert( double v ){ verticalUncertainty = v; }

// Get the vertical uncertainty
double Ping::getVerticalUncert() const{ return verticalUncertainty; }

// Set the uncertainty
void Ping::setUncert( double u ){ uncertainty = u; }

// Get the uncertainty
double Ping::getUncert() const{ return uncertainty; }

// getLong returns the current longitude
double Ping::getLong() const { return longitude; }

// getLat returns the current latitude
double Ping::getLat() const { return latitude; }

// getDepth returns the current depth
double Ping::getDepth() const { return depth; }

// Overloads the = operator for assignment operations
// allows for the Ping1 = Ping2;
Ping& Ping::operator=(const Ping& p)
{
	longitude = p.getLong();
	latitude = p.getLat();
	depth = p.getDepth();
	uncertainty = p.getUncert();
	horizontalUncertainty = p.getHorizontalUncert();
	verticalUncertainty = p.getVerticalUncert();
	return *this;  // REO Inserted this.
}

// Overloads the == operator for Boolean operations
// if p1 = p2 then it returns true
// otherwise returns false
bool Ping::operator==(const Ping& p) const
{
	if( longitude == p.getLong() &&
		latitude == p.getLat() &&
		depth == p.getDepth() )
	{ return true; }
	return false;
}

// Overloads the != operator for Boolean operations
bool Ping::operator!=(const Ping& p) const
{
	if( longitude != p.getLong() ||
		latitude != p.getLat() ||
		depth != p.getDepth() )
	{	return true; }
	return false;
}

// Overloads the < operator for use with list sort
// function in order to make lists unique
bool Ping::operator<(const Ping& p) const
{
	if(p.getLong() < longitude)
		return false;
	else if(p.getLong() > longitude)
		return true;

	if(p.getLat() < latitude)
		return false;
	else if(p.getLat() > latitude)
		return true;

	if(p.getDepth() < depth)
		return false;
	return true;
}

void Ping::print() const
{
	cout << setprecision(15) << "Long: " << longitude << " Lat: " << latitude << " Depth: " << depth << endl;
}

/*
	End Class Ping
*/
#pragma endregion Ping Class

#pragma region -- Point Class
/*
	Class Point
*/
Point::Point(const double& x, const double& y, const double& z) : x(x), y(y), z(z) {}

Point::Point(const Point& p)
{
	x = p.x;
	y = p.y;
	z = p.z;
	u = p.u;
	hU = p.hU;
	vU = p.vU;
	id=p.id;
	//Store additional uncertainty calculations
	u2=p.u2;
	u3=p.u3;
	u4=p.u4;
	u5=p.u5;

	//raster******************************************
	//e = p.e;
	nei = p.nei;
	rei = p.rei;
	z0 = p.z0;
	e0 = p.e0;
	zK = p.zK;
	eK = p.eK;
}

Point Point::operator=(const Point& p)
{
	x = p.x;
	y = p.y;
	z = p.z;
	u = p.u;
	hU = p.hU;
	vU = p.vU;
	id=p.id;
	//Store additional uncertainty computations
	u2=p.u2;
	u3=p.u3;
	u4=p.u4;
	u5=p.u5;

	//raster******************************************
	//e = p.e;
	nei = p.nei;
	rei = p.rei;
	z0 = p.z0;
	e0 = p.e0;
	zK = p.zK;
	eK = p.eK;

	return p;
}

// Overloads the == operator for Boolean operations
// if p1 = p2 then it returns true
// otherwise returns false
bool Point::operator==(const Point& p) const
{
	if( x == p.x &&
		y == p.y &&
		z == p.z )
	{ return true; }
	return false;
}

// Overloads the != operator for Boolean operations
bool Point::operator!=(const Point& p) const
{
	if( x != p.x ||
		y != p.y ||
		z != p.z )
	{	return true; }
	return false;
}

// Overloads the < operator for use with list sort
// function in order to make lists unique
bool Point::operator<(const Point& p) const
{
	if(p.x < x)
		return false;
	else if(p.x > x)
		return true;

	if(p.y < y)
		return false;
	else if(p.y > y)
		return true;

	if(p.z < z)
		return false;
	return true;
}

double Point::distance2d(const Point& p) const
{
	return sqrt((((x - p.x) * (x - p.x)) +
		((y - p.y) * (y - p.y))));
}

Vector2d Point::operator-(const Point& p) const
{
	return Vector2d(x - p.x, y - p.y);
}

bool Point::eq2d(const Point& p) const
{
	if(x == p.x && y == p.y)
		return true;
	return false;
}

bool Point::eqErr(const Point& p, const double err) const
{
	if(fabs(x - p.x) < err && fabs(y - p.y) < err)
		return true;
	return false;
}

void Point::print() const
{
	cout << setprecision(15) << "X: " << x << " Y: " << y << " Z: " << z << endl;
}
/*
	End Class Point
*/
#pragma endregion Point Class

#pragma region -- Triangle Class
/*
	Class Triangle
*/
// Construct a triangle with 3 Pings
Triangle::Triangle(const Point& V1, const Point& V2, const Point& V3) : v1(V1), v2(V2), v3(V3) {}

// Construct a triangle at all zeros
Triangle::Triangle() {}

// Returns the midpoint of the triangle's hypotenuse (Used with rTin)
Point Triangle::midHyp() const
{
	return Point(((v2.x + v3.x) / 2), ((v2.y + v3.y) / 2), ((v2.z + v3.z) / 2));
}

double Triangle::edgeLen() const
{
	return sqrt((v1.x - v2.x) * (v1.x - v2.x) + (v1.y - v2.y) * (v1.y - v2.y));
}

// Sets all pings in triangle
void Triangle::set(Point p1, Point p2, Point p3)
{
	v1 = p1;
	v2 = p2;
	v3 = p3;
}

Triangle Triangle::operator=(const Triangle& t)
{
	v1 = t.v1;
	v2 = t.v2;
	v3 = t.v3;
	return t;
}

// Returns true if point is inside triangle
// ignores depth
// might have some error with
// 	points that lie really close to an edge
bool Triangle::isInside(const Point& p) const
{
	return pointInTriangle(p, v1, v2, v3);
}

// Overloads < operator for sorting, simply adds values of all vertices
// and finds the triangle with the lowest total
bool Triangle::operator<(const Triangle& t) const
{
	double selfAvgLong = (v1.x + v2.x + v3.x) / 3;
	double selfAvgLat = (v1.y + v2.y + v3.y) / 3;
	double selfAvgDepth = (v1.z + v2.z + v3.z) / 3;

	double tAvgLong = (t.v1.x + t.v2.x + t.v3.x) / 3;
	double tAvgLat = (t.v1.y + t.v2.y + t.v3.y) / 3;
	double tAvgDepth = (t.v1.z + t.v2.z + t.v3.z) / 3;

	if(tAvgLong < selfAvgLong)
		return false;
	else if(tAvgLong > selfAvgLong)
		return true;

	if(tAvgLat < selfAvgLat)
		return false;
	else if(tAvgLat > selfAvgLat)
		return true;

	if(tAvgDepth < selfAvgDepth)
		return false;
	return true;
}

// Overloads == operator for Boolean comparison
bool Triangle::operator==(const Triangle& t) const
{
	if(v1 == t.v1)
	{
		if(v2 == t.v2 && v3 == t.v3)
			return true;
	}
	else if(v1 == t.v2)
	{
		if(v2 == t.v3 && v3 == t.v1)
			return true;
	}
	else if(v1 == t.v3)
	{
		if(v2 == t.v1 && v3 == t.v2)
			return true;
	}
		return false;
}

double Triangle::sample(const Point& p) const
{
	return samplePlane(p, v1, v2, v3);
}

Gradient Triangle::getGradient() const
{
	return gradient(v1, v2, v3);
}

#pragma endregion Traingle Class

#pragma region -- Gradient Class
//This was original function for determining the gradient of
//a given triangle.  This was replaced this with gradientGrid(),
//a version of Paul's matlab gradientm_manual_v4.m for finding all gradients
//over a grid.
Gradient gradient(const Point& v1, const Point& v2, const Point& v3)
{
	double A = (v1.y * (v2.z - v3.z)) +
				(v2.y * (v3.z - v1.z)) +
				(v3.y * (v1.z - v2.z));
	double B = (v1.z * (v2.x - v3.x)) +
				(v2.z * (v3.x - v1.x)) +
				(v3.z * (v1.x - v2.x));
	double C = (v1.x * (v2.y - v3.y)) +
				(v2.x * (v3.y - v1.y)) +
				(v3.x * (v1.y - v2.y));
	double D = -1 * ((v1.x * ((v2.y * v3.z) -
				(v3.y * v2.z))) +
				(v2.x * ((v3.y * v1.z) -
				(v1.y * v3.z))) +
				(v3.x * ((v1.y * v2.z) -
				(v2.y * v1.z))));

	Gradient g(v1.x, v1.y, v1.z, v2.x, v2.y, v2.z,
				v3.x, v3.y, v3.z, 0, 0, 0, 0);

	if(C == 0)
		return g;

	double z1, z3, z7;
	z1 = -1 * (D / C);
	z3 = -1 * ((A + D) / C);
	z7 = -1 * ((B + D) / C);
	g.FX = (z3 - z1) / 2;
	g.FY = (z7 - z1) / 2;
	g.S = atan(sqrt(g.FX*g.FX + g.FY*g.FY));
	g.A = PI - atan(g.FY / g.FX) + PI / 2 * (g.FX / fabs(g.FX));

	return g;
}
#pragma endregion Gradient Class

#pragma region --Edge Class

Edge::Edge() { data = 0; }

// Return the dual of the current edge, directed from its right to its left
Edge* Edge::rot() { return (num < 3) ? this + 1 : this - 3; }

// Return the dual of the current edge, directed from its left to its right
Edge* Edge::invRot() { return (num > 0) ? this - 1 : this + 3; }

// Return the edge from the destination to the origin of the current edge
Edge* Edge::sym() { return (num < 2) ? this + 2: this - 2; }

// Return the next ccw edge around (from) the origin of the current edge
Edge* Edge::oNext() { return next; }

// Return the next cw edge around (from) the origin of the current edge
Edge* Edge::oPrev() { return rot()->oNext()->rot(); }

// Return the next ccw edge around (into) the destination of the current edge
Edge* Edge::dNext() { return sym()->oNext()->sym(); }

// Return the next cw edge around (into) the destination of the current edge
Edge* Edge::dPrev() { return invRot()->oNext()->invRot(); }

// Return the ccw edge around the left face following the current edge
Edge* Edge::lNext() { return invRot()->oNext()->rot(); }

// Return the ccw edge around the left face before the current edge
Edge* Edge::lPrev() { return oNext()->sym(); }

// Return the edeg around the right face ccw following the current edge
Edge* Edge::rNext() { return rot()->oNext()->invRot(); }

// Return the edge around the right face ccw before the current edge
Edge* Edge::rPrev() { return sym()->oNext(); }

// Access Pointers //
Point* Edge::org() { return data; }
Point* Edge::dest() { return sym()->data; }
const Point& Edge::org2d() const { return *data; }
const Point& Edge::dest2d() const { return (num < 2) ? *((this+2)->data) : *((this - 2)->data); }
void Edge::endPoints(Point* org, Point* dest)
{
	data = org;
	sym()->data = dest;
}
#pragma endregion Edge Class

#pragma region --QuadEdge Class

QuadEdge* Edge::qEdge(){ return (QuadEdge *)(this - num); }

QuadEdge::QuadEdge(int id)
{
	e[0].num = 0, e[1].num = 1, e[2].num = 2, e[3].num = 3;
	e[0].next = &(e[0]); e[1].next = &(e[3]);
	e[2].next = &(e[2]); e[3].next = &(e[1]);
	// My edit for keeping track with edges
	e[0].id = id;
	e[1].id = id;
	e[2].id = id;
	e[3].id = id;
	// End my edit
}

QuadEdge::QuadEdge()
{
	e[0].num = 0, e[1].num = 1, e[2].num = 2, e[3].num = 3;
	e[0].next = &(e[0]); e[1].next = &(e[3]);
	e[2].next = &(e[2]); e[3].next = &(e[1]);
}
#pragma endregion QuadEdge Class

#pragma region --Line Class

Line::Line(){}

Line::Line(const Point& p, const Point& q)
{
	Vector2d t = q - p;
	double len = t.norm();
	a =   t.y / len;
	b = - t.x / len;
	c = -(a*p.x + b*p.y);
}

double Line::eval(const Point& p) const
{
	return (a * p.x + b * p.y + c);
}

int Line::classify(const Point& p) const
{
	double d = eval(p);
	return (d < -EPS) ? -1 : (d > EPS ? 1 : 0);
}
#pragma endregion Line Class

#pragma region --Delaunay Manipulations
void splice(Edge* a, Edge* b)
// This operator affects the two edge rings around the origins of a and b,
// and, independently, the two edge rings around the left faces of a and b.
// In each case, (i) if the two rings are distinct, Splice will combine
// them into one; (ii) if the two are the same ring, Splice will break it
// into two separate pieces.
// Thus, Splice can be used both to attach the two edges together, and
// to break them apart. See Guibas and Stolfi (1985) p.96 for more details
// and illustrations.
{
	Edge* alpha = a->oNext()->rot();
	Edge* beta = b->oNext()->rot();
	Edge* t1 = b->oNext();
	Edge* t2 = a->oNext();
	Edge* t3 = beta->oNext();
	Edge* t4 = alpha->oNext();

	a->next = t1;
	b->next = t2;
	alpha->next = t3;
	beta->next = t4;
}

void deleteEdge(Edge* e)
{
	splice(e, e->oPrev());
	splice(e->sym(), e->sym()->oPrev());
	delete e->qEdge();
}

Edge* makeEdge(const int& id)
{
	QuadEdge *ql = new QuadEdge(id);
	return ql->e;
}

Edge* makeEdge()
{
	QuadEdge *ql = new QuadEdge();
	return ql->e;
}

Edge* connect(Edge* a, Edge* b, const int& id)
// Add a new edge e connecting the destination of a to the
// origin of b, in such a way that all three have the same
// left face after the connection is complete.
// Additionally, the data pointers of the new edge are set.
{
	Edge* e = makeEdge(id);
	splice(e, a->lNext());
	splice(e->sym(), b);
	e->endPoints(a->dest(), b->org());
	return e;
}

Edge* connect(Edge* a, Edge* b)
// Add a new edge e connecting the destination of a to the
// origin of b, in such a way that all three have the same
// left face after the connection is complete.
// Additionally, the data pointers of the new edge are set.
{
	Edge* e = makeEdge();
	splice(e, a->lNext());
	splice(e->sym(), b);
	e->endPoints(a->dest(), b->org());
	return e;
}

void swap(Edge* e)
// Essentially turns edge e counterclockwise inside its enclosing
// quadrilateral. The data pointers are modified accordingly.
{
	Edge* a = e->oPrev();
	Edge* b = e->sym()->oPrev();
	splice(e, a);
	splice(e->sym(), b);
	splice(e, a->lNext());
	splice(e->sym(), b->lNext());
	e->endPoints(a->dest(), b->dest());
}
#pragma endregion Delaunay Manipulations

#pragma region --Delaunay Geo Predicates
/*************** Geometric Predicates for Delaunay Diagrams *****************/
double triArea(const Point& a, const Point& b, const Point& c)
// Returns twice the area of the oriented triangle (a, b, c), i.e., the
// area is positive if the triangle is oriented counterclockwise.
{//SJZ added
	double val = (b.x - a.x)*(c.y - a.y) -
				 (b.y - a.y)*(c.x - a.x);

	return (b.x - a.x)*(c.y - a.y) -
			(b.y - a.y)*(c.x - a.x);
}

// Returns TRUE if the point d is inside the circle defined by the
// points a, b, c. See Guibas and Stolfi (1985) p.107.
int inCircle(const Point& a, const Point& b, const Point& c, const Point& d)
{
	double val2= (a.x-d.x)*(b.y-d.y)*(((c.x*c.x)-(d.x*d.x))+((c.y*c.y)-(d.y*d.y)))+
	(a.y-d.y)*(c.x-d.x)*(((b.x*b.x)-(d.x*d.x))+((b.y*b.y)-(d.y*d.y)))+
	(b.x-d.x)*(c.y-d.y)*(((a.x*a.x)-(d.x*d.x))+((a.y*a.y)-(d.y*d.y)))-
	(c.x-d.x)*(b.y-d.y)*(((a.x*a.x)-(d.x*d.x))+((a.y*a.y)-(d.y*d.y)))-
	(b.x-d.x)*(a.y-d.y)*(((c.x*c.x)-(d.x*d.x))+((c.y*c.y)-(d.y*d.y)))-
	(a.x-d.x)*(c.y-d.y)*(((b.x*b.x)-(d.x*d.x))+((b.y*b.y)-(d.y*d.y)));

	double val = (a.x*a.x + a.y*a.y) * triArea(b, c, d) -
				(b.x*b.x + b.y*b.y) * triArea(a, c, d) +
				(c.x*c.x + c.y*c.y) * triArea(a, b, d) -
				(d.x*d.x + d.y*d.y) * triArea(a, b, c);

	if (val==0)
		return -1;
	return (a.x*a.x + a.y*a.y) * triArea(b, c, d) -
	(b.x*b.x + b.y*b.y) * triArea(a, c, d) +
	(c.x*c.x + c.y*c.y) * triArea(a, b, d) -
	(d.x*d.x + d.y*d.y) * triArea(a, b, c) > 0;
}
int isConvexQuad(const Point& a, const Point& b, const Point& c, const Point& d)
{
	if ((ccw(a,b,c) && ccw(b,c,d)) && (ccw(c,d,a) && ccw(d,a,b)))
		return 1;
	//else if (triArea(a,b,c)==0 ||triArea(b,c,d)==0||triArea(c,d,a)==0||triArea(d,a,b))
	else if (triArea(b,c,d)==0||triArea(c,d,a)==0||triArea(d,a,b))
		return 2;
	else return 0;
}
// Returns TRUE if the points a, b, c are in a counterclockwise order
int ccw(const Point& a, const Point& b, const Point& c)
{
	return (triArea(a, b, c) > 0);
}

int rightOf(const Point& x, Edge* e)
{
	return ccw(x, e->dest2d(), e->org2d());
}

int LeftOf(const Point& x, Edge* e)
{
	return ccw(x, e->org2d(), e->dest2d());
}

int onEdge(const Point& x, Edge* e)
// A predicate that determines if the point x is on the edge e.
// The point is considered on if it is in the EPS-neighborhood
// of the edge.
{
	double t1, t2, t3;
	t1 = (x - e->org2d()).norm();
	t2 = (x - e->dest2d()).norm();
	if (t1 < EPS || t2 < EPS)
		return TRUE;
	t3 = (e->org2d() - e->dest2d()).norm();
	if (t1 > t3 || t2 > t3)
		return FALSE;
	Line line(e->org2d(), e->dest2d());
	return (fabs(line.eval(x)) < EPS);
}

/*
*	Returns the minimum angle difference in the two
*		planes formed by this triangle and the triangle t1
*	arg0 - Triangle t1 - Test triangle
*/
double differenceInPlanes(const Triangle& t1, const Triangle& t2)
{
	double A1 = (t1.v1.y * (t1.v2.z - t1.v3.z)) +
			(t1.v2.y * (t1.v3.z - t1.v1.z)) + (t1.v3.y * (t1.v1.z - t1.v2.z));
	double B1 = (t1.v1.z * (t1.v2.x - t1.v3.x)) +
			(t1.v2.z * (t1.v3.x - t1.v1.x)) + (t1.v3.z * (t1.v1.x - t1.v2.x));
	double C1 = (t1.v1.x * (t1.v2.y - t1.v3.y)) +
			(t1.v2.x * (t1.v3.y - t1.v1.y)) + (t1.v3.x * (t1.v1.y - t1.v2.y));
	double A2 = (t2.v1.y * (t2.v2.z - t2.v3.z)) +
			(t2.v2.y * (t2.v3.z - t2.v1.z)) + (t2.v3.y * (t2.v1.z - t2.v2.z));
	double B2 = (t2.v1.z * (t2.v2.x - t2.v3.x)) +
			(t2.v2.z * (t2.v3.x - t2.v1.x)) + (t2.v3.z * (t2.v1.x - t2.v2.x));
	double C2 = (t2.v1.x * (t2.v2.y - t2.v3.y)) +
			(t2.v2.x * (t2.v3.y - t2.v1.y)) + (t2.v3.x * (t2.v1.y - t2.v2.y));
	double num = fabs(A1 * A2 + B1 * B2 + C1 * C2);
	double den = sqrt(pow(A1, 2) + pow(B1, 2) + pow(C1, 2)) *
							sqrt(pow(A2, 2) + pow(B2, 2) + pow(C2, 2));
	double x = num / den;
	double theta = acos(x);
	return theta * 180 / 3.14159265;
}

// Returns positive number if the point P1 is on one side of the line
// and returns a negative number if it is on the other
double sign(const Point& p1, const Point& p2, const Point& p3)
{
	return (p3.x - p2.x) * (p1.y - p2.y) - (p3.y - p2.y) * (p1.x - p2.x);
}

bool pointInTriangle(const Point& p, const Point& v1, const Point& v2, const Point& v3)
{
	bool b1, b2, b3;
	b1 = sign(p, v1, v2) < 0.0f;
   	b2 = sign(p, v2, v3) < 0.0f;
	b3 = sign(p, v3, v1) < 0.0f;
    return ((b1 == b2) && (b2 == b3));
}

double samplePlane(const Point& p, const Point& v1, const Point& v2, const Point& v3)
{
	double A = (v1.y * (v2.z - v3.z)) +
				(v2.y * (v3.z - v1.z)) +
				(v3.y * (v1.z - v2.z));
	double B = (v1.z * (v2.x - v3.x)) +
				(v2.z * (v3.x - v1.x)) +
				(v3.z * (v1.x - v2.x));
	double C = (v1.x * (v2.y - v3.y)) +
				(v2.x * (v3.y - v1.y)) +
				(v3.x * (v1.y - v2.y));
	double D = -1 * ((v1.x * ((v2.y * v3.z) -
				(v3.y * v2.z))) +
				(v2.x * ((v3.y * v1.z) -
				(v1.y * v3.z))) +
				(v3.x * ((v1.y * v2.z) -
				(v2.y * v1.z))));

	if(C == 0)
		return (v1.z + v2.z + v3.z) / 3;

	return -1 * (((A * p.x) + (B * p.y) + D) / C);
}

Point centerOfCircumCircle(const Point& p1, const Point& p2, const Point& p3)
{
	double d = (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
	double x = (((p1.x - p3.x) * (p1.x + p3.x) + (p1.y - p3.y) * (p1.y + p3.y)) / 2 * (p2.y - p3.y)
			   - ((p2.x - p3.x) * (p2.x + p3.x) + (p2.y - p3.y) * (p2.y + p3.y)) / 2 * (p1.y - p3.y)) / d;

	double y = (((p2.x - p3.x) * (p2.x + p3.x) + (p2.y - p3.y) * (p2.y + p3.y)) / 2 * (p1.x - p3.x)
			- ((p1.x - p3.x) * (p1.x + p3.x) + (p1.y - p3.y) * (p1.y + p3.y)) / 2 * (p2.x - p3.x)) / d;
	Point p(x, y, 0);
	return p;
}

double radiusOfCircumCircle(const Point& p1, const Point& p2, const Point& p3)
{
	double a = CFPE(p1.distance2d(p2));
	double b = CFPE(p2.distance2d(p3));
	double c = CFPE(p3.distance2d(p1));
	double num = CFPE(a * b * c);
	double den = CFPE(sqrt((a + b + c) * (b + c - a ) * (c + a - b) * (a + b - c)));
	double ans = CFPE(num / den);

	if(den == 0)
		return 99999999.0;
	if(den >= 0 && den < 0)
		return 99999999.0;
	//if(isinf(ans) || isnan(ans) || ans <= 0)
	if(ans <= 0)
		return 99999999.0;
	return ans;
}
#pragma endregion Delaunay Geo Predicates

#pragma region --SJZ Added
int collinear(const Point& x, Edge* e)
{
	return triArea(x, e->dest2d(), e->org2d()) == 0;
}
int collinear(const Point& x, Point& b, Point& c)
{
	return triArea(x, b, c) == 0;
}
double CFPE(const double number)
{
	int precision = 5;
	//default returns (10000 * number) / 10000
	//should correct very small floating point errors

	double correction = pow(10.0, precision);
	return roundDouble(correction * number) / correction;
}
/**
 * Tests if two numbers are almost equal.
 */
bool fuzzyEquals(const Point& x, const Point& e)
{
	double number1=x.x;
	double number2=x.y;
	double number3=e.x;
	double number4=e.y;

	int precision = 5;
	double difference13= number1 - number3;
	double difference24= number2 - number4;

	double range= pow(10.0, precision);

	//default check:
	//0.00001 < difference > -0.00001
	if ((difference13 < range && difference13 > -range) && (difference24 < range && difference24 > -range)) return true;
	else return false;
	//return difference < range && difference > -range;
}

#pragma endregion SJZ
