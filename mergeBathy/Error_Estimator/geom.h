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
#include <sstream>
#include <vector>
#include "../grid.h"
#include "../constants.h"

#ifndef ABS
#define ABS(a)((a) >= 0 ? (a) : -(a))
#endif
#ifndef MAX
#define MAX(a, b)((a) >= (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a, b)((a) <= (b) ? (a) : (b))
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif
#ifndef EPS
#define EPS 1e-6
#endif
//#ifndef PI
//#define PI 3.14159265
//#endif
/**
* A 2d Vector Class for performing vector arithmetic and normalization.
*/
class Vector2d {
	public:
		/**
		* The x value of the vector.
		*/
		double x;

		/**
		* The y value of the vector.
		*/
		double y;

		/**
		* A constructor for the vector.
		* Takes no parameters and sets x and y to 0.
		*/
		Vector2d();

		/**
		* A constructor for the vector.
		* Sets the values of x and y to a and b respectively.
		* @param a - Argument that sets the value of x.
		* @param b - Argument that sets the value of y.
		*/
		Vector2d(const double& a, const double& b);

		/**
		* A member that returns the length of the vector.
		* @return length of the vector.
		*/
		double norm() const;

		/**
		* A member that divides the x and y of the vector
		* by the length of the vector.
		*/
		void normalize();

		/**
		* A member that overloads the + operator and allows for
		* the addition of two Vector2d objects, uses standard
		* vector addition.
		* @param v - Argument to be added with "this" object.
		* @return The new resulting Vector2d object.
		*/
		Vector2d operator+(const Vector2d& v) const;

		/**
		* A member that overloads the - operator and allows for
		* the subtraction of two Vector2d objects, uses standard
		* vector subtraction.
		* @param v - Argument to be subtracted with "this" object
		* @return The new resulting Vector2d object.
		*/
		Vector2d operator-(const Vector2d& v) const;

		/**
		* A friend member that overloads the * operator and allows for
		* the cross multiplication of two Vector2d objects, uses standard
		* vector cross multiplication.
		* @param c - Argument multiplied times vector2d v.
		* @param v - Argument to be multiplied by c.
		* @return The new resulting Vector2d object.
		*/
		friend Vector2d operator*(const double&, const Vector2d&);

		/**
		* A friend member that allows for the dot multiplication
		* of two Vector2d objects, uses standard vector dot multiplication.
		* @param u - Argument to be multiplied times vector2d v.
		* @param v - Argument to be multiplied times vector2d u.
		* @return The new resulting Vector2d object.
		*/
		friend double dot(const Vector2d&, const Vector2d&);
};

/**
* Class for string lat long and depth of a ping. Also used for storing uncertainty
*/
class Ping
{
	private:
		/** Longitude of the Ping. */
		double longitude;

		/** Latitude of the Ping. */
		double latitude;

		/** Depth of the Ping. */
		double depth;

		/** Uncertainty of the Ping */
		double uncertainty;

		/** Horizontal Uncertainty of the Ping */
		double horizontalUncertainty;

		/** Vertical Uncertainty of the Ping */
		double verticalUncertainty;

	public:

		/** Constructor that takes in and sets the value of lat long and depth.
		* @param lon - Argument for the longitude with default value of 0.
		* @param lat - Argument for the latitude with default value of 0.
		* @param dep - Argument for the depth with default value of 0.
		*/
		Ping(double lon = 0, double lat = 0, double dep = 0);

		/** Sets the vertical uncertainty of the ping.
		* @param v - double being stored as the new uncertainty.
		*/
		void setVerticalUncert( double v );

		/** Returns the vertical uncertainty of the Ping.
		* @return Uncertainty of the Ping.
		*/
		double getVerticalUncert() const;

		/** Sets the horizontal uncertainty of the ping.
		* @param h - double being stored as the new uncertainty.
		*/
		void setHorizontalUncert( double h );

		/** Returns the horizontal uncertainty of the Ping.
		* @return Uncertainty of the Ping.
		*/
		double getHorizontalUncert() const;

		/** Sets the uncertainty of the ping.
		* @param u - double being stored as the new uncertainty.
		*/
		void setUncert( double u );

		/** Returns the uncertainty of the Ping.
		* @return Uncertainty of the Ping.
		*/
		double getUncert() const;

		/** Copy Constructor
		* @param p - Ping being copied
		*/
		Ping( const Ping& p );

		/** Takes in and sets the longitude.
		* @param lon - Argument for the new longitude.
		*/
		void setLong(const double& lon);

		/** Takes in and sets the latitude.
		* @param lat - Argument for the new latitude.
		*/
		void setLat(const double& lat);

		/** Takes in and sets the depth.
		* @param dep - Argument for the new depth.
		*/
		void setDepth(const double& dep);

		/** Set longitude, latitude, and depth all at once.
		* @param lon - Argument for new longtiude to be set.
		* @param lat - Argument for new latitude to be set.
		* @param dep - Argument for new depth to be set.
		*/
		void set(const double& lon, const double& lat, const double& dep);

		/** Returns the value of longitude.
		* @return Value of longtiude.
		*/
		double getLong() const;

		/** Returns the value of latitude.
		* @return Value of latitude.
		*/
		double getLat() const;

		/** Returns the value of depth.
		* @return Value of depth.
		*/
		double getDepth() const;

		/** Overloads the = operator for assignment operations, allowing
		* for the usage of Ping1 = Ping2.
		* @param p - Argument Ping whose values will be stored into "this" Ping.
		* @return Results of Ping1 = Ping2 will be that Ping1 will now contain
		* the values of Ping2.
		*/
		Ping& operator=(const Ping& p);

		/** Overloads the == operator for boolean comparison operations.
		* @return The result of Ping1 == Ping2 will return true if
		* the values of Ping1 and Ping2 are equivalent and will return
		* false otherwise.
		*/
		bool operator==(const Ping& p) const;

		/** Overloads the != operator for boolean comparison operations.
		* @return The result of Ping1 != Ping2 will return false if
		* the values of Ping1 and Ping2 are equivalent and will return
		* false otherwise.
		*/
		bool operator!=(const Ping& p) const;

		/** Overloads the < operator for returning the Ping that is less by
		* longitude -> latitude -> depth.  Used for sorting to easilty remove duplicates.
		* @return True if this Ping is < p, False otherwise
		*/
		bool operator<(const Ping& p) const;

		/** Prints the values of the Ping */
		void print() const;
};
/** Class for storing the relevant values of a Gradient and also
* provides a toString ability for easy export.
* Use the Triangle class to calculate the gradient and get a gradient object.
*/
class Gradient//SJZ
{
	public:
		/** X for the first Vertex in the represented Plane */
		double X1;

		/** X for the second Vertex in the represented Plane */
		double X2;

		/** X for the third Vertex in the represented Plane */
		double X3;

		/** Y for the first Vertex in the represented Plane */
		double Y1;

		/** Y for the second Vertex in the represented Plane */
		double Y2;

		/** Y for the third Vertex in the represented Plane */
		double Y3;

		/** Z for the first Vertex in the represented Plane */
		double Z1;

		/** Z for the second Vertex in the represented Plane */
		double Z2;

		/** Z for the third Vertex in the represented Plane */
		double Z3;

		/** Slope of the Gradient */
		double S;

		/** Aspect of the Gradient */
		double A;

		/** X direction Gradient */
		double FX;

		/** Y direction Gradient */
		double FY;

		/** Constructor for the Gradient Object
		* @param x1 - X for vertex 1 of the plane
		* @param y1 - y for vertex 1 of the plane
		* @param z1 - z for vertex 1 of the plane
		* @param x2 - x for vertex 2 of the plane
		* @param y2 - y for vertex 2 of the plane
		* @param z2 - z for vertex 2 of the plane
		* @param x3 - x for vertex 3 of the plane
		* @param y3 - y for vertex 3 of the plane
		* @param z3 - z for vertex 3 of the plane
		* @param s - slope of the gradient
		* @param a - aspect of the gradient
		* @param fx - x direction gradient
		* @param fy - y direction gradient
		*/
		Gradient(double x1, double y1, double z1, double x2,
			double y2, double z2, double x3, double y3,
			double z3, double s, double a, double fx, double fy) :
				X1(x1), Y1(y1), Z1(z1), X2(x2), Y2(y2), Z2(z2),
				X3(x3), Y3(y3), Z3(z3), S(s), A(a), FX(fx), FY(fy) {}

		/** Generates a string for exporting relevant gradient values in format: Slope Aspect FX FY
		* @return String of relevant values.
		*/
		std::string toString()
		{
			std::ostringstream oss;
			oss << S << " " << A << " " << FX  << " " << FY;
			return oss.str();
		}
};

/**
* Class for storing a Point in a surface (Normalized version of Ping)
* By normalized I mean that it has been converted to meters using the
* Vincenty functions and using a relative Ping.
* Can also be used for synthetic data points that may or may not be
* represented as data in meters.
*/
class Point
{
	public:
		/** X position of the Point */
		double x;

		/** Y position of the Point */
		double y;

		/** Z Position of the Point */
		double z;

		/** Uncertainty of the Point */
		double u;
		double u2;
		double u3;
		double u4;
		double u5;

		/** Gradient Slope of the point */
		double s;

		/** Vertex id of the Point */
		double id; //SJZ added for print

		/** Horizontal Uncertainty of the Point */
		double hU;

		/** Vertical Uncertainty of the Point */
		double vU;
		
		//double e; //bilinear error for raster uses u for now
		//Raster values
		double z0;
		double e0;
		double zK;
		double eK;
		double nei;
		double rei;

		/** Constructor for the Point
		* @param x - value for x.
		* @param y - value for y.
		* @param z - value for z.
		*/
		Point(const double& x = 0, const double& y = 0, const double& z = 0);

		/** Copy constructor
		* @param p - Point being copied.
		*/
		Point(const Point& p);

		/** Overloads the = operator for assignment operations.
		* @param p - Argument Point whose values will be stored in this Point.
		*/
		Point operator=(const Point& p);

		/** Overloads the == operator for Boolean comparison operations.
		* @return True if the Points have equal values, false otherwise (ignores uncertainty)
		*/
		bool operator==(const Point& p) const;

		/** Overloads the != operator for Boolean comparison operations.
		* @return True if the Points have different values, false otherwise (ignores uncertainty)
		*/
		bool operator!=(const Point& p) const;

		/** Overloads the < operator for Boolean comparison operations by
		* x -> y -> z, used for sorting to help remove duplicate points easily.
		* @return true if this Point is less than p, false otherwise
		*/
		bool operator<(const Point& p) const;

		/** Calculates and returns the distance between "this" Point and the Point p (ignores z)
		* @param p - Argument Point for finding distance from "this" Point
		* @return Distance from "this" Point to Point p
		*/
		double distance2d(const Point& p) const;

		/** Overloads the - operator for minus operations that returns
		* a resulting 2d vector (Vector2d) obect of Point1 - Point2.
		* @param p - Argument Point to subtract from "this" ping
		* @return Vector2d result of the minus operation.
		*/
		Vector2d operator-(const Point& p) const;

		/** A Boolean function for testing with two pings are equivalent
		* while ignoring the depth.
		* @param p - Argument Point to compare to "this" Point.
		* @return True if the x and y values are equal, false otherwise
		*/
		bool eq2d(const Point& p) const;

		/** A Boolean function for testing equality with a provided amount of "give"
		* @param p - Point being compared to this Point for equality.
		* @param err - The amount of allowed error for the equality test.
		* @return True if the the Points are equal within the provided error, false otherwise.
		*/
		bool eqErr(const Point& p, const double err) const;

		/** Outputs the values of this Point. */
		void print() const;
};

/** A Triangle class for representing and manipulating the vertices of
* a triangle as Points. Each vertice of the triangle is represented as
* a Point.
* The class offers helpful functions for sorting, storing, and retrieving
* information to and from the triangle.
*/
class Triangle
{
	public:
		/** Vertex 1 of the Triangle */
		Point v1;

		/** Vertex 2 of the Triangle */
		Point v2;

		/** Third vertice of the Triangle */
		Point v3;

		/**Constructor that takes in 3 Points and forms the Triangle.
	 	* @param V1 - Vertex 1
	 	* @param V2 - Vertex 2
	 	* @param V3 - Vertex 3
		*/
		Triangle(const Point& V1, const Point& V2, const Point& V3);

		/** Constructor that generates a triangle at vertices containing all zeros. */
		Triangle();

		/** Returns the midpoint of the Triangle's hypotenuse, (Assumes Triangle was formed for an rTin).
		* @return Midpoint of hypotenuse
		*/
		Point midHyp() const;

		/** Returns the length of a leg in the Triangle (Assumes Triangle was formed for an rTin)
		* @return double representing the length of one of the Triangle's legs
		*/
		double edgeLen() const;

		/** Sets the values of all vertices at once.
		* @param p1 - New vertex 1.
		* @param p2 - New vertex 2.
		* @param p3 - New vertex 3.
		*/
		void set(Point p1, Point p2, Point p3);

		/** Overloads the = operator for assignment operations on Triangles.
		* @param t - Triangle to be stored.
		*/
		Triangle operator=(const Triangle& t);

		/** Returns whether the Point p is inside "this" Triangle, ignoring depth.
		* @param p - Argument Point to test inclusion in "this" Triangle.
		* @return True if Point p is inside "this" Triangle and False otherwise.
		*/
		bool isInside(const Point& p) const;

		/** Overloads the < operator for comparison, Used for sorting and removing duplicates
		* @param t - Triangle to be compared against "this" Triangle.
		* @return True if the added totals of "this" Triangle are
		* less than that of Triangle t and False otherwise.
		*/
		bool operator<(const Triangle& t) const;

		/** Overloads operator == for boolean comparison.
		* @param t - Argument Triangle for comparison with "this" Triangle.
		* @return Result of comparison will be True of all vertices of
		* "this" Triangle are equivalent to that of Triangle t.
		*/
		bool operator==(const Triangle& t) const;

		/** Returns the height at which a vertical line through Point p would cross the plane formed by "this" Triangle.
		* @param p - Point for testing height of plane.
		* @return Height sample in the plane formed by "this" Triangle.
		*/
		double sample(const Point& p) const;

		/** Returns a Gradient object with the calculated Gradient of this Triangle
		* @return Gradient of this Triangle
		*/
		Gradient getGradient() const;
};

class QuadEdge;

/** Class for representing an Edge structure, contains functions for performing
* Edge algebra in an Quad-Edge data structure as defined by Guibas and Stolfi.
*/
class Edge
{
	/** This operator affects the two edge rings around the origins of a and b,
	* and, independently, the two edge rings around the left faces of a and b.
	* In each case, (i) if the two rings are distinct, Splice will combine
	* them into one; (ii) if the two are the same ring, Splice will break it
	* into two separate pieces.
	* Thus, Splice can be used both to attach the two edges together, and
	* to break them apart. See Guibas and Stolfi (1985) p.96 for more details
	* and illustrations.
	* @param a - Argument Edge a to be split from or combined with b.
	* @param b - Argument Edge b to be split from or combined with a.
	*/
	friend void splice(Edge* a, Edge* b);

	public:

		/** Edge number for keeping up with which edge of the quad edge this is. */
		int num;

		/** Edge id for keeping all quadEdges unique. */
		int id;

		/** Pointer to the next edge in the Quad-Edge structure. */
		Edge *next;

		/** Pointer to the data contained within "this" Edge. */
		Point *data;

		/** Constructor for instantianting the Edge. */
		Edge();

		/** Return the dual of the current edge, directed from its right to its left.
		* @return Dual of the current edge, directed from its right to its left.
		*/
		Edge* rot();

		/** Return the dual of the current edge, directed from its left to its right.
		* @return Dual of the current edge, directed from its left to its right.
		*/
		Edge* invRot();

		/** Return the edge from the destination to the origin of the current edge.
		* @return Edge from the destination to the origin of the current edge.
		*/
		Edge* sym();

		/** Return the next ccw edge around (from) the origin of the current edge.
		* @return Next ccw edge around (from) the origin of the current edge.
		*/
		Edge* oNext();

		/** Return the next cw edge around (from) the origin of the current edge.
		* @return Next cw edge around (from) the origin of the current edge.
		*/
		Edge* oPrev();

		/** Return the next ccw edge around (into) the destination of the current edge.
		* @return Next ccw edge around (into) the destination of the current edge.
		*/
		Edge* dNext();

		/** Return the next cw edge around (into) the destination of the current edge.
		* @return Next cw edge around (into) the destination of the current edge.
		*/
		Edge* dPrev();

		/** Return the ccw edge around the left face following the current edge.
		* @return ccw edge around the left face following the current edge.
		*/
		Edge* lNext();

		/** Return the ccw edge around the left face before the current edge.
		* @return ccw edge around the left face before the current edge.
		*/
		Edge* lPrev();

		/** Return the edge around the right face ccw following the current edge.
		* @return Edge around the right face ccw following the current edge.
		*/
		Edge* rNext();

		/** Return the edge around the right face ccw before the current edge.
		* @return Edge around the right face ccw before the current edge.
		*/
		Edge* rPrev();

		// Access Pointers //

		/** Return the origin of the edge.
		* @return Origin of the edge.
		*/
		Point* org();

		/** Return the destination of the edge.
		* @return Destination of the edge.
		*/
		Point* dest();

		/** Return the origin of the edge.
		* @return Origin of the edge.
		*/
		const Point& org2d() const;

		/** Return the destination of the edge.
		* @return Destination of the edge.
		*/
		const Point& dest2d() const;

		/** Sets the endpoints of the edge.
		* @param org - Argument Point for the origin of the Edge.
		* @param dest - Argument Point for the destination of the Edge.
		*/
		void endPoints(Point* org, Point* dest);

		/** Returns the QuadEdge of this edge.
		* @returns QuadEdge of this edge.
		*/
		QuadEdge* qEdge();
};

/** Class for representing the Quad-Edge structure as an
* array of 4 edge structures.
*/
class QuadEdge
{
	/** A friend member that returns a new Edge
	* @param id - Argument for setting the id of the new edge
	* @return newly created edge
	*/
	friend Edge *makeEdge(const int& id);
	friend Edge *makeEdge();

	public:
		/** Array of 4 edges forming the QuadEdge */
		Edge e[4];

		/** Constructor that sets the id of the new edges.
		* @param id - Argument id for the new edges
		*/
		QuadEdge(int id);
		QuadEdge();
};

/** Class for representing a geometric line */
class Line {
	public:
		/** Constructor for the Line
		* Does nothing
		*/
		Line();

		/** Constructor for the Line that sets a b and c
		* of the line
		* @param p - Argument representing one point on the line
		* @param q - Argument representing another point on the line
		*/
		Line(const Point& p, const Point& q);

		/** Member function for evaluating a relationship
		* of a Point to the line.
		* @param p - Point to compare to the line
		* @return Non zero number means the point is not on the line
		*/
		double eval(const Point& p) const;

		/** Member function for classifying a point on or off a line
		* @param p - Point to compare to the line
		* @return 0 for point on the line, 1 for right, -1 for left
		*/
		int classify(const Point& p) const;

	private:
		double a, b, c;
};

// This operator affects the two edge rings around the origins of a and b,
// and, independently, the two edge rings around the left faces of a and b.
// In each case, (i) if the two rings are distinct, Splice will combine
// them into one; (ii) if the two are the same ring, Splice will break it
// into two separate pieces.
// Thus, Splice can be used both to attach the two edges together, and
// to break them apart. See Guibas and Stolfi (1985) p.96 for more details
// and illustrations.
void splice(Edge* a, Edge* b);

void deleteEdge(Edge* e);

Edge* makeEdge(const int& id);
Edge* makeEdge();

// Add a new edge e connecting the destination of a to the
// origin of b, in such a way that all three have the same
// left face after the connection is complete.
// Additionally, the data pointers of the new edge are set.
Edge* connect(Edge* a, Edge* b, const int& id);
Edge* connect(Edge* a, Edge* b);

// Essentially turns edge e counterclockwise inside its enclosing
// quadrilateral. The data pointers are modified accordingly.
void swap(Edge* e);

/*************** Geometric Predicates for Delaunay Diagrams *****************/
// Returns twice the area of the oriented triangle (a, b, c), i.e., the
// area is positive if the triangle is oriented counterclockwise.
double triArea(const Point& a, const Point& b, const Point& c);

// Returns TRUE if the point d is inside the circle defined by the
// points a, b, c. See Guibas and Stolfi (1985) p.107.
int inCircle(const Point& a, const Point& b, const Point& c, const Point& d);
int inCircle2(const Point& a, const Point& b, const Point& c, const Point& d);

// Returns TRUE if the points a, b, c are in a counterclockwise order
int ccw(const Point& a, const Point& b, const Point& c);

int rightOf(const Point& x, Edge* e);

int LeftOf(const Point& x, Edge* e);

int isConvexQuad(const Point& a, const Point& b, const Point& c, const Point& d);

// A predicate that determines if the point x is on the edge e.
// The point is considered on if it is in the EPS-neighborhood
// of the edge.
int onEdge(const Point& x, Edge* e);

double differenceInPlanes(const Triangle& t1, const Triangle& t2);

double sign(const Point& p1, const Point& p2, const Point& p3);

bool pointInTriangle(const Point& p, const Point& v1, const Point& v2, const Point& v3);

double samplePlane(const Point& p, const Point& v1, const Point& v2, const Point& v3);

Gradient gradient(const Point& v1, const Point& v2, const Point& v3);

Point centerOfCircumCircle(const Point& p1, const Point& p2, const Point& p3);
double radiusOfCircumCircle(const Point& p1, const Point& p2, const Point& p3);

struct MinDistance {
	Point c;
	bool operator() (Point p1, Point p2) { return (c.distance2d(p1) < c.distance2d(p2)); }
};

/*Added by SJZ. can delete*/
#pragma region SJZ Added
int collinear(const Point& x, Edge* e);
int collinear(const Point& x, Point& b, Point& c);
bool fuzzyEquals(const Point& x, const Point& e);
double CFPE(double number);
#pragma endregion
