#include <iomanip>
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
//#include "s_hull.h"
#include "geom.h"
#include "pointList.h"
#include "pingList.h"
#include "vincenty.h"

using namespace std;

PointList::PointList () { }

void PointList::clear()
{
	positions.clear();
	maxX = MIN_INT;
	maxY = MIN_INT;
	avgZ = 0;
	minZ = MAX_INT;
	maxZ = MIN_INT;
}

bool uniqueLL (Point first, Point second)
{ return ( first.x == second.x && first.y == second.y) ; }
bool sortLL (Point first, Point second)
{ return ( first.x < second.x );}

PointList::PointList (const PointList& p)
{
	positions = p.getPositions();
	maxY = p.getMaxY();
	maxX = p.getMaxX();
	avgZ = p.getAvgZ();
	minZ = p.getMinZ();
	maxZ = p.getMaxZ();
	offsetX = p.getOffsetX();
	offsetY = p.getOffsetY();

	/*avgZ0 = p.getAvgZ0();
	minZ0 = p.getMinZ0();
	maxZ0 = p.getMaxZ0();
	avgZK = p.getAvgZK();
	minZK = p.getMinZK();
	maxZK = p.getMaxZK();*/
	
}

void PointList::operator=(const PointList& p)
{
	positions = p.getPositions();
	maxY = p.getMaxY();
	maxX = p.getMaxX();
	avgZ = p.getAvgZ();
	minZ = p.getMinZ();
	maxZ = p.getMaxZ();
	offsetX = p.getOffsetX();
	offsetY = p.getOffsetY();
}

PointList::~PointList(){ clear(); }

vector<Point> PointList::getPositions() const { return positions; }

//Sam Not sure when this function is used or what it was for but it needs to be tested.
//The offset also needs to be examined as it was found in the case of the setFromVectors 
//functions, where one creates a pointList for the Delaunay Triangulation and the other 
//creates the pointList for querying the tin.
//For this reason, the offsets need to be the same so that both are translated the same.  
PointList::PointList (vector<Point>& p)
{
	if(!DEBUG_DISABLE_PTLOFFSET)
	{
		cout << endl;
		cout << "This function has not been tested, nor its purpose found therefore the following warning is exhibited." << endl;
		cout << "OFFSET TRANSLATION WARNING!!!!!" << endl;
		cout << "An offset will be calculated for all x values and all y values and then added to each x and y, respectively." << endl;  
		cout << "This translates the set of points to the origin." << endl;
		cout << "This is fine as long as these points will not be used in conjunction with another point list similarly translated with its own computed offsets." << endl;
		cout << "By using different offsets, both point lists will center around the origin individually instead of together." << endl;
		cout << "Therefore, both point lists should used the same offsets for their translations." << endl;
		cout << "If point list interaction is unknown, it is advised set DEBUG_DISABLE_PTLOFFSET = 1" << endl;
		cout << "Press any key to continue... " << endl;
		cin.get();
	}

	// Prepare bounding points to be found
	double minX = MAX_INT;
	double minY = MAX_INT;

	maxX = MIN_INT;
	maxY = MIN_INT;
	avgZ = 0;
	minZ = MAX_INT;
	maxZ = MIN_INT;
	offsetX = 0;
	offsetY = 0;

	for(int i = 0; i < (const int)p.size(); i++)
	{
		minX = min(p[i].x, minX);
		minY = min(p[i].y, minY);
		maxX = max(p[i].x, maxX);
		maxY = max(p[i].y, maxY);
		minZ = min(p[i].z, minZ);
		maxZ = max(p[i].z, maxZ);
		avgZ += p[i].z;
	}

	avgZ /= p.size();

	positions.reserve(p.size());

	offsetX -= minX;
	offsetY -= minY;

	minX += offsetX;
	minY += offsetY;
	maxX += offsetX;
	maxY += offsetY;
	for(int i = 0; i < (const int)p.size(); i++)
	{
		if(DEBUG_DISABLE_PTLOFFSET){ offsetX = 0; offsetY = 0;}
		p[i].x += offsetX;
		p[i].y += offsetY;
		positions.push_back(p[i]);
	}
}

Point PointList::operator[](const int index) const
{
	assert(index >= 0 && index < positions.size());
	return positions[index];
}

void PointList::readText(const string& fileName)
{
	if(!DEBUG_DISABLE_PTLOFFSET)
	{
		cout << endl;
		cout << "This function has not been tested, nor its purpose found therefore the following warning is exhibited." << endl;
		cout << "OFFSET TRANSLATION WARNING!!!!!" << endl;
		cout << "An offset will be calculated for all x values and all y values and then added to each x and y, respectively." << endl;  
		cout << "This translates the set of points to the origin." << endl;
		cout << "This is fine as long as these points will not be used in conjunction with another point list similarly translated with its own computed offsets." << endl;
		cout << "By using different offsets, both point lists will center around the origin individually instead of together." << endl;
		cout << "Therefore, both point lists should used the same offsets for their translations." << endl;
		cout << "If point list interaction is unknown, it is advised set DEBUG_DISABLE_PTLOFFSET = 1" << endl;
		cout << "Press any key to continue... " << endl;
		cin.get();
	}
	// Input file stream
	ifstream inData;

	// Open the file "fileName" with the file stream
	inData.open(fileName.c_str());

	// Holder lat, long, depth, and ping for read in
	double x, y, z;

	// Prepare bounding points to be found
	double minX = MAX_INT;
	double minY = MAX_INT;
	
	maxX = MIN_INT;
	maxY = MIN_INT;
	avgZ = 0;
	minZ = MAX_INT;
	maxZ = MIN_INT;
	offsetX = 0;
	offsetY = 0;

	//list<Point> pList;
	//list<Point>::iterator pIt;
	vector<Point> pList;
	vector<Point>::iterator pIt;
	Point prev(0,0,0);

	// While not end of file for stream
	while(!inData.eof())
	{
		// Read in lat long and depth
		//inData >> lat >> lon >> dep;
		inData >> x >> y >> z;

		// Invert Depth to be negative
		//z = fabs(z) * -1;

		// Continue search for bounding points
		minX = min(x, minX);
		minY = min(y, minY);
		maxX = max(x, maxX);
		maxY = max(y, maxY);
		minZ = min(z, minZ);
		maxZ = max(z, maxZ);
		avgZ += z;

		// Set the ping to new values
		Point p(x, y, z);

		// push the new ping to the end of the list
		if(p != prev)
			pList.push_back(p);
		prev = p;
	}
	//pList.sort(&sortLL);//SJZ
	std::sort(pList.begin(),pList.end(), &sortLL);
	pList.erase(std::unique(pList.begin(), pList.end(), &uniqueLL), pList.end());

	avgZ /= pList.size();

	// close the input stream
	inData.close();

	positions.reserve(pList.size());

	offsetX -= minX;
	offsetY -= minY;

	minX += offsetX;
	minY += offsetY;
	maxX += offsetX;
	maxY += offsetY;

	for(pIt = pList.begin(); pIt != pList.end(); pIt++)
	{
		if(DEBUG_DISABLE_PTLOFFSET){ offsetX = 0; offsetY = 0;}
		pIt->x += offsetX;
		pIt->y += offsetY;
		positions.push_back(*pIt);
	}
}


void PointList::setFromVectorsRaster(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn, const vector<double>& eIn/*, vector<double>& nEiIn, vector<double> rEiIn*//*, const vector<double>& z0In, const vector<double>& e0In, const vector<double>& zKIn, const vector<double>& eKIn*/)
{
	// Holder lat, long, depth, and ping for read in
	double x, y, z, e;//, z0, e0, zK, eK;

	// Prepare bounding points to be found
	int currentLocation = 0;
	double minX = (double)MAX_INT;
	double minY = (double)MAX_INT;

	maxX = (double)MIN_INT;
	maxY = (double)MIN_INT;
	avgZ = 0;
	minZ = (double)MAX_INT;
	maxZ = (double)MIN_INT;
	/*minZ0 = (double)MAX_INT;
	maxZ0 = (double)MIN_INT;
	minZK = (double)MAX_INT;
	maxZK = (double)MIN_INT;*/

	offsetX = 0;
	offsetY = 0;

	//list<Point> pList;
	//list<Point>::iterator pIt;
	vector<Point> pList;
	vector<Point>::iterator pIt;
	Point prev(0,0,0);

	while(currentLocation < (const int)xIn.size())
	{
		// Read in lat long and depth
		x = xIn[currentLocation];
		y = yIn[currentLocation];
		z = zIn[currentLocation];
		e = eIn[currentLocation];
		
		// Invert Depth to be negative
		//z = fabs(z) * -1;

		// Continue search for bounding points
		minX = min(x, minX);
		minY = min(y, minY);
		maxX = max(x, maxX);
		maxY = max(y, maxY);
		minZ = min(z, minZ);
		maxZ = max(z, maxZ);
		avgZ += z;

		/*minZ0 = min(z, minZ0);
		maxZ0 = max(z, maxZ0);
		avgZ0 += z0;

		minZK = min(zK, minZK);
		maxZK = max(zK, maxZK);
		avgZK += zk;*/

		// Set the ping to new values
		Point p(x, y, z);
		p.z = z;
		//**********************************************************
//		p.e = e;
		p.u = e;
		
		// push the new ping to the end of the list
		// rule needs to be placed/checked for which to use sort unique last first
		// the matlab version only cares about lat and lon, not depth
		if(p != prev)
			pList.push_back(p);
		prev = p;
		currentLocation += 1;
	}
	//pList.sort(&sortLL);
	std::sort(pList.begin(),pList.end(), &sortLL);//SJZ
	pList.erase(std::unique(pList.begin(), pList.end(), &uniqueLL), pList.end());

	avgZ /= (double)pList.size();
	//avgZ0 /= (double)pList.size();
	//avgZK /= (double)pList.size();

	positions.reserve(pList.size());

	offsetX -= minX;
	offsetY -= minY;

	minX += offsetX;
	minY += offsetY;
	maxX += offsetX;
	maxY += offsetY;

	int ind=0;
	for(pIt = pList.begin(); pIt != pList.end(); pIt++)
	{
		if(DEBUG_DISABLE_PTLOFFSET){ offsetX = 0; offsetY = 0;}
		pIt->x += offsetX;
		pIt->y += offsetY;
		pIt->id = ind;
		ind++;
		positions.push_back(*pIt);
	}
}


void PointList::setFromVectorsRaster(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn, const vector<double>& eIn, vector<double>& neiIn, vector<double> reiIn, const vector<double>& z0In, const vector<double>& e0In, const vector<double>& zKIn, const vector<double>& eKIn)
{
	// Holder lat, long, depth, and ping for read in
	double x, y, z, e, z0, e0, zK, eK, nei, rei;

	// Prepare bounding points to be found
	int currentLocation = 0;
	double minX = (double)MAX_INT;
	double minY = (double)MAX_INT;

	maxX = (double)MIN_INT;
	maxY = (double)MIN_INT;
	avgZ = 0;
	minZ = (double)MAX_INT;
	maxZ = (double)MIN_INT;
	minZ0 = (double)MAX_INT;
	maxZ0 = (double)MIN_INT;
	minZK = (double)MAX_INT;
	maxZK = (double)MIN_INT;

	offsetX = 0;
	offsetY = 0;

	//list<Point> pList;
	//list<Point>::iterator pIt;
	vector<Point> pList;
	vector<Point>::iterator pIt;
	Point prev(0,0,0);
	pList.reserve(xIn.size());
	while(currentLocation < (const int)xIn.size())
	{
		// Read in lat long and depth
		x = xIn[currentLocation];
		y = yIn[currentLocation];
		z = zIn[currentLocation];
		e = eIn[currentLocation];
		z0 = z0In[currentLocation];
		e0 = e0In[currentLocation];
		zK = zKIn[currentLocation];
		eK = eKIn[currentLocation];
		nei = neiIn[currentLocation];
		rei = reiIn[currentLocation];

		// Invert Depth to be negative
		//z = fabs(z) * -1;

		// Continue search for bounding points
		minX = min(x, minX);
		minY = min(y, minY);
		maxX = max(x, maxX);
		maxY = max(y, maxY);
		minZ = min(z, minZ);
		maxZ = max(z, maxZ);
		avgZ += z;

		minZ0 = min(z, minZ0);
		maxZ0 = max(z, maxZ0);
		avgZ0 += z0;

		minZK = min(zK, minZK);
		maxZK = max(zK, maxZK);
		avgZK += zK;

		// Set the ping to new values
		Point p(x, y, z);
		p.z = z;
		//**********************************************************
//		p.e = e;
		p.u = e;
		p.z0 = z0;
		p.e0 = e0;
		p.zK = zK;
		p.eK = eK;
		p.nei = nei;
		p.rei = rei;

		// push the new ping to the end of the list
		// rule needs to be placed/checked for which to use sort unique last first
		// the matlab version only cares about lat and lon, not depth
		if(p != prev)
			pList.push_back(p);
		prev = p;
		currentLocation += 1;
	}
	//pList.sort(&sortLL);
	std::sort(pList.begin(),pList.end(), &sortLL);
	pList.erase(std::unique(pList.begin(), pList.end(), &uniqueLL), pList.end());

	avgZ /= (double)pList.size();
	avgZ0 /= (double)pList.size();
	avgZK /= (double)pList.size();

	positions.reserve(pList.size());

	offsetX -= minX;
	offsetY -= minY;

	minX += offsetX;
	minY += offsetY;
	maxX += offsetX;
	maxY += offsetY;

	int ind=0;
	for(pIt = pList.begin(); pIt != pList.end(); pIt++)
	{
		if(DEBUG_DISABLE_PTLOFFSET){ offsetX = 0; offsetY = 0;}
		pIt->x += offsetX;
		pIt->y += offsetY;
		pIt->id = ind;
		ind++;
		positions.push_back(*pIt);
	}
}


void PointList::setFromVectors(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn, const vector<double>& hIn, const vector<double>& vIn)
{
	// Holder lat, long, depth, and ping for read in
	double x, y, z, h, v;

	// Prepare bounding points to be found
	int currentLocation = 0;
	double minX = (double)MAX_INT;
	double minY = (double)MAX_INT;

	maxX = (double)MIN_INT;
	maxY = (double)MIN_INT;
	avgZ = 0;
	minZ = (double)MAX_INT;
	maxZ = (double)MIN_INT;
	offsetX = 0;
	offsetY = 0;

	//list<Point> pList;
	//list<Point>::iterator pIt;
	vector<Point> pList;
	vector<Point>::iterator pIt;
	Point prev(0,0,0);

	while(currentLocation < (const int)xIn.size())
	{
		// Read in lat long and depth
		x = xIn[currentLocation];
		y = yIn[currentLocation];
		z = zIn[currentLocation];
		h = hIn[currentLocation];
		v = vIn[currentLocation];

		// Invert Depth to be negative
		//z = fabs(z) * -1;

		// Continue search for bounding points
		minX = min(x, minX);
		minY = min(y, minY);
		maxX = max(x, maxX);
		maxY = max(y, maxY);
		minZ = min(z, minZ);
		maxZ = max(z, maxZ);
		avgZ += z;

		// Set the ping to new values
		Point p(x, y, z);
		p.hU = h;
		p.vU = v;

		// push the new ping to the end of the list
		// rule needs to be placed/checked for which to use sort unique last first
		// the matlab version only cares about lat and lon, not depth
		if(p != prev)
			pList.push_back(p);
		prev = p;
		currentLocation += 1;
	}
	//pList.sort(&sortLL);
	std::sort(pList.begin(),pList.end(), &sortLL); //SJZ
	pList.erase(std::unique(pList.begin(), pList.end(), &uniqueLL), pList.end());

	avgZ /= (double)pList.size();

	positions.reserve(pList.size());

	offsetX -= minX;
	offsetY -= minY;

	minX += offsetX;
	minY += offsetY;
	maxX += offsetX;
	maxY += offsetY;

	int ind=0;
	for(pIt = pList.begin(); pIt != pList.end(); pIt++)
	{
		if(DEBUG_DISABLE_PTLOFFSET){ offsetX = 0; offsetY = 0;}
		pIt->x += offsetX;
		pIt->y += offsetY;
		pIt->id = ind;
		ind++;
		positions.push_back(*pIt);
	}
}

//Create PointList for querying Delaunay Triangulation using the same offsetX and offsetY used for the triangulation's PointList.
void PointList::setFromVectors(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn, double offsetX, double offsetY)
{
	// Holder lat, long, depth, and ping for read in
	double x, y, z;

	// Prepare bounding points to be found
	int currentLocation = 0;
	double minX = (double)MAX_INT;
	double minY = (double)MAX_INT;

	maxX = (double)MIN_INT;
	maxY = (double)MIN_INT;
	avgZ = 0;
	minZ = (double)MAX_INT;
	maxZ = (double)MIN_INT;
	offsetX = 0;
	offsetY = 0;

	//list<Point> pList;
	//list<Point>::iterator pIt;
	vector<Point> pList;
	vector<Point>::iterator pIt;
	Point prev(0,0,0);

	// While not end of file for stream
	while(currentLocation < (const int)xIn.size())
	{
		// Read in lat long and depth
		x = xIn[currentLocation];
		y = yIn[currentLocation];
		z = zIn[currentLocation];

		// Invert Depth to be negative
		//z = fabs(z) * -1;

		// Continue search for bounding points
		minX = min(x, minX);
		minY = min(y, minY);
		maxX = max(x, maxX);
		maxY = max(y, maxY);
		minZ = min(z, minZ);
		maxZ = max(z, maxZ);
		avgZ += z;

		// Set the ping to new values
		Point p(x, y, z);
		p.hU = 0.00;
		p.vU = 0.00;

		// push the new ping to the end of the list
		if(p != prev)
			pList.push_back(p);
		prev = p;
		currentLocation += 1;
	}

	avgZ /= (double)pList.size();

	positions.reserve(pList.size());

	minX += offsetX;
	minY += offsetY;
	maxX += offsetX;
	maxY += offsetY;

	for(pIt = pList.begin(); pIt != pList.end(); pIt++)
	{
		if(DEBUG_DISABLE_PTLOFFSET){ offsetX = 0; offsetY = 0;}
		pIt->x += offsetX;
		pIt->y += offsetY;
		positions.push_back(*pIt);
	}
}

int PointList::size() const { return (const int)positions.size(); }

void PointList::decimate(double x)
{
	assert(x > 0 && x < 100);
	rand();
	double p = x/100;
	double num = p * positions.size();
	for(int i = 0; i < num; i++)
		positions.pop_back();
}

double PointList::min (const double& x, const double& y) const { if(x < y) { return x; } return y; }

double PointList::max (const double& x, const double& y) const { if(x > y) { return x; } return y; }

void PointList::rand(){ std::random_shuffle(positions.begin(), positions.end()); }

void PointList::sort(){ std::sort(positions.begin(), positions.end()); }

// Sorts radially by distance from p
void PointList::sort(const Point& p)
{
	MinDistance md;
	md.c = p;
	std::sort(positions.begin(), positions.end(), md);
}

void PointList::unique(){ std::unique(positions.begin(), positions.end()); }

//void PointList::reserve(const int size){ positions.reserve(size); }

/*void PointList::push_back(const Point p)
{
	int newOffsetX = 0;
	int newOffsetY = 0;

	maxX = max(p.x, maxX);
	maxY = max(p.y, maxY);

	minZ = min(p.z, minZ);
	maxZ = max(p.z, maxZ);
	avgZ *= positions.size();
	avgZ += p.z;

	positions.push_back(p);
	avgZ /= positions.size();
}*/

//void PointList::pop_back() { positions.pop_back(); }

/*PingList* PointList::getPingList(const Ping p) const
{
	PingList *pList;
	pList->reserve(positions.size());
	Vincenty v(p);

	for(int i = 0; i < positions.size(); i++)
	{
		Point pt(positions[i].x - offsetX,
				positions[i].y - offsetY,
				positions[i].z)
		pList->push_back(v.normalizePoint(pt));
	}

	return pList;
}*/

/*void PointList::calcStats()
{
	minX = MAX_INT;
	maxX = MIN_INT;
	minY = MAX_INT;
	maxY = MIN_INT;
	avgZ = 0;
	minZ = MAX_INT;
	maxZ = MIN_INT;

	Point p;

	for(int i = 0; i < positions.size(); i++)
	{
		p = positions[i];
		minX = min(p.x, minX);
		minY = min(p.y, minY);
		maxX = max(p.x, maxX);
		maxY = max(p.y, maxY);
		minZ = min(p.z, minZ);
		maxZ = max(p.z, maxZ);
		avgZ += p.z;
	}

	avgZ /= positions.size();
}*/

void PointList::printStats()
{
	cout.precision(8);
	cout << fixed;
	cout << "Size:      " << size() << endl;
	cout << "Maximum x: " << maxX << endl;
	cout << "Maximum y: " << maxY << endl;
	cout << "Minimum z: " << minZ << endl;
	cout << "Maximum z: " << maxZ << endl;
	cout << "Average z: " << avgZ << endl;
}


