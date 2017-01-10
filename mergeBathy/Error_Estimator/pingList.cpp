
#include "pointList.h"
#include "pingList.h"
#include "geom.h"
#include "vincenty.h"
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>

#ifndef MAX_INT
#define MAX_INT 9999999
#define MIN_INT -9999999
#endif
#define MAX_LAT 361
#define MIN_LAT -361
#define MAX_LONG 361
#define MIN_LONG -361

using namespace std;

PingList::PingList () {}

PingList::PingList (const PingList& p)
{
	positions = p.getPositions();
	maxLat = p.getMaxLat();
	minLat = p.getMinLat();
	maxLong = p.getMaxLong();
	minLong = p.getMinLong();
	avgDepth = p.getAvgDepth();
	minDepth = p.getMinDepth();
	maxDepth = p.getMaxDepth();
}

void PingList::operator=(const PingList& p)
{
	positions = p.getPositions();
	maxLat = p.getMaxLat();
	minLat = p.getMinLat();
	maxLong = p.getMaxLong();
	minLong = p.getMinLong();
	avgDepth = p.getAvgDepth();
	minDepth = p.getMinDepth();
	maxDepth = p.getMaxDepth();
}

PingList::~PingList() { clear(); }

vector<Ping> PingList::getPositions() const { return positions; }

Ping PingList::operator[](const int index) const
{
	assert(index >= 0 && index < positions.size());
	return positions[index];
}

void PingList::read(const string& fileName)
{
	// Input file stream
	ifstream inData;

	// Open the file "fileName" with the filestream
	inData.open(fileName.c_str());

	// holder lat, long, depth, and ping for read in
	double lon, lat, dep;
	
	// Prepare bounding points to be found
	minLat = MAX_LAT;
	maxLat = MIN_LAT;
	minLong = MAX_LONG;
	maxLong = MIN_LONG;
	avgDepth = 0;
	minDepth = MAX_INT;
	maxDepth = MIN_INT;

	list<Ping> pList;
	list<Ping>::iterator pIt;

	// While not end of file for stream
	while(!inData.eof())
	{
		// Read in lat long and depth
		//inData >> lat >> lon >> dep;
		inData >> lon >> lat >> dep;
		
		// Invert Depth to be negative
		dep = fabs(dep) * -1;

		// Continue search for bounding points
		minLat = min(lat, minLat);
		minLong = min(lon, minLong);
		maxLat = max(lat, maxLat);
		maxLong = max(lon, maxLong);
		minDepth = min(dep, minDepth);
		maxDepth = max(dep, maxDepth);
		avgDepth += dep;

		// Set the ping to new values
		Ping p(lon, lat, dep);
		// push the new ping to the end of the list
		pList.push_back(p);
	}
	
	avgDepth /= positions.size();
	
	// close the input stream
	inData.close();
	
	positions.reserve(pList.size());
	for(pIt = pList.begin(); pIt != pList.end(); pIt++)
		positions.push_back(*pIt);
}

int PingList::size() const { return (const int)positions.size(); }

void PingList::decimate(double x)
{
	assert(x > 0 && x < 100);
	rand();
	double p = x/100;
	double num = p * positions.size();
	for(int i = 0; i < num; i++)
		positions.pop_back();
}

double PingList::min (const double& x, const double& y) const { if(x < y) { return x; } return y; }

double PingList::max (const double& x, const double& y) const { if(x > y) { return x; } return y; }

void PingList::sort(){ std::sort(positions.begin(), positions.end()); }

void PingList::unique(){ std::unique(positions.begin(), positions.end()); }

void PingList::rand(){ std::random_shuffle(positions.begin(), positions.end()); }

//void PingList::reserve(const int size){ positions.reserve(size); }

/*void PingList::push_back(const Ping p) 
{ 
	minLat = min(p.getLat(), minLat);
	minLong = min(p.getLong(), minLong);
	maxLat = max(p.getLat(), maxLat);
	maxLong = max(p.getLong(), maxLong);
	minDepth = min(p.getDepth(), minDepth);
	maxDepth = max(p.getDepth(), maxDepth);
	avgDepth *= positions.size();
	avgDepth += p.getDepth();
	positions.push_back(p); 
	avgDepth /= positions.size();
}*/

//void PingList::pop_back() { positions.pop_back(); }

void PingList::clear() 
{ 
	positions.clear(); 
	minLat = MAX_LAT;
	maxLat = MIN_LAT;
	minLong = MAX_LONG;
	maxLong = MIN_LONG;
	avgDepth = 0;
	minDepth = MAX_INT;
	maxDepth = MIN_INT;
}

/*PointList* PingList::getPointList() const
{
	PointList *pList;
	pList->reserve(positions.size());
	Ping p(minLong, minLat);
	Vincenty v(p);

	for(int i = 0; i < positions.size(); i++)
		pList->push_back(v.normalizePing(positions[i]));

	return pList;
}*/

PointList* PingList::getPointList() const
{
	vector<Point> pv;

	Ping p(minLong, minLat);
	Vincenty v(p);
	Point pt;
	
	for(int i = 0; i < (const int)positions.size(); i++)
	{
		pt = v.normalizePing(positions[i]);
		if((pt.x == pt.x) && (pt.y == pt.y)) 
			pv.push_back(pt);
	}
	
	PointList *pList;
	pList = new PointList(pv);
	return pList;
}

/*void PingList::calcStats()
{
	minLat = MAX_LAT;
	maxLat = MIN_LAT;
	minLong = MAX_LONG;
	maxLong = MIN_LONG;
	avgDepth = 0;
	minDepth = MAX_INT;
	maxDepth = MIN_INT;

	Ping p;
	
	for(int i = 0; i < positions.size(); i++)
	{
		p = positions[i];
		minLat = min(p.getLat(), minLat);
		minLong = min(p.getLong(), minLong);
		maxLat = max(p.getLat(), maxLat);
		maxLong = max(p.getLong(), maxLong);
		minDepth = min(p.getDepth(), minDepth);
		maxDepth = max(p.getDepth(), maxDepth);
		avgDepth += p.getDepth();
	}
	
	avgDepth /= positions.size();
}*/

void PingList::printStats()
{
	cout.precision(8);
	cout << fixed;
	cout << "Size:              " << size() << endl;
	cout << "Minimum Longitude: " << minLong << endl;
	cout << "Maximum Longitude: " << maxLong << endl;
	cout << "Minimum Latitude:  " << minLat << endl;
	cout << "Maximum Latitude:  " << maxLat << endl;
	cout << "Minimum Depth:     " << minDepth << endl;
	cout << "Maximum Depth:     " << maxDepth << endl;
	cout << "Average Depth:     " << avgDepth << endl;
}
