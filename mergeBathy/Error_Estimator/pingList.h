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
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include "geom.h"
#include "pointList.h"

//using namespace std;

class PointList;

/** Class for reading in, storing, and accessing Pings from text files */
class PingList
{
	private:
		/** Vector of Pings (the list) */
		std::vector<Ping> positions;
		
		/** Maximum Latitude */
		double maxLat;
		
		/** Minimum Latitude */
		double minLat;
		
		/** Maximum Longitude */
		double maxLong;
		
		/** Minumum Longitude */
		double minLong;
		
		/** Average Depth */
		double avgDepth;

		/** Minimum Depth */
		double minDepth;
		
		/** Maximum Depth */
		double maxDepth;
		
		/** Simple min function
		* @param x - Double to be compared.
		* @param y - Double to be compared.
		* @return - lowest value of x and y.
		*/
		double min(const double& x, const double& y) const;
		
		/** Simple max function
		* @param x - Double to be compared.
		* @param y - Double to be compared.
		* @return - Maximum value of x and y.
		*/
		double max(const double& x, const double& y) const;

	public:

		/** Constructor, simply initializes variables */
		PingList();

		/** Copy Constructor, stores contents of p in the new PingList
		* @param p - PingList being copied.
		*/
		PingList(const PingList& p);

		/** Copies contents of PingList p into this PingList
		* @param p - PingList being copied.
		*/
		void operator=(const PingList& p);

		/** Destructor clears PingList contents */
		~PingList();
		
		/** Returns the vector of Pings in the structure.
		* @return vector of Pings in the structure.
		*/
		std::vector<Ping> getPositions() const;
		
		PointList* getPointList() const;

		/** Reads a list of Pings from a text file into the list, formated: Long Lat Depth.
		* @param fileName - File to be read.
		*/
		void read(const std::string& fileName);
		
		/** Returns the number of Pings in the list.
		* @return Number of Pings in the list.
		*/
		int size() const;

		/** Returns the Average Depth of the Pings in the list.
		* @return Average Depth of the Pings in the list.
		*/
		double getAvgDepth() const { return avgDepth; }
		
		/** Returns the Minimum Latitude of the Pings in the list.
		* @return Minimum Latitude of the Pings in the list.
		*/
		double getMinLat() const { return minLat; }
		
		/** Returns the Maximum Latitude of the Pings in the list.
		* @return Maximum Latitude of the Pings in the list.
		*/
		double getMaxLat() const { return maxLat; }

		/** Returns the Minimum Longitude of the Pings in the list.
		* @return Minimum Longitude of the Pings in the list.
		*/
		double getMinLong() const { return minLong; }
		
		/** Returns the Maximum Longitude of the Pings in the list.
		* @return Maximum Longitude of the Pings in the list.
		*/
		double getMaxLong() const { return maxLong; }

		/** Returns the Minimum Depth of the Pings in the list.
		* @return Minimum depth of the Pings in the list.
		*/
		double getMinDepth() const { return minDepth; }
		
		/** Returns the Maximum Depth of the Pings in the list.
		* @return Maximum Depth of the Pings in the list.
		*/
		double getMaxDepth() const { return maxDepth; }

		/** Overloads the [] operator to allow for indexed access to list elements.
		* @param index - Integer 0 - this.size()-1 .
		* @return Ping in lits at position index.
		*/
		Ping operator[](const int index) const;

		/** Randomizes the list */
		void rand();

		/** Sorts the list */
		void sort();

		/** Removes duplicate elements (should be sorted first) */
		void unique();

		/** Drops the last x percent of Pings in the list.
		* @param x - Percentage of Pings to be dropped.
		*/
		void decimate(double x);

		/** Clears the list */
		void clear();

		/** Prints the stats of the list */
		void printStats();
		
		//void reserve(const int size);

		//void push_back(const Ping p);

		//void pop_back();
		
		//void calcStats();
};

