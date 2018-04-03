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
#include <string>
#include <vector>
#include "geom.h"
#include "pingList.h"
#include "../constants.h"
#include "math.h"

class PingList;

/** Class for reading in, storing, and accessing Points from text files and images, also
* converts PingLists to projected Points and back if needed.
*/
class PointList
{
	private:
		/** Vector of Points (the List) */
		std::vector<Point> positions;

		/** Maximum X dimension */
		double maxX;

		/** Maximum Y dimension */
		double maxY;

		/** Average Z value (depth) */
		double avgZ;

		/** Minimum Z value (depth) */
		double minZ;

		/** Maximum Z value (depth) */
		double maxZ;

		
		//raster**********************************************************
		/** average z value (depth) */
		double avgZ0;

		/** minimum z value (depth) */
		double minZ0;

		/** maximum z value (depth) */
		double maxZ0;
		
		/** average z value (depth) */
		double avgZK;

		/** minimum z value (depth) */
		double minZK;

		/** maximum z value (depth) */
		double maxZK;
		//***************************************************************





		/** Offset of the data to set the lowest X to 0 */
		double offsetX;

		/** Offset of the data to set the lowest Y to 0 */
		double offsetY;

		/** A simple min function for returning the minimum double
		* @param x - double to be compared.
		* @param y - double to be compared.
		* @return - lowest value of x and y.
		*/
		double min(const double& x, const double& y) const;

		/** A simple max function for returning the maximum double
		* @param x - double to be compared.
		* @param y - double to be compared.
		* @return - maximum value of x and y.
		*/
		double max(const double& x, const double& y) const;

	public:

		/** Constructor, simply initializes variables */
		PointList ();

		/** Copy Constructor, stores contents of p in the new PointList
		* @param p - PointList being copied.
		*/
		PointList (const PointList& p);

		/** Constructor from a vector of Points
		* @param p - std::vector<Point> to be stored as the values for this PointList
		*/
		PointList (std::vector<Point>& p);

		/** Copies contents of PointList p into this PointList
		* @param p - PointList being copied.
		*/
		void operator=(const PointList& p);

		/** Destructor clears PointList contents */
		~PointList();

		/** Returns the vector of Points in the structure.
		* @return vector of Points in the structure.
		*/
		std::vector<Point> getPositions() const;

		//PingList* getPingList(const Ping p) const;

		/** Reads a list of Points from a text file into the list, formatted: X Y Z
		* @param fileName - File to be read.
		*/
		void readText(const std::string& fileName);

		//raster**********************************************************
		void setFromVectorsRaster(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn, const vector<double>& eIn);
		
		void setFromVectorsRaster(const vector<double>& xIn, const vector<double>& yIn, const vector<double>& zIn, const vector<double>& eIn, vector<double>& neiIn, vector<double> reiIn, const vector<double>& z0In, const vector<double>& e0In, const vector<double>& zKIn, const vector<double>& eKIn);
		//raster**********************************************************

		/**  Places points from a text file into the list and computed the offset for other pointLists to use, formatted: X Y Z H V
		* @param xIn - vector of data points x.
		* @param yIn - vector of data points y.
		* @param zIn - vector of data points z.
		* @param hIn - vector of data points h.
		* @param vIn - vector of data points v.
		*/
		void setFromVectors(const std::vector<double>& xIn, const std::vector<double>& yIn, const std::vector<double>& zIn, const std::vector<double>& hIn, const std::vector<double>& vIn);

		/** Places points from a text file into the list and applies a previously found offset, formatted: X Y Z 
		* @param xIn - vector of data points x.
		* @param yIn - vector of data points y.
		* @param zIn - vector of data points z.
		* @param offsetX - x offset value.
		* @param offsetY - y offset value.
		*/
		void setFromVectors(const std::vector<double>& xIn, const std::vector<double>& yIn, const std::vector<double>& zIn, double offsetX, double offsetY);

		/** Returns the number of Points in the list.
		* @return Number of Points in the list currently.
		*/
		int size() const;

		/** Returns the Average Z (depth) of the Points in the list.
		* @return Average Z (depth) of the Points in the list.
		*/
		double getAvgZ() const { return avgZ; }

		/** Returns the Maximum X of the Points in the list.
		* @return Maximum X of the Points in the list.
		*/
		double getMaxX() const { return maxX; }

		/** Returns the Maximum Y of the Points in the list.
		* @return the Maximum Y of the Points in the list.
		*/
		double getMaxY() const { return maxY; }

		/** Returns the Minimum Z (depth) of the Points in the list.
		* @return Minimum Z (depth) of the Points in the list.
		*/
		double getMinZ() const { return minZ; }

		/** Returns the Maximum Z (depth) of the Points in the list.
		* @return Maximum Z (depth) of the Points in the list.
		*/
		double getMaxZ() const { return maxZ; }

		/** Returns the X Offset to set the lowest X value to 0.
		* @return X Offset to set the lowest X value to 0.
		*/
		double getOffsetX() const { return offsetX; }

		/** Returns the Y Offset to set the lowest Y value to 0.
		* @return Y Offset to set the lowest Y value to 0.
		*/
		double getOffsetY() const { return offsetY; }

		/** Sets the Average Z (depth) of the Points in the list.
		* @param temp - value to set.
		*/
		void setAvgZ(double temp){ avgZ = temp; }

		/** Sets the Maximum Z (depth) of the Points in the list.
		* @param temp - Double value to set.
		*/
		void setMaxZ(double temp){ maxZ = temp; }

		/** Sets the Minimum Z (depth) of the Points in the list.
		* @param temp - Double value to set.
		*/
		void setMinZ(double temp){ minZ = temp; }

		/** Sets the Z (depth) of a Points in the list.
		* @param i - Integer 0 - this.size()-1 .
		* @param temp - Double value to set.
		*/
		void setZ(int i, double temp){ positions[i].z = temp;}
		void setE(int i, double temp){ positions[i].u = temp;}
//		void setE(int i, double temp){ positions[i].e = temp;}



		//raster**********************************************************
		void setZ0(int i, double temp){ positions[i].z0 = temp;}
		void setE0(int i, double temp){ positions[i].e0 = temp;}
		void setZK(int i, double temp){ positions[i].zK = temp;}
		void setEK(int i, double temp){ positions[i].eK = temp;}
		void setNEI(int i, double temp){ positions[i].nei = temp;}
		void setREI(int i, double temp){ positions[i].rei = temp;}

		//Raster Propagated Uncertainty and Kalman Getters
		/** Returns the Average Z (depth) of the Points in the list.
		* @return Average Z (depth) of the Points in the list.
		*/
		double getAvgZ0() const { return avgZ0; }

		/** Returns the Minimum Z (depth) of the Points in the list.
		* @return Minimum Z (depth) of the Points in the list.
		*/
		double getMinZ0() const { return minZ; }

		/** Returns the Maximum Z (depth) of the Points in the list.
		* @return Maximum Z (depth) of the Points in the list.
		*/
		double getMaxZ0() const { return maxZ; }

		/** Returns the Average Z (depth) of the Points in the list.
		* @return Average Z (depth) of the Points in the list.
		*/
		double getAvgZK() const { return avgZK; }

		/** Returns the Minimum Z (depth) of the Points in the list.
		* @return Minimum Z (depth) of the Points in the list.
		*/
		double getMinZK() const { return minZ; }

		/** Returns the Maximum Z (depth) of the Points in the list.
		* @return Maximum Z (depth) of the Points in the list.
		*/
		double getMaxZK() const { return maxZ; }

		//Raster Propagated Uncertainty and Kalman Setters
		/** Sets the Average Z (depth) of the Points in the list.
		* @param temp - value to set.
		*/
		void setAvgZ0(double temp){ avgZ0 = temp; }

		/** Sets the Maximum Z (depth) of the Points in the list.
		* @param temp - Double value to set.
		*/
		void setMaxZ0(double temp){ maxZ0 = temp; }

		/** Sets the Minimum Z (depth) of the Points in the list.
		* @param temp - Double value to set.
		*/
		void setMinZ0(double temp){ minZ0 = temp; }

		/** Sets the Average Z (depth) of the Points in the list.
		* @param temp - value to set.
		*/
		void setAvgZK(double temp){ avgZK = temp; }

		/** Sets the Maximum Z (depth) of the Points in the list.
		* @param temp - Double value to set.
		*/
		void setMaxZK(double temp){ maxZK = temp; }

		/** Sets the Minimum Z (depth) of the Points in the list.
		* @param temp - Double value to set.
		*/
		void setMinZK(double temp){ minZK = temp; }
		//raster**********************************************************






		/** Overloads the [] operator to allow for indexed access to list elements.
		* @param index - Integer 0 - this.size()-1 .
		* @return Point in list at position index.
		*/
		Point operator[](const int index) const;

		/** Randomizes the list */
		void rand();

		/** Sorts the list */
		void sort();

		/** Sorts radially by distance from p */
		void sort(const Point& p);

		/** Removes duplicate elements (should be sorted first) */
		void unique();

		/** Drops the last x percent of Points in the list.
		* @param x - Percentage of the Points to be dropped.
		*/
		void decimate(double x);

		/** Clears the list */
		void clear();

		/** Prints the stats of the list */
		void printStats();

		//void reserve(const int size);

		//void push_back(const Point p);

		//void pop_back();

		//void calcStats();

		//double getMinX() const { return MinX; }

		//double getMinY() const { return MinY; }

};
