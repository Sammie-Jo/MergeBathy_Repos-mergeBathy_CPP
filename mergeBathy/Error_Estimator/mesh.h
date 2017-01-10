#pragma once
#include "fstream"
#include <vector>
#include "geom.h"
#include "pointList.h"

/** Class for converting Triangles from other data structures to a generic Mesh which
* can be exported and viewed in software such as Matlab.
*/
class Mesh
{
	private:	
		/** List of Positions. */
		std::vector<Point> positions;		

		/** List of Indices. */
		std::vector<int> indices;

	public:
		/** Constructor for Mesh object, does nothing. */
		Mesh(){}

		/** Returns the entire list of positions.
		* @return List of positions.
		*/
		std::vector<Point> getPositions() const;
		
		/** Returns the entire list of indices.
		* @return List of indices.
		*/
		std::vector<int> getIndices() const;
		
		/** Inserts all the indices for the Triangle t.
		* Assumes positions have already been inserted.
		* @param t - Triangle to be inserted into indices.
		*/
		void insertIndices(const Triangle& t);

		/** Inserts indices for the Point p.
		* Assumes position for p has already been inserted.
		* @param p - Point to be inserted into indices.
		*/
		void insertIndices(const Point& p);

		/** Inserts all the Points from p as positions in the Mesh.
		* @param p - PointList to be inserted into the Mesh.
		*/
		void insertPoints(PointList& p);
		
		/** Inserts all the Points from p as positions in the Mesh.
		* @param p - list of Points to be inserted into the Mesh.
		*/
		void insertPoints(std::list<Point>& p);
		
		/** Inserts all the Points from p as positions in the Mesh.
		* @param p - vector of Points to be inserted into the Mesh.
		*/
		void insertPoints(std::vector<Point>& p);

		/** Clears the mesh completely. */
		void clear();
		
		/** Returns the size of the mesh in number of triangles.
		* @return Number of triangles in the Mesh.
		*/
		int size() const;
		
		/** Overloads the [] operator to return the Triangle of the 
		* index as long as it is a valid value from 0 - mesh.size().
		* @return Returns Triangle of index.
		*/
		Triangle operator[](const int& index);

		/** Writes the Mesh Positions and indices to two separate files named fileName.positions and fileName.indices.
		* @param fileName - base name for both position and indicies files to be written too.
		*/
		void write(std::string fileName) const;
};

