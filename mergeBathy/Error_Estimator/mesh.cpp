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
/**
*
*
*/
#pragma once

#ifndef MESH_CPP
#define MESH_CPP

#include "fstream"
#include <vector>
#include <list>
#include <algorithm>
#include <cassert>
#include <iostream>
#include "mesh.h"
#include "geom.h"
#include "pointList.h"


using namespace std;

// Returns the entire list of positions
vector<Point> Mesh::getPositions() const { return positions; }
		
// Returns the entire list of indices
vector<int> Mesh::getIndices() const { return indices; }
		
// calls createIndicesFromPoint for all pings in the triangle
void Mesh::insertIndices(const Triangle& t)
{
	insertIndices(t.v1);
	insertIndices(t.v2);
	insertIndices(t.v3);
}

// Requires that positions for all triangles being inserted already be present
// creates the indices for the triangle and inserts them
void Mesh::insertIndices(const Point& p)
{
	int index, indexMove;
	indexMove = (const int)positions.size() / 4;
	index = (const int)positions.size() / 2;
	while(true)
	{
		if(positions[index] == p)
			break;
		else if(positions[index] < p)
			index += indexMove;
		else
			index -= indexMove;
		if(indexMove != 1)
			indexMove = indexMove / 2;
		if(indexMove == 0)
			indexMove++;
	}

	indices.push_back(index);
}

void Mesh::insertPoints(list<Point>& p)
{
	p.sort();
	p.unique();
	positions.reserve(p.size());
	positions.assign(p.begin(), p.end());
}

// inserts positions from a pre-sorted list (MUCH FASTER)
void Mesh::insertPoints(PointList& p)
{
	p.sort();
	p.unique();
	positions.reserve(p.size());
	for(int i = 0; i < p.size(); i++)
		positions.push_back(p[i]);
}
		
// inserts positions from a pre-sorted list (MUCH FASTER)
void Mesh::insertPoints(vector<Point>& p)
{
	sort(p.begin(), p.end());
	unique(p.begin(), p.end());
	positions = p;
}

// clears the mesh completely
void Mesh::clear() 
{ 
	positions.clear(); 
	indices.clear();
}
		
// Returns the size of the mesh in number of  triangles
int Mesh::size() const { return (const int)indices.size() / 3; }
		
// Operator for indexing specific triangles using []
Triangle Mesh::operator[](const int& index)
{
	Point p(positions[indices[index*3]]);
	Point v1(p);
	p = positions[indices[index*3+1]];
	Point v2(p);
	p = positions[indices[index*3+2]];
	Point v3(p);
	Triangle t(v1, v2, v3);
	return t;
}

void Mesh::write(string fileName) const
{
	Point p;
	ofstream outData;
	string fn = fileName + ".positions";
	outData.open(fn.c_str());
	int k = 0;
	for(int i = 0; i < (const int)positions.size(); i++)
	{
			p = positions[i];
			outData << p.x << " " << p.y << " " << p.z << endl;
	}
	outData.close();
	fn = fileName + ".indices";
	outData.open(fn.c_str());
	for(int i = 0; i < (const int)indices.size(); i++)
	{
		outData << indices[i++] << " ";
		outData << indices[i++] << " ";
		outData << indices[i] << endl;
	}
}

#endif
