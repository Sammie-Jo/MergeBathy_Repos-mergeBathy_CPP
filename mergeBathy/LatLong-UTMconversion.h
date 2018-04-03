/**********************************************************************
* CC0 License
**********************************************************************
* MergeBathy - Tool to combine one or more bathymetric data files onto a single input grid.
* Modified in 2015 by Samantha J.Zambo(samantha.zambo@gmail.com) while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Todd Holland while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Nathaniel Plant while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Kevin Duvieilh while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Paul Elmore while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Will Avera while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by Brian Bourgeois while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by A.Louise Perkins while employed by the U.S.Naval Research Laboratory.
* Modified in 2015 by David Lalejini while employed by the U.S.Naval Research Laboratory.
* To the extent possible under law, the author(s) and the U.S.Naval Research Laboratory have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide.This software is distributed without any warranty.
* You should have received a copy of the CC0 Public Domain Dedication along with this software.If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
**********************************************************************/
/**
* @file			LatLong-UTMconversion.h
* @brief		Definitions for lat/long to UTM and UTM to lat/lng conversions.
*
*/

//LatLong- UTM conversion..h
//definitions for lat/long to UTM and UTM to lat/lng conversions
//#ifndef LATLONGCONV_H_
//#define LATLONGCONV_H_
#pragma once 
#include <string>
#include <cstring> //for UNIX sam
/*
converts from lat lon to UTM coordinates
converts from UTM coordinates to lat lon
*/

//Convert from Lat/Lon coordinates to UTM coordinates
void LLtoUTM(int ReferenceEllipsoid, const double Lat, const double Long, 
			 double &UTMNorthing, double &UTMEasting, char* UTMZone, const double meanLong);

//Convert from UTM coordinates to Lat/Lon coordinates 
void UTMtoLL(int ReferenceEllipsoid, const double UTMNorthing, const double UTMEasting, const char* UTMZone,
			  double& Lat,  double& Long );

char UTMLetterDesignator(double Lat);

//Define the basic ellipsoid class for use in converting between coordinate systems.
class Ellipsoid
{
public:
	Ellipsoid(){};
	Ellipsoid(int Id, char* name, double radius, double ecc)
	{
		id = Id; ellipsoidName = name; 
		EquatorialRadius = radius; eccentricitySquared = ecc;
	}

	int id;
	char* ellipsoidName;
	double EquatorialRadius; 
	double eccentricitySquared;  

};



//#endif
