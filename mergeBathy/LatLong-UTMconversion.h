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
