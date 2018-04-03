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
//LatLong- UTM conversion.cpp
//Lat Long - UTM, UTM - Lat Long conversions

//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	#include "WarningStates.h"		//Disable all Warnings!!!
#endif

#include "LatLong-UTMconversion.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "constants.h"
#include <iostream>

/*Reference ellipsoids derived from Peter H. Dana's website-
http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
Department of Geography, University of Texas at Austin
Internet: pdana@mail.utexas.edu
3/22/95

Source
Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to Department of Defense World Geodetic System
1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency
*/

void LLtoUTM(int ReferenceEllipsoid, const double Lat, const double Long, double &UTMNorthing, double &UTMEasting, char* UTMZone, const double meanLong)
{
	//converts lat/long to UTM coords.  Equations from USGS Bulletin 1532
	//East Longitudes are positive, West longitudes are negative.
	//North latitudes are positive, South latitudes are negative
	//Lat and Long are in decimal degrees
	//Written by Chuck Gantz- chuck.gantz@globalstar.com
	Ellipsoid ellipsoid[] =
	{//  id, Ellipsoid name, Equatorial Radius, square of eccentricity
	Ellipsoid( -1, "Placeholder", 0, 0),//placeholder only, To allow array indices to match id numbers
	Ellipsoid( 1, "Airy", 6377563, 0.00667054),
	Ellipsoid( 2, "Australian National", 6378160, 0.006694542),
	Ellipsoid( 3, "Bessel 1841", 6377397, 0.006674372),
	Ellipsoid( 4, "Bessel 1841 (Nambia) ", 6377484, 0.006674372),
	Ellipsoid( 5, "Clarke 1866", 6378206, 0.006768658),
	Ellipsoid( 6, "Clarke 1880", 6378249, 0.006803511),
	Ellipsoid( 7, "Everest", 6377276, 0.006637847),
	Ellipsoid( 8, "Fischer 1960 (Mercury) ", 6378166, 0.006693422),
	Ellipsoid( 9, "Fischer 1968", 6378150, 0.006693422),
	Ellipsoid( 10, "GRS 1967", 6378160, 0.006694605),
	Ellipsoid( 11, "GRS 1980", 6378137, 0.00669438),
	Ellipsoid( 12, "Helmert 1906", 6378200, 0.006693422),
	Ellipsoid( 13, "Hough", 6378270, 0.00672267),
	Ellipsoid( 14, "International", 6378388, 0.00672267),
	Ellipsoid( 15, "Krassovsky", 6378245, 0.006693422),
	Ellipsoid( 16, "Modified Airy", 6377340, 0.00667054),
	Ellipsoid( 17, "Modified Everest", 6377304, 0.006637847),
	Ellipsoid( 18, "Modified Fischer 1960", 6378155, 0.006693422),
	Ellipsoid( 19, "South American 1969", 6378160, 0.006694542),
	Ellipsoid( 20, "WGS 60", 6378165, 0.006693422),
	Ellipsoid( 21, "WGS 66", 6378145, 0.006694542),
	Ellipsoid( 22, "WGS-72", 6378135, 0.006694318),
	Ellipsoid( 23, "WGS-84", 6378137, 0.00669437999014)
	};

	double a = ellipsoid[ReferenceEllipsoid].EquatorialRadius;
	double eccSquared = ellipsoid[ReferenceEllipsoid].eccentricitySquared;
	double k0 = 0.9996;

	double LongOrigin;
	double eccPrimeSquared;
	double N, T, C, A, M;


	//Make sure the longitude is between -180.00 .. 179.9
	double LongTemp = (Long);
//	double LongTemp = (Long+180.00)-int((Long+180.00)/360.00)*360.00-180.00; // -180.00 .. 179.9; // SJZ
	//std::cout << "LT: " << LongTemp << std::endl;

	double  LatRad = Lat*deg2rad;
	double  LongRad = LongTemp*deg2rad;
	double  LongOriginRad;
	int     ZoneNumber	= 0;
	int		tempLoop	= 0;
	char zChar;
	//Find UTM Zone
	if(*UTMZone!=NULL) // SJZ
	{
		char ZoneNum[4];
		int end=strlen(UTMZone)-1;
		strncpy(ZoneNum,UTMZone,end);
		ZoneNum[end] = '\0';
		ZoneNumber=(int)strtol(ZoneNum,(char **)NULL, 10);
		zChar=(UTMZone)[end];
	}
	else
	{// SJZ
		
	for (tempLoop = -180; tempLoop < 180; tempLoop += 6) // SJZ <=
	//for (tempLoop = 0; tempLoop <= 360; tempLoop += 6)
	{
		if (tempLoop <= LongTemp)//meanLong) // SJZ
			ZoneNumber = ZoneNumber + 1;
		else
			break;
	}

	if( Lat >= 56.00 && Lat < 64.00 && LongTemp >= 3.00 && LongTemp < 12.00 )
		ZoneNumber = 32;

	// Special zones for Svalbard
	if( Lat >= 72.00 && Lat < 84.00 )
	{
		if	   ( LongTemp >= 0.00  && LongTemp <  9.00 ) ZoneNumber = 31;
		else if( LongTemp >= 9.00  && LongTemp < 21.00 ) ZoneNumber = 33;
		else if( LongTemp >= 21.00 && LongTemp < 33.00 ) ZoneNumber = 35;
		else if( LongTemp >= 33.00 && LongTemp < 42.00 ) ZoneNumber = 37;
	}
	zChar = UTMLetterDesignator(Lat);
	}// SJZ
	LongOrigin = (ZoneNumber - 1)*6.00 - 180.00 + 3.00;  //+3 puts origin in middle of zone
	//printf("ZN ECS: %d %f\n",ZoneNumber, eccSquared);
	//printf("ORIG: %f\n",LongOrigin);
	LongOriginRad = LongOrigin * deg2rad;

	//compute the UTM Zone from the latitude and longitude
	sprintf(UTMZone, "%d%c", ZoneNumber, UTMLetterDesignator(Lat));
	sprintf(UTMZone, "%d%c", ZoneNumber, zChar);
	//end UTM Zone

	eccPrimeSquared = (eccSquared)/(1.00-eccSquared);

	N = a/sqrt(1.00-eccSquared*std::sin(LatRad)*std::sin(LatRad));
	T = std::tan(LatRad)*std::tan(LatRad);
	C = eccPrimeSquared*cos(LatRad)*std::cos(LatRad);
	A = std::cos(LatRad)*(LongRad-LongOriginRad);

	M = a*((1.00	- eccSquared/4.00		- 3.00*eccSquared*eccSquared/64.00	- 5.00*eccSquared*eccSquared*eccSquared/256.00)*LatRad
		- (3.00*eccSquared/8.00	+ 3.00*eccSquared*eccSquared/32.00	+ 45.00*eccSquared*eccSquared*eccSquared/1024.00)*std::sin(2.00*LatRad)
									+ (15.00*eccSquared*eccSquared/256.00 + 45.00*eccSquared*eccSquared*eccSquared/1024.00)*std::sin(4.00*LatRad)
									- (35.00*eccSquared*eccSquared*eccSquared/3072.00)*std::sin(6.00*LatRad));

	UTMEasting = (double)(k0*N*(A+(1.00-T+C)*A*A*A/6.00
					+ (5.00-18.00*T+T*T+72.00*C-58.00*eccPrimeSquared)*A*A*A*A*A/120.00)
					+ 500000.00); // 500000 False Easting at Origin

	UTMNorthing = (double)(k0*(M+N*tan(LatRad)*(A*A/2.00+(5.00-T+9.00*C+4.00*C*C)*A*A*A*A/24.00
				 + (61.00-58.00*T+T*T+600.00*C-330.00*eccPrimeSquared)*A*A*A*A*A*A/720.00)));
	
	/* Added 10/24/14 SJZ
	* Changed so we properly handle cases where the data spans the equator but
	* the user wants data south of the equator projected to a zone north of the
	* equator. (DML) 
	*/
	
	/* Old code: */
	//if(Lat < 0)
	//	UTMNorthing += 10000000.0; //10000000 meter offset for southern hemisphere

	/* New code (DML): */
//	char zChar = UTMLetterDesignator(Lat);
	if (zChar=='9')
		std::cout << "Assuming Northern Hemisphere." << std::endl;
	else if (zChar=='Z')//ignore this if using -inputinMeters
		std::cout << "Error: Latitude is outside the UTM limits. Ignore if already in UTM." << std::endl;
	else
	{
		if(zChar <= 'M') //southern hemisphere
			UTMNorthing += 10000000.0; // 10000000 meter offset for southern hemisphere - aka False Northing at Origin
	}
}

char UTMLetterDesignator(double Lat)
{
	//This routine determines the correct UTM letter designator for the given latitude
	//returns 'Z' if latitude is outside the UTM limits of 84N to 80S
	//Written by Chuck Gantz- chuck.gantz@globalstar.com
	char LetterDesignator;

	if((84 >= Lat) && (Lat >= 72)) LetterDesignator = 'X';
	else if((72 > Lat) && (Lat >= 64)) LetterDesignator = 'W';
	else if((64 > Lat) && (Lat >= 56)) LetterDesignator = 'V';
	else if((56 > Lat) && (Lat >= 48)) LetterDesignator = 'U';
	else if((48 > Lat) && (Lat >= 40)) LetterDesignator = 'T';
	else if((40 > Lat) && (Lat >= 32)) LetterDesignator = 'S';
	else if((32 > Lat) && (Lat >= 24)) LetterDesignator = 'R';
	else if((24 > Lat) && (Lat >= 16)) LetterDesignator = 'Q';
	else if((16 > Lat) && (Lat >= 8)) LetterDesignator = 'P';
	else if(( 8 > Lat) && (Lat >= 0)) LetterDesignator = 'N';
	else if(( 0 > Lat) && (Lat >= -8)) LetterDesignator = 'M';
	else if((-8> Lat) && (Lat >= -16)) LetterDesignator = 'L';
	else if((-16 > Lat) && (Lat >= -24)) LetterDesignator = 'K';
	else if((-24 > Lat) && (Lat >= -32)) LetterDesignator = 'J';
	else if((-32 > Lat) && (Lat >= -40)) LetterDesignator = 'H';
	else if((-40 > Lat) && (Lat >= -48)) LetterDesignator = 'G';
	else if((-48 > Lat) && (Lat >= -56)) LetterDesignator = 'F';
	else if((-56 > Lat) && (Lat >= -64)) LetterDesignator = 'E';
	else if((-64 > Lat) && (Lat >= -72)) LetterDesignator = 'D';
	else if((-72 > Lat) && (Lat >= -80)) LetterDesignator = 'C';
	else LetterDesignator = 'Z'; //This is here as an error flag to show that the Latitude is outside the UTM limits

	return LetterDesignator;
}

void UTMtoLL(int ReferenceEllipsoid, const double UTMNorthing, const double UTMEasting, const char* UTMZone,
			  double& Lat,  double& Long )
{
	//converts UTM coords to lat/long.  Equations from USGS Bulletin 1532
	//East Longitudes are positive, West longitudes are negative.
	//North latitudes are positive, South latitudes are negative
	//Lat and Long are in decimal degrees.
	//Written by Chuck Gantz- chuck.gantz@globalstar.com
	Ellipsoid ellipsoid[] =
	{//  id, Ellipsoid name, Equatorial Radius, square of eccentricity
	Ellipsoid( -1, "Placeholder", 0, 0),//placeholder only, To allow array indices to match id numbers
	Ellipsoid( 1, "Airy", 6377563, 0.00667054),
	Ellipsoid( 2, "Australian National", 6378160, 0.006694542),
	Ellipsoid( 3, "Bessel 1841", 6377397, 0.006674372),
	Ellipsoid( 4, "Bessel 1841 (Nambia) ", 6377484, 0.006674372),
	Ellipsoid( 5, "Clarke 1866", 6378206, 0.006768658),
	Ellipsoid( 6, "Clarke 1880", 6378249, 0.006803511),
	Ellipsoid( 7, "Everest", 6377276, 0.006637847),
	Ellipsoid( 8, "Fischer 1960 (Mercury) ", 6378166, 0.006693422),
	Ellipsoid( 9, "Fischer 1968", 6378150, 0.006693422),
	Ellipsoid( 10, "GRS 1967", 6378160, 0.006694605),
	Ellipsoid( 11, "GRS 1980", 6378137, 0.00669438),
	Ellipsoid( 12, "Helmert 1906", 6378200, 0.006693422),
	Ellipsoid( 13, "Hough", 6378270, 0.00672267),
	Ellipsoid( 14, "International", 6378388, 0.00672267),
	Ellipsoid( 15, "Krassovsky", 6378245, 0.006693422),
	Ellipsoid( 16, "Modified Airy", 6377340, 0.00667054),
	Ellipsoid( 17, "Modified Everest", 6377304, 0.006637847),
	Ellipsoid( 18, "Modified Fischer 1960", 6378155, 0.006693422),
	Ellipsoid( 19, "South American 1969", 6378160, 0.006694542),
	Ellipsoid( 20, "WGS 60", 6378165, 0.006693422),
	Ellipsoid( 21, "WGS 66", 6378145, 0.006694542),
	Ellipsoid( 22, "WGS-72", 6378135, 0.006694318),
	Ellipsoid( 23, "WGS-84", 6378137, 0.00669437999014)
	};
	double k0 = 0.9996;
	double a = ellipsoid[ReferenceEllipsoid].EquatorialRadius;
	double eccSquared = ellipsoid[ReferenceEllipsoid].eccentricitySquared;
	double eccPrimeSquared;
	double e1 = (1.00-sqrt(1.00-eccSquared))/(1.00+sqrt(1.00-eccSquared));
	double N1, T1, C1, R1, D, M;
	double LongOrigin;
	double mu, phi1, phi1Rad;
	double x, y;
	int ZoneNumber;
	char* ZoneLetter;
	int NorthernHemisphere; //1 for northern hemisphere, 0 for southern

	x = UTMEasting - 500000.00; //remove 500,000 meter offset for longitude
	y = UTMNorthing;

	ZoneNumber = strtoul(UTMZone, &ZoneLetter, 10);
	if((*ZoneLetter - 'N') >= 0)
		NorthernHemisphere = 1;//point is in northern hemisphere
	else
	{
		NorthernHemisphere = 0;//point is in southern hemisphere
		y -= 10000000.0;//remove 10,000,000 meter offset used for southern hemisphere
	}

	LongOrigin = (ZoneNumber - 1)*6.00 - 180.00 + 3.00;  //+3 puts origin in middle of zone

	eccPrimeSquared = (eccSquared)/(1.00-eccSquared);

	M = y / k0;
	mu = M/(a*(1.00-eccSquared/4.00-3.00*eccSquared*eccSquared/64.00-5.00*eccSquared*eccSquared*eccSquared/256.00));

	phi1Rad = mu	+ (3.00*e1/2.00-27.00*e1*e1*e1/32.00)*std::sin(2.00*mu)
				+ (21.00*e1*e1/16.00-55.00*e1*e1*e1*e1/32.00)*std::sin(4.00*mu)
				+(151.00*e1*e1*e1/96.00)*std::sin(6.00*mu);
	phi1 = phi1Rad*rad2deg;

	N1 = a/sqrt(1.00-eccSquared*std::sin(phi1Rad)*std::sin(phi1Rad));
	T1 = std::tan(phi1Rad)*std::tan(phi1Rad);
	C1 = eccPrimeSquared*std::cos(phi1Rad)*std::cos(phi1Rad);
	R1 = a*(1.00-eccSquared)/pow(1.00-eccSquared*std::sin(phi1Rad)*std::sin(phi1Rad), 1.50);
	D = x/(N1*k0);

	Lat = phi1Rad - (N1*std::tan(phi1Rad)/R1)*(D*D/2-(5.00+3.00*T1+10.00*C1-4.00*C1*C1-9.00*eccPrimeSquared)*D*D*D*D/24.00
					+(61.00+90.00*T1+298.00*C1+45.00*T1*T1-252.00*eccPrimeSquared-3.00*C1*C1)*D*D*D*D*D*D/720.00);
	Lat = Lat * rad2deg;

	Long = (D-(1+2*T1+C1)*D*D*D/6.00+(5.00-2.00*C1+28.00*T1-3.00*C1*C1+8.00*eccPrimeSquared+24.00*T1*T1)
					*D*D*D*D*D/120.00)/std::cos(phi1Rad);
	Long = LongOrigin + Long * rad2deg;
}

#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning( pop )			//Restore warning state
#endif
