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
#ifndef _ELLIPSOID_H
#define _ELLIPSOID_H

typedef enum
{
AA, AM, 
AN, 
BN, BR, 
CC, CD, 
EA, EB, EC, ED, EE, EF, 
FA, FB, FI, 
HA, 
HE, 
HO, 
ID, 
IN, 
KA, 
RF, 
SA, 
WA, WB, WD, WE 
} ELLIPSOID;

struct XELLIP_PARAMS {	const char *name, *abbrv; double axis, rflat;  }; //Added const for UNIX sam

#define NAMEP(k)		ELLIPSOID_LIST[k].name
#define ABBRVP(k)		ELLIPSOID_LIST[k].abbrv
#define AXIS(k)			ELLIPSOID_LIST[k].axis
#define RFLAT(k)		ELLIPSOID_LIST[k].rflat

static struct XELLIP_PARAMS ELLIPSOID_LIST[] = 
{

"Airy 1830",			"AA", 	 6377563.396,	299.3249646, 
"Modified Airy",		"AM", 	 6377340.189,	299.3249646, 
"Australian National",		"AN",	 6378160.000,	298.25,      
"Bessel 1841 (Namibia)",	"BN",	 6377483.865,	299.1528128, 
"Bessel 1841",			"BR",	 6377397.155,	299.1528128, 
"Clarke 1866",			"CC",	 6378206.4,	294.9786982, 
"Clarke 1880",			"CD",	 6378249.145,	293.465,     
"Everest 1830",			"EA",	 6377276.345,	300.8017,    
"Everest (Sabah & Sarawak)",	"EB",	 6377298.556,	300.8017,    
"Everest 1956",			"EC",	 6377301.243,	300.8017, 
"Everest 1969",			"ED",	 6377295.664,	300.8017,    
"Everest 1948",			"EE",	 6377304.063,	300.8017,    
"Everest (Pakistan)",		"EF",	 6377309.613,	300.8017,    
"Modified Fischer 1960",	"FA",	 6378155.000,	298.3,       
"Fischer 1960",			"FB",	 6378166.000,	298.3,       
"Fischer 1968",			"FI",	 6378150.000,	298.3,       
"Hayford",			"HA",	 6378388.000,	297.0,       
"Helmert 1906",			"HE",	 6378200.000,	298.3,       
"Hough 1960",			"HO",	 6378270.000,	297.0,       
"Indonesian 1974",		"ID",	 6378160.000,	298.247,	 
"International",		"IN",	 6378388.000,	297.0,       
"Krassovsky",			"KA",	 6378245.000,	298.3,       
"GRS 1980",			"RF",	 6378137.000,	298.257222101,
"South American 1969",		"SA",	 6378160.000,	298.25,      
"WGS 60",			"WA",	 6378165.000,	298.3,       
"WGS 66",			"WB",	 6378145.000,	298.25,      
"WGS 72",			"WD",	 6378135.000,	298.26,      
"WGS 84",			"WE",	 6378137.000,	298.257223563

}; 

#endif /* _ELLIPSOID_H */
