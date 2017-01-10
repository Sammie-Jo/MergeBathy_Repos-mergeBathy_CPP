/**
* @file			gsfReader.h
* @brief		Act as the front end to the GMT File Reader routines to make them easier to use.
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/

#ifndef GSFREADER_H
#define GSFREADER_H
#include <stdio.h>
#include <math.h>
#include "gsf.h"
#include "ellipsoid.h"	/* for AXIS(WE) and RFLAT(WE) */
#include "geod.h"	/* for geo_direct prototype */
#include "geodesic.h"	/* for newgp prototype */

#define NONSENSE	1.0e10

//Open a gsf file and get the file handle 
int open_gsffile(char *filename);

//Read an independent record from the file
int read_rec(int gsfHandle, gsfRecords *gsfRec, gsfDataID *id);

//Locate a beam from within a record
void place_beam2(int prime_meridian, int beamon,
			double *lonmin, double *latmin, double *depth,
			double *h, double *v, gsfRecords gsfRec, gsfDataID id);

//Parse independent beams to locate viable data
int gsfReader(double *lonmin, double *latmin, double *depth, double *h, double *v,
		int prime_meridian, int gsfHandle, int *beamon, int *stat, gsfRecords *gsfRec, gsfDataID *id);

//Close a gsf file
void close_gsffile(int *handle);

#endif

