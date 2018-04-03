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

