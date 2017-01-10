/*      file: gsfReader.c
        function: gsfReader
        purpose: to read in lonmin-latmin-depth triples from gmf files.
        date: 22 Nov 94
	history:  originally "readgmf" converted to "readgsf" for
		HMPS97 Improvements, Nov 1997
        author: David H. Fabre
	additional history: reader and associated functions were taken
		from "readgsf" to be used in the mergeBathy software package.  
		Some of the fucntions were changed to allow integration with 
		mergeBathy software
		modified by: Kevin M. Duvieilh


                       The application programmer must supply a function
                       'reader' that reads the data in the input files.
                       The argument list for the function is as follows :

   int reader (float *lon, float *lat, float *zvalue, int prime_meridian,
       char *input_filenames[], int numfiles)

       lon             -   longitude value (in minutes, west negative) to
                           pass to caller
       lat             -   latitude value (in minutes, south negative) to
                           pass to caller
       zvalue          -   z value (user units) to pass to caller
       prime_meridian  -   passed to reader, 1 if area covered crosses the
                           prime_meridian
       input_filenames -   passed to reader, array of pointers to input
                           filenames to be read
       numfiles        -   passed to reader, number of input filenames

       int returned    -   0 on a good read, 1 if all data has been read

        (the last few lines of the reader_gen.c (general reader) are:)

    *lon *= 60.0;       (these convert deg to mins)
    *lat *= 60.0;

    if (prime_meridian) *lon += 21600.0; (this adds 360 deg to the lon)

    return (0);

*/

//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	#include "../WarningStates.h"	//Disable all Warnings!!!
#endif

#include "gsfReader.h"

/******************************************************************************/
int open_gsffile(char *filename)
/* open a file */
{
	int gsfHandle;
	if ( gsfOpen(filename, GSF_READONLY, &gsfHandle) == -1 )
	{
		printf("pblm opening <%s>\n", filename);
		return -1;
	}	
	else
	{
		gsfSeek(gsfHandle, GSF_REWIND); /* just in case */
		return gsfHandle;
	}

} /* open_gsffile */


/******************************************************************************/
int read_rec(int gsfHandle, gsfRecords *gsfRec, gsfDataID *id)
{
    int stat;
	int i;

	while ( (stat = gsfRead(gsfHandle, GSF_NEXT_RECORD,
				id, gsfRec, NULL, 0)) != -1 
		&& (*id).recordID != GSF_RECORD_SWATH_BATHYMETRY_PING );
	return (stat != -1 && (*id).recordID == GSF_RECORD_SWATH_BATHYMETRY_PING);

} /* read_rec */

/******************************************************************************/
void place_beam2(int prime_meridian, int beamon,
                        double *lonmin, double *latmin, double *depth,
			double *h, double *v, gsfRecords gsfRec, gsfDataID id)
{
	static ELLIPSOID ellip=WE;
	short stat=1;
        double lat, lon, az2;
	int special;

	if ( (gsfRec.mb_ping.ping_flags & GSF_IGNORE_PING)
		|| (gsfRec.mb_ping.beam_flags[beamon] & GSF_IGNORE_BEAM) )
	{
		*depth = NONSENSE;
	}
	else
	{

		special = 0;
		if ( gsfRec.mb_ping.along_track[beamon]  == 0.0
			&& gsfRec.mb_ping.across_track[beamon] == 0.0)
		{
			/* very special case - probably single beam */
			lat = gsfRec.mb_ping.latitude;
			lon = gsfRec.mb_ping.longitude;
			special = 1;
		}

		if (fabs(gsfRec.mb_ping.along_track[beamon]) > 0.0)
			newgp( gsfRec.mb_ping.latitude,
					gsfRec.mb_ping.longitude,
					gsfRec.mb_ping.heading,
					gsfRec.mb_ping.along_track[beamon],
					&lat, &lon );
		stat = 0;
		if (stat == 0 && fabs(gsfRec.mb_ping.across_track[beamon])>0.0)
			newgp( lat, lon, gsfRec.mb_ping.heading + 90.0,
					gsfRec.mb_ping.across_track[beamon],
					&lat, &lon );
		stat = 0;

		if (stat == 0 || special) /* it worked or special */
		{
			*lonmin = 60.0*lon;	/* convert to min */
			if (prime_meridian) *lonmin += 21600.0;/* +360 deg min*/
			*latmin = 60.0*lat;	/* convert to min */
			*depth = gsfRec.mb_ping.depth[beamon];
		}
		else
		{
			*depth = NONSENSE;
		}
	}

} /* place_beam2 */

/******************************************************************************/
int gsfReader(double *lonmin, double *latmin, double *depth, double *h, double *v,
		int prime_meridian, int gsfHandle, int *beamon, int *stat, gsfRecords *gsfRec, gsfDataID *id)
/*
	notes:
	1. this routine is a bit convoluted because i had to use the single
	reader for a list of filenames.  it would fit better if the higher level
	routine would have done the file opening and passed in a file ptr,
	maybe.
	
	2. this routine is keyed upon the beamon flag.  it should only be -2 the
	first time through.  thereafter it'll be -1, 0, ..., number_beams-1.
	the -1 case indicates that we need to read another record cause we just
	opened the next file on the previous attempt to return.  the 0 thru
	number_beams-1 set indicates just to take the next beamon.
	when it gets to be number_beams we try to read another record.  if
	the stat of the read is invalid we try to open another file.  if the
	file opening doesn't work return a 1 indicating that we're through.

	3. by design we'll be returning bad lat,lon pairs whenever we open a 
	new file after the first one.  the main routine has a gross check
	that should see these values as nonsensical and consider it a bad
	point.

	4. an invalid depth value is set to NONSENSE, which should be non-
	sensical enough to be ignored by the main routine .
*/
{
	int ans=0;

	if ( (*beamon) == -2 && ((*stat) = read_rec(gsfHandle, gsfRec, id)) )
	{
		(*beamon) = 0;
	}
	else if ( (*beamon) == -1 &&
		((*stat) = read_rec(gsfHandle, gsfRec, id)) )
	{
		(*beamon) = 0;
	}
	else if ((*beamon) < (*gsfRec).mb_ping.number_beams-1)
	{
		++(*beamon);
	}	
	else if ( ((*stat) = read_rec(gsfHandle, gsfRec, id)) ) 
	{
		(*beamon) = 0;
	}	

	if (!(*stat))
	{
		ans=1;
	}	
	else
	{
		place_beam2(prime_meridian, *beamon, lonmin, latmin, depth, h, v, *gsfRec, *id);	
	}

	if (*depth == NONSENSE)
		*lonmin = *latmin = *h = *v = NONSENSE;

	return ans;

} /* reader */

void close_gsffile(int *handle)
{
	gsfClose((*handle));
}

//Restore warning state -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( pop )
#endif 