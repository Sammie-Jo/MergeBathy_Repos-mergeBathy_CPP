/*

    libsurf.a - SURF gridding library

    Author : Tom Gray and Jan Depner

    Date : 10/10/08

    This library contains an implementation of the GMT surface program as a function.  The GMT blockmean
    program has been replaced since it required all of the data to reside in memory.  Some of our grids
    will contain more than 1,000,000,000 input points so, in some cases, we would run out of memory.
    In the following comments we'll refer to our replacement for GMT's blockmean as "blockmean".
    After we "blockmean" the data we pass it to the modified surface function as 2D arrays of data
    containing X, Y, Z, and CNT values for each grid.  The following code consists of four main functions:

        surf_init - initializes all variables and structures needed for the gridding operations
        surf_load - loads all of the user's input data into a sparse grid
        surf_proc - finalizes the "blockmean" operation and, optionally, calls the surface function
        surf_rtrv - retrieves the final gridded data from surface as 2D arrays of Z and CNT
        surf_cleanup - Frees memory after the caller is finished with it.

    If you do not have enough memory to grid the area you have specified you have a few options.  Move
    to another system with more memory, buy more memory, re-write this to do subgridding and edge matching,
    or go find something else to do.

    Happy gridding!
    Tom Gray and Jan Depner

*/


/*Experimental implementation for median (instead of mean) was done because everything I've read points to blockmedian being the more robust, and the way to go, for bathymetry.  But in my quick attempt of a day, there was not enough memory (among other problems I'm sure) to complete a successful run.  Therefore, it was aborted since implementation was tried only to see if it would provide a better result and account for some of the visual artifacts seen. For future reference, they were left, denoted as experimental, and commented out.  GEBCO cookbook states that GMT's blockmedian may not be able to compute the median for large files due to memory.-- SJZ 2/12/16.*/


//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	#include "../WarningStates.h"	//Disable Warnings!!!
#endif

#include "surf.h"
#include "version.h"
#define DEBUG 0

/*  Static variables and arrays.  */

static NV_INT32      rows, cols;
static NV_FLOAT64    x_grid, y_grid;
static NV_FLOAT64    x_orig, y_orig;
static NV_FLOAT64    x_max, y_max;
static NV_FLOAT32    min_z, max_z;
static NV_FLOAT64    tension;
static NV_FLOAT64    *x = NULL, *y = NULL, *z = NULL;
static NV_FLOAT64    *x_clean = NULL, *y_clean = NULL, *z_clean = NULL;
/*
//Experimenting with median for blockmedian
static NV_FLOAT64    **x = NULL, **y = NULL, **z = NULL;
static NV_FLOAT64    **x_clean = NULL, **y_clean = NULL, **z_clean = NULL;
*/
static NV_INT32      *cnt = NULL;
static NV_INT32      *cnt_clean = NULL;
static NV_INT32      num_pts = 0;
static NV_INT32      surf_err = 0;
static NV_INT32      surf_errno = 0;
static NV_F64_COORD3 coord3;
static NV_INT32      out_of_bounds = 0;
static NV_INT32      process_state = 0;

/*

  Function:         surf_init         -  Initialize variables needed for "blockmean" and surface.  Also allocate
                                         and initialize the X, Y, Z, and CNT arrays.

  Arguments:        mbr               -  Minimum bounding rectangle for area to be gridded
                    x_interval        -  X grid interval
                    y_interval        -  Y grid interval
		    minz              -  Minimum allowable Z vaue
		    maxz              -  Maximum allowable Z vaue
		    tension_value     -  Tension setting for surface (0.0 - 1.0, if set to 0.0 default is 0.35)

  Returns:          NV_INT32          -  0 if successful, on failure check error by calling get_surf_error
                                         and checking the surf_err_str, and surf_errno values.  See get_surf_error
					 for an explanation of the error values.  Important note - negative return
					 values are catastrophic errors (i.e. you need to display the error and
					 then GTHOOD [Get The Hell Out Of Dodge]).

  Caveats:          The MBR must have maximums that are larger than the minimums.  If you're trying to grid
                    actual lat/lon values and you want to cross the dateline you must use 0-360 world
		    coordinates.  For example, if you want to grid from 160E to 160W you must specify
		    160 to 200 not 160 to -160.

*/

NV_INT32 surf_init (NV_F64_XYMBR mbr, NV_FLOAT64 x_interval, NV_FLOAT64 y_interval, NV_FLOAT32 minz, NV_FLOAT32 maxz,
		    NV_FLOAT64 tension_value)

{
	NV_FLOAT32        temp;
	//Experimenting with median for blockmedian
	//NV_INT32 i;

	/*  Save the important stuff  */

	if(DEBUG)
		printf("Entering surf_init.  mbr.max_x = %lf, mbr_maxy = %lf\n", mbr.max_x, mbr.max_y);

	if (x_interval < 0.0 || y_interval < 0.0)
	{
		surf_err = -1;
		surf_errno = 0;
		return (surf_err);
	}

	x_grid = x_interval;
	y_grid = y_interval;

	if (mbr.max_x < mbr.min_x || mbr.max_y < mbr.min_y)
	{
		surf_err = -2;
		surf_errno = 0;
		return (surf_err);
	}

	x_orig = mbr.min_x;
	y_orig = mbr.min_y;

	if (maxz < minz)
	{
		temp = minz;
		minz = maxz;
		maxz = temp;
	}

	min_z = minz;
	max_z = maxz;

	if (tension_value <= 0.0 || tension_value > 1.0)
	{
		tension = 0.35;
	}
	else
	{
		tension = tension_value;
	}

	/*  Compute rows and columns  */
	cols = NINT(((mbr.max_x - mbr.min_x) / x_interval));
	rows = NINT(((mbr.max_y - mbr.min_y) / y_interval));
	/*
	cols = ceil((mbr.max_x - mbr.min_x) / x_interval);
	rows = ceil((mbr.max_y - mbr.min_y) / y_interval);
	*/

	if(DEBUG)
	{
		printf("surf_init: (mbr.max_x - mbr.min_x) = %f\n", (mbr.max_x - mbr.min_x));
		printf("surf_init: ((mbr.max_x - mbr.min_x) / x_interval) = %f\n",((mbr.max_x - mbr.min_x) / x_interval));
		printf("surf_init:  NINT ((mbr.max_x - mbr.min_x) / x_interval) = %d\n", NINT ((mbr.max_x - mbr.min_x) / x_interval));
		printf("surf_init:  ceil((mbr.max_x - mbr.min_x) / x_interval) = %lf\n", ceil ((mbr.max_x - mbr.min_x) / x_interval));
		printf("surf_init:  NINT ((mbr.max_y - mbr.min_y) / y_interval) = %d\n", NINT ((mbr.max_y - mbr.min_y) / y_interval));
		printf("surf_init:  ceil ((mbr.max_y - mbr.min_y) / y_interval) = %lf\n", ceil ((mbr.max_y - mbr.min_y) / y_interval));
		printf("surf_init:  (int) ((a) + 0.5 = %d\n", (int) ((mbr.max_x - mbr.min_x) / x_interval));
		printf("surf_init: x_interval = %f\n", x_interval);
		printf("surf_init: y_interval = %f\n", x_interval);
		printf("surf_init: cols = %d\n", cols);
		printf("surf_init: rows = %d\n", rows);
	}

	/*  Compute the other side of the "ROI" for surface.  */
	/*
	x_max = x_orig + ((NV_FLOAT64) cols * x_grid);
	y_max = y_orig + ((NV_FLOAT64) rows * y_grid);
	*/
	x_max = mbr.max_x;
	y_max = mbr.max_y;

	if(DEBUG)
	{
		printf("surf_init: cols = %d\n", cols);
		printf("surf_init: rows = %d\n", rows);
		printf("x_max is %lf\n", x_max);
		printf("y_max is %lf\n", y_max);
	}

	/*  Allocate the "row" arrays.  */
	/*
	//Experimenting with median for blockmedian
	x = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	y = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	z = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	//x_clean = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	//y_clean = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	//z_clean = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	for( i=0; i<(rows*cols); i++){
		x[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		y[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		z[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		//x_clean[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		//y_clean[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		//z_clean[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
	}
	*/
	x = (NV_FLOAT64 *) calloc (rows * cols, sizeof (NV_FLOAT64));
	if (x == NULL)
	{
		surf_err = -3;
		surf_errno = errno;
		surf_cleanup();
		return (surf_err);
	}

	x_clean = (NV_FLOAT64 *) calloc (rows * cols, sizeof (NV_FLOAT64));
	if (x_clean == NULL)
	{
		surf_err = -3;
		surf_errno = errno;
		surf_cleanup();
		return (surf_err);
	}

	y = (NV_FLOAT64 *) calloc (rows * cols, sizeof (NV_FLOAT64));
	if (y == NULL)
	{
		surf_err = -4;
		surf_errno = errno;
		surf_cleanup();
		return (surf_err);
	}

	y_clean = (NV_FLOAT64 *) calloc (rows * cols, sizeof (NV_FLOAT64));
	if (y_clean == NULL)
	{
		surf_err = -4;
		surf_errno = errno;
		surf_cleanup();
		return (surf_err);
	}

	z = (NV_FLOAT64 *) calloc ((rows + 1) * (cols + 1), sizeof (NV_FLOAT64));
	if (z == NULL)
	{
		surf_err = -5;
		surf_errno = errno;
		surf_cleanup();
		return (surf_err);
	}

	z_clean = (NV_FLOAT64 *) calloc ((rows + 1) * (cols + 1), sizeof (NV_FLOAT64));
	if (z_clean == NULL)
	{
		surf_err = -5;
		surf_errno = errno;
		surf_cleanup();
		return (surf_err);
	}

	cnt = (NV_INT32 *) calloc (rows * cols, sizeof (NV_INT32));
	if (cnt == NULL)
	{
		surf_err = -6;
		surf_errno = errno;
		surf_cleanup();
		return (surf_err);
	}

	cnt_clean = (NV_INT32 *) calloc (rows * cols, sizeof (NV_INT32));
	if (cnt_clean == NULL)
	{
		surf_err = -6;
		surf_errno = errno;
		surf_cleanup();
		return (surf_err);
	}

	/*  Keep track of where we are in the process.  */
	process_state = 1;
	return (0);
}

/*

  Function:         surf_load         -  Add single points to the sum values in the X, Y, Z, and CNT arrays.

  Arguments:        xyz               -  Coordinate triplet to be added to the X, Y, Z, and CNT arrays.

  Returns:          NV_INT32          -  0 if successful, on failure check error by calling get_surf_error
                                         and checking the surf_err_str, and surf_errno values.  See get_surf_error
					 for an explanation of the error values.  Important note - negative return
					 values are catastrophic errors (i.e. you need to display the error and
					 then GTHOOD [Get The Hell Out Of Dodge]).

  Caveats:          There are a number of entry points for this function.  These allow you to use either separate X,
                    Y, and Z values or one of the NV_ structures.  The calls are as follows:

		    surf_load_XYZ (NV_FLOAT64 x, NV_FLOAT64 y, NV_FLOAT64 z);
		    surf_load_XYz (NV_FLOAT64 x, NV_FLOAT64 y, NV_FLOAT32 z);
		    surf_load_xyz (NV_FLOAT32 x, NV_FLOAT32 y, NV_FLOAT32 z);
		    surf_load_NV_F64_COORD2_Z (NV_F64_COORD2 xy, NV_FLOAT64 z);
		    surf_load_NV_F64_COORD2_z (NV_F64_COORD2 xy, NV_FLOAT32 z);
		    surf_load_NV_F32_COORD2_z (NV_F32_COORD2 xy, NV_FLOAT32 z);
		    surf_load_NV_F32_COORD3 (NV_F32_COORD3 xyz);

		    These are just convenience functions in case your application has the values in these types
		    of variables instead of the normal NV_F64_COORD3.  Hopefully the arguments to the above are
		    intuitively obvious to the most casual observer ;-)

*/

NV_INT32 surf_load (NV_F64_COORD3 xyz)
{
	NV_INT32          row, col, row_ndx;
	NV_INT32 jval,ival;
	/*  Make sure we already called surf_init.  */

	if (!process_state)
	{
		surf_err = -7;
		surf_errno = 0;
		surf_cleanup();
		return (surf_err);
	}

	/*  Compute the row and column for this data point.  */
	if(xyz.y == y_max) /* Want to use points that fall on the border */
		row = rows - 1 ;
	else
		row = (xyz.y - y_orig) / y_grid;

	if(xyz.x == x_max) /* Want to use points that fall on the border */
		col = cols - 1;
	else
		col = (xyz.x - x_orig) / x_grid;

	if(DEBUG)
	{
		printf("surf_load: x value = %lf\n", xyz.x);
		printf("surf_load: y value = %lf\n", xyz.y);
		printf("surf_load: z value = %lf\n", xyz.z);
		printf("surf_load: row = %d\n", row);
		printf("surf_load: rows = %d\n", rows);
		printf("surf_load: col = %d\n", col);
		printf("surf_load: cols = %d\n", cols);
		fflush(stdout);
	}

	/* Make sure z values are in the user-specified range */
	if((xyz.z < min_z)  || (xyz.z > max_z))
		return(0);

	/*  Check for the point being in bounds.  */
	if (row < 0 || row >= rows || col < 0 || col >= cols)
	{
		/* TG - Don't throw an error here. */
		/* surf_err = 1; */
		out_of_bounds++;
		/* return (surf_err); */
		return(0);
	}

	row_ndx = row * cols;
	x[row_ndx + col] += xyz.x;
	y[row_ndx + col] += xyz.y;
	z[row_ndx + col] += xyz.z;
	/*
	//Experimenting with median for blockmedian
	ival=row_ndx + col;
	jval=cnt[row_ndx+col];
	x[ival][jval] = xyz.x;
	y[ival][jval] = xyz.y;
	z[ival][jval] = xyz.z;*/
	cnt[row_ndx + col]++;
	num_pts++;

	if(DEBUG)
	{
		printf("surf_load: row_ndx = %d\n", row_ndx);
		printf("surf_load: x[%d + %d] = %f\n", row_ndx, col, x[row_ndx + col]);
		printf("surf_load: y[%d + %d] = %f\n", row_ndx, col, y[row_ndx + col]);
		printf("surf_load: z[%d + %d] = %f\n", row_ndx, col, z[row_ndx + col]);
		printf("surf_load: cnt[%d + %d] = %d\n", row_ndx, col, cnt[row_ndx + col]);
		printf("cnt[0] = %d\n", cnt[0]);
		printf("cnt[1] = %d\n", cnt[1]);
		printf("cnt[2] = %d\n", cnt[2]);
		printf("cnt[3] = %d\n", cnt[3]);
		fflush(stdout);
	}

	/*  Keep track of where we are in the process.  */
	process_state = 2;
	return (0);
}

NV_INT32 surf_load_XYZ (NV_FLOAT64 x, NV_FLOAT64 y, NV_FLOAT64 z)
{
	coord3.x = x;
	coord3.y = y;
	coord3.z = z;

	return (surf_load (coord3));
}

NV_INT32 surf_load_XYz (NV_FLOAT64 x, NV_FLOAT64 y, NV_FLOAT32 z)
{
	coord3.x = x;
	coord3.y = y;
	coord3.z = (NV_FLOAT64) z;

	return (surf_load (coord3));
}

NV_INT32 surf_load_xyz (NV_FLOAT32 x, NV_FLOAT32 y, NV_FLOAT32 z)
{
	coord3.x = (NV_FLOAT64) x;
	coord3.y = (NV_FLOAT64) y;
	coord3.z = (NV_FLOAT64) z;

	return (surf_load (coord3));
}

NV_INT32 surf_load_NV_F64_COORD2_Z (NV_F64_COORD2 xy, NV_FLOAT64 z)
{
	coord3.x = xy.x;
	coord3.y = xy.y;
	coord3.z = z;

	return (surf_load (coord3));
}

NV_INT32 surf_load_NV_F64_COORD2_z (NV_F64_COORD2 xy, NV_FLOAT32 z)
{
	coord3.x = xy.x;
	coord3.y = xy.y;
	coord3.z = (NV_FLOAT64) z;

	return (surf_load (coord3));
}

NV_INT32 surf_load_NV_F32_COORD2_z (NV_F32_COORD2 xy, NV_FLOAT32 z)
{
	coord3.x = (NV_FLOAT64) xy.x;
	coord3.y = (NV_FLOAT64) xy.y;
	coord3.z = (NV_FLOAT64) z;

	return (surf_load (coord3));
}

NV_INT32 surf_load_NV_F32_COORD3 (NV_F32_COORD3 xyz)
{
	coord3.x = (NV_FLOAT64) xyz.x;
	coord3.y = (NV_FLOAT64) xyz.y;
	coord3.z = (NV_FLOAT64) xyz.z;

	return (surf_load (coord3));
}

/*

  Function:         surf_proc         -  Compute the "blockmean" values and run the surface process.  You may
                                         optionally bypass the surface processing and use surf_rtrv to retrieve
					 the "blockmean" values without doing the spline interpolation (thanks
					 and a tip of the hat to Bill Rankin ;-)

  Arguments:        surface           -  Set to NVTrue if you want to run the surface process.

  Returns:          NV_INT32          -  0 if successful, on failure check error by calling get_surf_error
                                         and checking the surf_err_str, and surf_errno values.  See get_surf_error
					 for an explanation of the error values.  Important note - negative return
					 values are catastrophic errors (i.e. you need to display the error and
					 then GTHOOD [Get The Hell Out Of Dodge]).

*/

NV_INT32 surf_proc (NV_BOOL surface)
{
	NV_INT32          row, col, row_ndx, col_ndx, percent = 0, prev_percent = -1;
	NV_INT32 rc;
	long i;
	long j = 0;
	//Experimenting with median for blockmedian
	//NV_FLOAT64 tempx, tempy, tempz;

	/*  Make sure we successfully called surf_load at least once.  */
	if (process_state != 2)
	{
		surf_err = -12;
		surf_errno = 0;
		surf_cleanup();
		return (surf_err);
	}

	/*  Compute the "blockmean" values.  */
	for (row = 0 ; row < rows ; row++)
	{
		row_ndx = row * cols;
		for (col = 0 ; col < cols ; col++)
		{
			if (cnt[row_ndx + col])
			{
				x[row_ndx + col] /= (NV_FLOAT64) cnt[row_ndx + col];
				y[row_ndx + col] /= (NV_FLOAT64) cnt[row_ndx + col];
				z[row_ndx + col] /= (NV_FLOAT64) cnt[row_ndx + col];
				/*
				//Experimenting with median for blockmedian
				tempx = median(cnt[row_ndx + col], x[row_ndx + col]);// /= (NV_FLOAT64) cnt[row_ndx + col];
				tempy = median(cnt[row_ndx + col], y[row_ndx + col]);// /= (NV_FLOAT64) cnt[row_ndx + col];
				tempz = median(cnt[row_ndx + col], z[row_ndx + col]);// /= (NV_FLOAT64) cnt[row_ndx + col];
				*/
			}
			else
			{
				/* Flag these cells as containing no data */
				x[row_ndx + col] = -999999.999;
				y[row_ndx + col] = -999999.999;
				z[row_ndx + col] = -999999.999;
			}
			/*
			//Experimenting with median for blockmedian
			free(x[row_ndx+col]);
			free(y[row_ndx+col]);
			free(z[row_ndx+col]);
			*x[row_ndx+col] = tempx;
			*y[row_ndx+col] = tempy;
			*z[row_ndx+col] = tempz;
			*/
		}

		/*  If the caller has registered a callback function, give them some info.  */
		if (surf_progress_callback_registered ())
		{
			percent = NINT (((NV_FLOAT32) row / (NV_FLOAT32) rows) * 100.0);
			if (percent != prev_percent) surf_progress (NAVO_BLOCKMEAN, percent);
		}
	}

    if(DEBUG)
    {
       for (row = 0 ; row < rows ; row++)
       {
           row_ndx = row * cols;
           for (col = 0 ; col < cols ; col++)
           {
                printf("x[%d] = %f\n", (row_ndx + col), x[row_ndx + col]);
                printf("y[%d] = %f\n", (row_ndx + col), y[row_ndx + col]);
                printf("z[%d] = %f\n", (row_ndx + col), z[row_ndx + col]);
                printf("cnt[%d] = %d\n", (row_ndx + col), cnt[row_ndx + col]);
           }
       }
    }

	/*
	//Experimenting with median for blockmedian
	x_clean = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	y_clean = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	z_clean = (NV_FLOAT64 **) calloc (rows*cols, sizeof (NV_FLOAT64*));
	for( i=0; i<(rows*cols); i++){
		//x[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		//y[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		//z[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		x_clean[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		y_clean[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
		z_clean[i] = (NV_FLOAT64 *) calloc (rows*cols, sizeof (NV_FLOAT64));
	}
*/
	/*    for (row = 0 ; row < rows ; row++) */
	for (row = rows - 1 ; row >= 0; row--)
	{
		row_ndx = row * cols;
		for (col = 0 ; col < cols ; col++)
		{
			if(!((x[row_ndx + col] == -999999.999) && (y[row_ndx + col] == -999999.999)))
			{
				x_clean[j] = x[row_ndx + col];
				y_clean[j] = y[row_ndx + col];
				z_clean[j] = z[row_ndx + col];
				if(DEBUG)
				{
					printf("x_clean[%ld] = %f\n", j, x_clean[j]);
					printf("y_clean[%ld] = %f\n", j, y_clean[j]);
					printf("z_clean[%ld] = %f\n", j, z_clean[j]);
					printf("%f ", x_clean[j]);
					printf("%f ", y_clean[j]);
					printf("%f\n", z_clean[j]);
  				}
				j++;
			}
		}
	}

	if(DEBUG)
	{
		for(i = 0; i < j; i++)
		{
			printf("%6.3f ", x_clean[i]);
			printf("%7.4f ", y_clean[i]);
			printf("%5.0f\n", z_clean[i]);
		}
	}

	if(DEBUG)
	{
		for (col = 0; col < cols; col++)
		{
			col_ndx = col * rows;
			for(row = 0; row < rows; row++)
			{
			if(!((x[col_ndx + row] == -999999.999) && (y[col_ndx + row] == -999999.999)))
			{
				printf("%f ", x[col_ndx + row]);
				printf("%f ", y[col_ndx + row]);
				printf("%f\n", z[col_ndx + row]);
			}
			}
		}
	}

	/*  If the caller has registered a callback function, let them know that "blockmean" is done.  */
	if (surf_progress_callback_registered ()) surf_progress (NAVO_BLOCKMEAN, 100);

	/*  If we don't want to run surface, return now (to the days of yesteryear, the Lone Ranger rides... but I digress).  */
	if (!surface)
	{
		/*  Keep track of where we are in the process.  */
		process_state = 3;
		return (0);
    }

  /* Some of the cells may not have data so will need to be removed. */

  /*  Here's where we do the surface thingy Tom ;-)  */

  /*
      The following arrays and variables will be needed by the surface function:

      x       -  array of X values
      y       -  array of Y values
  */
/*
  x_orig = x_orig + (0.5 * x_grid);
  y_orig = y_orig + (0.5 * y_grid);
*/
/*
  x_max = x_orig + (cols * x_grid) - x_grid;
  y_max = y_orig + (rows * y_grid) - y_grid;
*/

	num_pts = j;

	if(DEBUG)
	{
		printf("Calling create_surface_from_array.\n");
		printf("x_grid = %lf\ty_grid = %lf\tx_orig = %lf\ty_orig = %lf\n", x_grid, y_grid, x_orig, y_orig);
		printf("x_max = %lf\ty_max = %lf\n", x_max, y_max);
		printf("num points = %d\n", num_pts);
		printf("rows = %ld\n", (long) rows);
		printf("cols = %ld\n", (long) cols);
		printf("tension = %lf\n", tension);
		fflush(stdout);
	
		printf("surf_proc: Calling create_surface_from_array. num_pts = %ld\n", (long) num_pts);
		for(i=0; i<num_pts; i++)
		{
			printf("surf_proc: x_clean[%ld] = %lf\ty_clean[%ld] = %lf\tz_clean[%ld] = %lf\n", i, x_clean[i], i, y_clean[i], i, z_clean[i]);
			fflush(stdout);
		}
	}

	/* The GMT Surface library will crash if less than 4 points are used */
	if(num_pts < 4)
	{
		surf_err = -14;
		surf_errno = 0;
		surf_cleanup();
		return (surf_err);
	}

	rc = create_surface_from_array(x_clean, y_clean, z_clean, num_pts, x_grid, y_grid, x_orig, y_orig, x_max, y_max, tension, &rows, &cols);
	if(rc < 0)
	{
		surf_err = -15;
		surf_errno = 0;
		surf_cleanup();
		return (surf_err);
	}

	/*  Keep track of where we are in the process.  */
	process_state = 3;
	return (0);
}

/*

  Function:         surf_rtrv         -  Returns the Z and CNT arrays to the caller.  The CNT array can be used to
                                         determine if there were actual input values in any grid cell.

  Arguments:        z_array           -  Pointer to Z array.
                    cnt_array         -  Pointer to CNT array.
		    final_rows        -  Number of rows gridded.
		    final_cols        -  Number of columns gridded.

  Returns:          NV_INT32          -  0 if successful, on failure check error by calling get_surf_error
                                         and checking the surf_err_str, and surf_errno values.  See get_surf_error
					 for an explanation of the error values.  Important note - negative return
					 values are catastrophic errors (i.e. you need to display the error and
					 then GTHOOD [Get The Hell Out Of Dodge]).

*/

NV_INT32 surf_rtrv (NV_FLOAT64 **z_array, NV_INT32 **cnt_array, NV_INT32 *final_rows, NV_INT32 *final_cols)
{
	NV_INT32 i;
	/*  Make sure we called surf_proc prior to calling surf_rtrv.  */
	if (process_state != 3)
	{
		surf_err = -13;
		surf_errno = 0;
		surf_cleanup();
		return (surf_err);
	}

	/*  Free the unused memory.  */
  	if (x)
	{
	    /*for( i=0; i<(rows*cols); i++)
		{
			free (x[i]);
			x[i] = NULL;
		}*/
		free (x);
		x = NULL;
	}
	if (y)
	{
		/*for( i=0; i<(rows*cols); i++)
		{
			free (y[i]);
			y[i] = NULL;
		}*/
		free (y);
		y = NULL;
	}

  *final_cols = cols;
  *final_rows = rows;

  *z_array = z_clean;
  *cnt_array = cnt;
  return (0);
}

/*

  Function:         surf_cleanup      -  Free the remaining memory after the caller has finished with it.

  Returns:          NV_INT32          -  0 if successful, on failure check error by calling get_surf_error
                                         and checking the surf_err_str, and surf_errno values.  See get_surf_error
					 for an explanation of the error values.  Important note - negative return
					 values are catastrophic errors (i.e. you need to display the error and
					 then GTHOOD [Get The Hell Out Of Dodge]).

*/

NV_INT32 surf_cleanup ()
{
	//NV_INT32 i;
	/*  Just in case they quit early.  */
	if (x)
	{
		/*for( i=0; i<(rows*cols); i++)
		{
			free (x[i]);
			x[i] = NULL;
		}*/
		free (x);
		x = NULL;
	}
	if (y)
	{
		/*for( i=0; i<(rows*cols); i++)
		{
			free (y[i]);
			y[i] = NULL;
		}*/
		free (y);
		y = NULL;
	}
	 /*  Free what was left after surf_rtrv.  */
	if (z)
    {
		/*for( i=0; i<(rows*cols); i++)
		{
			free (z[i]);
			z[i] = NULL;
		}*/
		free (z);
		z = NULL;
    }
	if (cnt)
    {
		free (cnt);
		cnt = NULL;
    }
	
	/*  Might as well check all just in case to avoid memory leaks.  */
	if (x_clean)
	{
		/*for( i=0; i<(rows*cols); i++)
		{
			free (x_clean[i]);
			x_clean[i] = NULL;
		}*/
		free (x_clean);
		x_clean = NULL;
	}
	if (y_clean)
	{
		/*for( i=0; i<(rows*cols); i++)
		{
			free (y_clean[i]);
			y_clean[i] = NULL;
		}*/
		free (y_clean);
		y_clean = NULL;
	}
	if (z_clean)
	{
		/*for( i=0; i<(rows*cols); i++)
		{
			free (z_clean[i]);
			z_clean[i] = NULL;
		}*/
		free (z_clean);
		z_clean = NULL;
	}
	if (cnt_clean)
	{
		free (cnt_clean);
		cnt_clean = NULL;
	}
	return (0);
}

/*

  Function:         get_surf_error    -  Returns error information about the last encountered error.

  return (0);
}

*/

/*

  Function:         get_surf_oob      -  Returns the surf out of bounds point count.

*/

NV_INT32 get_surf_oob ()
{
  return (out_of_bounds);
}

/*

  Function:         get_surf_error    -  Returns error information about the last encountered error.

  Returns:          error_string      -  Explanation of error (less than 128 characters).
                    system_errno      -  The system errno value from the last error.

*/

void get_surf_error (NV_CHAR *error_string, NV_INT32 *system_errno)
{
  *system_errno = surf_errno;

  switch (surf_err)
    {
    case -1:
      strcpy (error_string, "Negative X or Y grid interval is invalid");
      break;

    case -2:
      strcpy (error_string, "MBR maximum value is less than minimum value.  Use 0-360 world for dateline crossing.");
      break;

    case -3:
      strcpy (error_string, "Error allocating X array in surf_init");
      break;

    case -4:
      strcpy (error_string, "Error allocating Y array in surf_init");
      break;

    case -5:
      strcpy (error_string, "Error allocating Z array in surf_init");
      break;

    case -6:
      strcpy (error_string, "Error allocating CNT array in surf_init");
      break;

    case -7:
      strcpy (error_string, "surf_init must be called prior to calling surf_load");
      break;

    case -8:
      strcpy (error_string, "surf_load either wasn't called or there were no valid data points loaded");
      break;

    case -9:
      strcpy (error_string, "surf_proc must be called prior to calling surf_rtrv");
      break;

    case -14:
      strcpy (error_string, "attempting to create surface with less than four points (after blockmean)");
      break;

    case -15:
      strcpy (error_string, "error in GMT surface");
      break;

    case 0:
      strcpy (error_string, "Success!");
      break;

    case 1:
      strcpy (error_string, "Input point outside of MBR bounds");
      break;

    default:
      strcpy (error_string, "Unknown surf error value");
      break;
    }

  surf_err = 0;
  surf_errno = 0;
}

/*

  The following functions are pretty straight forward.  The calling application may want some progress information
  while "blockmean" and surface are running.  The caller can simply register a callback function like this:

  surf_register_progress_callback (surf_progress_callback);

  And then they need a function to display the information in some manner (I normally use a QProgressDialog for
  Qt graphics applications or fprintf for command line applications).  For example:

  static void surf_progress_callback (NV_INT32 state, NV_INT32 percent)
    {
      static NV_BOOL first = NVTrue;
      static NV_BOOL surface_processing = NVFalse;
      static NV_BOOL indeterminate = NVFalse;

      if (first)
        {
          surfProgressDialog->setLabelText ("Blockmean processing");
          surfProgressDialog->setRange (0, 100);
          first = NVFalse;
        }

      if (!surface_processing && state == NAVO_SURFACE)
        {
          surfProgressDialog->setLabelText ("Surface processing");
	  surface_processing = NVTrue;
        }

      if (percent < 0 && !indeterminate)
        {
    If you have a process that has an indeterminate time frame (that is, there is no way to
    compute the percentage done) then just call surf progress before starting with a state value
    and a percent value of -1.  This will tell the caller that something is going on.  When the
    time indeteminate section is complete, call surf_progress with the same state and percent
    set to 100.

*/

static  SURF_PROGRESS_CALLBACK  surf_progress_callback = NULL;
void surf_register_progress_callback (SURF_PROGRESS_CALLBACK progressCB)
{
    surf_progress_callback = progressCB;
}

NV_BOOL surf_progress_callback_registered ()
{
	if (surf_progress_callback)
	{
		return (NVTrue);
	}
	else
    {
      return (NVFalse);
    }
}

void surf_progress (NV_INT32 state, NV_INT32 percent)
{
  (*surf_progress_callback) (state, percent);
}
//Experimenting with median for blockmedian
//NV_FLOAT64 median(int n, NV_FLOAT64 x[]) {
//    NV_FLOAT64 temp;
//    int i, j;
//    // the following two loops sort the array x in ascending order
//    for(i=0; i<n-1; i++) {
//        for(j=i+1; j<n; j++) {
//            if(x[j] < x[i]) {
//                // swap elements
//                temp = x[i];
//                x[i] = x[j];
//                x[j] = temp;
//            }
//        }
//    }
//
//    if(n%2==0) {
//        // if there is an even number of elements, return mean of the two elements in the middle
//        return((x[n/2] + x[n/2 - 1]) / 2.0);
//    } else {
//        // else return the element in the middle
//        return x[n/2];
//    }
//}
//Restore warning state -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( pop )
#endif 