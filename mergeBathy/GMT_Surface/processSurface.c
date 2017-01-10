/* 
* File:			processSurface.c
* Purpose:		Act as the front end to the GMT Surface routines to make them easier to use
* Author:		Kevin Duvieilh
* Date:			06 June 2011
*
*/

#include <stdio.h>
#include <stdlib.h>
//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	#include "../WarningStates.h"	//Disable all Warnings!!!
#endif

#include "processSurface.h"
#include "surf.h"
#define DEBUG 0


extern NV_INT32 surf_init (NV_F64_XYMBR mbr, 
  			   NV_FLOAT64 x_interval, 
                           NV_FLOAT64 y_interval, 
                           NV_FLOAT32 minz, 
                           NV_FLOAT32 maxz,
                           NV_FLOAT64 tension_value);

extern NV_INT32 surf_load (NV_F64_COORD3 xyz);

extern  NV_INT32 surf_proc (NV_BOOL surface);

extern NV_INT32 surf_rtrv (NV_FLOAT64 **z_array, 
                 	   NV_INT32 **cnt_array, 
         		   NV_INT32 *final_rows, 
			   NV_INT32 *final_cols);

extern NV_INT32 surf_cleanup ();

/************************************************************************************
*
*       Module: processSurface
*       Parameters:
*		argv[1]:  input file name for ASCII xyz data file.
*		argv[2]:  x grid increment
*		argv[3]:  y grid increment
*		argv[4]:  tension factor 
*               argv[5]:  minz
*		argv[6]:  maxz
*               argv[7]:  mbr min x
*               argv[8]: mbr max x
*               argv[9]: mbr min y
*               argv[10]: mbr max y
*       Description: Driver program used to test the GMT-based surface library.
*               
*		This is essentially a re-coded version of surf_tst.c from the original GMT Surface library that plays nice with mergeBathy
************************************************************************************/

int processSurface(double *xData, double *yData, double *zData, int inputDataSize, double x0, double y0, double z0, double x1, double y1, double z1, double spacingX, double spacingY, double tension, double **xPostSurface, double **yPostSurface, double **zPostSurface, int *postSurfaceSize)
{
	NV_INT32 rc = 0; /* return code */
	long i, j, k, l;
	char buffer[256];
	double * coord_x, * coord_y, *coord_z; /* x y z coordinate arrays */
    NV_F64_XYMBR mbr;
	NV_FLOAT64 x_interval;
	NV_FLOAT64 y_interval;
	NV_FLOAT32 minz;
	NV_FLOAT32 maxz;
	NV_FLOAT64 tension_value;
	NV_F64_COORD3 xyz;
	NV_BOOL surface_flg;
    NV_FLOAT64 * z_final;
    NV_INT32 * cnt_array;
    NV_INT32 final_rows;
    NV_INT32 final_cols;
	NV_CHAR error_str[128];

	coord_x = coord_y = coord_z = NULL;
	i = j = k = l = 0;
        z_final = NULL;
	cnt_array = NULL;
        final_rows = final_cols = 0;

	buffer[0] = '\0';

	mbr.min_x = x0;
	mbr.max_x = x1;
	mbr.min_y = y0;
	mbr.max_y = y1;
	x_interval = spacingX;
	y_interval = spacingY;

	minz = z0;
	maxz = z1;
	tension_value = tension;
	surface_flg = 1;


 	if(DEBUG)
	{
		printf("Calling surf_init.\n");
		printf("\tmin_x = %13.10lf\n", mbr.min_x);
		printf("\tmax_x = %13.10lf\n", mbr.max_x);
		printf("\tmin_y = %13.10lf\n", mbr.min_y);
		printf("\tmax_y = %13.10lf\n", mbr.max_y);
		printf("\tx_interval = %13.10lf\n", x_interval);
		printf("\ty_interval = %13.10lf\n", y_interval);
		printf("\tminz = %f\n", minz);
		printf("\tmaxz = %f\n", maxz);
		printf("\ttension = %f\n", tension_value);
		
	}
	rc = surf_init(mbr, x_interval, y_interval, minz, maxz, tension_value);
	if(rc)
	{	
		error_str[0] = '\0';
		get_surf_error(error_str, &rc);
		if(error_str)
			printf("surf_init Error: %s\n", error_str);
		exit(1);
	}

	coord_x = (double *) malloc(inputDataSize * sizeof(double)); 	
	if (coord_x == NULL)
	{
		printf("Error allocating memory for x array\n");
		exit(1);
	}

	coord_y = (double *) malloc(inputDataSize * sizeof(double)); 	
	if (coord_y == NULL)
	{
		if(coord_x)
		{
			free(coord_x);
			coord_x = NULL;
		}
		printf("Error allocating memory for y array\n");
		exit(1);
	}

	coord_z = (double *) malloc(inputDataSize * sizeof(double)); 	
	if (coord_z == NULL)
	{
		if(coord_x)
		{
			free(coord_x);
			coord_x = NULL;
		}
		if(coord_y)
		{
			free(coord_y);
			coord_y = NULL;
		}
		printf("Error allocating memory for z array\n");
		exit(1);
	}

	/* Load the input xyz data. */
	for(j = 0; j < inputDataSize; j++)
	{
		coord_x[j] = xData[j];
		coord_y[j] = yData[j];
		coord_z[j] = zData[j];
		 
		if(DEBUG)
			printf("%ld. x=%lf y=%lf z=%lf\n", j, coord_x[j], coord_y[j], coord_z[j]);    
		xyz.x = coord_x[j];
		xyz.y = coord_y[j];
		xyz.z = coord_z[j];
        if(DEBUG)
			printf("%ld. x=%lf y=%lf z=%lf\n", j, xyz.x, xyz.y, xyz.z);    
		rc = surf_load(xyz);
		if(rc)
		{
			error_str[0] = '\0';
			get_surf_error(error_str, &rc);
			if(error_str)
				printf("surf_load Error: %s\n", error_str);
			exit(1);
		}
	}

	//printf("About to call surf_proc with surface flag set to %d\n", surface_flg);fflush(stdout);
    rc = surf_proc(surface_flg);	
	//printf("Rc after surf_proc is %d\n", rc);fflush(stdout);
	if(rc)
	{
		error_str[0] = '\0';
		get_surf_error(error_str, &rc);
		if(error_str)
			printf("surf_proc Error: %s\n", error_str);
		exit(1);
	}
	
	rc = surf_rtrv(&z_final, &cnt_array, &final_rows, &final_cols);

	printf("\tNumber of computed rows: %d\n", final_rows);
	printf("\tNumber of computed columns: %d\n", final_cols);

	//printf("After surf_rtrv rc is %d\n", rc);
	if(rc)
	{
		error_str[0] = '\0';
		get_surf_error(error_str, &rc);
		if(error_str)
			printf("surf_rtrv Error: %s\n", error_str);
		exit(1);
	}

	(*postSurfaceSize) = (final_rows*final_cols);
	(*xPostSurface) = (double *) malloc((*postSurfaceSize) * sizeof(double)); 
	(*yPostSurface) = (double *) malloc((*postSurfaceSize) * sizeof(double)); 
	(*zPostSurface) = (double *) malloc((*postSurfaceSize) * sizeof(double)); 

	i = 0;
	for(k = (final_rows - 1); k >= 0; k--)
	{
		for(l = 0; l < final_cols; l++)
		{
			(*xPostSurface)[i] = mbr.min_x + (l*x_interval);
			(*yPostSurface)[i] = mbr.min_y + (k*y_interval);
			(*zPostSurface)[i] = z_final[l+(k*final_cols)];
			i = i + 1;
		}
	}

	if(DEBUG)
  	{
               // printf("After surf_rtrv: \nx\ty\tz\n");
		for(k = (final_rows - 1); k >= 0; k--)
		{
		  for(l = 0; l < final_cols; l++)
		  {
                    printf("%f  %f  ", (mbr.min_x + (l*x_interval)), (mbr.min_y + (k*y_interval))); 
		    printf(" %f\n", z_final[l+(k*final_cols)]);
		  }
		}
	}

	rc = surf_cleanup();
	if(rc)
	{
		error_str[0] = '\0';
		get_surf_error(error_str, &rc);
		if(error_str)
			printf("surf_cleanup Error: %s\n", error_str);
		exit(1);
	}

	if(coord_x)
	{
		free(coord_x);
		coord_x = NULL;
	}
	if(coord_y)
	{
		free(coord_y);
		coord_y = NULL;
	}
	if(coord_z)
	{
		free(coord_z);
		coord_z = NULL;
	}

	return (int) rc;
}

//Restore warning state -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( pop )
#endif 