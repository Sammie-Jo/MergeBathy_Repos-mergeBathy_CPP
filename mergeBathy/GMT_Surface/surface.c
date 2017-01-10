/*--------------------------------------------------------------------
 *	$Id: surface.c,v 1.4 2008/08/28 21:33:26 tgray Exp tgray $
 *
 *	Copyright (c) 1991-2008 by P. Wessel and W. H. F. Smith
 *	See COPYING file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; version 2 of the License.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * surface.c:  a gridding program.
 * reads xyz triples and fits a surface to the data.
 * surface satisfies (1 - T) D4 z - T D2 z = 0,
 * where D4 is the biharmonic operator,
 * D2 is the Laplacian,
 * and T is a "tension factor" between 0 and 1.
 * End member T = 0 is the classical minimum curvature
 * surface.  T = 1 gives a harmonic surface.  Use T = 0.25
 * or so for potential data; something more for topography.
 *
 * Program includes overrelaxation for fast convergence and
 * automatic optimal grid factorization.
 *
 * See Smith & Wessel (Geophysics, 3, 293-305, 1990) for details.
 *
 * Authors: Walter H. F. Smith and Paul Wessel
 * Date: April, 1988.
 *
 * Version 3.0 1 Nov 1988:  New argument switches standardized.
 *	W. H. F. Smith
 *
 * Version 3.2 5 Dec 1988:	Search radius default = 0.0
 *	to skip this step - unnecessary now that planes are
 *	removed, unless grid starts at 1 for diabolical prime lattice.
 *	Ver 3.1 removed best fit plane, undocumented changes.
 *	W. H. F. Smith
 *
 * Version 4.0 -- 5 Dec 1988:	The scratch file has been eliminated,
 *	resulting in a C->small gain in speed.   It is replaced by the 
 *	in core struct SURFACE_BRIGGS.  Division by xinc or grid_xinc has been
 *	replaced with multiplies by reciprocals, which also made a C->small
 *	improvement in speed.  The computation of i, j to index a point
 *	was changed from rint(x) to floor(x + 0.5) to make it follow the
 *	convention of blockmean.  This is actually significant at times.
 *	The significance was discovered when the routine "throw_away_unusables"
 *	was written.  This routine and struct SURFACE_BRIGGS are the major differences
 *	between V3.2 and V4.0.  Modifications by w.h.f. smith
 *
 * V 4.1 -- 26 July, 1989:	W.H.F. Smith fixed up the arg loop and usage
 *	message, and added automatic scaling of the range of z values as an
 *	experiment.  Often this improves the accuracy/convergence of numerical
 *	work; simple testing today suggests it makes no difference, but I will
 *	leave it in.
 *
 * V 4.2 -- 4 June, 1991-2000:	P. Wessel upgraded to GMT v2.0 grid file i/o.
 *	Added feature to constrain solution to be within lower and/or
 *	upper bounds.  Bounds can be min/max input data value, another
 *	value outside the input data range, or be provided by a grid file
 *	whose values all are outside the input data range.
 *
 * V 4.3 -- 26 February, 1992:	W. H. F. Smith added option to suggest better
 *	dimensions, and fixed bug in -L option to -Lu -Ll so that full paths
 *	to lower/upper bound files may be specified.
 *
 * 21may98 whfs changed warning about prime dimensions to make it use -V
 *
 * 21-JUN-1998 PW: Upgraded to GMT 3.1 with -b option
 *
 * 24-JAN-2000 PW/WHFS: Fixed bug related to use of -L option for "near-node" constraints.
 * 10-JUL-2000 3.3.5 PW: Added -L plain.
 * Version:	4
 * Version:	4.1: 14-SEP-2005 by PW, Added enhanced -I
 * Version:	4.1.2 PW: Minimized (but not eliminated) global variables
 * 64-bit Compliant: Yes
 *
 *
 */


#include "gmt.h"
#include "surf.h"


struct SURFACE_CTRL {
	struct A {	/* -A<aspect_ratio> */
		BOOLEAN active;
		double value;
	} A;
	struct C {	/* -C<converge_limit> */
		BOOLEAN active;
		double value;
	} C;
	struct G {	/* -G<file> */
		BOOLEAN active;
		char *file;
	} G;
	struct I {	/* -Idx[/dy] */
		BOOLEAN active;
		double xinc, yinc;
	} I;
	struct L {	/* -Ll|u<limit> */
		BOOLEAN active;
		char *low, *high;
		double min, max;
		int lmode, hmode;
	} L;
	struct N {	/* -N<max_iterations> */
		BOOLEAN active;
		int value;
	} N;
	struct Q {	/* -Q */
		BOOLEAN active;
	} Q;
	struct S {	/* -S<radius>[m|c] */
		BOOLEAN active;
		double radius;
		char unit;
	} S;
	struct T {	/* -T<tension>[i][b] */
		BOOLEAN active;
		double b_tension, i_tension;
	} T;
	struct Z {	/* -Z<over_relaxation_parameter> */
		BOOLEAN active;
		double value;
	} Z;
};

#define SURFACE_OUTSIDE LONG_MAX	/* Index number indicating data is outside usable area */

struct SURFACE_DATA {	/* Data point and index to node it currently constrains  */
	float x;
	float y;
	float z;
	GMT_LONG index;
};

struct SURFACE_BRIGGS {		/* Coefficients in Taylor series for Laplacian(z) a la I. C. Briggs (1974)  */
	double b[6];
};

struct SURFACE_GLOBAL {		/* Things needed inside compare function must be global for now */
	GMT_LONG block_ny;			/* Number of nodes in y-dir for a given grid factor */
	double grid_xinc, grid_yinc;	/* size of each grid cell for a given grid factor */
	double x_min, y_min;		/* Lower left corner of grid */
} GMT_Surface_Global;

struct SURFACE_INFO {	/* Control structure for surface setup and execution */
	char *iu;			/* Pointer to grid info array */
	char mode_type[2];		/* D means include data points when iterating
					 * I means just interpolate from larger grid */
	char format[BUFSIZ];
	char *low_file, *high_file;	/* Pointers to grids with low and high limits, if selected */
	int grid, old_grid;		/* Node spacings  */
	int n_fact;			/* Number of factors in common (ny-1, nx-1) */
	int factors[32];		/* Array of common factors */
	int long_verbose;
	int set_low;			/* 0 unconstrained,1 = by min data value, 2 = by user value */
	int set_high;			/* 0 unconstrained,1 = by max data value, 2 = by user value */
	size_t n_alloc;
	GMT_LONG npoints;			/* Number of data points */
	GMT_LONG ij_sw_corner, ij_se_corner,ij_nw_corner, ij_ne_corner;
	GMT_LONG n_empty;			/* No of unconstrained nodes at initialization  */
	GMT_LONG nx;				/* Number of nodes in x-dir. */
	GMT_LONG ny;				/* Number of nodes in y-dir. (Final grid) */
	GMT_LONG nxny;			/* Total number of grid nodes without boundaries  */
	GMT_LONG mx;
	GMT_LONG my;
	GMT_LONG mxmy;			/* Total number of grid nodes with boundaries  */
	GMT_LONG block_nx;			/* Number of nodes in x-dir for a given grid factor */
	GMT_LONG block_ny;			/* Number of nodes in y-dir for a given grid factor */
	GMT_LONG max_iterations;		/* Max iter per call to iterate */
	GMT_LONG total_iterations;
	GMT_LONG grid_east;
	GMT_LONG offset[25][12];		/* Indices of 12 nearby points in 25 cases of edge conditions  */
	BOOLEAN constrained;		/* TRUE if set_low or set_high is TRUE */
	float *lower, *upper;		/* arrays for minmax values, if set */
	float *u;			/* Pointer to grid array */
	double low_limit, high_limit;	/* Constrains on range of solution */
	double grid_xinc, grid_yinc;	/* size of each grid cell for a given grid factor */
	double r_xinc, r_yinc, r_grid_xinc, r_grid_yinc;	/* Reciprocals  */
	double converge_limit;		/* Convergence limit */
	double radius;			/* Search radius for initializing grid  */
	double tension;			/* Tension parameter on the surface  */
	double boundary_tension;
	double interior_tension;
	double a0_const_1, a0_const_2;	/* Constants for off grid point equation  */
	double e_2, e_m2, one_plus_e2;
	double eps_p2, eps_m2, two_plus_ep2, two_plus_em2;
	double x_edge_const, y_edge_const;
	double l_epsilon;
	double z_mean;
	double z_scale;			/* Root mean square range of z after removing planar trend  */
	double r_z_scale;		/* reciprocal of z_scale  */
	double plane_c0, plane_c1, plane_c2;	/* Coefficients of best fitting plane to data  */
	double small;			/* Let data point coincide with node if distance < C->small */
	double coeff[2][12];		/* Coefficients for 12 nearby points, constrained and unconstrained  */
	double relax_old, relax_new;	/* Coefficients for relaxation factor to speed up convergence */
	struct SURFACE_DATA  *data;
	struct SURFACE_BRIGGS *briggs;
	struct GRD_HEADER h;
};

struct SURFACE_SUGGESTION {	/* Used to find top ten list of faster grid dimensions  */
	int nx;
	int ny;
	double factor;	/* Speed up by a factor of factor  */
};

int compare_points (const void *point_1v, const void *point_2v);
GMT_LONG gcd_euclid (GMT_LONG a, GMT_LONG b);	/* Finds the greatest common divisor  */
int get_prime_factors (GMT_LONG n, int *f);
int iterate (struct SURFACE_INFO *C, int mode);
void load_constraints (struct SURFACE_INFO *C, BOOLEAN transform);

/*
int create_surface_from_array(char *outfile_name,
			     char *increment_str,
			     char *roi_str,
			     double tension_factor,
			     NV_FLOAT64 *x,
			     NV_FLOAT64 *y,
                       	     NV_FLOAT64 *z,
			     long num_pts);
*/

int  create_surface_from_array(NV_FLOAT64 *x,
                            NV_FLOAT64 *y,
                            NV_FLOAT64 *z,
                            NV_INT32 num_pts,
                            NV_FLOAT64 x_grid,
                            NV_FLOAT64 y_grid,
                            NV_FLOAT64 x_orig,
                            NV_FLOAT64 y_orig,
                            NV_FLOAT64 x_max,
                            NV_FLOAT64 y_max,
                            NV_FLOAT64 tension,
                            NV_INT32 * rows,
                            NV_INT32 * cols);

/************************************************************************************************
*
*	Module: create_surface_from_array
*	Parameters:
*		outfile_name - The directory and path of the output surface file (netcdf format).
*		increment_str - The grid spacing.  Suffix modifiers include m (arcminutes), 
*				c (arcseconds), e (meters), k (kilometers), i (miiles), and
*				n (nautical miles).  See GMT documentation for the '-I' parameter
*				to the Surface program for further details.
*		roi_str - The region of interest in the format 'xmin/xmax/ymin/ymax'.  See GMT
*				documentation for the '-R' paramter to the Surface program
*				for further details.
*		tension_factor - Laplacian operator between 0 and 1.  See GMT documentation
*				for the '-T' parameter to the Surface program for further
*				details.	
*		x - Array of x coordinate values.
*		y - Array of y coordinate values.
*		z - Array of z coordinate values.
*		num_pts - The number of x,y,z coordinate values.
*	Return Values:
*		0 - Success
*		1 - Error
*	Description:
*		This function is a modification of the main() function of the GMT 'surface'
*		program that allows it to be compiled as a library.  In addition, modifications
*		were made to allow an array of x,y,z values to be passed in instead of reading
*		from an input file.  A gridded surface is created from the input coordinates
*		using the increment value and tension factor, and the results are output to
*		file in NetCDF format.
*
************************************************************************************************/
NV_INT32  create_surface_from_array(NV_FLOAT64 *x,
                            NV_FLOAT64 *y,
                            NV_FLOAT64 *z,
                            NV_INT32 num_pts,
                            NV_FLOAT64 x_grid,
                            NV_FLOAT64 y_grid,
                            NV_FLOAT64 x_orig,
                            NV_FLOAT64 y_orig,
                            NV_FLOAT64 x_max,
                            NV_FLOAT64 y_max,
                            NV_FLOAT64 tension,
                            NV_INT32 * rows,
                            NV_INT32 * cols)
{
	#define DEBUG 0	

	int argc = 5;
	char **argv = NULL;
	int rc = 0;
	char	modifier;

	int	i, j, error = FALSE, n_files = 0;

	FILE	*fp_in = NULL;	/* File pointer  */
	
	struct SURFACE_INFO C;
	struct SURFACE_CTRL *Ctrl;
	
	void suggest_sizes_for_surface (int factors[], int nx, int ny);
	void set_grid_parameters (struct SURFACE_INFO *C);
	void read_data (FILE *fp, struct SURFACE_INFO *C);
	void read_array (double *x, double *y, double *z, long num_pts, struct SURFACE_INFO *C);
	void throw_away_unusables (struct SURFACE_INFO *C);
	void remove_planar_trend (struct SURFACE_INFO *C);
	void rescale_z_values (struct SURFACE_INFO *C);
	void smart_divide (struct SURFACE_INFO *C);
	void set_offset (struct SURFACE_INFO *C);
	void set_index (struct SURFACE_INFO *C);
	void initialize_grid (struct SURFACE_INFO *C);
	void set_coefficients (struct SURFACE_INFO *C);
	void find_nearest_point (struct SURFACE_INFO *C);
	void fill_in_forecast (struct SURFACE_INFO *C);
	void check_errors (struct SURFACE_INFO *C);
	void replace_planar_trend (struct SURFACE_INFO *C);
	void write_output (struct SURFACE_INFO *C, char *grdfile, NV_FLOAT64 * z);
	void load_parameters (struct SURFACE_INFO *C, struct SURFACE_CTRL *Ctrl);
	void *New_Surface_Ctrl (), Free_Surface_Ctrl (struct SURFACE_CTRL *C);
	
/* TG	argc = GMT_begin (argc, argv); */

	Ctrl = (struct SURFACE_CTRL *)New_Surface_Ctrl ();	/* Allocate and initialize a new control structure */

	/* Check input parameters */
	if((x == NULL) || (y == NULL) || (z == NULL))
	{
		rc = 1;
		return rc;
	}

	if(num_pts == 0)
        {
                rc = 1;
                return rc;

        }

	if((tension < 0.0) || (tension > 1.0))
	{
		rc = 1;
		return rc;
	}
/*
	if(roi_str == NULL)
	{
		rc = 1;
		return rc;
	}

	if(increment_str == NULL)
	{
		rc = 1;
		return rc;
	}

	if(outfile_name == NULL)
	{
		rc = 1;
		return rc;
	}
*/
	argv = (char **) malloc(5*sizeof(char*));	
	if(argv == NULL)
	{
		rc = 1;
		return rc;
	}
		
 	for(i = 0; i<5; i++)
	{
		argv[i] = (char *) malloc(256*sizeof(char));
		if(argv[i] == NULL)
		{
			if(argv)
				free(argv);
			for(j = 0; j < i; j++)
			{
				if(argv[j])
					free(argv[j]);
			}
			rc = 1;
			return rc;
		}
	}

	strcpy(argv[0], "surface");
	strcpy(argv[1], "-Gout.tmp");
	sprintf(argv[2], "-I%13.10lf/%13.10lf", x_grid, y_grid);
	sprintf(argv[3], "-R%13.10lf/%13.10lf/%13.10lf/%13.10lf",x_orig, x_max, y_orig, y_max);
	sprintf(argv[4], "-T%f", tension);
	
	if(DEBUG)
	{	
		printf("argv[1] is %s\n", argv[1]);
		fflush(stdout);
		printf("argv[2] is %s\n", argv[2]);
		fflush(stdout);
		printf("argv[3] is %s\n", argv[3]);
		fflush(stdout);
		printf("argv[4] is %s\n", argv[4]);
		fflush(stdout);
	}

       argc = GMT_begin (argc, argv); 
	
	memset ((void *)&C, 0, sizeof (struct SURFACE_INFO));
	memset ((void *)&GMT_Surface_Global, 0, sizeof (struct SURFACE_GLOBAL));
	C.n_alloc = GMT_CHUNK;
	C.z_scale = C.r_z_scale = 1.0;
	C.mode_type[0] = 'I';
	C.mode_type[1] = 'D';	/* D means include data points when iterating */

	GMT_grd_init (&C.h, argc, argv, FALSE);

	if(DEBUG)
	{
		printf("After GMT_grd_init. argc = %d\n", argc);
		fflush(stdout);
	}

	for (i = 1; i < argc; i++) {
		if(DEBUG)
		{ 
			printf("parsing parameters..argv[%d] is %s\n", i, argv[i]);
		 	fflush(stdout);
		}

		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
              
				/* Common parameters */
                      
				case 'V':
					if (argv[i][2] == 'L' || argv[i][2] == 'l') C.long_verbose = TRUE;
				case 'H':
				case 'R':
				case ':':
				case 'b':	/* Input triplets are binary, not ascii */
				case 'f':
				case '\0':
				if(DEBUG)
				{
					printf("Calling GMT_parse_common_options. argv[%d] is %s\n",i,argv[i]);
					fflush(stdout);
				}
                                      error += GMT_parse_common_options (argv[i], &C.h.x_min, &C.h.x_max, &C.h.y_min, &C.h.y_max);
                                      break;
                              
				/* Supplemental parameters */
                              
				case 'A':
					Ctrl->A.active = TRUE;
					Ctrl->A.value = atof (&argv[i][2]);
					break;
				case 'C':
					Ctrl->C.active = TRUE;
					Ctrl->C.value = atof (&argv[i][2]);
					break;
				case 'G':
					Ctrl->G.active = TRUE;
					Ctrl->G.file = strdup (&argv[i][2]);
					break;
				case 'I':
					Ctrl->I.active = TRUE;
					if (GMT_getinc (&argv[i][2], &Ctrl->I.xinc, &Ctrl->I.yinc)) {
						GMT_inc_syntax ('I', 1);
						error = TRUE;
					}
					break;
				case 'L':	/* Set limits */
					Ctrl->L.active = TRUE;
					switch (argv[i][2]) {
						case 'l':	/* Lower limit  */
							if (argv[i][3] == 0) {
								/*fprintf (stderr, "%s: GMT SYNTAX ERROR -Ll option: No argument given\n", GMT_program); */
								error++;
							}
							Ctrl->L.low = strdup (&argv[i][3]);
							if (!GMT_access (Ctrl->L.low, R_OK))	/* file exists */
								Ctrl->L.lmode = 3;
							else if (Ctrl->L.low[0] == 'd')		/* Use data minimum */
								Ctrl->L.lmode = 1;
							else {
								Ctrl->L.lmode = 2;		/* Use given value */
								Ctrl->L.min = atof (&argv[i][3]);
							}
							break;
						case 'u':	/* Upper limit  */
							if (argv[i][3] == 0) {
								/*fprintf (stderr, "%s: GMT SYNTAX ERROR -Lu option: No argument given\n", GMT_program); */
								error++;
							}
							Ctrl->L.high = strdup (&argv[i][3]);
							if (!GMT_access (Ctrl->L.high, R_OK))	/* file exists */
								Ctrl->L.hmode = 3;
							else if (Ctrl->L.high[0] == 'd')	/* Use data maximum */
								Ctrl->L.hmode = 1;
							else {
								Ctrl->L.hmode = 2;		/* Use given value */
								Ctrl->L.max = atof (&argv[i][3]);
							}
							break;
						default:	/* 360-periodicity option */
							GMT_io.in_col_type[0] = GMT_io.out_col_type[0] = GMT_IS_LON;
							GMT_io.in_col_type[1] = GMT_io.out_col_type[1] = GMT_IS_LAT;
/*							fprintf (stderr, "%s: Option -L is obsolete (but is processed correctly).  Please use -f instead\n", GMT_program);*/
							break;
					}
					break;
				case 'N':
					Ctrl->N.active = TRUE;
					Ctrl->N.value = atoi (&argv[i][2]);
					break;
				case 'S':
					Ctrl->S.active = TRUE;
					Ctrl->S.radius = atof (&argv[i][2]);
					Ctrl->S.unit = argv[i][strlen(argv[i])-1];
					break;
				case 'T':
					Ctrl->T.active = TRUE;
					modifier = argv[i][strlen(argv[i])-1];
					if (modifier == 'b' || modifier == 'B') {
						Ctrl->T.b_tension = atof (&argv[i][2]);
					}
					else if (modifier == 'i' || modifier == 'I') {
						Ctrl->T.i_tension = atof (&argv[i][2]);
					}
					else if (modifier >= '0' && modifier <= '9') {
						Ctrl->T.i_tension = Ctrl->T.b_tension = atof (&argv[i][2]);
					}
					else {
						/*fprintf (stderr, "%s: GMT SYNTAX ERROR -T option: Unrecognized modifier %c\n", GMT_program, modifier); */
						error = TRUE;
					}
					break;
				case 'Q':
					Ctrl->Q.active = TRUE;
					break;
				case 'Z':
					Ctrl->Z.active = TRUE;
					Ctrl->Z.value = atof (&argv[i][2]);
					break;
				default:
					error = TRUE;
					GMT_default_error (argv[i][1]);
					break;
			}
		}
/*
		else {
			if ((fp_in = GMT_fopen (argv[i], GMT_io.r_mode)) == NULL) {
				fprintf (stderr, "surface: cannot open input data file %s\n", argv[i]);
				exit (EXIT_FAILURE);
			}
			n_files++;
		}
*/
	}
	if (argc == 1 || GMT_give_synopsis_and_exit) {	/* Display usage */
		fprintf (stderr, "surface %s - Adjustable tension continuous curvature surface gridding\n\n", GMT_VERSION);
		fprintf (stderr, "usage: surface [xyz-file] -G<output_grdfile_name> %s\n", GMT_I_OPT);
		fprintf (stderr, "\t%s [-A<aspect_ratio>] [-C<convergence_limit>] [%s]\n", GMT_Rgeo_OPT, GMT_H_OPT);
		fprintf (stderr, "\t[-Ll<limit>] [-Lu<limit>] [-N<n_iterations>] ] [-S<search_radius>[m|c]] [-T<tension>[i][b]]\n");
		fprintf (stderr, "\t[-Q] [-V[l]] [-Z<over_relaxation_parameter>] [%s] [%s] [%s]\n\n", GMT_t_OPT, GMT_bi_OPT, GMT_f_OPT);

		if (GMT_give_synopsis_and_exit) exit (EXIT_FAILURE);

		fprintf (stderr, "\tsurface will read from standard input or a single <xyz-file>.\n\n");
		fprintf (stderr, "\tRequired arguments to surface:\n");
		fprintf (stderr, "\t-G sets output grid file name\n");
		GMT_inc_syntax ('I', 0);
		fprintf (stderr, "\t\tNote that only gridline registration can be used.\n");
		GMT_explain_option ('R');
		fprintf (stderr, "\n\tOPTIONS:\n");
		fprintf (stderr, "\t-A<aspect_ratio>  = 1.0  by default which gives an isotropic solution.\n");
		fprintf (stderr, "\t\ti.e. xinc and yinc assumed to give derivatives of equal weight; if not, specify\n");
		fprintf (stderr, "\t\t<aspect_ratio> such that yinc = xinc / <aspect_ratio>.\n");
		fprintf (stderr, "\t\te.g. if gridding lon,lat use <aspect_ratio> = cosine(middle of lat range).\n");
		fprintf (stderr, "\t-C<convergence_limit> iteration stops when max abs change is less than <c.l.>\n");
		fprintf (stderr, "\t\tdefault will choose 0.001 of the range of your z data (1 ppt precision).\n");
		fprintf (stderr, "\t\tEnter your own convergence limit in same units as z data.\n");
		GMT_explain_option ('H');
		fprintf (stderr, "\t-L constrain the range of output values:\n");
		fprintf (stderr, "\t\t-Ll<limit> specifies lower limit; forces solution to be >= <limit>.\n");
		fprintf (stderr, "\t\t-Lu<limit> specifies upper limit; forces solution to be <= <limit>.\n");
		fprintf (stderr, "\t\t<limit> can be any number, or the letter d for min (or max) input data value,\n");
		fprintf (stderr, "\t\tor the filename of a grid with bounding values.  [Default solution unconstrained].\n");
		fprintf (stderr, "\t\tExample:  -Ll0 gives a non-negative solution.\n");
		fprintf (stderr, "\t-N sets max <n_iterations> in each cycle; default = 250.\n");
		fprintf (stderr, "\t-S sets <search_radius> to initialize grid; default = 0 will skip this step.\n");
		fprintf (stderr, "\t\tThis step is slow and not needed unless grid dimensions are pathological;\n");
		fprintf (stderr, "\t\ti.e., have few or no common factors.\n");
		fprintf (stderr, "\t\tAppend m or c to give <search_radius> in minutes or seconds.\n");
		fprintf (stderr, "\t-T adds Tension to the gridding equation; use a value between 0 and 1.\n");
		fprintf (stderr, "\t\tdefault = 0 gives minimum curvature (smoothest; bicubic) solution.\n");
		fprintf (stderr, "\t\t1 gives a harmonic spline solution (local max/min occur only at data points).\n");
		fprintf (stderr, "\t\ttypically 0.25 or more is good for potential field (smooth) data;\n");
		fprintf (stderr, "\t\t0.75 or so for topography.  Experiment.\n");
		fprintf (stderr, "\t\tAppend B or b to set tension in boundary conditions only;\n");
		fprintf (stderr, "\t\tAppend I or i to set tension in interior equations only;\n");
		fprintf (stderr, "\t\tNo appended letter sets tension for both to same value.\n");
		fprintf (stderr, "\t-Q Query for grid sizes that might run faster than your -R -I give.\n");
		GMT_explain_option ('V');
		fprintf (stderr, "\t\tAppend l for long verbose\n");
		fprintf (stderr, "\t-Z sets <over_relaxation parameter>.  Default = 1.4\n");
		fprintf (stderr, "\t\tUse a value between 1 and 2.  Larger number accelerates convergence but can be unstable.\n");
		fprintf (stderr, "\t\tUse 1 if you want to be sure to have (slow) stable convergence.\n\n");
		GMT_explain_option (':');
		GMT_explain_option ('i');
		GMT_explain_option ('n');
		fprintf (stderr, "\t\tDefault is 3 input columns.\n\n");
		GMT_explain_option ('f');
		GMT_explain_option ('.');
		fprintf (stderr, "\t(For additional details, see Smith & Wessel, Geophysics, 55, 293-305, 1990.)\n");
		return(-1);
		/*exit (EXIT_FAILURE); */
	}

	if (!project_info.region_supplied) {
	/*	fprintf (stderr, "%s: GMT SYNTAX ERROR:  Must specify -R option\n", GMT_program); */
		error++;
	}
	if (Ctrl->I.xinc <= 0.0 || Ctrl->I.yinc <= 0.0) {
	/*	fprintf (stderr, "%s: GMT SYNTAX ERROR -I option.  Must specify positive increment(s)\n", GMT_program); */
		error++;
	}
	if (Ctrl->N.value < 1) {
	/*	fprintf (stderr, "%s: GMT SYNTAX ERROR -N option.  Max iterations must be nonzero\n", GMT_program); */
		error++;
	}
	if (Ctrl->Z.value < 1.0 || Ctrl->Z.value > 2.0) {
	/*	fprintf (stderr, "%s: GMT SYNTAX ERROR -Z option.  Relaxation value must be 1 <= z <= 2\n", GMT_program); */
		error++;
	}
	if (!Ctrl->G.file) {
	/*	fprintf (stderr, "%s: GMT SYNTAX ERROR option -G:  Must specify output file\n", GMT_program); */
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) {
	/*	fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data cannot have header -H\n", GMT_program); */
		error++;
	}
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] == 0) GMT_io.ncol[GMT_IN] = 3;
	if (GMT_io.binary[GMT_IN] && GMT_io.ncol[GMT_IN] < 3) {
	/*	fprintf (stderr, "%s: GMT SYNTAX ERROR.  Binary input data (-bi) must have at least 3 columns\n", GMT_program); */
		error++;
	}
	if (n_files > 1) {
        /*	fprintf (stderr, "%s: GMT SYNTAX ERROR.  Only one input file allowed.  Preprocess you files with blockmean/median/mode first\n", GMT_program); */
		error++;
	}

/*	if (error) exit (EXIT_FAILURE); */	
	if(error)
	{
		rc = 1;
		return rc;
	}

	if (GMT_io.binary[GMT_IN] && gmtdefs.verbose) {
		/* char *type[2] = {"double", "single"}; */
		/*fprintf (stderr, "%s: Expects %d-column %s-precision binary data\n", GMT_program, GMT_io.ncol[GMT_IN], type[GMT_io.single_precision[GMT_IN]]); */
	}

	load_parameters (&C, Ctrl);	/* Pass parameters from parsing control to surface INFO structure */
	
	GMT_RI_prepare (&C.h);	/* Ensure -R -I consistency and set nx, ny */
	GMT_err_fail (GMT_grd_RI_verify (&C.h, 1), Ctrl->G.file);

	C.relax_old = 1.0 - C.relax_new;

	C.nx = (GMT_LONG)C.h.nx;
	C.ny = (GMT_LONG)C.h.ny;
	C.nxny = C.nx * C.ny;
	C.mx = C.nx + 4;
	C.my = C.ny + 4;
	C.mxmy = C.mx * C.my;
	C.r_xinc = 1.0 / C.h.x_inc;
	C.r_yinc = 1.0 / C.h.y_inc;
	GMT_Surface_Global.x_min = C.h.x_min;
	GMT_Surface_Global.y_min = C.h.y_min;

	*rows = C.ny;
	*cols = C.nx;
	/* New stuff here for v4.3:  Check out the grid dimensions:  */
	C.grid = gcd_euclid (C.nx-1, C.ny-1);

	if (gmtdefs.verbose || Ctrl->Q.active) {
		sprintf (C.format, "%s: W: %s E: %s S: %s N: %s nx: %%ld ny: %%ld\n", GMT_program, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
	/*	fprintf (stderr, C.format, C.h.x_min, C.h.x_max, C.h.y_min, C.h.y_max, C.nx-1, C.ny-1); */
	}
	if (C.grid == 1 && gmtdefs.verbose) fprintf (stderr, "%s:  WARNING:  Your grid dimensions are mutually prime.\n", GMT_program);
	if ((C.grid == 1 && gmtdefs.verbose) || Ctrl->Q.active) suggest_sizes_for_surface (C.factors, (int)C.nx-1, (int)C.ny-1);
	if (Ctrl->Q.active) 
	{
		return(1);
		/*exit (EXIT_SUCCESS);*/
	}

	/* New idea: set grid = 1, read data, setting index.  Then throw
		away data that can't be used in end game, constraining
		size of briggs->b[6] structure.  */

	C.grid = 1;
	set_grid_parameters (&C);
	if (fp_in == NULL) {
		fp_in = GMT_stdin;
#ifdef SET_IO_MODE
		GMT_setmode (GMT_IN);
#endif
	}
/*	read_data (fp_in, &C); */
	if(DEBUG)
		printf("Calling read_array.  x[0] is %f\ty[0] is %f\tz[0] is %f\tnum_pts is %ld\n", 
                       x[0], y[0], z[0], (long) num_pts);

	read_array ((double *)x, (double *)y, (double *)z, (long)num_pts, &C);
/*TG****************

	FILE * OUT1 = NULL;
 	long p;
	OUT1 = fopen("outtest.txt", "w");
	for (p=0; p<C.npoints;p++)
	{
		fprintf(OUT1, "%ld %f %f %f\n", p, C.data[p].x, C.data[p].y, C.data[p].z);
	}
	fclose(OUT1);
}
*****************TG*/
	throw_away_unusables (&C);
	remove_planar_trend (&C);
	rescale_z_values (&C);
	load_constraints (&C, TRUE);
	/* Set up factors and reset grid to first value  */

	C.grid = gcd_euclid (C.nx-1, C.ny-1);
	C.n_fact = get_prime_factors ((GMT_LONG)C.grid, C.factors);
	set_grid_parameters (&C);
	while (C.block_nx < 4 || C.block_ny < 4) 
        {
		if(DEBUG)
		{
			printf("Calling smart_divide\tC.block_nx = %ld\tC.block_ny = %ld\n", C.block_nx,C.block_ny);
			fflush(stdout);
		}
		smart_divide (&C);
		set_grid_parameters (&C);
	}

	set_offset (&C);
	set_index (&C);
	/* Now the data are ready to go for the first iteration.  */
	/* Allocate more space  */

	C.briggs = (struct SURFACE_BRIGGS *) GMT_memory (VNULL, (size_t)C.npoints, sizeof(struct SURFACE_BRIGGS), GMT_program);
	C.iu = (char *) GMT_memory (VNULL, (size_t)C.mxmy, sizeof(char), GMT_program);
	C.u = (float *) GMT_memory (VNULL, (size_t)C.mxmy, sizeof(float), GMT_program);


	if (C.radius > 0) initialize_grid (&C); /* Fill in nodes with a weighted avg in a search radius  */

	if (gmtdefs.verbose) fprintf (stderr, "%s: Grid\tMode\tIteration\tMax Change\tConv Limit\tTotal Iterations\n", GMT_program);

	set_coefficients (&C);


	C.old_grid = C.grid;
	find_nearest_point (&C);
	iterate (&C, 1);
	 
	while (C.grid > 1) {
		smart_divide (&C);
		set_grid_parameters (&C);
		set_offset (&C);
		set_index (&C);
		fill_in_forecast (&C);
		iterate(&C, 0);
		C.old_grid = C.grid;
		find_nearest_point (&C);
		iterate (&C, 1);
	}

	if (gmtdefs.verbose) check_errors (&C);


	replace_planar_trend (&C);
/*
{
FILE * grid_file = NULL;
long i, j, index = C.ij_sw_corner;

grid_file = fopen("./surf_lib.z", "w");
fprintf(grid_file, "C.ij_sw_corner = %ld\n", C.ij_sw_corner);
fprintf(grid_file, "C.nx = %ld\n", C.nx);
fprintf(grid_file, "C.ny = %ld\n", C.ny);
	for (i = 0; i < C.nx; i++, index += C.my) {
		for (j = 0; j < C.ny; j++) {
			fprintf(grid_file, "C.u[%ld] = %f\n", (index + C.ny - j - 1), C.u[index + C.ny - j - 1]);
		}
	}
fclose(grid_file);
}
*/
	GMT_free ((void *) C.data);
	GMT_free ((void *) C.briggs);
	GMT_free ((void *)C.iu);
	if (C.set_low) GMT_free ((void *)C.lower);
	if (C.set_high) GMT_free ((void *)C.upper);

	write_output (&C, Ctrl->G.file, z); 

	GMT_free ((void *) C.u);

	Free_Surface_Ctrl (Ctrl);	/* Deallocate control structure */

	GMT_end (argc, argv); 

	return rc;
/*	exit (EXIT_SUCCESS); */
}



void set_coefficients (struct SURFACE_INFO *C)
{
	double	e_4, loose, a0;

	loose = 1.0 - C->interior_tension;
	C->e_2 = C->l_epsilon * C->l_epsilon;
	e_4 = C->e_2 * C->e_2;
	C->eps_p2 = C->e_2;
	C->eps_m2 = 1.0/C->e_2;
	C->one_plus_e2 = 1.0 + C->e_2;
	C->two_plus_ep2 = 2.0 + 2.0*C->eps_p2;
	C->two_plus_em2 = 2.0 + 2.0*C->eps_m2;

	C->x_edge_const = 4 * C->one_plus_e2 - 2 * (C->interior_tension / loose);
	C->e_m2 = 1.0 / C->e_2;
	C->y_edge_const = 4 * (1.0 + C->e_m2) - 2 * (C->interior_tension * C->e_m2 / loose);


	a0 = 1.0 / ( (6 * e_4 * loose + 10 * C->e_2 * loose + 8 * loose - 2 * C->one_plus_e2) + 4*C->interior_tension*C->one_plus_e2);
	C->a0_const_1 = 2 * loose * (1.0 + e_4);
	C->a0_const_2 = 2.0 - C->interior_tension + 2 * loose * C->e_2;

	C->coeff[1][4] = C->coeff[1][7] = -loose;
	C->coeff[1][0] = C->coeff[1][11] = -loose * e_4;
	C->coeff[0][4] = C->coeff[0][7] = -loose * a0;
	C->coeff[0][0] = C->coeff[0][11] = -loose * e_4 * a0;
	C->coeff[1][5] = C->coeff[1][6] = 2 * loose * C->one_plus_e2;
	C->coeff[0][5] = C->coeff[0][6] = (2 * C->coeff[1][5] + C->interior_tension) * a0;
	C->coeff[1][2] = C->coeff[1][9] = C->coeff[1][5] * C->e_2;
	C->coeff[0][2] = C->coeff[0][9] = C->coeff[0][5] * C->e_2;
	C->coeff[1][1] = C->coeff[1][3] = C->coeff[1][8] = C->coeff[1][10] = -2 * loose * C->e_2;
	C->coeff[0][1] = C->coeff[0][3] = C->coeff[0][8] = C->coeff[0][10] = C->coeff[1][1] * a0;

	C->e_2 *= 2;		/* We will need these in boundary conditions  */
	C->e_m2 *= 2;

	C->ij_sw_corner = 2 * C->my + 2;			/*  Corners of array of actual data  */
	C->ij_se_corner = C->ij_sw_corner + (C->nx - 1) * C->my;
	C->ij_nw_corner = C->ij_sw_corner + C->ny - 1;
	C->ij_ne_corner = C->ij_se_corner + C->ny - 1;

}

void set_offset (struct SURFACE_INFO *C)
{
	GMT_LONG	add_w[5], add_e[5], add_s[5], add_n[5], add_w2[5], add_e2[5], add_s2[5], add_n2[5];
	int	i, j, kase;

	add_w[0] = -C->my; add_w[1] = add_w[2] = add_w[3] = add_w[4] = -C->grid_east;
	add_w2[0] = -2 * C->my;  add_w2[1] = -C->my - C->grid_east;  add_w2[2] = add_w2[3] = add_w2[4] = -2 * C->grid_east;
	add_e[4] = C->my; add_e[0] = add_e[1] = add_e[2] = add_e[3] = C->grid_east;
	add_e2[4] = 2 * C->my;  add_e2[3] = C->my + C->grid_east;  add_e2[2] = add_e2[1] = add_e2[0] = 2 * C->grid_east;

	add_n[4] = 1; add_n[3] = add_n[2] = add_n[1] = add_n[0] = C->grid;
	add_n2[4] = 2;  add_n2[3] = C->grid + 1;  add_n2[2] = add_n2[1] = add_n2[0] = 2 * C->grid;
	add_s[0] = -1; add_s[1] = add_s[2] = add_s[3] = add_s[4] = -C->grid;
	add_s2[0] = -2;  add_s2[1] = -C->grid - 1;  add_s2[2] = add_s2[3] = add_s2[4] = -2 * C->grid;

	for (i = 0, kase = 0; i < 5; i++) {
		for (j = 0; j < 5; j++, kase++) {
			C->offset[kase][0] = add_n2[j];
			C->offset[kase][1] = add_n[j] + add_w[i];
			C->offset[kase][2] = add_n[j];
			C->offset[kase][3] = add_n[j] + add_e[i];
			C->offset[kase][4] = add_w2[i];
			C->offset[kase][5] = add_w[i];
			C->offset[kase][6] = add_e[i];
			C->offset[kase][7] = add_e2[i];
			C->offset[kase][8] = add_s[j] + add_w[i];
			C->offset[kase][9] = add_s[j];
			C->offset[kase][10] = add_s[j] + add_e[i];
			C->offset[kase][11] = add_s2[j];
		}
	}
}



void fill_in_forecast (struct SURFACE_INFO *C) {

	/* Fills in bilinear estimates into new node locations
	   after grid is divided.   
	 */

	GMT_LONG i, j, ii, jj, index_0, index_1, index_2, index_3;
	GMT_LONG index_new;
	char *iu;
	double delta_x, delta_y, a0, a1, a2, a3;
	double old_size;
	float *u;
	u = C->u;
	iu = C->iu;

	old_size = 1.0 / (double)C->old_grid;

	/* first do from southwest corner */

	for (i = 0; i < C->nx-1; i += C->old_grid) {

		for (j = 0; j < C->ny-1; j += C->old_grid) {

			/* get indices of bilinear square */
			index_0 = C->ij_sw_corner + i * C->my + j;
			index_1 = index_0 + C->old_grid * C->my;
			index_2 = index_1 + C->old_grid;
			index_3 = index_0 + C->old_grid;

			/* get coefficients */
			a0 = u[index_0];
			a1 = u[index_1] - a0;
			a2 = u[index_3] - a0;
			a3 = u[index_2] - a0 - a1 - a2;

			/* find all possible new fill ins */

			for (ii = i;  ii < i + C->old_grid; ii += C->grid) {
				delta_x = (ii - i) * old_size;
				for (jj = j;  jj < j + C->old_grid; jj += C->grid) {
					index_new = C->ij_sw_corner + ii * C->my + jj;
					if (index_new == index_0) continue;
					delta_y = (jj - j) * old_size;
					u[index_new] = (float)(a0 + a1 * delta_x + delta_y * ( a2 + a3 * delta_x));
					iu[index_new] = 0;
				}
			}
			iu[index_0] = 5;
		}
	}

	/* now do linear guess along east edge */

	for (j = 0; j < (C->ny-1); j += C->old_grid) {
		index_0 = C->ij_se_corner + j;
		index_3 = index_0 + C->old_grid;
		for (jj = j;  jj < j + C->old_grid; jj += C->grid) {
			index_new = C->ij_se_corner + jj;
			delta_y = (jj - j) * old_size;
			u[index_new] = u[index_0] + (float)(delta_y * (u[index_3] - u[index_0]));
			iu[index_new] = 0;
		}
		iu[index_0] = 5;
	}
	/* now do linear guess along north edge */
	for (i = 0; i < (C->nx-1); i += C->old_grid) {
		index_0 = C->ij_nw_corner + i * C->my;
		index_1 = index_0 + C->old_grid * C->my;
		for (ii = i;  ii < i + C->old_grid; ii += C->grid) {
			index_new = C->ij_nw_corner + ii * C->my;
			delta_x = (ii - i) * old_size;
			u[index_new] = u[index_0] + (float)(delta_x * (u[index_1] - u[index_0]));
			iu[index_new] = 0;
		}
		iu[index_0] = 5;
	}
	/* now set northeast corner to fixed and we're done */
	iu[C->ij_ne_corner] = 5;
}

int compare_points (const void *point_1v, const void *point_2v)
{
		/*  Routine for qsort to sort data structure for fast access to data by node location.
		    Sorts on index first, then on radius to node corresponding to index, so that index
		    goes from low to high, and so does radius.
		*/
	GMT_LONG block_i, block_j, index_1, index_2;
	double x0, y0, dist_1, dist_2;
	struct SURFACE_DATA *point_1, *point_2;
	
	point_1 = (struct SURFACE_DATA *)point_1v;
	point_2 = (struct SURFACE_DATA *)point_2v;
	index_1 = point_1->index;
	index_2 = point_2->index;
	if (index_1 < index_2)
		return (-1);
	else if (index_1 > index_2)
		return (1);
	else if (index_1 == SURFACE_OUTSIDE)
		return (0);
	else {	/* Points are in same grid cell, find the one who is nearest to grid point */
		block_i = point_1->index/GMT_Surface_Global.block_ny;
		block_j = point_1->index%GMT_Surface_Global.block_ny;
		x0 = GMT_Surface_Global.x_min + block_i * GMT_Surface_Global.grid_xinc;
		y0 = GMT_Surface_Global.y_min + block_j * GMT_Surface_Global.grid_yinc;
		dist_1 = (point_1->x - x0) * (point_1->x - x0) + (point_1->y - y0) * (point_1->y - y0);
		dist_2 = (point_2->x - x0) * (point_2->x - x0) + (point_2->y - y0) * (point_2->y - y0);
		if (dist_1 < dist_2)
			return (-1);
		if (dist_1 > dist_2)
			return (1);
		else
			return (0);
	}
}

void smart_divide (struct SURFACE_INFO *C) {
		/* Divide grid by its largest prime factor */
	if(DEBUG)
	{ 
		printf("smart_divide: C->n_fact is %d\n", C->n_fact);
		fflush(stdout);
		printf("smart_divide: C->factors[%d - 1] is %d\n", C->n_fact, C->factors[C->n_fact - 1]);
		fflush(stdout);
	}
	C->grid /= C->factors[C->n_fact - 1];
	C->n_fact--;
	if(DEBUG)
	{
		printf("Exiting smart_divide\n");
		fflush(stdout);
	}
}

void set_index (struct SURFACE_INFO *C) {
		/* recomputes data[k].index for new value of grid,
		   sorts data on index and radii, and throws away
		   data which are now outside the usable limits. */
	GMT_LONG i, j, k, k_skipped = 0;

	for (k = 0; k < C->npoints; k++) {
		i = (GMT_LONG)floor(((C->data[k].x-C->h.x_min)*C->r_grid_xinc) + 0.5);
		j = (GMT_LONG)floor(((C->data[k].y-C->h.y_min)*C->r_grid_yinc) + 0.5);
		if (i < 0 || i >= C->block_nx || j < 0 || j >= C->block_ny) {
			C->data[k].index = SURFACE_OUTSIDE;
			k_skipped++;
		}
		else
			C->data[k].index = i * C->block_ny + j;
	}

	qsort ((void *)C->data, (size_t)C->npoints, sizeof (struct SURFACE_DATA), compare_points);

	C->npoints -= k_skipped;

}

void find_nearest_point (struct SURFACE_INFO *C) {
	GMT_LONG i, j, ij_v2, k, last_index, block_i, block_j, iu_index, briggs_index;
	double x0, y0, dx, dy, xys, xy1, btemp;
	double b0, b1, b2, b3, b4, b5;
	float z_at_node, *u;
	char *iu;
	u = C->u;
	iu = C->iu;

	last_index = -1;
	C->small = 0.05 * ((C->grid_xinc < C->grid_yinc) ? C->grid_xinc : C->grid_yinc);

	for (i = 0; i < C->nx; i += C->grid)	/* Reset grid info */
		for (j = 0; j < C->ny; j += C->grid)
			iu[C->ij_sw_corner + i*C->my + j] = 0;

	briggs_index = 0;
	for (k = 0; k < C->npoints; k++) {	/* Find constraining value  */
		if (C->data[k].index != last_index) {
			block_i = C->data[k].index/C->block_ny;
			block_j = C->data[k].index%C->block_ny;
			last_index = C->data[k].index;
	 		iu_index = C->ij_sw_corner + (block_i * C->my + block_j) * C->grid;
	 		x0 = C->h.x_min + block_i*C->grid_xinc;
	 		y0 = C->h.y_min + block_j*C->grid_yinc;
	 		dx = (C->data[k].x - x0)*C->r_grid_xinc;
	 		dy = (C->data[k].y - y0)*C->r_grid_yinc;
	 		if (fabs(dx) < C->small && fabs(dy) < C->small) {	/* Close enough to assign value to node */
	 			iu[iu_index] = 5;
	 			/* v3.3.4: NEW CODE
	 			 * Since point is basically moved from (dx, dy) to (0,0) we must adjust for
	 			 * the C->small change in the planar trend between the two locations, and then
	 			 * possibly clip the range if constraining surfaces were given.  Note that
	 			 * dx, dy is in -1/1 range normalized by (grid * x|y_inc) so to recover the
	 			 * dx,dy in final grid fractions we must scale by grid */
	 			 
	 			z_at_node = C->data[k].z + (float) (C->r_z_scale * C->grid * (C->plane_c1 * dx + C->plane_c2 * dy));
	 			if (C->constrained) {
					ij_v2 = (C->ny - block_j * C->grid - 1) * C->nx + block_i * C->grid;
					if (C->set_low  && !GMT_is_fnan (C->lower[ij_v2]) && z_at_node < C->lower[ij_v2])
						z_at_node = C->lower[ij_v2];
					else if (C->set_high && !GMT_is_fnan (C->upper[ij_v2]) && z_at_node > C->upper[ij_v2])
						z_at_node = C->upper[ij_v2];
	 			}
	 			u[iu_index] = z_at_node;
	 		}
	 		else {
	 			if (dx >= 0.0) {
	 				if (dy >= 0.0)
	 					iu[iu_index] = 1;
	 				else
	 					iu[iu_index] = 4;
	 			}
	 			else {
	 				if (dy >= 0.0)
	 					iu[iu_index] = 2;
	 				else
	 					iu[iu_index] = 3;
	 			}
	 			dx = fabs(dx);
	 			dy = fabs(dy);
	 			btemp = 2 * C->one_plus_e2 / ( (dx + dy) * (1.0 + dx + dy) );
	 			b0 = 1.0 - 0.5 * (dx + (dx * dx)) * btemp;
	 			b3 = 0.5 * (C->e_2 - (dy + (dy * dy)) * btemp);
	 			xys = 1.0 + dx + dy;
	 			xy1 = 1.0 / xys;
	 			b1 = (C->e_2 * xys - 4 * dy) * xy1;
	 			b2 = 2 * (dy - dx + 1.0) * xy1;
	 			b4 = b0 + b1 + b2 + b3 + btemp;
	 			b5 = btemp * C->data[k].z;
	 			C->briggs[briggs_index].b[0] = b0;
	 			C->briggs[briggs_index].b[1] = b1;
	 			C->briggs[briggs_index].b[2] = b2;
	 			C->briggs[briggs_index].b[3] = b3;
	 			C->briggs[briggs_index].b[4] = b4;
	 			C->briggs[briggs_index].b[5] = b5;
	 			briggs_index++;
	 		}
	 	}
	 }
}


void set_grid_parameters (struct SURFACE_INFO *C)
{
	GMT_Surface_Global.block_ny = C->block_ny = (C->ny - 1) / C->grid + 1;
	C->block_nx = (C->nx - 1) / C->grid + 1;
	GMT_Surface_Global.grid_xinc = C->grid_xinc = C->grid * C->h.x_inc;
	GMT_Surface_Global.grid_yinc = C->grid_yinc = C->grid * C->h.y_inc;
	C->grid_east = C->grid * C->my;
	C->r_grid_xinc = 1.0 / C->grid_xinc;
	C->r_grid_yinc = 1.0 / C->grid_yinc;
}

void initialize_grid (struct SURFACE_INFO *C)
{	/*
	 * For the initial gridsize, compute weighted averages of data inside the search radius
	 * and assign the values to u[i,j] where i,j are multiples of gridsize.
	 */
	GMT_LONG	irad, jrad, i, j, imin, imax, jmin, jmax, index_1, index_2, k, ki, kj, k_index;
	double	r, rfact, sum_w, sum_zw, weight, x0, y0;
	float *u;
	u = C->u;

	 irad = (GMT_LONG)ceil(C->radius/C->grid_xinc);
	 jrad = (GMT_LONG)ceil(C->radius/C->grid_yinc);
	 rfact = -4.5/(C->radius*C->radius);
	 for (i = 0; i < C->block_nx; i ++ ) {
	 	x0 = C->h.x_min + i*C->grid_xinc;
	 	for (j = 0; j < C->block_ny; j ++ ) {
	 		y0 = C->h.y_min + j*C->grid_yinc;
	 		imin = i - irad;
	 		if (imin < 0) imin = 0;
	 		imax = i + irad;
	 		if (imax >= C->block_nx) imax = C->block_nx - 1;
	 		jmin = j - jrad;
	 		if (jmin < 0) jmin = 0;
	 		jmax = j + jrad;
	 		if (jmax >= C->block_ny) jmax = C->block_ny - 1;
	 		index_1 = imin*C->block_ny + jmin;
	 		index_2 = imax*C->block_ny + jmax + 1;
	 		sum_w = sum_zw = 0.0;
	 		k = 0;
	 		while (k < C->npoints && C->data[k].index < index_1) k++;
	 		for (ki = imin; k < C->npoints && ki <= imax && C->data[k].index < index_2; ki++) {
	 			for (kj = jmin; k < C->npoints && kj <= jmax && C->data[k].index < index_2; kj++) {
	 				k_index = ki*C->block_ny + kj;
	 				while (k < C->npoints && C->data[k].index < k_index) k++;
	 				while (k < C->npoints && C->data[k].index == k_index) {
	 					r = (C->data[k].x-x0)*(C->data[k].x-x0) + (C->data[k].y-y0)*(C->data[k].y-y0);
	 					weight = exp (rfact*r);
	 					sum_w += weight;
	 					sum_zw += weight*C->data[k].z;
	 					k++;
	 				}
	 			}
	 		}
	 		if (sum_w == 0.0) {
	 			sprintf (C->format, "%%s: Warning: no data inside search radius at: %s %s\n", gmtdefs.d_format, gmtdefs.d_format);
	 			fprintf (stderr, C->format, GMT_program, x0, y0);
	 			u[C->ij_sw_corner + (i * C->my + j) * C->grid] = (float)C->z_mean;
	 		}
	 		else {
	 			u[C->ij_sw_corner + (i*C->my+j)*C->grid] = (float)(sum_zw/sum_w);
	 		}
		}
	}
}


void new_initialize_grid (struct SURFACE_INFO *C)
{	/*
	 * For the initial gridsize, load constrained nodes with weighted avg of their data;
	 * and then do something with the unconstrained ones.
	 */
	GMT_LONG	k, k_index, u_index, block_i, block_j;
	char *iu;
	double	sum_w, sum_zw, weight, x0, y0, dx, dy, dx_scale, dy_scale;
	float *u;
	u = C->u;
	iu = C->iu;

	dx_scale = 4.0 / C->grid_xinc;
	dy_scale = 4.0 / C->grid_yinc;
	C->n_empty = C->block_ny * C->block_nx;
	k = 0;
	while (k < C->npoints) {
		block_i = C->data[k].index / C->block_ny;
		block_j = C->data[k].index % C->block_ny;
		x0 = C->h.x_min + block_i*C->grid_xinc;
		y0 = C->h.y_min + block_j*C->grid_yinc;
		u_index = C->ij_sw_corner + (block_i*C->my + block_j) * C->grid;
		k_index = C->data[k].index;

		dy = (C->data[k].y - y0) * dy_scale;
		dx = (C->data[k].x - x0) * dx_scale;
		sum_w = 1.0 / (1.0 + dx*dx + dy*dy);
		sum_zw = C->data[k].z * sum_w;
		k++;

		while (k < C->npoints && C->data[k].index == k_index) {

			dy = (C->data[k].y - y0) * dy_scale;
			dx = (C->data[k].x - x0) * dx_scale;
			weight = 1.0 / (1.0 + dx*dx + dy*dy);
			sum_zw += C->data[k].z * weight;
			sum_w += weight;
			sum_zw += weight*C->data[k].z;
			k++;
	 	}
	 	u[u_index] = (float)(sum_zw/sum_w);
	 	iu[u_index] = 5;
	 	C->n_empty--;
	 }
}

void read_data (FILE *fp_in, struct SURFACE_INFO *C)
{
	GMT_LONG	i, j, k, kmax = 0, kmin = 0;
	int n_fields, n_expected_fields;
	double	*in, zmin = DBL_MAX, zmax = -DBL_MAX;
	char	buffer[BUFSIZ];

	C->data = (struct SURFACE_DATA *) GMT_memory (VNULL, (size_t)C->n_alloc, sizeof(struct SURFACE_DATA), GMT_program);

	/* Read in xyz data and computes index no and store it in a structure */

	k = 0;
	C->z_mean = 0;
	n_expected_fields = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3;
	if (!GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (buffer, BUFSIZ, fp_in);

	while ((n_fields = GMT_input (fp_in, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) {	/* Not yet EOF */

		if (GMT_io.status & GMT_IO_MISMATCH) {
			/*fprintf (stderr, "%s: Mismatch between actual (%d) and expected (%d) fields near line %ld\n", GMT_program, n_fields, n_expected_fields, k);
			exit (EXIT_FAILURE); */
		}

		if (GMT_is_dnan (in[2])) continue;

		i = (GMT_LONG)floor(((in[0]-C->h.x_min)*C->r_grid_xinc) + 0.5);
		if (i < 0 || i >= C->block_nx) continue;
		j = (GMT_LONG)floor(((in[1]-C->h.y_min)*C->r_grid_yinc) + 0.5);
		if (j < 0 || j >= C->block_ny) continue;

		C->data[k].index = i * C->block_ny + j;
		C->data[k].x = (float)in[0];
		C->data[k].y = (float)in[1];
		C->data[k].z = (float)in[2];
		if (zmin > in[2]) zmin = in[2], kmin = k;
		if (zmax < in[2]) zmax = in[2], kmax = k;
		k++;
		C->z_mean += in[2];
		if (k == (GMT_LONG)C->n_alloc) {
			C->n_alloc <<= 1;
			C->data = (struct SURFACE_DATA *) GMT_memory ((void *)C->data, (size_t)C->n_alloc, sizeof(struct SURFACE_DATA), GMT_program);
		}
	}

	if (fp_in != GMT_stdin) GMT_fclose (fp_in);

	C->npoints = k;

	if (C->npoints == 0) {
		/*fprintf (stderr, "%s:  No datapoints inside region, aborts\n", GMT_program);
		exit (EXIT_FAILURE);*/
	}

	C->z_mean /= k;
	if (gmtdefs.verbose) {
		sprintf(C->format, "%s: %s %s %s\n", GMT_program, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, "%s: Minimum value of your dataset x,y,z at: ", GMT_program);
		fprintf (stderr, C->format, (double)C->data[kmin].x, (double)C->data[kmin].y, (double)C->data[kmin].z);
		fprintf (stderr, "%s: Maximum value of your dataset x,y,z at: ", GMT_program);
		fprintf (stderr, C->format, (double)C->data[kmax].x, (double)C->data[kmax].y, (double)C->data[kmax].z);
	}
	C->data = (struct SURFACE_DATA *) GMT_memory ((void *)C->data, (size_t)C->npoints, sizeof(struct SURFACE_DATA), GMT_program);

	if (C->set_low == 1)
		C->low_limit = C->data[kmin].z;
	else if (C->set_low == 2 && C->low_limit > C->data[kmin].z) {
	/*	C->low_limit = data[kmin].z;	*/
		fprintf (stderr, "%s: Warning:  Your lower value is > than min data value.\n", GMT_program);
	}
	if (C->set_high == 1)
		C->high_limit = C->data[kmax].z;
	else if (C->set_high == 2 && C->high_limit < C->data[kmax].z) {
	/*	C->high_limit = data[kmax].z;	*/
		fprintf (stderr, "%s: Warning:  Your upper value is < than max data value.\n", GMT_program);
	}

}

void read_array (double *x, double *y, double *z, long num_pts, struct SURFACE_INFO *C)
{
	GMT_LONG	i, j, k, kmax = 0, kmin = 0;
	int n_expected_fields;
	double	zmin = DBL_MAX, zmax = -DBL_MAX;
	long l = 0;

	C->data = (struct SURFACE_DATA *) GMT_memory (VNULL, (size_t)C->n_alloc, sizeof(struct SURFACE_DATA), GMT_program);

	/* Read in xyz data and computes index no and store it in a structure */

	k = 0;
	C->z_mean = 0;
	n_expected_fields = (GMT_io.binary[GMT_IN]) ? GMT_io.ncol[GMT_IN] : 3;
/*	if (!GMT_io.binary[GMT_IN] && GMT_io.io_header[GMT_IN]) for (i = 0; i < GMT_io.n_header_recs; i++) GMT_fgets (buffer, BUFSIZ, fp_in);
*/
/*	while ((n_fields = GMT_input (fp_in, &n_expected_fields, &in)) >= 0 && !(GMT_io.status & GMT_IO_EOF)) { */	/* Not yet EOF */
	for(l = 0; l < num_pts;  l++)
	{
/*
		if (GMT_io.status & GMT_IO_MISMATCH) {
			fprintf (stderr, "%s: Mismatch between actual (%d) and expected (%d) fields near line %ld\n", GMT_program, n_fields, n_expected_fields, k);
			exit (EXIT_FAILURE);
		}
*/
		if (GMT_is_dnan (z[l])) continue;
		i = (GMT_LONG)floor(((x[l]-C->h.x_min)*C->r_grid_xinc) + 0.5);
		if (i < 0 || i >= C->block_nx) continue;
		j = (GMT_LONG)floor(((y[l]-C->h.y_min)*C->r_grid_yinc) + 0.5);
		if (j < 0 || j >= C->block_ny) continue;

		C->data[k].index = i * C->block_ny + j;
		C->data[k].x = (float)x[l];
		C->data[k].y = (float)y[l];
		C->data[k].z = (float)z[l];
		if (zmin > z[l]) zmin = z[l], kmin = k;
		if (zmax < z[l]) zmax = z[l], kmax = k;
		k++;
		C->z_mean += z[l];
		if (k == (GMT_LONG)C->n_alloc) {
			C->n_alloc <<= 1;
			C->data = (struct SURFACE_DATA *) GMT_memory ((void *)C->data, (size_t)C->n_alloc, sizeof(struct SURFACE_DATA), GMT_program);
		}
	}

/*	if (fp_in != GMT_stdin) GMT_fclose (fp_in); */

	C->npoints = k;
	if (C->npoints == 0) {
		/*fprintf (stderr, "%s:  No datapoints inside region, aborts\n", GMT_program);
		exit (EXIT_FAILURE);*/
	}

	C->z_mean /= k;
	if (gmtdefs.verbose) {
		sprintf(C->format, "%s: %s %s %s\n", GMT_program, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
		fprintf (stderr, "%s: Minimum value of your dataset x,y,z at: ", GMT_program);
		fprintf (stderr, C->format, (double)C->data[kmin].x, (double)C->data[kmin].y, (double)C->data[kmin].z);
		fprintf (stderr, "%s: Maximum value of your dataset x,y,z at: ", GMT_program);
		fprintf (stderr, C->format, (double)C->data[kmax].x, (double)C->data[kmax].y, (double)C->data[kmax].z);
	}
	C->data = (struct SURFACE_DATA *) GMT_memory ((void *)C->data, (size_t)C->npoints, sizeof(struct SURFACE_DATA), GMT_program);

	if (C->set_low == 1)
		C->low_limit = C->data[kmin].z;
	else if (C->set_low == 2 && C->low_limit > C->data[kmin].z) {
	/*	C->low_limit = data[kmin].z;	*/
		fprintf (stderr, "%s: Warning:  Your lower value is > than min data value.\n", GMT_program);
	}
	if (C->set_high == 1)
		C->high_limit = C->data[kmax].z;
	else if (C->set_high == 2 && C->high_limit < C->data[kmax].z) {
	/*	C->high_limit = data[kmax].z;	*/
		fprintf (stderr, "%s: Warning:  Your upper value is < than max data value.\n", GMT_program);
	}
}

void write_output (struct SURFACE_INFO *C, char *grdfile, NV_FLOAT64 *z)
{	/* Uses v.2.0 netCDF grd format - hence need to transpose original grid to be GMT compatible.  This will be rewritten, eventually */
	GMT_LONG	index, i, j, k;
	float *u, *v2;
	u = C->u;

	load_constraints (C, FALSE);	/* Reload constraints but this time do not transform data */
		
	strcpy (C->h.title, GMT_program);

	v2 = (float *) GMT_memory (VNULL, (size_t)C->nxny, sizeof (float), GMT_program);
	index = C->ij_sw_corner;
	for (i = 0; i < C->nx; i++, index += C->my) {
		for (j = 0; j < C->ny; j++) {
			k = j * C->nx + i;
			v2[k] = u[index + C->ny - j - 1];
/*printf("write_output: v2[%ld] = %f\n", k, u[index + C->ny - j - 1]);   */
			if (C->set_low  && !GMT_is_fnan (C->lower[k]) && v2[k] < C->lower[k]) v2[k] = C->lower[k];
			if (C->set_high && !GMT_is_fnan (C->upper[k]) && v2[k] > C->upper[k]) v2[k] = C->upper[k];
		}
	}
	
	GMT_err_fail (GMT_write_grd (grdfile, &(C->h), v2, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE, z), grdfile);
	GMT_free ((void *)v2);
	if (C->set_low) GMT_free ((void *)C->lower);
	if (C->set_high) GMT_free ((void *)C->upper);
}

int iterate (struct SURFACE_INFO *C, int mode)
{

	GMT_LONG	i, j, k, ij, kase, briggs_index, ij_v2;
	GMT_LONG	x_case, y_case, x_w_case, x_e_case, y_s_case, y_n_case;
	GMT_LONG	iteration_count = 0;
	char *iu;

	double	current_limit = C->converge_limit / C->grid;
	double	change, max_change = 0.0, busum, sum_ij;
	double	b0, b1, b2, b3, b4, b5;
	float *u;

	double	x_0_const = 4.0 * (1.0 - C->boundary_tension) / (2.0 - C->boundary_tension);
	double	x_1_const = (3 * C->boundary_tension - 2.0) / (2.0 - C->boundary_tension);
	double	y_denom = 2 * C->l_epsilon * (1.0 - C->boundary_tension) + C->boundary_tension;
	double	y_0_const = 4 * C->l_epsilon * (1.0 - C->boundary_tension) / y_denom;
	double	y_1_const = (C->boundary_tension - 2 * C->l_epsilon * (1.0 - C->boundary_tension) ) / y_denom;

	sprintf(C->format,"%s: %%4ld\t%%c\t%%8ld\t%s\t%s\t%%10ld\n", GMT_program, gmtdefs.d_format, gmtdefs.d_format);

	u = C->u;
	iu = C->iu;
	do {
		briggs_index = 0;	/* Reset the constraint table stack pointer  */

		max_change = -1.0;

		/* Fill in auxiliary boundary values (in new way) */

		/* First set d2[]/dn2 = 0 along edges:  */
		/* New experiment : (1-T)d2[]/dn2 + Td[]/dn = 0  */

		for (i = 0; i < C->nx; i += C->grid) {
			/* set d2[]/dy2 = 0 on south side:  */
			ij = C->ij_sw_corner + i * C->my;
			/* u[ij - 1] = 2 * u[ij] - u[ij + grid];  */
			u[ij - 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij + C->grid]);
			/* set d2[]/dy2 = 0 on north side:  */
			ij = C->ij_nw_corner + i * C->my;
			/* u[ij + 1] = 2 * u[ij] - u[ij - grid];  */
			u[ij + 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij - C->grid]);

		}

		for (j = 0; j < C->ny; j += C->grid) {
			/* set d2[]/dx2 = 0 on west side:  */
			ij = C->ij_sw_corner + j;
			/* u[ij - my] = 2 * u[ij] - u[ij + grid_east];  */
			u[ij - C->my] = (float)(x_1_const * u[ij + C->grid_east] + x_0_const * u[ij]);
			/* set d2[]/dx2 = 0 on east side:  */
			ij = C->ij_se_corner + j;
			/* u[ij + my] = 2 * u[ij] - u[ij - grid_east];  */
			u[ij + C->my] = (float)(x_1_const * u[ij - C->grid_east] + x_0_const * u[ij]);
		}

		/* Now set d2[]/dxdy = 0 at each corner:  */

		ij = C->ij_sw_corner;
		u[ij - C->my - 1] = u[ij + C->grid_east - 1] + u[ij - C->my + C->grid] - u[ij + C->grid_east + C->grid];

		ij = C->ij_nw_corner;
		u[ij - C->my + 1] = u[ij + C->grid_east + 1] + u[ij - C->my - C->grid] - u[ij + C->grid_east - C->grid];

		ij = C->ij_se_corner;
		u[ij + C->my - 1] = u[ij - C->grid_east - 1] + u[ij + C->my + C->grid] - u[ij - C->grid_east + C->grid];

		ij = C->ij_ne_corner;
		u[ij + C->my + 1] = u[ij - C->grid_east + 1] + u[ij + C->my - C->grid] - u[ij - C->grid_east - C->grid];

		/* Now set (1-T)dC/dn + Tdu/dn = 0 at each edge :  */
		/* New experiment:  only dC/dn = 0  */

		x_w_case = 0;
		x_e_case = C->block_nx - 1;
		for (i = 0; i < C->nx; i += C->grid, x_w_case++, x_e_case--) {

			if(x_w_case < 2)
				x_case = x_w_case;
			else if(x_e_case < 2)
				x_case = 4 - x_e_case;
			else
				x_case = 2;

			/* South side :  */
			kase = x_case * 5;
			ij = C->ij_sw_corner + i * C->my;
			u[ij + C->offset[kase][11]] = 
				(float)(u[ij + C->offset[kase][0]] + C->eps_m2*(u[ij + C->offset[kase][1]] + u[ij + C->offset[kase][3]]
					- u[ij + C->offset[kase][8]] - u[ij + C->offset[kase][10]])
					+ C->two_plus_em2 * (u[ij + C->offset[kase][9]] - u[ij + C->offset[kase][2]]) );
				/*  + tense * C->eps_m2 * (u[ij + C->offset[kase][2]] - u[ij + C->offset[kase][9]]) / (1.0 - tense);  */
			/* North side :  */
			kase = x_case * 5 + 4;
			ij = C->ij_nw_corner + i * C->my;
			u[ij + C->offset[kase][0]] = 
				-(float)(-u[ij + C->offset[kase][11]] + C->eps_m2 * (u[ij + C->offset[kase][1]] + u[ij + C->offset[kase][3]]
					- u[ij + C->offset[kase][8]] - u[ij + C->offset[kase][10]])
					+ C->two_plus_em2 * (u[ij + C->offset[kase][9]] - u[ij + C->offset[kase][2]]) );
				/*  - tense * C->eps_m2 * (u[ij + C->offset[kase][2]] - u[ij + C->offset[kase][9]]) / (1.0 - tense);  */
		}

		y_s_case = 0;
		y_n_case = C->block_ny - 1;
		for (j = 0; j < C->ny; j += C->grid, y_s_case++, y_n_case--) {

			if(y_s_case < 2)
				y_case = y_s_case;
			else if(y_n_case < 2)
				y_case = 4 - y_n_case;
			else
				y_case = 2;

			/* West side :  */
			kase = y_case;
			ij = C->ij_sw_corner + j;
			u[ij+C->offset[kase][4]] = 
				u[ij + C->offset[kase][7]] + (float)(C->eps_p2 * (u[ij + C->offset[kase][3]] + u[ij + C->offset[kase][10]]
				-u[ij + C->offset[kase][1]] - u[ij + C->offset[kase][8]])
				+ C->two_plus_ep2 * (u[ij + C->offset[kase][5]] - u[ij + C->offset[kase][6]]));
				/*  + tense * (u[ij + C->offset[kase][6]] - u[ij + C->offset[kase][5]]) / (1.0 - tense);  */
			/* East side :  */
			kase = 20 + y_case;
			ij = C->ij_se_corner + j;
			u[ij + C->offset[kase][7]] = 
				- (float)(-u[ij + C->offset[kase][4]] + C->eps_p2 * (u[ij + C->offset[kase][3]] + u[ij + C->offset[kase][10]]
				- u[ij + C->offset[kase][1]] - u[ij + C->offset[kase][8]])
				+ C->two_plus_ep2 * (u[ij + C->offset[kase][5]] - u[ij + C->offset[kase][6]]) );
				/*  - tense * (u[ij + C->offset[kase][6]] - u[ij + C->offset[kase][5]]) / (1.0 - tense);  */
		}



		/* That's it for the boundary points.  Now loop over all data  */

		x_w_case = 0;
		x_e_case = C->block_nx - 1;
		for (i = 0; i < C->nx; i += C->grid, x_w_case++, x_e_case--) {

			if(x_w_case < 2)
				x_case = x_w_case;
			else if(x_e_case < 2)
				x_case = 4 - x_e_case;
			else
				x_case = 2;

			y_s_case = 0;
			y_n_case = C->block_ny - 1;

			ij = C->ij_sw_corner + i * C->my;

			for (j = 0; j < C->ny; j += C->grid, ij += C->grid, y_s_case++, y_n_case--) {

				if (iu[ij] == 5) continue;	/* Point is fixed  */

				if(y_s_case < 2)
					y_case = y_s_case;
				else if(y_n_case < 2)
					y_case = 4 - y_n_case;
				else
					y_case = 2;

				kase = x_case * 5 + y_case;
				sum_ij = 0.0;

				if (iu[ij] == 0) {		/* Point is unconstrained  */
					for (k = 0; k < 12; k++) {
						sum_ij += (u[ij + C->offset[kase][k]] * C->coeff[0][k]);
					}
				}
				else {				/* Point is constrained  */

					b0 = C->briggs[briggs_index].b[0];
					b1 = C->briggs[briggs_index].b[1];
					b2 = C->briggs[briggs_index].b[2];
					b3 = C->briggs[briggs_index].b[3];
					b4 = C->briggs[briggs_index].b[4];
					b5 = C->briggs[briggs_index].b[5];
					briggs_index++;
					if (iu[ij] < 3) {
						if (iu[ij] == 1) {	/* Point is in quadrant 1  */
							busum = b0 * u[ij + C->offset[kase][10]]
								+ b1 * u[ij + C->offset[kase][9]]
								+ b2 * u[ij + C->offset[kase][5]]
								+ b3 * u[ij + C->offset[kase][1]];
						}
						else {			/* Point is in quadrant 2  */
							busum = b0 * u[ij + C->offset[kase][8]]
								+ b1 * u[ij + C->offset[kase][9]]
								+ b2 * u[ij + C->offset[kase][6]]
								+ b3 * u[ij + C->offset[kase][3]];
						}
					}
					else {
						if (iu[ij] == 3) {	/* Point is in quadrant 3  */
							busum = b0 * u[ij + C->offset[kase][1]]
								+ b1 * u[ij + C->offset[kase][2]]
								+ b2 * u[ij + C->offset[kase][6]]
								+ b3 * u[ij + C->offset[kase][10]];
						}
						else {		/* Point is in quadrant 4  */
							busum = b0 * u[ij + C->offset[kase][3]]
								+ b1 * u[ij + C->offset[kase][2]]
								+ b2 * u[ij + C->offset[kase][5]]
								+ b3 * u[ij + C->offset[kase][8]];
						}
					}
					for (k = 0; k < 12; k++) {
						sum_ij += (u[ij + C->offset[kase][k]] * C->coeff[1][k]);
					}
					sum_ij = (sum_ij + C->a0_const_2 * (busum + b5))
						/ (C->a0_const_1 + C->a0_const_2 * b4);
				}

				/* New relaxation here  */
				sum_ij = u[ij] * C->relax_old + sum_ij * C->relax_new;

				if (C->constrained) {	/* Must check limits.  Note lower/upper is v2 format and need ij_v2! */
					ij_v2 = (C->ny - j - 1) * C->nx + i;
					if (C->set_low && !GMT_is_fnan (C->lower[ij_v2]) && sum_ij < C->lower[ij_v2])
						sum_ij = C->lower[ij_v2];
					else if (C->set_high && !GMT_is_fnan (C->upper[ij_v2]) && sum_ij > C->upper[ij_v2])
						sum_ij = C->upper[ij_v2];
				}

				change = fabs(sum_ij - u[ij]);
				u[ij] = (float)sum_ij;
				if (change > max_change) max_change = change;
			}
		}
		iteration_count++;
		C->total_iterations++;
		max_change *= C->z_scale;	/* Put max_change into z units  */
		if (C->long_verbose) fprintf (stderr, C->format,
			C->grid, C->mode_type[mode], iteration_count, max_change, current_limit, C->total_iterations);

	} while (max_change > current_limit && iteration_count < C->max_iterations);

	if (gmtdefs.verbose && !C->long_verbose) fprintf (stderr, C->format,
		C->grid, C->mode_type[mode], iteration_count, max_change, current_limit, C->total_iterations);

	return (iteration_count);
}

void check_errors (struct SURFACE_INFO *C) {

	GMT_LONG	i, j, k, ij, move_over[12];
	char *iu;	/* move_over = C->offset[kase][12], but grid = 1 so move_over is easy  */

	double	x0, y0, dx, dy, mean_error, mean_squared_error, z_est, z_err, curvature, c;
	double	du_dx, du_dy, d2u_dx2, d2u_dxdy, d2u_dy2, d3u_dx3, d3u_dx2dy, d3u_dxdy2, d3u_dy3;

	double	x_0_const = 4.0 * (1.0 - C->boundary_tension) / (2.0 - C->boundary_tension);
	double	x_1_const = (3 * C->boundary_tension - 2.0) / (2.0 - C->boundary_tension);
	double	y_denom = 2 * C->l_epsilon * (1.0 - C->boundary_tension) + C->boundary_tension;
	double	y_0_const = 4 * C->l_epsilon * (1.0 - C->boundary_tension) / y_denom;
	double	y_1_const = (C->boundary_tension - 2 * C->l_epsilon * (1.0 - C->boundary_tension) ) / y_denom;
	float *u;

	u = C->u;
	iu = C->iu;
	
	move_over[0] = 2;
	move_over[1] = 1 - C->my;
	move_over[2] = 1;
	move_over[3] = 1 + C->my;
	move_over[4] = -2 * C->my;
	move_over[5] = -C->my;
	move_over[6] = C->my;
	move_over[7] = 2 * C->my;
	move_over[8] = -1 - C->my;
	move_over[9] = -1;
	move_over[10] = -1 + C->my;
	move_over[11] = -2;

	mean_error = 0;
	mean_squared_error = 0;

	/* First update the boundary values  */

	for (i = 0; i < C->nx; i ++) {
		ij = C->ij_sw_corner + i * C->my;
		u[ij - 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij + 1]);
		ij = C->ij_nw_corner + i * C->my;
		u[ij + 1] = (float)(y_0_const * u[ij] + y_1_const * u[ij - 1]);
	}

	for (j = 0; j < C->ny; j ++) {
		ij = C->ij_sw_corner + j;
		u[ij - C->my] = (float)(x_1_const * u[ij + C->my] + x_0_const * u[ij]);
		ij = C->ij_se_corner + j;
		u[ij + C->my] = (float)(x_1_const * u[ij - C->my] + x_0_const * u[ij]);
	}

	ij = C->ij_sw_corner;
	u[ij - C->my - 1] = u[ij + C->my - 1] + u[ij - C->my + 1] - u[ij + C->my + 1];
	ij = C->ij_nw_corner;
	u[ij - C->my + 1] = u[ij + C->my + 1] + u[ij - C->my - 1] - u[ij + C->my - 1];
	ij = C->ij_se_corner;
	u[ij + C->my - 1] = u[ij - C->my - 1] + u[ij + C->my + 1] - u[ij - C->my + 1];
	ij = C->ij_ne_corner;
	u[ij + C->my + 1] = u[ij - C->my + 1] + u[ij + C->my - 1] - u[ij - C->my - 1];

	for (i = 0; i < C->nx; i ++) {

		ij = C->ij_sw_corner + i * C->my;
		u[ij + move_over[11]] = 
			(float)(u[ij + move_over[0]] + C->eps_m2*(u[ij + move_over[1]] + u[ij + move_over[3]]
				- u[ij + move_over[8]] - u[ij + move_over[10]])
				+ C->two_plus_em2 * (u[ij + move_over[9]] - u[ij + move_over[2]]) );

		ij = C->ij_nw_corner + i * C->my;
		u[ij + move_over[0]] = 
			-(float)(-u[ij + move_over[11]] + C->eps_m2 * (u[ij + move_over[1]] + u[ij + move_over[3]]
				- u[ij + move_over[8]] - u[ij + move_over[10]])
				+ C->two_plus_em2 * (u[ij + move_over[9]] - u[ij + move_over[2]]) );
	}

	for (j = 0; j < C->ny; j ++) {

		ij = C->ij_sw_corner + j;
		u[ij+move_over[4]] = 
			u[ij + move_over[7]] + (float)(C->eps_p2 * (u[ij + move_over[3]] + u[ij + move_over[10]]
			-u[ij + move_over[1]] - u[ij + move_over[8]])
			+ C->two_plus_ep2 * (u[ij + move_over[5]] - u[ij + move_over[6]]));

		ij = C->ij_se_corner + j;
		u[ij + move_over[7]] = 
			- (float)(-u[ij + move_over[4]] + C->eps_p2 * (u[ij + move_over[3]] + u[ij + move_over[10]]
			- u[ij + move_over[1]] - u[ij + move_over[8]])
			+ C->two_plus_ep2 * (u[ij + move_over[5]] - u[ij + move_over[6]]) );
	}

	/* That resets the boundary values.  Now we can test all data.  
		Note that this loop checks all values, even though only nearest were used.  */

	for (k = 0; k < C->npoints; k++) {
		i = C->data[k].index/C->ny;
		j = C->data[k].index%C->ny;
	 	ij = C->ij_sw_corner + i * C->my + j;
	 	if ( iu[ij] == 5 ) continue;
	 	x0 = C->h.x_min + i*C->h.x_inc;
	 	y0 = C->h.y_min + j*C->h.y_inc;
	 	dx = (C->data[k].x - x0)*C->r_xinc;
	 	dy = (C->data[k].y - y0)*C->r_yinc;
 
	 	du_dx = 0.5 * (u[ij + move_over[6]] - u[ij + move_over[5]]);
	 	du_dy = 0.5 * (u[ij + move_over[2]] - u[ij + move_over[9]]);
	 	d2u_dx2 = u[ij + move_over[6]] + u[ij + move_over[5]] - 2 * u[ij];
	 	d2u_dy2 = u[ij + move_over[2]] + u[ij + move_over[9]] - 2 * u[ij];
	 	d2u_dxdy = 0.25 * (u[ij + move_over[3]] - u[ij + move_over[1]]
	 			- u[ij + move_over[10]] + u[ij + move_over[8]]);
	 	d3u_dx3 = 0.5 * ( u[ij + move_over[7]] - 2 * u[ij + move_over[6]]
	 				+ 2 * u[ij + move_over[5]] - u[ij + move_over[4]]);
	 	d3u_dy3 = 0.5 * ( u[ij + move_over[0]] - 2 * u[ij + move_over[2]]
	 				+ 2 * u[ij + move_over[9]] - u[ij + move_over[11]]);
	 	d3u_dx2dy = 0.5 * ( ( u[ij + move_over[3]] + u[ij + move_over[1]] - 2 * u[ij + move_over[2]] )
	 				- ( u[ij + move_over[10]] + u[ij + move_over[8]] - 2 * u[ij + move_over[9]] ) );
	 	d3u_dxdy2 = 0.5 * ( ( u[ij + move_over[3]] + u[ij + move_over[10]] - 2 * u[ij + move_over[6]] )
	 				- ( u[ij + move_over[1]] + u[ij + move_over[8]] - 2 * u[ij + move_over[5]] ) );

	 	/* 3rd order Taylor approx:  */
	 
	 	z_est = u[ij] + dx * (du_dx +  dx * ( (0.5 * d2u_dx2) + dx * (d3u_dx3 / 6.0) ) )
				+ dy * (du_dy +  dy * ( (0.5 * d2u_dy2) + dy * (d3u_dy3 / 6.0) ) )
	 			+ dx * dy * (d2u_dxdy) + (0.5 * dx * d3u_dx2dy) + (0.5 * dy * d3u_dxdy2);
	 
	 	z_err = z_est - C->data[k].z;
	 	mean_error += z_err;
	 	mean_squared_error += (z_err * z_err);
	 }
	 mean_error /= C->npoints;
	 mean_squared_error = sqrt( mean_squared_error / C->npoints);
	 
	 curvature = 0.0;
	 
	 for (i = 0; i < C->nx; i++) {
	 	for (j = 0; j < C->ny; j++) {
	 		ij = C->ij_sw_corner + i * C->my + j;
	 		c = u[ij + move_over[6]] + u[ij + move_over[5]]
	 			+ u[ij + move_over[2]] + u[ij + move_over[9]] - 4.0 * u[ij + move_over[6]];
			curvature += (c * c);
		}
	}

	 fprintf (stderr, "%s: Fit info: N data points  N nodes\tmean error\trms error\tcurvature\n", GMT_program);
	 sprintf (C->format,"%s:\t%%8ld\t%%8ld\t%s\t%s\t%s\n", GMT_program, gmtdefs.d_format, gmtdefs.d_format, gmtdefs.d_format);
	 fprintf (stderr, C->format, C->npoints, C->nxny, mean_error, mean_squared_error, curvature);
 }

void remove_planar_trend (struct SURFACE_INFO *C)
{

	GMT_LONG	i;
	double	a, b, c, d, xx, yy, zz;
	double	sx, sy, sz, sxx, sxy, sxz, syy, syz;

	sx = sy = sz = sxx = sxy = sxz = syy = syz = 0.0;

	for (i = 0; i < C->npoints; i++) {

		xx = (C->data[i].x - C->h.x_min) * C->r_xinc;
		yy = (C->data[i].y - C->h.y_min) * C->r_yinc;
		zz = C->data[i].z;

		sx += xx;
		sy += yy;
		sz += zz;
		sxx +=(xx * xx);
		sxy +=(xx * yy);
		sxz +=(xx * zz);
		syy +=(yy * yy);
		syz +=(yy * zz);
	}

	d = C->npoints*sxx*syy + 2*sx*sy*sxy - C->npoints*sxy*sxy - sx*sx*syy - sy*sy*sxx;

	if (d == 0.0) {
		C->plane_c0 = C->plane_c1 = C->plane_c2 = 0.0;
		return;
	}

	a = sz*sxx*syy + sx*sxy*syz + sy*sxy*sxz - sz*sxy*sxy - sx*sxz*syy - sy*syz*sxx;
	b = C->npoints*sxz*syy + sz*sy*sxy + sy*sx*syz - C->npoints*sxy*syz - sz*sx*syy - sy*sy*sxz;
	c = C->npoints*sxx*syz + sx*sy*sxz + sz*sx*sxy - C->npoints*sxy*sxz - sx*sx*syz - sz*sy*sxx;

	C->plane_c0 = a / d;
	C->plane_c1 = b / d;
	C->plane_c2 = c / d;

	for (i = 0; i < C->npoints; i++) {

		xx = (C->data[i].x - C->h.x_min) * C->r_xinc;
		yy = (C->data[i].y - C->h.y_min) * C->r_yinc;

		C->data[i].z -= (float)(C->plane_c0 + C->plane_c1 * xx + C->plane_c2 * yy);
	}
}

void replace_planar_trend (struct SURFACE_INFO *C)
{
	GMT_LONG	i, j, ij;
	float *u;
	u = C->u;

	 for (i = 0; i < C->nx; i++) {
	 	for (j = 0; j < C->ny; j++) {
	 		ij = C->ij_sw_corner + i * C->my + j;
	 		u[ij] = (float)((u[ij] * C->z_scale) + (C->plane_c0 + C->plane_c1 * i + C->plane_c2 * j));
		}
	}
}

void throw_away_unusables (struct SURFACE_INFO *C)
{
	/* This is a new routine to eliminate data which will become
		unusable on the final iteration, when grid = 1.
		It assumes grid = 1 and set_grid_parameters has been
		called.  We sort, mark redundant data as SURFACE_OUTSIDE, and
		sort again, chopping off the excess.

		Experimental modification 5 Dec 1988 by Smith, as part
		of a new implementation using core memory for b[6]
		coefficients, eliminating calls to temp file.
	*/

	GMT_LONG	last_index, n_outside, k;

	/* Sort the data  */

	qsort ((void *)C->data, (size_t)C->npoints, sizeof (struct SURFACE_DATA), compare_points);

	/* If more than one datum is indexed to same node, only the first should be kept.
		Mark the additional ones as SURFACE_OUTSIDE
	*/
	last_index = -1;
	n_outside = 0;
	for (k = 0; k < C->npoints; k++) {
		if (C->data[k].index == last_index) {
			C->data[k].index = SURFACE_OUTSIDE;
			n_outside++;
		}
		else {
			last_index = C->data[k].index;
		}
	}
	/* Sort again; this time the SURFACE_OUTSIDE points will be thrown away  */

	qsort ((void *)C->data, (size_t)C->npoints, sizeof (struct SURFACE_DATA), compare_points);
	C->npoints -= n_outside;
	C->data = (struct SURFACE_DATA *) GMT_memory ((void *)C->data, (size_t)C->npoints, sizeof(struct SURFACE_DATA), GMT_program);
	if (gmtdefs.verbose && (n_outside)) {
		fprintf (stderr, "%s: %ld unusable points were supplied; these will be ignored.\n", GMT_program, n_outside);
		fprintf (stderr, "%s: You should have pre-processed the data with block-mean, -median, or -mode.\n", GMT_program);
	}
}

void rescale_z_values (struct SURFACE_INFO *C)
{
	GMT_LONG	i;
	double	ssz = 0.0;

	for (i = 0; i < C->npoints; i++) ssz += (C->data[i].z * C->data[i].z);

	/* Set z_scale = rms(z):  */
	C->z_scale = sqrt (ssz / C->npoints);

	if (C->z_scale == 0.0) {
/*		fprintf (stderr, "%s: WARNING: Input data lie exactly on a plane - no solution (aborting).\n", GMT_program); */
		C->r_z_scale = C->z_scale = 1.0;
	}
	else
		C->r_z_scale = 1.0 / C->z_scale;

	for (i = 0; i < C->npoints; i++) C->data[i].z *= (float)C->r_z_scale;

	if (C->converge_limit == 0.0) C->converge_limit = 0.001 * C->z_scale; /* i.e., 1 ppt of L2 scale */
}

void load_constraints (struct SURFACE_INFO *C, BOOLEAN transform)
{
	GMT_LONG i, j, ij;
	double yy;
	struct GRD_HEADER hdr;

	/* Load lower/upper limits, verify range, deplane, and rescale */

	if (C->set_low > 0) {
		C->lower = (float *) GMT_memory (VNULL, (size_t)C->nxny, sizeof (float), GMT_program);
		if (C->set_low < 3)
			for (i = 0; i < C->nxny; i++) C->lower[i] = (float)C->low_limit;
		else {
			GMT_err_fail (GMT_read_grd_info (C->low_file, &hdr), C->low_file);
			if (hdr.nx != C->nx || hdr.ny != C->ny) {
				/*fprintf (stderr, "%s: lower limit file not of proper dimension!\n", GMT_program);*/
				exit (EXIT_FAILURE);
			}
			GMT_err_fail (GMT_read_grd (C->low_file, &hdr, C->lower, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), C->low_file);
		}
		if (transform) {
			for (j = ij = 0; j < C->ny; j++) {
				yy = C->ny - j - 1;
				for (i = 0; i < C->nx; i++, ij++) {
					if (GMT_is_fnan (C->lower[ij])) continue;
					C->lower[ij] -= (float)(C->plane_c0 + C->plane_c1 * i + C->plane_c2 * yy);
					C->lower[ij] *= (float)C->r_z_scale;
				}
			}
		}
		C->constrained = TRUE;
	}
	if (C->set_high > 0) {
		C->upper = (float *) GMT_memory (VNULL, (size_t)C->nxny, sizeof (float), GMT_program);
		if (C->set_high < 3)
			for (i = 0; i < C->nxny; i++) C->upper[i] = (float)C->high_limit;
		else {
			GMT_err_fail (GMT_read_grd_info (C->high_file, &hdr), C->high_file);
			if (hdr.nx != C->nx || hdr.ny != C->ny) {
				/*fprintf (stderr, "%s: upper limit file not of proper dimension!\n", GMT_program); */
				exit (EXIT_FAILURE);
			}
			GMT_err_fail (GMT_read_grd (C->high_file, &hdr, C->upper, 0.0, 0.0, 0.0, 0.0, GMT_pad, FALSE), C->high_file);
		}
		if (transform) {
			for (j = ij = 0; j < C->ny; j++) {
				yy = C->ny - j - 1;
				for (i = 0; i < C->nx; i++, ij++) {
					if (GMT_is_fnan (C->upper[ij])) continue;
					C->upper[ij] -= (float)(C->plane_c0 + C->plane_c1 * i + C->plane_c2 * yy);
					C->upper[ij] *= (float)C->r_z_scale;
				}
			}
		}
		C->constrained = TRUE;
	}
}

double guess_surface_time (int factors[], int nx, int ny)
{
	/* Routine to guess a number proportional to the operations
	 * required by surface working on a user-desired grid of
	 * size nx by ny, where nx = (x_max - x_min)/dx, and same for
	 * ny.  (That is, one less than actually used in routine.)
	 *
	 * This is based on the following untested conjecture:
	 * 	The operations are proportional to T = nxg*nyg*L,
	 *	where L is a measure of the distance that data
	 *	constraints must propagate, and nxg, nyg are the
	 * 	current size of the grid.
	 *	For nx,ny relatively prime, we will go through only
	 * 	one grid cycle, L = max(nx,ny), and T = nx*ny*L.
	 *	But for nx,ny whose greatest common divisor is a highly
	 * 	composite number, we will have L equal to the division
	 * 	step made at each new grid cycle, and nxg,nyg will
	 * 	also be smaller than nx,ny.  Thus we can hope to find
	 *	some nx,ny for which the total value of T is C->small.
	 *
	 * The above is pure speculation and has not been derived
	 * empirically.  In actual practice, the distribution of the
	 * data, both spatially and in terms of their values, will
	 * have a strong effect on convergence.
	 *
	 * W. H. F. Smith, 26 Feb 1992.  */

	GMT_LONG	gcd;		/* Current value of the gcd  */
	int	nxg, nyg;	/* Current value of the grid dimensions  */
	int	nfactors = 0;	/* Number of prime factors of current gcd  */
	int	factor;		/* Currently used factor  */
	/* Doubles are used below, even though the values will be integers,
		because the multiplications might reach sizes of O(n**3)  */
	double	t_sum;		/* Sum of values of T at each grid cycle  */
	double	length;		/* Current propagation distance.  */


	gcd = gcd_euclid((GMT_LONG)nx, (GMT_LONG)ny);
	if (gcd > 1) {
		nfactors = get_prime_factors(gcd, factors);
		nxg = nx/gcd;
		nyg = ny/gcd;
		if (nxg < 3 || nyg < 3) {
			factor = factors[nfactors - 1];
			nfactors--;
			gcd /= factor;
			nxg *= factor;
			nyg *= factor;
		}
	}
	else {
		nxg = nx;
		nyg = ny;
	}
	length = MAX(nxg, nyg);
	t_sum = nxg * (nyg * length);	/* Make it double at each multiply  */

	/* Are there more grid cycles ?  */
	while (gcd > 1) {
		factor = factors[nfactors - 1];
		nfactors--;
		gcd /= factor;
		nxg *= factor;
		nyg *= factor;
		length = factor;
		t_sum += nxg * (nyg * length);
	}
	return(t_sum);
}

void suggest_sizes_for_surface (int factors[], int nx, int ny)
{
	/* Calls guess_surface_time for a variety of trial grid
	 * sizes, where the trials are highly composite numbers
	 * with lots of factors of 2, 3, and 5.  The sizes are
	 * within the range (nx,ny) - (2*nx, 2*ny).  Prints to
	 * stderr the values which are an improvement over the
	 * user's original nx,ny.
	 * Should be called with nx=(x_max-x_min)/dx, and ditto
	 * for ny; that is, one smaller than the lattice used
	 * in surface.c
	 *
	 * W. H. F. Smith, 26 Feb 1992.  */

	double	guess_surface_time(int factors[], int nx, int ny);
	double	users_time;	/* Time for user's nx, ny  */
	double	current_time;	/* Time for current nxg, nyg  */
	int	i;
	int	nxg, nyg;	/* Guessed by this routine  */
	int	nx2, ny2, nx3, ny3, nx5, ny5;	/* For powers  */
	int	xstop, ystop;	/* Set to 2*nx, 2*ny  */
	int	n_sug = 0;	/* N of suggestions found  */
	int	compare_sugs(const void *point_1, const void *point_2);	/* Sort suggestions decreasing  */
	struct SURFACE_SUGGESTION *sug = NULL;

	users_time = guess_surface_time (factors, nx, ny);
	xstop = 2*nx;
	ystop = 2*ny;

	for (nx2 = 2; nx2 <= xstop; nx2 *= 2) {
	  for (nx3 = 1; nx3 <= xstop; nx3 *= 3) {
	    for (nx5 = 1; nx5 <= xstop; nx5 *= 5) {
		nxg = nx2 * nx3 * nx5;
		if (nxg < nx || nxg > xstop) continue;

		for (ny2 = 2; ny2 <= ystop; ny2 *= 2) {
		  for (ny3 = 1; ny3 <= ystop; ny3 *= 3) {
		    for (ny5 = 1; ny5 <= ystop; ny5 *= 5) {
			nyg = ny2 * ny3 * ny5;
			if (nyg < ny || nyg > ystop) continue;

			current_time = guess_surface_time (factors, nxg, nyg);
			if (current_time < users_time) {
				n_sug++;
				sug = (struct SURFACE_SUGGESTION *)GMT_memory ((void *)sug, (size_t)n_sug, sizeof(struct SURFACE_SUGGESTION), GMT_program);
				sug[n_sug-1].nx = nxg;
				sug[n_sug-1].ny = nyg;
				sug[n_sug-1].factor = users_time/current_time;
			}

		    }
		  }
		}
	    }
	  }
	}

	if (n_sug) {
		qsort((void *)sug, (size_t)n_sug, sizeof(struct SURFACE_SUGGESTION), compare_sugs);
		for (i = 0; i < n_sug && i < 10; i++) {
			fprintf (stderr, "%s:  HINT:  Choosing nx = %d, ny = %d might cut run time by a factor of %.8g\n",
				GMT_program, sug[i].nx, sug[i].ny, sug[i].factor);
		}
		GMT_free ((void *)sug);
	}
	else {
		/*fprintf (stderr, "%s: Cannot suggest any nx,ny better than your -R -I define.\n", GMT_program);*/

	}
	return;
}

int compare_sugs (const void *point_1, const void *point_2)
{
	/* Sorts sugs into DESCENDING order!  */
	if ( ((struct SURFACE_SUGGESTION *)point_1)->factor < ((struct SURFACE_SUGGESTION *)point_2)->factor)
		return(1);
	else if ( ((struct SURFACE_SUGGESTION *)point_1)->factor > ((struct SURFACE_SUGGESTION *)point_2)->factor)
		return(-1);
	else
		return(0);
}

int get_prime_factors (GMT_LONG n, int *f)
{
	/* Fills the integer array f with the prime factors of n.
	 * Returns the number of locations filled in f, which is
	 * one if n is prime.
	 *
	 * f[] should have been malloc'ed to enough space before
	 * calling prime_factors().  We can be certain that f[32]
	 * is enough space, for if n fits in a long, then n < 2**32,
	 * and so it must have fewer than 32 prime factors.  I think
	 * that in general, ceil(log2((double)n)) is enough storage
	 * space for f[].
	 *
	 * Tries 2,3,5 explicitly; then alternately adds 2 or 4
	 * to the previously tried factor to obtain the next trial
	 * factor.  This is done with the variable two_four_toggle.
	 * With this method we try 7,11,13,17,19,23,25,29,31,35,...
	 * up to a maximum of sqrt(n).  This shortened list results
	 * in 1/3 fewer divisions than if we simply tried all integers
	 * between 5 and sqrt(n).  We can reduce the size of the list
	 * of trials by an additional 20% by removing the multiples
	 * of 5, which are equal to 30m +/- 5, where m >= 1.  Starting
	 * from 25, these are found by alternately adding 10 or 20.
	 * To do this, we use the variable ten_twenty_toggle.
	 *
	 * W. H. F. Smith, 26 Feb 1992, after D.E. Knuth, vol. II  */

	GMT_LONG	current_factor;	/* The factor currently being tried  */
	GMT_LONG	max_factor;	/* Don't try any factors bigger than this  */
	int	n_factors = 0;	/* Returned; one if n is prime  */
	GMT_LONG	two_four_toggle = 0;	/* Used to add 2 or 4 to get next trial factor  */
	GMT_LONG	ten_twenty_toggle = 0;	/* Used to add 10 or 20 to skip_five  */
	GMT_LONG	skip_five = 25;	/* Used to skip multiples of 5 in the list  */
	GMT_LONG	m;	/* Used to keep a working copy of n  */


	/* Initialize m and max_factor  */
#ifdef __LP64__
	m = labs(n);
#else
	m = abs(n);
#endif
	if (m < 2) return(0);
	max_factor = (int)floor(sqrt((double)m));

	/* First find the 2s  */
	current_factor = 2;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Next find the 3s  */
	current_factor = 3;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Next find the 5s  */
	current_factor = 5;
	while(!(m%current_factor)) {
		m /= current_factor;
		f[n_factors] = current_factor;
		n_factors++;
	}
	if (m == 1) return(n_factors);

	/* Now try all the rest  */

	while (m > 1 && current_factor <= max_factor) {

		/* Current factor is either 2 or 4 more than previous value  */

		if (two_four_toggle) {
			current_factor += 4;
			two_four_toggle = 0;
		}
		else {
			current_factor += 2;
			two_four_toggle = 1;
		}

		/* If current factor is a multiple of 5, skip it.  But first,
			set next value of skip_five according to 10/20 toggle:  */

		if (current_factor == skip_five) {
			if (ten_twenty_toggle) {
				skip_five += 20;
				ten_twenty_toggle = 0;
			}
			else {
				skip_five += 10;
				ten_twenty_toggle = 1;
			}
			continue;
		}

		/* Get here when current_factor is not a multiple of 2,3 or 5:  */

		while(!(m%current_factor)) {
			m /= current_factor;
			f[n_factors] = current_factor;
			n_factors++;
		}
	}

	/* Get here when all factors up to floor(sqrt(n)) have been tried.  */

	if (m > 1) {
		/* m is an additional prime factor of n  */
		f[n_factors] = m;
		n_factors++;
	}
	return (n_factors);
}
/* gcd_euclid.c  Greatest common divisor routine  */

#define IABS(i)	(((i) < 0) ? -(i) : (i))

GMT_LONG gcd_euclid (GMT_LONG a, GMT_LONG b)
{
	/* Returns the greatest common divisor of u and v by Euclid's method.
	 * I have experimented also with Stein's method, which involves only
	 * subtraction and left/right shifting; Euclid is faster, both for
	 * integers of size 0 - 1024 and also for random integers of a size
	 * which fits in a long integer.  Stein's algorithm might be better
	 * when the integers are HUGE, but for our purposes, Euclid is fine.
	 *
	 * Walter H. F. Smith, 25 Feb 1992, after D. E. Knuth, vol. II  */

	int	u,v,r;

	u = MAX(IABS(a), IABS(b));
	v = MIN(IABS(a), IABS(b));

	while (v > 0) {
		r = u%v;	/* Knuth notes that u < 2v 40% of the time;  */
		u = v;		/* thus we could have tried a subtraction  */
		v = r;		/* followed by an if test to do r = u%v  */
	}
	return(u);
}

void load_parameters (struct SURFACE_INFO *C, struct SURFACE_CTRL *Ctrl)
{
	if (Ctrl->S.active) {
		if (Ctrl->S.unit == 'm' || Ctrl->S.unit == 'M') Ctrl->S.radius /= 60.0;
		if (Ctrl->S.unit == 'c' || Ctrl->S.unit == 'M') Ctrl->S.radius /= 3600.0;
	}
	C->radius = Ctrl->S.radius;
	C->h.x_inc = Ctrl->I.xinc;
	C->h.y_inc = Ctrl->I.yinc;
	C->relax_new = Ctrl->Z.value;
	C->max_iterations = Ctrl->N.value;
	C->radius = Ctrl->S.radius;
	C->low_file = Ctrl->L.low;
	C->high_file = Ctrl->L.high;
	C->set_low = Ctrl->L.lmode;
	C->low_limit = Ctrl->L.min;
	C->set_high = Ctrl->L.hmode;
	C->high_limit = Ctrl->L.max;
	C->boundary_tension = Ctrl->T.b_tension;
	C->interior_tension = Ctrl->T.i_tension;
	C->l_epsilon = Ctrl->A.value;
	C->converge_limit = Ctrl->C.value;
}

void *New_Surface_Ctrl () {	/* Allocate and initialize a new control structure */
	struct SURFACE_CTRL *C;
	
	C = (struct SURFACE_CTRL *) GMT_memory (VNULL, (size_t)1, sizeof (struct SURFACE_CTRL), "New_Surface_Ctrl");
	
	/* Initialize values whose defaults are not 0/FALSE/NULL */
	C->N.value = 250;
	C->A.value = 1.0;
	C->Z.value = 1.4;
		
	return ((void *)C);
}

void Free_Surface_Ctrl (struct SURFACE_CTRL *C) {	/* Deallocate control structure */
	if (C->G.file) 
		free ((void *)C->G.file);	 
	if (C->L.low) 
		free ((void *)C->L.low);	
	if (C->L.high) 
		free ((void *)C->L.high);	 
	GMT_free ((void *)C);	  /* Causing segmentation fault */
}

