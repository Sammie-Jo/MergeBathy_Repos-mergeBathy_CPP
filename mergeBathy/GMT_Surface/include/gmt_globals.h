/*--------------------------------------------------------------------
 *	$Id: gmt_globals.h,v 1.33 2008/04/15 15:55:29 remko Exp $
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
 *  This include file contains initializations of global GMT variables.
 *  MUST BE INCUDED AFTER gmt.h
 *
 * Author:	Paul Wessel
 * Date:	21-AUG-1995
 * Revised:	22-MAR-2006
 * Version:	4.1
 *
 */

#ifdef DEBUG
struct MEMORY_TRACKER *GMT_mem_keeper;
#endif

struct GMT_CTRL *GMT;			/* Pointer to struct with common option settings */

struct GMT_DEFAULTS gmtdefs = {		/* Initial default values */
	20.0,			/* ANNOT_MIN_ANGLE */
	0.0,			/* ANNOT_MIN_SPACING */
	{0, 0},			/* ANNOT_FONT_PRIMARY, ANNOT_FONT_SECONDARY*/
	{14, 16},		/* ANNOT_FONT_SIZE_PRIMARY, ANNOT_FONT_SIZE_SECONDARY */
	{0.075,	0.075},		/* ANNOT_OFFSET_PRIMARY, ANNOT_OFFSET_SECONDARY */
	"WESN",			/* BASEMAP_AXES */
	{0, 0, 0},		/* BASEMAP_FRAME_RGB */
	0,			/* BASEMAP_TYPE */
	{0, 0, 0},		/* COLOR_BACKGROUND */
	{255, 255, 255},	/* COLOR_FOREGROUND */
	{128, 128, 128},	/* COLOR_NAN */
	0,			/* COLOR_IMAGE */
	3,			/* COLOR_MODEL */
	"%g",			/* D_FORMAT */
	0,			/* DEGREE_FORMAT */
	300,			/* DOTS_PR_INCH */
	0,			/* ELLIPSOID */
	{1.25, 0.0, {0, 0, 0}, ""},	/* FRAME_PEN */
	0.075,			/* FRAME_WIDTH */
	1.0,			/* GLOBAL_X_SCALE */
	1.0,			/* GLOBAL_Y_SCALE */
	{0.0, 0.0},		/* GRID_CROSS_SIZE P/S */
	"nf",			/* GRID_FORMAT */
	{
	  {0.25, 0.0, {0, 0, 0}, ""},	/* GRID_PEN_PRIMARY */
	  {0.5,  0.0, {0, 0, 0}, ""}	/* GRID_PEN_SECONDARY */
	},
	FALSE,			/* GRIDFILE_SHORTHAND */
	0,			/* HEADER_FONT */
	36,			/* HEADER_FONT_SIZE */
	0.1875,			/* HEADER_OFFSET */
	1.0,			/* HSV_MIN_SATURATION */
	0.1,			/* HSV_MAX_SATURATION */
	0.3,			/* HSV_MIN_VALUE */
	1.0,			/* HSV_MAX_VALUE */
	1,			/* INTERPOLANT */
	{ FALSE, FALSE },	/* IO_HEADER [0|1] */
	1,			/* N_HEADER_RECS (if -H is set) */
	0,			/* LABEL_FONT */
	24,			/* LABEL_FONT_SIZE */
	0.1125,			/* LABEL_OFFSET */
	0.01,
	-1.0,			/* MAP_SCALE_FACTOR */
	0.075,			/* MAP_SCALE_HEIGHT */
	1,			/* MEASURE_UNIT set to inch */
	25,			/* MEDIA set to Letter */
	1,			/* N_COPIES */
	1,			/* OBLIQUE_ANNOTATION */
	{255, 255, 255},	/* PAGE_COLOR */
	FALSE,			/* PORTRAIT */
	{612, 792},		/* PAPER_MEDIA (US Letter) */
	{85.0, 90.0},		/* POLAR_CAP (controls gridlines near poles) */
	0,			/* PS_COLOR (2 = HSV, 1 = CMYK, 0 = RGB) */
	0,			/* PS_IMAGE_COMPRESS (0 = NONE, 1 = RLE, 2 = LZW) */
	TRUE,			/* PS_IMAGE_FORMAT (TRUE = HEX, FALSE = BIN) */
	0,			/* PS_LINE_CAP (0 = butt, 1 = round, 2 = square) */
	0,			/* PS_LINE_JOIN (0 = miter, 1 = arc, 2 = bevel) */
	0,			/* PS_MITER_LIMIT (0 = Default, or 1-180) */
	FALSE,			/* PS_VERBOSE (TRUE = write comments, FALSE = no comments) */
	0.075,			/* TICK_LENGTH */
	{0.5, 0.0, {0, 0, 0}, ""},	/* TICK_PEN */
	FALSE,			/* UNIX_TIME */
	1, {-0.75, -0.75},	/* UNIX_TIME_POS */
	"%Y %b %d %H:%M:%S",	/* UNIX_TIME_FORMAT */
	0.0,			/* VECTOR_SHAPE */
	FALSE,			/* VERBOSE */
	FALSE,			/* WANT_EURO_FONT */
	9.0,			/* X_AXIS_LENGTH */
	6.0,			/* Y_AXIS_LENGTH */
	1.0,			/* X_ORIGIN */
	1.0,			/* Y_ORIGIN */
	{FALSE,	FALSE},		/* XY_TOGGLE */
	0,			/* Y_AXIS_TYPE */
	{	/* Ellipsoid structure and its members */
#include "gmt_ellipsoids.h"	/* This is created by GNUmakefile - do not edit it manually */
	},
	{	/* Datum structure and its members */
#include "gmt_datums.h"		/* This is created by GNUmakefile - do not edit it manually */
	},
	"hh:mm:ss",		/* INPUT_CLOCK_FORMAT */
	"yyyy-mm-dd",		/* INPUT_DATE_FORMAT */
	"hh:mm:ss",		/* OUTPUT_CLOCK_FORMAT */
	"yyyy-mm-dd",		/* OUTPUT_DATE_FORMAT */
	"+D",			/* OUTPUT_DEGREE_FORMAT */
	"hh:mm:ss",		/* PLOT_CLOCK_FORMAT */
	"yyyy-mm-dd",		/* PLOT_DATE_FORMAT */
	"+ddd:mm:ss",		/* PLOT_DEGREE_FORMAT */
	{ "full", "full" },	/* TIME_FORMAT_PRIMARY, TIME_FORMAT_SECONDARY */
	FALSE,			/* TIME_IS_INTERVAL */
	0.5,			/* TIME_INTERVAL_FRACTION */
	FALSE,			/* WANT_LEAP_SECONDS */
	{ "2000-01-01T12:00:00", 'd',		/* TIME_EPOCH, TIME_UNIT (J2000) */
	730120, 0.5, 86400.0, 1.0/86400.0 },	/* rata_die, epoch_t0, scale, i_scale based on the above */
	0,			/* TIME_WEEK_START */
	"us",			/* TIME_LANGUAGE */
	1950,			/* Y2K_OFFSET_YEAR */
	"\t",			/* FIELD_DELIMITER */
	gmt_ring,		/* DEGREE_SYMBOL */
	{ "Standard", 		/* CHAR_ENCODING */
	  { 32, 32, 32, 32, 32 } }, /* PostScript codes for degree, ring, colon, singlequote, and doublequote [Initialized to space] */
	TRUE,			/* HISTORY */
};

struct GMT_FONT *GMT_font;		/* Name and height of fonts recognized by GMT */
int GMT_N_FONTS;				/* Total number of fonts returned by GMT_init_fonts */

struct GMT_HASH GMT_rgb_hashnode[GMT_N_COLOR_NAMES];/* Used to translate colornames to r/g/b */

struct GMT_HASH GMT_month_hashnode[12];		/* Used to translate months to 1-12 */

/*--------------------------------------------------------------------*/
/*	For plotting purposes */
/*--------------------------------------------------------------------*/

struct GMT_PS GMT_ps;				/* Hold parameters related to PS setup */
struct GMT_PLOT_FRAME frame_info;
struct GMT_PLOT_CALCLOCK GMT_plot_calclock;
struct GMT_TRUNCATE_TIME GMT_truncate_time;
double *GMT_x_plot = 0;				/* Holds the x/y (inches) of a line to be plotted */
double *GMT_y_plot = 0;
int *GMT_pen = 0;				/* Pen (3 = up, 2 = down) for these points */
GMT_LONG GMT_n_plot = 0;				/* Number of such points */
GMT_LONG GMT_n_alloc = 0;				/* Size of allocated plot arrays */
int GMT_x_status_new;				/* Tells us what quadrant old and new points are in */
int GMT_y_status_new;
int GMT_x_status_old;
int GMT_y_status_old;
int GMT_corner = 0;
BOOLEAN GMT_on_border_is_outside = FALSE;	/* TRUE if a point exactly on the map border shoud be considered outside the map */
BOOLEAN GMT_world_map = FALSE;			/* TRUE if map has 360 degrees of longitude range */
BOOLEAN GMT_world_map_tm = FALSE;		/* TRUE if GMT_TM map is global? */
double GMT_map_width;				/* Full width in inches of this world map */
double GMT_map_height;				/* Full height in inches of this world map */
double GMT_half_map_size;			/* Half width in inches of this world map */
double GMT_half_map_height;			/* Half height of this world map */
PFI GMT_outside;				/* pointer to function checking if a lon/lat point is outside map */
PFI GMT_crossing;				/* pointer to functions returning crossover point at boundary */
PFI GMT_overlap;				/* pointer to function checking for overlap between 2 regions */
PFI GMT_map_clip;				/* pointer to functions that clip a polygon to fit inside map */
PFD GMT_left_edge, GMT_right_edge;		/* pointers to functions that return left/right edge of map */
PFD GMT_distance_func;				/* pointer to function returning distance between two points points */
BOOLEAN GMT_z_periodic = FALSE;			/* TRUE if grid values are 0-360 degrees (phases etc) */
PFI GMT_wrap_around_check;			/* Does x or y wrap checks */
PFI GMT_map_jump;				/* TRUE if we jump in x or y */
PFB GMT_will_it_wrap;				/* TRUE if consecutive points indicate wrap */
PFB GMT_this_point_wraps;			/* Used in above */
PFV GMT_get_crossings;				/* Returns map crossings in x or y */
PFI GMT_truncate;				/* Truncate polygons agains boundaries */
BOOLEAN GMT_meridian_straight = FALSE;		/* TRUE if meridians plot as straight lines */
BOOLEAN GMT_parallel_straight = FALSE;		/* TRUE if parallels plot as straight lines */
int GMT_3D_mode = 3;				/* Determines if we draw fore and/or back 3-D box lines [Default is both] */
char *GMT_plot_format[3][2];			/* Keeps the 6 formats for dd:mm:ss plot output */
GMT_LONG GMT_n_lon_nodes = 360;			/* Somewhat arbitrary # of nodes for lines in longitude (may be reset in gmt_map.c) */
GMT_LONG GMT_n_lat_nodes = 180;			/* Somewhat arbitrary # of nodes for lines in latitude (may be reset in gmt_map.c) */
double GMT_dlon = 0.0;				/* Steps taken in longitude along gridlines (gets reset in gmt_init.c) */
double GMT_dlat = 0.0;				/* Steps taken in latitude along gridlines (gets reset in gmt_init.c) */
/*--------------------------------------------------------------------*/
/*	For color lookup purposes */
/*--------------------------------------------------------------------*/

struct GMT_LUT *GMT_lut;		/* CPT lookup table read by GMT_read_cpt */
struct GMT_BFN_COLOR GMT_bfn[3];	/* Structures with back/fore/nan colors */
int GMT_n_colors = 0;			/* Number of colors in CPT lookup table */
int GMT_cpt_flags = 0;			/* Flags controling use of BFN colors */
BOOLEAN GMT_gray;			/* TRUE if only grayshades are needed */
BOOLEAN GMT_b_and_w;			/* TRUE if only black and white are needed */
BOOLEAN GMT_continuous;			/* TRUE if continuous color tables have been given */
BOOLEAN GMT_cpt_pattern = FALSE;	/* TRUE if cpt file contains any patterns */
BOOLEAN GMT_cpt_skip = FALSE;		/* TRUE if current z-slice is to be skipped */
#ifdef GMT_CPT2	
BOOLEAN GMT_categorical = FALSE;	/* TRUE if CPT applies to categorical data */
#endif
/*--------------------------------------------------------------------*/
/*	For projection purposes */
/*--------------------------------------------------------------------*/

struct GMT_MAP_PROJECTIONS project_info;
struct GMT_THREE_D z_project;
struct GMT_DATUM_CONV GMT_datum;	/* For datum conversions */
PFI GMT_forward, GMT_inverse;		/* Pointers to the selected mapping functions */
PFI GMT_x_forward, GMT_x_inverse;	/* Pointers to the selected linear functions */
PFI GMT_y_forward, GMT_y_inverse;	/* Pointers to the selected linear functions */
PFI GMT_z_forward, GMT_z_inverse;	/* Pointers to the selected linear functions */

/*--------------------------------------------------------------------*/
/*	For i/o purposes */
/*--------------------------------------------------------------------*/

FILE *GMT_stdin;			/* Pointer for standard input */
FILE *GMT_stdout;			/* Pointer for standard output */
PFI GMT_input;				/* Pointer to function reading ascii or binary tables */
PFI GMT_input_ascii;			/* Pointer to function reading ascii tables only */
PFI GMT_output;				/* Pointer to function writing ascii or binary tables */
PFI GMT_output_ascii;			/* Pointer to function writing ascii tables only */
struct GMT_IO GMT_io = {
	{ FALSE, FALSE },
	{ FALSE, FALSE },
	{ FALSE, FALSE },
	{ FALSE, FALSE },
	{ FALSE, FALSE },
	{ FALSE, FALSE },
	FALSE,
	FALSE,
	0,
	{ 0, 0 },
	-1,
	-1,
	0,
	0,
	0,
	0,
	{ '>', '>'},
	"",
	"",
	"",
	"r",
	"w",
	"a+",
	(BOOLEAN *)NULL,
	(int *)NULL,
	(int *)NULL,
	0,
	0,
	0,
	0,
	(int *)NULL,
	(double *)NULL,
	(double *)NULL,
	(double *)NULL,
	{ {-1, -1, -1, -1}, {-1, -1, -1, -1}, FALSE, FALSE, "", FALSE, FALSE, FALSE, FALSE, { "", ""} },
	{ {-1, -1, -1, -1}, {-1, -1, -1, -1}, FALSE, FALSE, "", FALSE, FALSE, FALSE, FALSE, { "", ""} },
	{ {-1, -1, -1}, 0, 0.0, FALSE, FALSE, { "am", "pm" }, "", { "", ""} },
	{ {-1, -1, -1}, 0, 0.0, FALSE, FALSE, { "am", "pm" }, "", { "", ""} },
	{ {-1, -1, -1}, 0, FALSE, FALSE, FALSE, 0, 0.0, "", "", { "", ""} }
};
struct GMT_Y2K_FIX GMT_Y2K_fix;		/* Used to convert 2-digit years to 4-digit years */
int GMT_n_file_suffix;
PFI GMT_io_readinfo[GMT_N_GRD_FORMATS];
PFI GMT_io_updateinfo[GMT_N_GRD_FORMATS];
PFI GMT_io_writeinfo[GMT_N_GRD_FORMATS];
PFI GMT_io_readgrd[GMT_N_GRD_FORMATS];
PFI GMT_io_writegrd[GMT_N_GRD_FORMATS];
double *GMT_file_scale, *GMT_file_offset, *GMT_file_nan;
int *GMT_file_id;
char **GMT_file_suffix;
double GMT_data[BUFSIZ];
int GMT_pad[4] = {0, 0, 0, 0};
int GMT_inc_code[2] = {0, 0};

/*--------------------------------------------------------------------*/
/*	For misc purposes */
/*--------------------------------------------------------------------*/


char *GMT_SHAREDIR = CNULL;
char *GMT_HOMEDIR = CNULL;
char *GMT_USERDIR = CNULL;
char *GMT_DATADIR = CNULL;
char *GMT_GRIDDIR = CNULL;
char *GMT_IMGDIR = CNULL;
char *GMT_TMPDIR = CNULL;
char *GMT_program;	/* Name of current GMT program */

float GMT_f_NaN;
double GMT_d_NaN;

/* GMT_u2u is the conversion matrix for cm, inch, m, pt: */

double GMT_u2u[4][4] = {
	{   1.00,    1.0/2.54,    0.01,         72.0/2.54 },
	{   2.54,    1.0,         0.0254,       72.0 },
	{ 100.00,    100.0/2.54,  1.0,          72.0/0.0254 },
	{ 2.54/72.0, 1.0/72.0,    0.0254/72.0,  1.0 }
};
char *GMT_unit_names[4] = {"cm", "inch", "m", "point"};
int GMT_no_rgb[3] = {-1, -1, -1};
int GMT_oldargc;
char *GMT_oldargv[GMT_N_UNIQUE];	/* Pointers to old common arguments */

BOOLEAN GMT_give_synopsis_and_exit = FALSE;

struct GMT_TIME_LANGUAGE GMT_time_language;	/* For time axis */

/* For custom symbol plotting in psxy[z] */

int GMT_n_custom_symbols = 0;
struct GMT_CUSTOM_SYMBOL **GMT_custom_symbol;
