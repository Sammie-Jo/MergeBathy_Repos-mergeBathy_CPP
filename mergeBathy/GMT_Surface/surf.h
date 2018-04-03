
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
#ifndef __SURF_H__
#define __SURF_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef NVWIN3X
    #ifdef __MINGW32__
        #include <sys/types.h>
        #include <sys/stat.h>
        #include <unistd.h>
        #include <process.h>
    #endif
#endif
#include <string.h>
#include <errno.h>
#include "nvtypes.h"
#include "nvdef.h"


#ifdef  __cplusplus
extern "C" {
#endif


#define NAVO_SURFNULL      10000000000000000.0
#define NAVO_BLOCKMEAN     1
#define NAVO_SURFACE       2

//Experimenting with median for blockmedian
//NV_FLOAT64 median(int n, NV_FLOAT64 x[]);

  NV_INT32 surf_init (NV_F64_XYMBR mbr, NV_FLOAT64 x_interval, NV_FLOAT64 y_interval, NV_FLOAT32 minz, NV_FLOAT32 maxz,
		      NV_FLOAT64 tension_value);

  NV_INT32 surf_load (NV_F64_COORD3 xyz);
  NV_INT32 surf_load_XYZ (NV_FLOAT64 x, NV_FLOAT64 y, NV_FLOAT64 z);
  NV_INT32 surf_load_XYz (NV_FLOAT64 x, NV_FLOAT64 y, NV_FLOAT32 z);
  NV_INT32 surf_load_xyz (NV_FLOAT32 x, NV_FLOAT32 y, NV_FLOAT32 z);
  NV_INT32 surf_load_NV_F64_COORD2_Z (NV_F64_COORD2 xy, NV_FLOAT64 z);
  NV_INT32 surf_load_NV_F64_COORD2_z (NV_F64_COORD2 xy, NV_FLOAT32 z);
  NV_INT32 surf_load_NV_F32_COORD2_z (NV_F32_COORD2 xy, NV_FLOAT32 z);
  NV_INT32 surf_load_NV_F32_COORD3 (NV_F32_COORD3 xyz);

  NV_INT32 surf_proc (NV_BOOL surface);

  NV_INT32 surf_rtrv (NV_FLOAT64 **z_array, NV_INT32 **cnt_array, NV_INT32 *final_rows, NV_INT32 *final_cols);

  NV_INT32 surf_cleanup ();

  NV_INT32 get_surf_oob ();

  void get_surf_error (NV_CHAR *error_string, NV_INT32 *system_errno);

  NV_BOOL surf_progress_callback_registered ();

  void surf_progress (NV_INT32 state, NV_INT32 percent);

  typedef void (*SURF_PROGRESS_CALLBACK) (NV_INT32 state, NV_INT32 percent);

  void surf_register_progress_callback (SURF_PROGRESS_CALLBACK progressCB);

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
			    NV_INT32 *nrows,
			    NV_INT32 *ncols);					

#ifdef  __cplusplus
}
#endif

#endif
