/** 
* @file			mb_zgrid.h
* @brief		Act as the front end to the MB_ZGrid routine to make it easier to use.
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int mb_zgrid(float *z, int *nx, int *ny, 
		float *x1, float *y1, float *dx, float *dy, float *xyz, 
		int *n, float *zpij, int *knxt, int *imnew, 
		float *cay, int *nrng);


#ifdef __cplusplus
}
#endif






