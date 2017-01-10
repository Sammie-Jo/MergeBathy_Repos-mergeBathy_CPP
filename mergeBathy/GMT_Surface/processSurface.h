/** 
* @file			processSurface.h
* @brief		Act as the front end to the GMT Surface routines to make them easier to use.
* @author		Kevin Duvieilh
* @date			06 June 2011
*
*/


#ifndef PROCESSSURFACE_H_
#define PROCESSSURFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

/**
* This function calls the GMT Surface routine.
* @param xData - 1 dimensional double array of initial known X coordinates.
* @param yData - 1 dimensional double array of initial known Y coordinates.
* @param zData - 1 dimensional double array of initial known depth values.
* @param inputDataSize - Size of the xData, yData, and zData vectors.
* @param x0 - Minimum value of xData.
* @param y0 - Minimum value of yData.
* @param z0 - Minimum value of zData.
* @param x1 - Maximum value of xData.
* @param y1 - Maximum value of yData.
* @param z1 - Maximum value of zData.
* @param spacingX - Computational grid spacing in the X direction.
* @param spacingY - Computational grid spacing in the Y direction.
* @param tension - Computational grid tension factor.
* @param xPostSurface - 1 dimensional double array of interpolated X coordinates. (Returned).
* @param yPostSurface - 1 dimensional double array of interpolated Y coordinates. (Returned).
* @param zPostSurface - 1 dimensional double array of interpolated depth values. (Returned).
* @param postSurfaceSize - Size of the xPostSurface, yPostSurface, and zPostSurface vectors.
* @return Success or failure vaule.
*/
int processSurface(double *xData, double *yData, double *zData, int inputDataSize, double x0, double y0, double z0, double x1, double y1, double z1, double spacingX, double spacingY, double tension, double **xPostSurface, double **yPostSurface, double **zPostSurface, int *postSurfaceSize);

#ifdef __cplusplus
}
#endif

#endif
