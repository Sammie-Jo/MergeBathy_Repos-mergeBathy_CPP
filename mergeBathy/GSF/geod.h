/* geod.h
       author:  Robert J. Lacey, NIMA (laceyb@nima.mil)

        Reference:
        Proceedings of the 7th International Symposium on Geodetic Computations,
        1985
        "The Nested Coefficients Method for Accurate Solutions of Direct
        and Inverse Geodetic Problems with Any Length"
        Zhang Xue-Lian
        p747-763.
*/


#ifndef _GEOD_H
#define	_GEOD_H

#ifdef  __cplusplus 
extern "C"{
#endif

extern short geo_direct
	(double a, double rf, double lat1, double lon1, double az1, double s, 
	    		double *lat2, double *lon2,  double *az2 );

extern short geo_inverse
	(double a, double rf, double lat1, double lon1, double lat2,
       			double lon2, double *az1, double *az2, double *s );

#ifdef  __cplusplus 
}
#endif

#endif	  /* _GEOD_H */
