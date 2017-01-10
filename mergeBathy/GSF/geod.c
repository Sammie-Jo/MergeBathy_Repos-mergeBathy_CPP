/* geod.c

	author:  Robert J. Lacey, NIMA (laceyb@nima.mil)

	Notes (by David H. Fabre, 5 Sep 97):

	This software is part of the MUSE 2.0 collection of software.
	It was taken from the MDTCC (MUSE Datum Transformation and
		Coordinate Conversion) Library's source file
		"dstazfns.c" (distance-azimuth functions).
	At this time, these are administrated by NIMA's employee
	Bob Lacey.  I converted NIMA's set of ellipsoids to the 
	companion header "ellipsoid.h".  The reference for the
	method of solving the geodetic forward and inverse problem is
	listed below.  These routines return 0 if successful.

	Reference:
	Proceedings of the 7th International Symposium on Geodetic Computations,
	1985
	"The Nested Coefficients Method for Accurate Solutions of Direct
	and Inverse Geodetic Problems with Any Length"
	Zhang Xue-Lian
	p747-763.
*/

#include <stdio.h>
#include <math.h>
#ifndef PI
#	define PI	(3.14159265358979323846) 
#endif

static double M0( double e2 )
	{	double e4 = e2*e2;
		return PI*(1.0 - e2*( 1.0/4.0 + e2*( 3.0/64.0 + e2*(5.0/256.0) )))/2.0;
	}
/* s == distance */
short geo_direct
	(	double a, double rf, double lat1, double lon1, double az1, double s, 
	    		double *lat2, double *lon2,  double *az2 )
{ 	double RADDEG = (PI)/180.0, testv = 1.0E-10;
	double f = ( rf > 0.0 ? 1.0/rf : 0.0 );
	double b = a*(1.0-f), e2 = f*(2.0-f);
	double phi1 = lat1*RADDEG, lam1 = lon1*RADDEG;
	double sinphi1 = sin(phi1), cosphi1 = cos(phi1);
	double azm1 = az1*RADDEG;
	double sinaz1 = sin(azm1), cosaz1 = cos(azm1);

/* dhf printf("a = %lf, b = %lf\n", a, b); */
	
	if( fabs(s) < 0.01 )  /* distance < centimeter => congruency */
	{	*lat2 = lat1;
		*lon2 = lon1;
		*az2 = 180.0 + az1;
		if( *az2 > 360.0 ) *az2 -= 360.0;
		return 0;
	}
	else
    if( cosphi1 )	 /* non-polar origin */
    {	/* u1 is reduced latitude */
    	double tanu1 = sqrt(1.0-e2)*sinphi1/cosphi1;
   		double sig1 = atan2(tanu1,cosaz1);
		double cosu1 = 1.0/sqrt( 1.0 + tanu1*tanu1 ), sinu1 = tanu1*cosu1;
		double sinaz =  cosu1*sinaz1, cos2saz = 1.0-sinaz*sinaz;
		double us = cos2saz*e2/(1.0-e2);
		/*	Terms */
		double	ta = 1.0+us*(4096.0+us*(-768.0+us*(320.0-175.0*us)))/16384.0,
				tb = us*(256.0+us*(-128.0+us*(74.0-47.0*us)))/1024.0,
				tc = 0;
		/*	FIRST ESTIMATE OF SIGMA (SIG) */		
		double first = s/(b*ta); /* !!*/
        double sig = first;
		double c2sigm, sinsig,cossig, temp,denom,rnumer, dlams, dlam;
		do
		{	c2sigm = cos(2.0*sig1+sig);
			sinsig = sin(sig); cossig = cos(sig);
	        temp = sig;
	        sig = first + 
	        		tb*sinsig*(c2sigm+tb*(cossig*(-1.0+2.0*pow(c2sigm,2.0)) - 
	        		tb*c2sigm*(-3.0+4.0*pow(sinsig,2.0))*(-3.0+4.0*pow(c2sigm,2.0))/6.0)/4.0);
        }
        while( fabs(sig-temp) > testv);
		/* 	LATITUDE OF POINT 2 */
      	/*	DENOMINATOR IN 2 PARTS (TEMP ALSO USED LATER) */
        temp = sinu1*sinsig-cosu1*cossig*cosaz1;
        denom = (1.0-f)*sqrt(sinaz*sinaz+temp*temp);
		/* NUMERATOR */
        rnumer = sinu1*cossig+cosu1*sinsig*cosaz1;
		*lat2 = atan2(rnumer,denom)/RADDEG;
		/* DIFFERENCE IN LONGITUDE ON AUXILARY SPHERE (DLAMS ) */
        rnumer = sinsig*sinaz1;
        denom = cosu1*cossig-sinu1*sinsig*cosaz1;
        dlams = atan2(rnumer,denom);
		/* TERM C */
        tc = f*cos2saz*(4.0+f*(4.0-3.0*cos2saz))/16.0;
		/* DIFFERENCE IN LONGITUDE */
        dlam = dlams-(1.0-tc)*f*sinaz*(sig+tc*sinsig*(c2sigm+tc*cossig*(-1.0+2.0*
          pow(c2sigm,2.0))));
		*lon2 = (lam1+dlam)/RADDEG;
        if(*lon2 > 180.0  ) *lon2 -= 360.0;
        if(*lon2 < -180.0 ) *lon2 += 360.0;
		/* AZIMUTH - FROM NORTH */
        *az2 = atan2(-sinaz,temp)/RADDEG;
		if( fabs(*az2) < testv ) *az2 = 0.0;
        if( *az2 < 0.0) *az2 += 360.0;
		return 0;
    }
    else /* phi1 == 90 degrees, polar origin  */
    {	double dM = a*M0(e2) - s;
		double paz = ( phi1 < 0.0 ? 180.0 : 0.0 );
    	return geo_direct( a,rf, 0.0, lon1, paz, dM,lat2,lon2,az2 );
    } 
}

short geo_inverse(	double a, double rf, double lat1, double lon1, double lat2,
       			double lon2, double *az1, double *az2, double *s )
{	short iter=0;
	/* dhf double RADDEG = (PI)/180.0, testv = 1.0E-10; */
	double RADDEG = (PI)/180.0, testv = 1.0E-11;
	double f = ( rf > 0.0 ? 1.0/rf : 0.0 );
	double b = a*(1.0-f), e2 = f*(2.0-f);
	double phi1 = lat1*RADDEG, lam1 = lon1*RADDEG;
	double sinphi1 = sin(phi1), cosphi1 = cos(phi1);
	double phi2 = lat2*RADDEG, lam2 = lon2*RADDEG;
	double sinphi2 = sin(phi2), cosphi2 = cos(phi2);
	
    if( (fabs(lat1-lat2) < testv && 
    	( fabs(lon1-lon2) < testv) || fabs(lat1-90.0) < testv ) )
    {	/* TWO STATIONS ARE IDENTICAL : SET DISTANCE & AZIMUTHS TO ZERO */
        *az1 = 0.0; *az2 = 0.0; *s = 0.0;
        return 0;
    }
	else
	if(  fabs(cosphi1) < testv ) /* initial point is polar */
	{	short k = geo_inverse( a,rf, lat2,lon2,lat1,lon1, az1,az2,s );
		b = *az1; *az1 = *az2; *az2 = b;
		return 0;
	}
	else
	if( fabs(cosphi2) < testv ) /* terminal point is polar */
	{	short k = geo_inverse( a,rf, lat1,lon1,lat1,lon1+180.0, az1,az2,s );
		*s /= 2.0;
		*az2 = *az1 + 180.0;
		if( *az2 > 360.0 ) *az2 -= 360.0; 
		return 0;
	}
	else  	/* Geodesic passes through the pole (antipodal) */
    if( (fabs( fabs(lon1-lon2) - 180 ) < testv) && (fabs(lat1+lat2) < testv) ) 
    {	double s1,s2;
    	geo_inverse( a,rf, lat1,lon1, lat1,lon2, az1,az2, &s1 );
    	geo_inverse( a,rf, lat2,lon2, lat1,lon2, az1,az2, &s2 );
    	*az2 = *az1;
    	*s = s1 + s2;
    	return 0;
    }
	else  /* antipodal and polar points don't get here */
	{	double	dlam = lam2 - lam1,
    			dlams = dlam;
		double sdlams,cdlams, sig,sinsig,cossig, sinaz,cos2saz, c2sigm;
		double tc,temp, us,rnumer,denom, ta,tb;
		double cosu1,sinu1, sinu2,cosu2;
		/* Reduced latitudes */
		temp = (1.0-f)*sinphi1/cosphi1;
		cosu1 = 1.0/sqrt(1.0+temp*temp);
		sinu1 = temp*cosu1;
		temp = (1.0-f)*sinphi2/cosphi2;
		cosu2 = 1.0/sqrt(1.0+temp*temp);
		sinu2 = temp*cosu2;
    
    	do
		{	sdlams = sin(dlams), cdlams = cos(dlams);
			sinsig = sqrt(pow(cosu2*sdlams,2.0)+pow(cosu1*sinu2-sinu1*cosu2*cdlams,2.0));
        	cossig = sinu1*sinu2+cosu1*cosu2*cdlams;
			
			sig = atan2(sinsig,cossig);
			sinaz = cosu1*cosu2*sdlams/sinsig;
			cos2saz = 1.0-sinaz*sinaz;
        	c2sigm = (sinu1 == 0.0 || sinu2 == 0.0 ? cossig : 
        					cossig-2.0*sinu1*sinu2/cos2saz);
			tc = f*cos2saz*(4.0+f*(4.0-3.0*cos2saz))/16.0;
			temp = dlams;
			dlams = dlam+(1.0-tc)*f*sinaz*(sig+tc*sinsig*(c2sigm+tc*cossig*(-1.0+2.0*
          				pow(c2sigm,2.0))));
        	if (fabs(dlams) > PI && iter++ > 50) 
        		return iter;
     	}
        while ( fabs(temp-dlams) > testv);
   		us = cos2saz*(pow(a,2.0)-pow(b,2.0))/pow(b,2.0);	 /* !! */
		/* BACK AZIMUTH FROM NORTH */
		rnumer = -(cosu1*sdlams);
    	denom = sinu1*cosu2-cosu1*sinu2*cdlams;
    	*az2 = atan2(rnumer,denom)/RADDEG;
		if( fabs(*az2) < testv ) *az2 = 0.0;
    	if(*az2 < 0.0) *az2 += 360.0;
		/* FORWARD AZIMUTH FROM NORTH */
    	rnumer = cosu2*sdlams;
    	denom = cosu1*sinu2-sinu1*cosu2*cdlams;
    	*az1 = atan2(rnumer,denom)/RADDEG;
		if( fabs(*az1) < testv ) *az1 = 0.0;
    	if(*az1 < 0.0) *az1 += 360.0;
		/* Terms a & b */
		ta = 1.0+us*(4096.0+us*(-768.0+us*(320.0-175.0*us)))/16384.0;
		tb = us*(256.0+us*(-128.0+us*(74.0-47.0*us)))/1024.0;
		/* GEODETIC DISTANCE */
    	*s = b*ta*(sig-tb*sinsig*(c2sigm+tb*(cossig*(-1.0+2.0*pow(c2sigm,2.0))-tb*
      		c2sigm*(-3.0+4.0*pow(sinsig,2.0))*(-3.0+4.0*pow(c2sigm,2.0))/6.0)/4.0));
		return 0;
	}
}

#ifdef TEST_GEOD

/* c89 -o geod geod.c -DTEST_GEOD -lm */

#include <stdio.h>
#include "geod.h"
#include "ellipsoid.h"

/******************************************************************************/
int main(int argc, char **argv)
{
	double lat, lon, lat1, lon1, azf, azb, s;
	ELLIPSOID choice=WE; /* WGS 84 */
	int i, which=0;

for(i = 0; i < 28; i++)
printf("<%-30.30s> %s=%2.2d %lf %.9lf\n",
NAMEP(i), ABBRVP(i), i, AXIS(i), RFLAT(i));
fflush(stdin);

if (argc > 1) 
{
sscanf(argv[1], "%d", &choice);
printf("argv[1] = %d\n", choice);
}

printf("using:  ");
printf("<%s> %s=%d %lf %.9lf\n",
NAMEP(choice), ABBRVP(choice), choice, AXIS(choice), RFLAT(choice));

printf("enter 0 (direct) or 1 (inverse): ");
scanf("%d", &which);

	
if (which == 0)
{
	while (1)
	{
		printf("enter lat, lon, azf, s: ");
		scanf("%lf %lf %lf %lf", &lat, &lon, &azf, &s);
		printf("stat = %d, ", 
			geo_direct(AXIS(choice), RFLAT(choice), 
				lat, lon, azf, s, &lat1, &lon1, &azb) );
		printf("lat1, lon1, azb:  %.11lf %.11lf %.11lf\n",
			lat1, lon1, azb);
	}
} /* direct */
else
{
	while (1)
	{
		printf("enter lat, lon, lat1, lon1: ");
		scanf("%lf %lf %lf %lf", &lat, &lon, &lat1, &lon1);
		printf("stat = %d, ", 
			geo_inverse(AXIS(choice), RFLAT(choice), 
				lat, lon, lat1, lon1, &azf, &azb, &s) );
		printf("azf, azb, s:  %.11lf %.11lf %.11lf\n",
			azf, azb, s);
	}
} /* inverse */

	return 0;
} /* main */
#endif
