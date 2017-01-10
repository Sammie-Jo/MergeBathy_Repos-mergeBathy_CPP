/* geodesic.c
	purpose:  use traditional analytic answer to the geodeisc on 
			wgs84.
*/

#include <math.h>
#include "geodesic.h"

static double azNtoS(double azFromNorth);
static int direct(double phi, double alam, double fazi, double s,
		double *phipri, double *alampr);

#define direct_angle(omg, sn, cs, c4, css, c3, c2)	( (omg) + (sn)*(cs)*(1.0e-6)*( (c4) + (css)*((c3) + (c2)*(css)) ) )

/* WGS-84 semi-major & semi-minor axes, reciprocal of flattening, resp. */
#define         WGS84_A0           6378137.0L
#define         WGS84_B0           6356752.314245L

/******************************************************************************/
static double azNtoS(double azFromNorth)
/* azFromNorth must be in [0, 360)
*/
{
	double azFromSouth;

	azFromSouth = azFromNorth + 180.0;

	if (azFromSouth > 360.0) azFromSouth -= 360.0;
	
	return azFromSouth;

} /* azNtoS */

/******************************************************************************/
static int direct(double phi, double alam, double fazi, double s,
		double *phipri, double *alampr)
/*
c
c     ..................................................................
c
c	computes the geodetic position (latitude,longitude) and   
c	azimuth of an observed station from a station of known      
c	geodetic position, with azimuth and distance to the observed
c	station given. 
c
c	the spheroid of reference is the clarke spheroid of 1866.
c
c	evaluation is based on equations for the forward position   
c	computation developed by the u.s. coast and geodetic survey.
c	this method is valid for distances up to 600 miles.         
c
c	reference:
c	simmons,lansing g., natural tables for the calculation of   
c	geodetic positions, u.s. coast and geodetic survey special  
c	publication no.241, 1949.  
c
c	description of parameters
c	phi    - geodetic latitude of known station, 
c                in seconds of arc (input)
c	alam   - geodetic longitude of known station,
c                in seconds of arc (input)
c	fazi   - geodetic forward azimuth of known station,
c                in seconds of arc (input)
c	s      - geodetic distance, in meters (input)
c	phipri - geodetic latitude of observed station,  
c                in seconds of arc (output)
c	alampr - geodetic longitude of observed station,  
c                in seconds of arc (output)
c
c	authors:
c	original author   - rosemary e. riordan /noaa/nos/c54x2.
c	subroutine author - milton stein /noaa/nos/c542. 
c
c	modified by yann king ,12/30/86 to include datum option.
c
c	modified 10/21/87 by tom stepka.
c
c	moved to vax/vms and fortran 77.
c	there are now integers or do loops to cause f77 conversion problems.
c
c	put all horizontal datum constants into common block hdatum.
c	call routine hdinit to load these constants with the proper values.
c
c	the original verison of this routine contained a discrepency in
c	its definitions of the constants z2 and z4.  in gpxy and marc these
c	are called w1 and w3, respectively, and are defined having the
c	same magnitude but have the opposite sign.  therefore, in order to
c	keep all common block definitions consistent for this version, z2 and
c	z4 will change signs in this verison of direct.  fortunately they
c	are only used once, as arguments to the inline funcion angle, so
c	to keep direct running the same, i changed the signs of z2 and z4
c	in this function call.
c
c	modified by tim rulon 22-mar-89
c
c	corrected underflow problem when faz very close to multiple
c	of 90 degrees and s is small. did not cause a problem under
c	ordinary circumstances. caused floating point trap.
c
c     ..................................................................
c
	Converted to C by David Fabre, Planning Systems Inc., 15 Mar 94
	for NRL for NAVO

        Had to change the macro "angle" to "direct_angle" due to conflict with 
        Qt 4 libraries.  Jan Depner, 10/19/07
*/
{
	/* wgs 84 / nad 83 parameters taken from hdinit routine */
	static double axis = 6378137.0,
		esq = 6.69438002292318e-3,
		z2 = 0.1104825480468667,
		z3 = 21.25880271589918,
		z4 = 5048.250737106752,
		z1 = 6367449.145771138,
		w2 = 0.23832988623869,
		w3 = 29.36926254762117,
		w4 = 5022.894084128594,
		w1 = 0.1570487611454482e-6;
	/* local constant */
	static double arc1 = 0.484813681110e-5;
	/* local variables */
	double fazj, x, y, b, phj, sina, ef, xcor, xpri, a, ycor, ypri, cosa,
		cssq, y0, y1, omega, sin1, cos1, css1, phi1, faca, v, va, y2,
		facb, facc, cay, y3, omegb, sin2, cos2, css2, h, applam,
		asinco, dellam;

/****	determined these unnecessary DHF
	double delphi, phimid, applap;
****/

	/*computation of geodetic position of observed station */

	fazj = fazi*arc1;
	x = s * sin(fazj);
	if (fabs(x) < 1.0e-3) x = 0.001;
	y = -s * cos(fazj);
	if (fabs(y) < 1.0e-3) y = 0.001;
	b = pow(y/10000.0, 2.0);
	phj = phi*arc1;
	sina = sin(phj);
	ef = (pow(1.0 - esq*sina*sina, 2.0)*1.0e15)/(3.0*axis*axis*(1.0 - esq));
	xcor = b*ef/2.0;
	xpri = x - xcor*x*1.0e-7;  
	a = pow(xpri/10000.0, 2.0);  
	ycor = ef*a;
	ypri = y + ycor*y*1.0e-7;  
	cosa = cos(phj);
	cssq = cosa*cosa;
	y0 = z1*direct_angle(phj, sina, cosa, -z4, cssq, z3, -z2);
	y1 = y0 + ypri;
	omega = w1*y1;
	sin1 = sin(omega);
	cos1 = cos(omega);
	css1 = cos1*cos1;
	phi1 = direct_angle(omega, sin1, cos1, w4, css1, w3, w2);
	sin1 = sin(phi1);
	cos1 = cos(phi1);
	faca = (sqrt(1.0 - esq*sin1*sin1)) / (2.0*axis);
	v = sin1/cos1*faca*1.0e8;
	va = v*a;
	y2 = y1 - va;
	facb = (1.0 + 3.0*( (sin1*sin1)/(cos1*cos1) )) / (3.0*(sin1/cos1));
	facc = (3.0*esq*sin1*cos1) / (1.0 - esq);
	cay = (faca*(facb - facc))*1.0e6;
	y3 = y2 + cay*pow(va/1000.0, 2.0);
	omegb = w1*y3;
	sin2 = sin(omegb);
	cos2 = cos(omegb);
	css2 = cos2*cos2;
	*phipri = direct_angle(omegb, sin2, cos2, w4, css2, w3, w2);

/****	delphi is related only to the appalp and phimid assignments below DHF
	delphi = phj - *phipri;
****/
	h = sqrt(1.0 - esq*pow(sin(*phipri), 2.0)) / (axis*cos(*phipri)*arc1);
	applam = h*xpri;
	asinco = (v*va)/15.0;
	dellam = applam + applam*asinco*1.0e-7;
	*alampr = alam - dellam;

/****	the 2 assignments below are made for no apparent reason so
	I took them out, perhaps they were output at one time but
	not necessary here DHF
	appalp = dellam * (sin(phj) + sin(*phipri)) / (1.0 + cos(delphi));
	phimid = (phj + *phipri)/2.0e0;
****/

	*phipri /= arc1;

	if (fabs(*alampr) > 648000.0)
	{
		if (*alampr > 0.0)

			*alampr -= 1296000.0;

		else if (*alampr < 0.0)

			*alampr += 1296000.0;
	}

  	return 1;

} /* direct */

/******************************************************************************/
void newgp(double latobs, double lonobs, double az, double dist,
		double *lat, double *lon)
{

/*	notes:
	1. This function was written to replace newgp of "when dinosaurs
		roamed".  It is simply a toplevel runner of the direct routine
		that was converted to C from NOAA's FORTRAN code.
	2. latobs, lonobs are the coordinates of the observed position (degrees)
	3. az is the azimuth from North and must be in [0, 360)
	4. dist is in meters
	5. *lat, *lon will contain the coordinate values of the directed
		position
	6. function direct works in arc seconds, azimuth from south, and meters
	7. function azNtoS takes the azimuth from North and gives you azimuth
		from South
	8. the ellipsoid constants for newgp are hard-coded to be those of
		wgs 84
*/
	double phipri, alampr;
	int stat;

	/* the conditional below was put in because direct bombed out when
		on the equator going east or west
	*/
	if ( latobs == 0.0 && (az == 90.0 || az == 270.0) )
		latobs += 1.0e-37;

	stat = direct(3600.0*latobs, 3600.0*lonobs, 3600.0*azNtoS(az), dist,
		&phipri, &alampr);

	*lat = phipri/3600.0;
	*lon = alampr/3600.0;

} /* newgp */

/***************************************************************************\
*                                                                           *
*   Module Name:        invgp                                               *
*                                                                           *
*   Programmer:         Unknown                                             *
*                                                                           *
*   Date Written:       When dinosaurs roamed the earth.                    *
*                                                                           *
*   Modified:           Jan C. Depner, converted to C and changed the       *
*                       arguments to degrees instead of radians.            *
*                                                                           *
*   Date:               September, 1992                                     *
*                                                                           *
*   Module Security                                                         *
*   Classification:     Unclassified                                        *
*                                                                           *
*   Data Security                                                           *
*   Classification:     Unknown                                             *
*                                                                           *
*   Purpose:            Given the semi-major axis, semi-minor axis, and     *
*                       two geographic positions this routine will compute  *
*                       the distance between the two gp's and the azimuth   *
*                       (clockwise from north) from the first gp to the     *
*                       second gp.                                          *
*                                                                           *
*   Inputs:             a0                  -   semi-major axis in meters   *
*                       b0                  -   semi-minor axis in meters   *
*                       rlat1               -   start latitude in degrees   *
*                       rlon1               -   start longitude in degrees  *
*                       rlat2               -   end latitude in degrees     *
*                       rlon2               -   end longitude in degrees    *
*                       dist                -   distance in meters          *
*                       az                  -   azimuth in degrees          *
*                                                                           *
*   Outputs:            none                                                *
*                                                                           *
\***************************************************************************/

void invgp (double rlat1, double rlon1, double rlat2, double rlon2,
		double *dist, double *az)
{
    double          drlat1, drlat2, drlon1, drlon2, dell, beta1, sbeta1,
                        cbeta1, beta2, sbeta2, cbeta2, adell, sidel, codel, a,
                        b, siphi, cophi, q1, q2, c, em, phi, phisq, csphi,
                        ctphi, psyco, term1, term2, term3, term4, term5, term6,
                        xlam1, tan;
    int            n;
    static int     first = 1;
    static double   pi = 3.141592653589793, twopi = 6.283185307179586,
                        rad_to_deg = 57.2957795147195,
                        tiny = .0000000000000000000000000000001, flat, flat2, f1,
                        f2, f3, f4, f5, f6, f7, f8;

    if (first)
    {
        flat = 1.0l - (WGS84_B0 / WGS84_A0);
        flat2 = flat * flat;
        f1 = flat2 * 1.25;
        f2 = flat2 * 0.5;
        f3 = flat2 * 0.25;
        f4 = flat2 * 0.125;
        f5 = flat2 * 0.0625;
        f6 = flat + flat2;
        f7 = f6 + 1.0;
        f8 = f6 * 0.5;

        first = 0;
    }

    drlat1 = rlat1 / rad_to_deg;
    drlat2 = rlat2 / rad_to_deg;
    drlon1 = rlon1 / rad_to_deg;
    drlon2 = rlon2 / rad_to_deg;

    beta1 = atan ((1.0 - flat) * sin (drlat1) / cos (drlat1));
    sbeta1 = sin (beta1);
    cbeta1 = cos (beta1);
    beta2 = atan ((1.0 - flat) * sin (drlat2) / cos (drlat2));
    sbeta2 = sin (beta2);
    cbeta2 = cos (beta2);

    dell = drlon1 - drlon2;
    adell = fabs (dell);

    if (drlon1 * drlon2 < 0.0)
    {
        adell = fabs (drlon1) + fabs (drlon2);
        dell = adell;
        if (drlon1 < 0.0) dell = - adell;

        if (adell > pi)
        {
            adell = twopi - adell;
            dell = adell;
            if (drlon1 > 0.0) dell = - adell;
        }
    }

    adell = twopi - adell;
    sidel = sin (adell);
    codel = cos (adell);
    a = sbeta1 * sbeta2;
    b = cbeta1 * cbeta2;
    cophi = a + b * codel;
    q1 = sidel * cbeta2;
    q1 *= q1;
    q2 = sbeta2 * cbeta1 - sbeta1 * cbeta2 * codel;
    q2 *= q2;
    siphi = sqrt (q1 + q2);
    c = b * sidel / siphi;
    em = 1.0 - c * c;

    phi = atan (siphi / (sqrt (1.0 - siphi * siphi) + tiny));

    if (cophi < 0.0) phi = pi - phi;
    phisq = phi * phi;
    csphi = 1.0 / siphi;
    ctphi = cophi / siphi;
    psyco = siphi / cophi;

    /*  Compute distance.                                               */

    term1 = f7 * phi;
    term2 = a * (f6 * siphi - f2 * phisq * csphi);
    term3 = em * (f2 * phisq * ctphi - f8 * (phi + psyco));
    term4 = a * a * f2 * psyco;
    term5 = em * em * (f5 * (phi + psyco) - f2 * phisq * ctphi - f4 * psyco *
        cophi * cophi);
    term6 = a * em * f2 * (phisq * csphi + psyco * cophi);
    *dist =  WGS84_B0 * (term1 + term2 + term3 - term4 + term5 + term6);

    /*  Compute azimuth.                                                */

    term1 = f6 * phi;
    term2 = a * (f2 * siphi + flat2 * phisq * csphi);
    term3 = em * (f3 * psyco + flat2 * phisq * ctphi - f1 * phi);
    xlam1 = c * (term1 - term2 + term3) + adell;
    q1 = sbeta2 * cbeta1 - cos (xlam1) * sbeta1 * cbeta2;
    q2 = sin (xlam1) * cbeta2;
    if (q1 == 0.0) q1 = tiny;
    tan = q2 / q1;
    *az = atan (tan);

    /*  Put azimuth in proper quadrant.                                 */

    n = 3;
    if (dell * tan < 0.0) n = 4;
    if (dell < 0.0) n = n - 2;
    if (q1 > 0.0 && dell == 0.0) n = 1;
    q2 = n;
    q1 = q2 * pi - pi - *az;
    if (n >= 3) q1 = (q2 - 2.0) * pi + *az;
    *az = q1 * rad_to_deg;

} /* invgp */

#ifdef TEST_GEODESIC

/* cc -o geodesic geodesic.c -DTEST_GEODESIC -lm */

#include <stdio.h>
#include "geodesic.h"

/******************************************************************************/
int main(int argc, char **argv)
{
	double lat, lon, lat1, lon1, azf, azb, s;
	int which=0;

printf("enter 0 (direct) or 1 (inverse): ");
scanf("%d", &which);

	
if (which == 0)
{
	while (1)
	{
		printf("enter lat, lon, azf, s: ");
		scanf("%lf %lf %lf %lf", &lat, &lon, &azf, &s);
		newgp( lat, lon, azf, s, &lat1, &lon1 );
		printf("lat1, lon1, azb:  %.11lf %.11lf\n",
			lat1, lon1);
	}
} /* direct */
else
{
	while (1)
	{
		printf("enter lat, lon, lat1, lon1: ");
		scanf("%lf %lf %lf %lf", &lat, &lon, &lat1, &lon1);
		invgp( lat, lon, lat1, lon1, &azf, &s );
		printf("azf, s:  %.11lf %.11lf\n",
			azf, s);
	}
} /* inverse */

	return 0;
} /* main */
#endif
