//	Author: Eric Mixon		Date:7/15/10
//	Description:
//	- Takes in Pings and converts them
//		to be closer to 0
//	- Must be instantiated using the
//		ping to be normalized around
//	- The ping that is used on instantiation
//		will be in position 0, 0, 0
//
//

//#ifndef VINCENTY_CPP
//#define VINCENTY_CPP

#include"vincenty.h"
#include"geom.h"
#include<math.h>
#include<iostream>

using namespace std;

//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	//#include "../WarningStates.h"	//Disable all Warnings!!!
	#pragma warning ( disable : 4700 )	//Deprecated call
#endif

// Must be instantiated with ping to be normalized around
Vincenty::Vincenty(const Ping& p) : cPing(p) 
{
	a = 6378137;
	b = 6356752.314245;
	f = (a - b) / a;
	pi = 3.14159265;
}

// Returns a Ping normalized as the distance between
// P1 and P2 along an oblateSpheroid
Point Vincenty::distance(const Ping& p) const
{
	if(p == cPing)
	{
		Point p1;;
		return p1;
	}
	double lat1 = cPing.getLat() * (pi / 180);
	double lat2 = p.getLat() * (pi / 180);
	double U1 = atan((1-f) * tan(lat1));
	double U2 = atan((1-f) * tan(lat2));
	double L = (p.getLong() - cPing.getLong()) * (pi / 180);
	double a1, a2, s, sinSigma, cosSigma, sigma,
			sinAlpha, cos2SigmaM, C, u2, A, B, deltaSigma,
			cos2Alpha;

	double lambda = L;
	int iterLimit = 100;
	double lambdaP;

	double sLam = sin(lambda);
	double cLam = cos(lambda);
	double sU1 = sin(U1);
	double cU1 = cos(U1);
	double sU2 = sin(U2);
	double cU2 = cos(U2);
	double accuracy = pow(10.0, -12);

	while(--iterLimit > 0 && fabs(lambda - lambdaP) > accuracy)
	{
		sLam = sin(lambda);
		cLam = cos(lambda);
		sinSigma = sqrt((cU2 * sLam) * (cU2 * sLam) + 
		(cU1 * sU2 - sU1 * cU2 * cLam) * (cU1 * sU2 - sU1 * cU2 * cLam));
		if(sinSigma == 0)
		{
			Point p1;
			return p1;
		}
		cosSigma = sU1 * sU2 + cU1 * cU2 * cLam;
		sigma = atan2(sinSigma, cosSigma);
		sinAlpha = cU1 * cU2 * sLam / sinSigma;
		cos2Alpha = 1 - sinAlpha * sinAlpha;
		cos2SigmaM = cosSigma - 2 * sU1 * sU2 / cos2Alpha;
		if(!(iterLimit > cos2SigmaM) && !(iterLimit < cos2SigmaM))
			cos2SigmaM = 0;
		C = f/16 * cos2Alpha * (4 + f * ( 4 - 3 * cos2Alpha));
		lambdaP = lambda;
		lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * 
		(cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
	}
	if (iterLimit == 0)
	{
		Point p1;
		return p1;
	}
	u2 = cos2Alpha * (a * a - b * b) / (b * b);
	A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
	B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
	deltaSigma = B * sinSigma * (cos2SigmaM + B/4 * (cosSigma * 
		(-1 + 2 * cos2SigmaM * cos2SigmaM) - B/6 * cos2SigmaM * 
		( -3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
	sLam = sin(lambda);
	cLam = cos(lambda);
	
	s = b * A * (sigma - deltaSigma);
	a1 = atan(cU2 * sLam / (cU1 * sU2 - sU1 * cU2 * cLam));
	a2 = atan(cU1 * sLam / (-1 * sU1 * cU2 + cU1 * sU2 * cLam));

	double x, y, z;
	y = s * cos(a1);
	x = -1 * s * sin(a1);
	z = p.getDepth();
	Point newPoint(fabs(x), fabs(y), z);
	return newPoint;
}

// Returns the Ping that is dist distance in az direction away
// from the given ping p
Ping Vincenty::destination(const Ping& p, const double& s, const double& azD) const
{
	double lat1 = p.getLat() * (pi / 180);
	double az1 = azD * (pi / 180);
	double U1 = atan((1-f) * tan(lat1));
	double tU1 = (1-f) * tan(lat1);
	double cU1 = 1 / sqrt((1 + tU1 * tU1));
	double sU1 = tU1 * cU1;
	double sigma1 = atan2(tU1, cos(az1));
	double sAz = cU1 * sin(az1);
	double c2Az = 1 - sAz * sAz;
	double u2 = c2Az * (a * a - b * b) / (b * b);
	double A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
	double B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
	double sigma = s / (b * A);
	double preSig = sigma;
	double accuracy = pow(10.0, -12);
	double dSig, sig2M;
	do
	{
		sig2M = 2 * sigma1 + sigma;
		dSig = B * sin(sigma) * ( cos(sig2M) + B/4 * (cos(sigma) * (-1 + 2 * cos(sig2M) * cos(sig2M)) - 
		B/6 * cos(sig2M) * (-3 + 4 * sin(sigma) * sin(sigma)) * (-3 + 4 * cos(sig2M) * cos(sig2M))));
		preSig = sigma;
		sigma = s / (b * A) + dSig;
	}while (fabs(preSig - sigma) > accuracy);
	double tmp = sU1 * sin(sigma) - cU1 * cos(sigma) * cos(az1);
	double lat2 = atan2(sU1 * cos(sigma) + cU1 * sin(sigma) * cos(az1), 
					(1-f) * sqrt(sAz * sAz + tmp * tmp));
	lat2 = lat2 * (180 / pi);
	double lamda = atan2(sin(sigma) * sin(az1), cU1 * cos(sigma) - sU1 * sin(sigma) * cos(az1));
	double C = f/16 * c2Az * (4 + f * (4 - 3 * c2Az));
	double L = lamda - (1 - C) * f * sAz * (sigma + C * sin(sigma) * (cos(sig2M) + C * cos(sigma) * 
		(-1 + 2 * cos(sig2M) * cos(sig2M))));
	double az2 = atan2(sAz, -tmp);
	double lon2 = p.getLong() + (L * 180 / pi);
	Ping p1(lon2, lat2, p.getDepth());
	return p1;
}
		
// Returns the Ping that is dist distance in az direction away
// from cPing
Ping Vincenty::destination(const double& s, const double& azD) const
{
	double lat1 = cPing.getLat() * (pi / 180);
	double az1 = azD * (pi / 180);
	double U1 = atan((1-f) * tan(lat1));
	double tU1 = (1-f) * tan(lat1);
	double cU1 = 1 / sqrt((1 + tU1 * tU1));
	double sU1 = tU1 * cU1;
	double sigma1 = atan2(tU1, cos(az1));
	double sAz = cU1 * sin(az1);
	double c2Az = 1 - sAz * sAz;
	double u2 = c2Az * (a * a - b * b) / (b * b);
	double A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
	double B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
	double sigma = s / (b * A);
	double preSig = sigma;
	double accuracy = pow(10.0, -12);
	double dSig, sig2M;
	do
	{
		sig2M = 2 * sigma1 + sigma;
		dSig = B * sin(sigma) * ( cos(sig2M) + B/4 * (cos(sigma) * (-1 + 2 * cos(sig2M) * cos(sig2M)) - 
		B/6 * cos(sig2M) * (-3 + 4 * sin(sigma) * sin(sigma)) * (-3 + 4 * cos(sig2M) * cos(sig2M))));
		preSig = sigma;
		sigma = s / (b * A) + dSig;
	}while (fabs(preSig - sigma) > accuracy);
	double tmp = sU1 * sin(sigma) - cU1 * cos(sigma) * cos(az1);
	double lat2 = atan2(sU1 * cos(sigma) + cU1 * sin(sigma) * cos(az1), 
				(1-f) * sqrt(sAz * sAz + tmp * tmp));
	lat2 = lat2 * (180 / pi);
	double lamda = atan2(sin(sigma) * sin(az1), cU1 * cos(sigma) - sU1 * sin(sigma) * cos(az1));
	double C = f/16 * c2Az * (4 + f * (4 - 3 * c2Az));
	double L = lamda - (1 - C) * f * sAz * (sigma + C * sin(sigma) * (cos(sig2M) + C * cos(sigma) * 
		(-1 + 2 * cos(sig2M) * cos(sig2M))));
	double az2 = atan2(sAz, -tmp);
	double lon2 = cPing.getLong() + (L * 180 / pi);
	Ping p1(lon2, lat2, cPing.getDepth());
	return p1;
}
		
double Vincenty::distance(const Ping& p1, const Ping& p2) const
{
	double lat1 = p1.getLat() * (pi / 180);
	double lat2 = p2.getLat() * (pi / 180);
	double U1 = atan((1-f) * tan(lat1));
	double U2 = atan((1-f) * tan(lat2));
	double L = (p2.getLong() - p1.getLong()) * (pi / 180);
	double /*a1, a2,*/ s, sinSigma, cosSigma, sigma,
			sinAlpha, cos2SigmaM, C, u2, A, B, deltaSigma,
			cos2Alpha;
	double lambda = L;
	int iterLimit = 100;
	double lambdaP;
	double sLam = sin(lambda);
	double cLam = cos(lambda);
	double sU1 = sin(U1);
	double cU1 = cos(U1);
	double sU2 = sin(U2);
	double cU2 = cos(U2);
	double accuracy = pow(10.0, -12);
	while(--iterLimit > 0 && fabs(lambda - lambdaP) > accuracy)
	{
		sLam = sin(lambda);
		cLam = cos(lambda);
		sinSigma = sqrt((cU2 * sLam) * (cU2 * sLam) + 
		(cU1 * sU2 - sU1 * cU2 * cLam) * (cU1 * sU2 - sU1 * cU2 * cLam));
		if(sinSigma == 0)
			return 0;
		cosSigma = sU1 * sU2 + cU1 * cU2 * cLam;
		sigma = atan2(sinSigma, cosSigma);
		sinAlpha = cU1 * cU2 * sLam / sinSigma;
		cos2Alpha = 1 - sinAlpha * sinAlpha;
		cos2SigmaM = cosSigma - 2 * sU1 * sU2 / cos2Alpha;
		if(!(iterLimit > cos2SigmaM) && !(iterLimit < cos2SigmaM))
			cos2SigmaM = 0;
		C = f/16 * cos2Alpha * (4 + f * ( 4 - 3 * cos2Alpha));
		lambdaP = lambda;
		lambda = L + (1 - C) * f * sinAlpha * (sigma + C * sinSigma * 
			(cos2SigmaM + C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));
	}

	if (iterLimit == 0)
		return -1;
	
	u2 = cos2Alpha * (a * a - b * b) / (b * b);
	A = 1 + u2 / 16384 * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)));
	B = u2 / 1024 * (256 + u2 * (-128 + u2 * (74 - 47 * u2)));
	deltaSigma = B * sinSigma * (cos2SigmaM + B/4 * (cosSigma * 
		(-1 + 2 * cos2SigmaM * cos2SigmaM) - B/6 * cos2SigmaM * 
		( -3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
	sLam = sin(lambda);
	cLam = cos(lambda);
			
	s = b * A * (sigma - deltaSigma);
	return s;
}

Point Vincenty::pointOnSphere(const Ping& p) const
{
	double lon = p.getLong() * (pi / 180);
	double lat = p.getLat() * (pi / 180);
	double Ex, Ey;
	if(p.getLat() == 90)
	{
		Ex = 0;
		Ey = b;
	}
	else
	{
		Ex = a * cos(lat);
		Ey = b * sin(lat);
	}
	
	double r = sqrt((Ex * Ex) + (Ey * Ey)) - p.getDepth();
	double x, y, z;
			
	y = r * sin(lat);
	double s;
	if(p.getLat() == 90)
		s = 0;
	else
		s = r * cos(lat);
	if(p.getLong() == 90)
		x = 0;
	else
		x = s * cos(lon);
	if(p.getLong() == 180)
		z = 0;
	else
		z = s * sin(lon);
	
	Point pN(z, y, x);

	return pN;
}

Point Vincenty::normalizePing(const Ping& p) const {return distance(p);}

Ping Vincenty::normalizePoint(const Point& p) const
{
	Ping p1;
	p1 = destination(p.x, 90);
	p1 = destination(p1, p.y, 0);
	p1.setDepth(p.z);
	return p1;
}

void Vincenty::setcPing(const Ping& p) { cPing = p; }

//#endif
//Restore warning state -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( pop )
#endif 
