#ifndef _GEODESIC_H
#define _GEODESIC_H

void newgp(double latobs, double lonobs, double az, double dist, double *lat, double *lon);
void invgp (double rlat1, double rlon1, double rlat2, double rlon2, double *dist, double *az);

#endif

