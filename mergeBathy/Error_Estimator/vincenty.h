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
#pragma once
#include"geom.h"
#include"math.h"

/** Class that implements vincenty formulae procedures to convert Pings to meters and Points in Meters to Pings */
class Vincenty
{
	private:
		/** Ping to subtract from all other pings when normalizing. */
		Ping cPing;	

		/** WGS-84 ellipsoid distance of semimajor axis. */
		double a;

		/** WGS-84 ellipsoid distance of semiminor axis. */
		double b;

		/** WGS-84 ellipsoid frequency. */
		double f;

		/** pi. */
		double pi;
	
	public:
		
		/** Constructor which sets the Ping to be Normalized around (generally bottom left Ping of a dataset)
		* @param p - Ping to be set as cPing.
		*/
		Vincenty(const Ping& p);

		/** Converts Ping p to a relative Point as meter distance from cPing.
		* @param p - Ping to be normalized.
		* @return Point with x and y set as meter distance of p from cPing and z = p.depth();
		*/
		Point normalizePing( const Ping& p) const;
		
		/** Converts Point p to a relative Ping as Longitude and Latitude.
		* @param p - Point to be normalized.
		* @return Ping position longitude and latitude destination when moving p.x / p.y distance from cPing.
		*/
		Ping normalizePoint( const Point& p) const;
		
		/** Sets the value of cPing to p.
		* @param p - New value to be stored to cPing.
		*/
		void setcPing(const Ping& p);

		/** Returns a Point with x and y representing the respective distances of p from cPing. (used in normalizePing)
		* @param p - Ping to be measured from cPing.
		* @return Point representing x and y equal to the distance in meters of p from cPing in the x and y direction
		*/
		Point distance(const Ping& p) const;

		/** Returns the resulting Ping when traveling from p at direction azD for distance s in meters.
		* @param p - Starting Ping.
		* @param s - Distance to travel in meters.
		* @param azD - Azimuth direction to travel.
		* @return - Resulting Ping destination when traveling from p.
		*/
		Ping destination(const Ping& p, const double& s, const double& azD) const;

		/** Returns the resulting Ping when traveling from cPing at direction azD for distance s in meters.
		* @param s - Distance to travel in meters.
		* @param azD - Azimuth direction to travel.
		* @return - Resulting Ping destination when traveling from cPing.
		*/
		Ping destination(const double& s, const double& azD) const;
		
		/** Returns the distance in meters between Pings p1, and p2.
		* @param p1 - First Ping.
		* @param p2 - Second Ping.
		* @return Distance in meters from p1 to p2 as a double.
		*/
		double distance(const Ping& p1, const Ping& p2) const;

		/** Returns absolute Point in 3d space in meters that p is positioned on the Earth with the
		* center of the Earth as the origin.
		* @param p - Ping being measured.
		* @return Point in 3d space that represents p with the Earth's center as the origin.
		*/
		Point pointOnSphere(const Ping& p) const;
};

