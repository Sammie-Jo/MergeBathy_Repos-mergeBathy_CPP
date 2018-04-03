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
//************************************************************************
// Created by Samantha Zambo	2014/08/28 
// WarningStates.h contains all warnings to disable for
// third-party files and libraries.
//************************************************************************


//Disable warnings since this is a Third-party file. -SJZ
#pragma warning ( disable : 4018 )	//Signed/unsigned mismatch
#pragma warning ( disable : 4101 )	//Unreferenced local variable
#pragma warning ( disable : 4244 )	//Type conversion; possible loss of data
#pragma warning ( disable : 4267 )  //Conversion from size_t; possible loss of data
#pragma warning ( disable : 4273 )  //Inconsistent dll linkage
#pragma warning ( disable : 4305 )	//Type truncation
#pragma warning ( disable : 4700 )	//Uninitialized local variable
#pragma warning ( disable : 4996 )	//Deprecated call



//#if defined __GNUC__
//#pragma GCC system_header
//#elif defined __SUNPRO_CC
//#pragma disable_warn
//#elif defined _MSC_VER
//#pragma warning(push, 1)
//#endif



//restore warning state
//#if defined __SUNPRO_CC
//#pragma enable_warn
//#elif defined _MSC_VER
//#pragma warning(pop)
//#endif