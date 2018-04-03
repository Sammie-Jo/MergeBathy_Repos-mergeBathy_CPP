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
/**
* @file			MB_Threads.h
* @brief		Thread class used to multi-thred mergeBathy.
* @author		Kevin Duvieilh
* @date			04 August 2011
*
*/
#pragma once
#include <stdlib.h>
#include <map>
#include <vector>
#include "outFileStructs.h"

#ifdef WIN32 
#include <Windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

/**
* A multi-threading class used in mergeBathy.
*/
class mbThreads
{

public:

	/**
	* A constructor for mbThread.
	* Initializes the arrays and sets the values number of total threads.
	* @param numThreads - Argument that sets the number of total threads.
	*/
	mbThreads(int numThreads);

	/**
	* Initializes the local data structures.
	* @param stdp - A pointer to the data structure to be passed to interp tile.
	*/
	void makeMBThread(SCALEC_TILE_DATA_POINTER stdp);

	/**
	* Initializes the local data structures.
	* @param sdp - A pointer to the data structure to be passed to interp tile.
	*/
	void makeMBThread(SCALEC_DATA_POINTER sdp);

	/**
	* Initializes the thread process.
	* This function branches of each thread for its own separate processing.  Used for scalecInterpTile.
	* @param runType - Run type to determine if kriging is to occur.
	*/
	void initMBThread_Tile(int runType);

	/**
	* Initializes the thread process.
	* This function branches of each thread for its own separate processing.  Used for scalecInterp.
	* @param runType - Run type to determine if kriging is to occur.
	*/
	void initMBThread(int runType);
	void initMBThread2(int runType);
	void initMBThread6(int runType);

	/**
	* Joins threads back to the main function.
	* Will wait for all threads to complete then return control back to main.
	*/
	void joinMBThread();

	/**
	* Kills al threads after completion.
	* Cleans up all variables and closes handles.
	*/
	void terminateMBThread();

private:

	/**
	* Data structure containing the arguments passed to the scalecInterpTile and scalecInterpTileKrig routines.
	*/
	std::vector<SCALEC_TILE_DATA> pDataArray;

	/**
	* Data structure containing the arguments passed to the scalecInterp routines.
	*/
	std::vector<SCALEC_DATA> pDataArray2;

#ifdef WIN32
	/**
	* Data structure containing the thread IDs. (Windows).
	*/
	 std::vector<DWORD> dwThreadIdArray;

	/**
	* Data structure containing the thread handles. (Windows).
	*/
	 std::vector<HANDLE> hThreadArray;
#else
	/**
	* Data structure containing the thread IDs. (Linux).
	*/
	 std::vector<int> dwThreadIdArray;

	/**
	* Data structure containing the thread handles. (Linux).
	*/
	 std::vector<pthread_t> hThreadArray;
#endif
};


