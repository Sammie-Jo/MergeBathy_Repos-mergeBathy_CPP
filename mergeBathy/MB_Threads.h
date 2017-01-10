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


