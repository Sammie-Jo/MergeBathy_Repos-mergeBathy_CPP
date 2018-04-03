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
#include "MB_Threads.h"
#include "scalecInterp.h"

// SJZ Added 1/22/15
//#include <tchar.h>
#include <errno.h>
#pragma comment(lib,"User32.lib")
#define BUF_SIZE 255
//************************************************************************************
// 0. Declare local functions that will call the scalecInterp* routines
//************************************************************************************
#ifdef WIN32
#include <strsafe.h>
void ErrorHandler(LPTSTR lpszFunction);
DWORD WINAPI threadInterp( LPVOID lpParam );
DWORD WINAPI threadInterp2( LPVOID lpParam );
DWORD WINAPI threadInterp6( LPVOID lpParam );
DWORD WINAPI threadInterpKrig( LPVOID lpParam );
DWORD WINAPI threadInterpTile( LPVOID lpParam );
DWORD WINAPI threadInterpTileKrig( LPVOID lpParam );
#else
void *threadInterp( void *lpParam );
void *threadInterp2( void *lpParam );
void *threadInterp6( void *lpParam );
void *threadInterpKrig( void *lpParam );
void *threadInterpTile( void *lpParam );
void *threadInterpTileKrig( void *lpParam );
void *SuspendThread();
void *ResumeThread();
timespec time_ns = {0, 50*1000*1000};
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  cond = PTHREAD_COND_INITIALIZER;
int play = 0;
#endif

//************************************************************************************
// I. Set up local variables that can be used for tracking
//************************************************************************************
map<int,int> locationThreadMap;
int numTotalThreads;

//************************************************************************************
// II. Constructor for mbThreads
//************************************************************************************
mbThreads::mbThreads(int numThreads)
{
	numTotalThreads = numThreads;
#ifdef WIN32
	dwThreadIdArray = vector<DWORD>(numThreads);
	hThreadArray = vector<HANDLE>(numThreads);
#else
	dwThreadIdArray = vector<int>(numThreads);
	hThreadArray = vector<pthread_t>(numThreads);
#endif
}

//************************************************************************************
// IIIa. Assign pointers to interpolation data to local data structures. Used for regular data.
//************************************************************************************
void mbThreads::makeMBThread(SCALEC_TILE_DATA_POINTER stdp)
{
	pDataArray = vector<SCALEC_TILE_DATA>(numTotalThreads);
	pDataArray2 = vector<SCALEC_DATA>();
	for (int i = 0; i < numTotalThreads; i++)
	{
		pDataArray[i] = *stdp;
	}
}

//************************************************************************************
// IIIb. Assign pointers to interpolation data to local data structures. Used for irregular data.
//************************************************************************************
void mbThreads::makeMBThread(SCALEC_DATA_POINTER sdp)
{
	pDataArray = vector<SCALEC_TILE_DATA>();
	pDataArray2 = vector<SCALEC_DATA>(numTotalThreads);
	for (int i = 0; i < numTotalThreads; i++)
	{
		pDataArray2[i] = *sdp;
	}
}

//************************************************************************************
// IVa. Create and spawn the independent threads that will process regular data
//************************************************************************************
void mbThreads::initMBThread_Tile(int runType)
{
#ifdef WIN32
	if (runType == 0)
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			hThreadArray[i] = CreateThread(
				NULL,					// default security attributes
				0,						// use default stack size
				threadInterpTile,	    // thread function name
				&pDataArray[i],			// argument to thread function
				CREATE_SUSPENDED,		// use default creation flags
				&dwThreadIdArray[i]);	// returns the thread identifier
			locationThreadMap[dwThreadIdArray[i]] = i;
			ResumeThread(hThreadArray[i]);

			// Check the return value for success.
			// If CreateThread fails, terminate execution. 
			// This will automatically clean up threads and memory. 
			if (hThreadArray[i] == NULL) 
			{
				ErrorHandler(TEXT("CreateThread"));
				ExitProcess(3);
			}
		}
	}else
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			hThreadArray[i] = CreateThread(
				NULL,				   // default security attributes
				0,					   // use default stack size
				threadInterpTileKrig,  // thread function name
				&pDataArray[i],		   // argument to thread function
				CREATE_SUSPENDED,	   // use default creation flags
				&dwThreadIdArray[i]);  // returns the thread identifier
			locationThreadMap[dwThreadIdArray[i]] = i;
			ResumeThread(hThreadArray[i]);

			// Check the return value for success.
			// If CreateThread fails, terminate execution. 
			// This will automatically clean up threads and memory. 
			if (hThreadArray[i] == NULL) 
			{
			   ErrorHandler(TEXT("CreateThread"));
			   ExitProcess(3);
			}
		}
	}
#else
	if (runType == 0)
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			SuspendThread();
			//fprintf(stderr,"i %d\n",i);
			//cout<<"i "<<i<<endl;

			dwThreadIdArray[i] = pthread_create( &hThreadArray[i], NULL, threadInterpTile, (void *) &pDataArray[i]);
			locationThreadMap[hThreadArray[i]] = i;
			
			//fprintf(stderr,"ID %d\n",hThreadArray[i]);
			//cout<<"ID "<<pthread_self()<<endl;
			ResumeThread();
			// Sleep to allow thread to run before suspending again.  This
			// is not necessary but added to see suspension and resume.
			//nanosleep(&time_ns, NULL); 

			if(dwThreadIdArray[i])
			{
				fprintf(stderr,"Error - pthread_create() return code: %d\n %s",dwThreadIdArray[i], strerror(dwThreadIdArray[i]));
				exit(EXIT_FAILURE);
			}
		}
	}else
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			SuspendThread();
			dwThreadIdArray[i] = pthread_create( &hThreadArray[i], NULL, threadInterpTileKrig, (void *) &pDataArray[i]);
			locationThreadMap[hThreadArray[i]] = i;
			ResumeThread();
			//nanosleep(&time_ns, NULL); 
			if(dwThreadIdArray[i])
			{
				fprintf(stderr,"Error - pthread_create() return code: %d\n %s",dwThreadIdArray[i], strerror(dwThreadIdArray[i]));
				exit(EXIT_FAILURE);
			}
		}
	}
#endif
}

//************************************************************************************
// IVb. Create and spawn the independent threads that will process irregular data
//************************************************************************************
void mbThreads::initMBThread(int runType)
{
#ifdef WIN32
	if (runType == 0)
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			hThreadArray[i] = CreateThread(
				NULL,				    // default security attributes
				0,					    // use default stack size
				threadInterp,			// thread function name
				&pDataArray2[i],		// argument to thread function
				CREATE_SUSPENDED,	    // use default creation flags
				&dwThreadIdArray[i]);   // returns the thread identifier
			locationThreadMap[dwThreadIdArray[i]] = i;
			ResumeThread(hThreadArray[i]);

			// Check the return value for success.
			// If CreateThread fails, terminate execution. 
			// This will automatically clean up threads and memory. 
			if (hThreadArray[i] == NULL) 
			{
			   ErrorHandler(TEXT("CreateThread"));
			   ExitProcess(3);
			}
		}
	}
	//This function go with the separate function scalecInterp_ProcessKrig.
	//Note that the function scalecInterp_ProcessKrig is broken and untested.
	/*else
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			hThreadArray[i] = CreateThread(
				NULL,				   // default security attributes
				0,					   // use default stack size
				threadInterpKrig,	   // thread function name
				&pDataArray2[i],	   // argument to thread function
				CREATE_SUSPENDED,	   // use default creation flags
				&dwThreadIdArray[i]);  // returns the thread identifier
			locationThreadMap[dwThreadIdArray[i]] = i;
			ResumeThread(hThreadArray[i]);

			// Check the return value for success.
			// If CreateThread fails, terminate execution. 
			// This will automatically clean up threads and memory. 
			if (hThreadArray[i] == NULL) 
			{
			   ErrorHandler(TEXT("CreateThread"));
			   ExitProcess(3);
			}
		}
	}*/
#else
	if (runType == 0)
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			SuspendThread();
			dwThreadIdArray[i] = pthread_create( &hThreadArray[i], NULL, threadInterp, (void *) &pDataArray2[i]);
			locationThreadMap[hThreadArray[i]] = i;
			ResumeThread();
			//nanosleep(&time_ns, NULL); 
			if(dwThreadIdArray[i])
			{
				fprintf(stderr,"Error - pthread_create() return code: %d\n %s",dwThreadIdArray[i], strerror(dwThreadIdArray[i]));
				exit(EXIT_FAILURE);
			}
		}
	}else
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			SuspendThread();
			dwThreadIdArray[i] = pthread_create( &hThreadArray[i], NULL, threadInterpKrig, (void *) &pDataArray2[i]);
			locationThreadMap[hThreadArray[i]] = i;
			ResumeThread();
			//nanosleep(&time_ns, NULL); 
			if(dwThreadIdArray[i])
			{
				fprintf(stderr,"Error - pthread_create() return code: %d\n %s",dwThreadIdArray[i], strerror(dwThreadIdArray[i]));
				exit(EXIT_FAILURE);
			}
		}
	}

#endif
}

void mbThreads::initMBThread2(int runType)
{
#ifdef WIN32
	if (runType == 0)
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			hThreadArray[i] = CreateThread(
				NULL,				    // default security attributes
				0,					    // use default stack size
				threadInterp2,			// thread function name
				&pDataArray2[i],		// argument to thread function
				CREATE_SUSPENDED,	    // use default creation flags
				&dwThreadIdArray[i]);   // returns the thread identifier
			locationThreadMap[dwThreadIdArray[i]] = i;
			ResumeThread(hThreadArray[i]);

			// Check the return value for success.
			// If CreateThread fails, terminate execution. 
			// This will automatically clean up threads and memory. 
			if (hThreadArray[i] == NULL) 
			{
			   ErrorHandler(TEXT("CreateThread"));
			   ExitProcess(3);
			}
		}
	}
#else
	if (runType == 0)
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			SuspendThread();
			dwThreadIdArray[i] = pthread_create( &hThreadArray[i], NULL, threadInterp2, (void *) &pDataArray2[i]);
			locationThreadMap[hThreadArray[i]] = i;
			ResumeThread();
			//nanosleep(&time_ns, NULL); 
			if(dwThreadIdArray[i])
			{
				fprintf(stderr,"Error - pthread_create() return code: %d\n %s",dwThreadIdArray[i], strerror(dwThreadIdArray[i]));
				exit(EXIT_FAILURE);
			}
		}
	}
#endif
}
void mbThreads::initMBThread6(int runType)
{
#ifdef WIN32
	if (runType == 0)
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			hThreadArray[i] = CreateThread(
				NULL,				    // default security attributes
				0,					    // use default stack size
				threadInterp6,			// thread function name
				&pDataArray2[i],		// argument to thread function
				CREATE_SUSPENDED,	    // use default creation flags
				&dwThreadIdArray[i]);   // returns the thread identifier
			locationThreadMap[dwThreadIdArray[i]] = i;
			ResumeThread(hThreadArray[i]);

			// Check the return value for success.
			// If CreateThread fails, terminate execution. 
			// This will automatically clean up threads and memory. 
			if (hThreadArray[i] == NULL) 
			{
			   ErrorHandler(TEXT("CreateThread"));
			   ExitProcess(3);
			}
		}
	}
#else
	if (runType == 0)
	{
		for (int i = 0; i < numTotalThreads; i++)
		{
			SuspendThread();
			dwThreadIdArray[i] = pthread_create( &hThreadArray[i], NULL, threadInterp6, (void *) &pDataArray2[i]);
			locationThreadMap[hThreadArray[i]] = i;
			ResumeThread();
			//nanosleep(&time_ns, NULL); 
			if(dwThreadIdArray[i])
			{
				fprintf(stderr,"Error - pthread_create() return code: %d\n %s",dwThreadIdArray[i], strerror(dwThreadIdArray[i]));
				exit(EXIT_FAILURE);
			}
		}
	}
#endif
}
//************************************************************************************
// V. Bring all threads back to a single original thread when they are each done with their computation
//************************************************************************************
void mbThreads::joinMBThread()
{
#ifdef WIN32
	for (int i = 0; i < numTotalThreads; i++)
	{
		WaitForSingleObject(hThreadArray[i], INFINITE);
	}
#else
	void *res;
	for (int i = 0; i < numTotalThreads; i++)
	{
		int s = pthread_join(hThreadArray[i], &res); // NULL):
		if (s == 0)
			free(res);
		else cerr << "Error: pthread_join" << endl;
	}
#endif
}

//************************************************************************************
// VI. Kill handles to created threads and clear up variables
//************************************************************************************
void mbThreads::terminateMBThread()
{
#ifdef WIN32
	for(int i = 0; i < numTotalThreads; i++)
	{
		CloseHandle(hThreadArray[i]);
	}
#else
	//The way pthreads work we don't actually need to close the handle it gets closed when it gets joined back to the main thread
#endif
	dwThreadIdArray.clear();
	hThreadArray.clear();
	pDataArray.clear();
	pDataArray2.clear();
	locationThreadMap.clear();
}

#ifdef WIN32

//************************************************************************************
// VII. Interpolate irregular data
//************************************************************************************
DWORD WINAPI threadInterp(LPVOID lpParam)
{
	SCALEC_DATA_POINTER pda;

	// Cast the parameter to the correct data type.
	// The pointer is known to be valid because
	// it was checked for NULL before the thread was created.

	pda = (SCALEC_DATA_POINTER)lpParam;
	int flag = scalecInterp_Process(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);

	return 0;
}
DWORD WINAPI threadInterp2(LPVOID lpParam)
{
	SCALEC_DATA_POINTER pda;

	// Cast the parameter to the correct data type.
	// The pointer is known to be valid because
	// it was checked for NULL before the thread was created.

	pda = (SCALEC_DATA_POINTER)lpParam;
	//int flag = scalecInterp_Process2A(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);//handles w and w/o kriging
	int flag = scalecInterp_Process4A(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);//handles w and w/o kriging
//	int flag = scalecInterp_Process2(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);//part 2 w/o kriging

	return 0;
}
DWORD WINAPI threadInterp6(LPVOID lpParam)
{
	SCALEC_DATA_POINTER pda;

	// Cast the parameter to the correct data type.
	// The pointer is known to be valid because
	// it was checked for NULL before the thread was created.

	pda = (SCALEC_DATA_POINTER)lpParam;
	int flag = scalecInterp_Process6A(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);//handles w and w/o kriging

	return 0;
}
//************************************************************************************
// VIII. Interpolate irregular data
//************************************************************************************
DWORD WINAPI threadInterpKrig(LPVOID lpParam)
{
	SCALEC_DATA_POINTER pda;

	// Cast the parameter to the correct data type.
	// The pointer is known to be valid because
	// it was checked for NULL before the thread was created.

	pda = (SCALEC_DATA_POINTER)lpParam;
	int flag = scalecInterp_ProcessKrig(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);

	return 0;
}

//************************************************************************************
// IX. Interpolate regular data
//************************************************************************************
DWORD WINAPI threadInterpTile(LPVOID lpParam)
{
	SCALEC_TILE_DATA_POINTER pda;

	// Cast the parameter to the correct data type.
	// The pointer is known to be valid because
	// it was checked for NULL before the thread was created.
	
	pda = (SCALEC_TILE_DATA_POINTER)lpParam;
	int flag = scalecInterpTile_ProcessA(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);
//	int flag = scalecInterpTile_Process(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);

	return 0;
}

//************************************************************************************p
// X. Interpolate regular data with kriging
//************************************************************************************
DWORD WINAPI threadInterpTileKrig(LPVOID lpParam)
{
	SCALEC_TILE_DATA_POINTER pda;

	// Cast the parameter to the correct data type.
	// The pointer is known to be valid because
	// it was checked for NULL before the thread was created.

	//Calling original function from mergeBathy v3.6 instead of the _Serial version for consistency.
	pda = (SCALEC_TILE_DATA_POINTER)lpParam;
	int flag = scalecInterpTile_ProcessKrig(pda, locationThreadMap[(int)GetCurrentThreadId()], numTotalThreads);

	return 0;
}

#else

void *threadInterp(void *lpParam)
{
	for(;;) 
	{
		pthread_mutex_lock(&lock);
		while(!play) { /* We're paused */
			pthread_cond_wait(&cond, &lock); /* Wait for play signal */
		}
		pthread_mutex_unlock(&lock);
		/* Continue */

		SCALEC_DATA_POINTER pda;

		// Cast the parameter to the correct data type.
		// The pointer is known to be valid because
		// it was checked for NULL before the thread was created.

		pda = (SCALEC_DATA_POINTER)lpParam;
		int flag = scalecInterp_Process(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);

		return 0;
	}
}

void *threadInterp2(void *lpParam)
{
	for(;;) 
	{
		pthread_mutex_lock(&lock);
		while(!play) { /* We're paused */
			pthread_cond_wait(&cond, &lock); /* Wait for play signal */
		}
		pthread_mutex_unlock(&lock);
		/* Continue */

		SCALEC_DATA_POINTER pda;

		// Cast the parameter to the correct data type.
		// The pointer is known to be valid because
		// it was checked for NULL before the thread was created.

		pda = (SCALEC_DATA_POINTER)lpParam;
		int flag = scalecInterp_Process4A(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);//handles w and w/o kriging
//		int flag = scalecInterp_Process2A(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);//handles w and w/o kriging
//		int flag = scalecInterp_Process2(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);

		return 0;
	}
}
void *threadInterp6(void *lpParam)
{
	for(;;) 
	{
		pthread_mutex_lock(&lock);
		while(!play) { /* We're paused */
			pthread_cond_wait(&cond, &lock); /* Wait for play signal */
		}
		pthread_mutex_unlock(&lock);
		/* Continue */

		SCALEC_DATA_POINTER pda;

		// Cast the parameter to the correct data type.
		// The pointer is known to be valid because
		// it was checked for NULL before the thread was created.

		pda = (SCALEC_DATA_POINTER)lpParam;
		int flag = scalecInterp_Process6A(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);//handles w and w/o kriging
//		int flag = scalecInterp_Process2A(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);//handles w and w/o kriging
//		int flag = scalecInterp_Process2(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);

		return 0;
	}
}
void *threadInterpKrig(void *lpParam)
{
	for(;;) 
	{
		pthread_mutex_lock(&lock);
		while(!play) { /* We're paused */
			pthread_cond_wait(&cond, &lock); /* Wait for play signal */
		}
		pthread_mutex_unlock(&lock);
		/* Continue */

		SCALEC_DATA_POINTER pda;

		// Cast the parameter to the correct data type.
		// The pointer is known to be valid because
		// it was checked for NULL before the thread was created.

		pda = (SCALEC_DATA_POINTER)lpParam;
		int flag = scalecInterp_ProcessKrig(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);

		return 0;
	}
}

void *threadInterpTile(void *lpParam)
{
	for(;;) 
	{
		pthread_mutex_lock(&lock);
		while(!play) { /* We're paused */
			pthread_cond_wait(&cond, &lock); /* Wait for play signal */
		}
		pthread_mutex_unlock(&lock);
		/* Continue */

		//cout<<"thread#"<<locationThreadMap[(unsigned int)pthread_self()]<<endl;
		//cout<<"threadID"<<(unsigned int)pthread_self()<<endl;
		//fprintf(stderr,"thread# %d\n",locationThreadMap[(unsigned int)pthread_self()]);
		//fprintf(stderr,"threadID %d\n",pthread_self());

		SCALEC_TILE_DATA_POINTER pda;

		// Cast the parameter to the correct data type.
		// The pointer is known to be valid because
		// it was checked for NULL before the thread was created.

		pda = (SCALEC_TILE_DATA_POINTER)lpParam;
	//	int flag = scalecInterpTile_Process(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);
		int flag = scalecInterpTile_ProcessA(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);

		return 0;
	}
}

void *threadInterpTileKrig(void *lpParam)
{
	for(;;) 
	{
		pthread_mutex_lock(&lock);
		while(!play) { /* We're paused */
			pthread_cond_wait(&cond, &lock); /* Wait for play signal */
		}
		pthread_mutex_unlock(&lock);
		/* Continue */

		SCALEC_TILE_DATA_POINTER pda;

		// Cast the parameter to the correct data type.
		// The pointer is known to be valid because
		// it was checked for NULL before the thread was created.

		pda = (SCALEC_TILE_DATA_POINTER)lpParam;
		int flag = scalecInterpTile_ProcessKrig(pda, locationThreadMap[(unsigned int)pthread_self()], numTotalThreads);

		return 0;
	}
}

#endif

//SJZ added to catch CreateThread Errors. 1/22/15
#ifdef WIN32
void ErrorHandler(LPTSTR lpszFunction) 
{ 
	// Retrieve the system error message for the last-error code.
	HANDLE hStdout;
	size_t cchStringSize;

	// Make sure there is a console to receive output results. 

	hStdout = GetStdHandle(STD_OUTPUT_HANDLE);
	if( hStdout == INVALID_HANDLE_VALUE )
		return;// 1;

	LPVOID lpMsgBuf;
	LPVOID lpDisplayBuf;
	DWORD dw = GetLastError(); 

	FormatMessage(
		FORMAT_MESSAGE_ALLOCATE_BUFFER | 
		FORMAT_MESSAGE_FROM_SYSTEM |
		FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		dw,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(LPTSTR) &lpMsgBuf,
		0, NULL );

	// Display the error message.

	lpDisplayBuf = (LPVOID)LocalAlloc(LMEM_ZEROINIT, (lstrlen((LPCTSTR) lpMsgBuf) + lstrlen((LPCTSTR) lpszFunction) + 40) * sizeof(TCHAR)); 
	StringCchPrintf((LPTSTR)lpDisplayBuf, LocalSize(lpDisplayBuf) / sizeof(TCHAR), TEXT("%s failed with error %d: %s"), lpszFunction, dw, lpMsgBuf); 
	
	// Print the parameter values using thread-safe functions.
		
	StringCchLength((LPCTSTR) lpDisplayBuf, BUF_SIZE, &cchStringSize);	
	WriteConsole(hStdout, (LPCTSTR) lpDisplayBuf, (DWORD)cchStringSize, &dw, NULL);
	
	// Free error-handling buffer allocations.

	LocalFree(lpMsgBuf);
	LocalFree(lpDisplayBuf);
}
#else
void *SuspendThread()
{
	pthread_mutex_lock(&lock);
	play = 0;
	pthread_mutex_unlock(&lock);
}
void *ResumeThread()
{
	pthread_mutex_lock(&lock);
	play = 1;
	pthread_cond_signal(&cond);
	pthread_mutex_unlock(&lock);
}
#endif
