:: Author: Samantha J Zambo
:: Date: 3/22/18
:: For OSP Publication 'MergeBathy 2015'
:: ::

@echo off
echo.
echo ***********************************************
echo * setenv.bat
echo * This script will set the usr environment variables to run MergeBathy.
echo * If set globally, environment variables are set permanently.
echo * If set locally, environment variables are set only for the current terminal/cmd window.
echo * 
echo * To work in both x86 and x64, and both Debug and Release,
echo * set environment variables locally for a single terminal/cmd instance.
echo * If set locally, you must run this script before calling MergeBathy for every new terminal/cmd window.
echo * Alternatively, make a script of MergeBathy calls and define the vars in the beginning of the script.
echo ***********************************************
echo.
::=====================================================================
:: Find platform and configuration 
echo Warning: GMT and MBZ pre-splining routines will crash in 32bit for large data sets.
set /p bitModeID="Select (1) 32bit (2) 64bit : " 
set /p configModeID="Select (1) Debug (2) Release : "
set instanceModeID=1
set /p instanceModeID="Select (1) Locally (current instance only) (2) Globally (permanently in usr envars) : "

if [%bitModeID%] EQU [1] ( 
	SET mPLATFORM=x86
	goto:config
)
if [%bitModeID%] EQU [2] ( 
	SET mPlATFORM=x64
	goto:config
) 
echo "Invalid Bit Selection" 
goto:endRun

:config
	if [%configModeID%] EQU [1] ( 
		SET mCONFIGURATION=Debug
		goto:setTemps
	)
	if [%configModeID%] EQU [2] ( 
		SET mCONFIGURATION=Release
		goto:setTemps
	) 
	echo "Invalid Configuration Selection" 
	goto:endRun

:setTemps
	SET BAG_HOME=%~dp0configdata
	SET MERGEBATHY_DLLS=%~dp0extlibs\libs_win\%mPLATFORM%\%mCONFIGURATION%
	SET MERGEBATHY_EXE=%~dp0%mPLATFORM%\%mCONFIGURATION%
	goto:instance
	
:instance
	if [%instanceModeID%] EQU [1] ( 
		goto:runTemps
	)
	if [%instanceModeID%] EQU [2] ( 
		goto:GetUserPath
	) 
	echo "Invalid Instance Selection: Local Instance Selected by Default." 
	goto:runTemps

:runTemps
	:: Assign usr environment variables
	if not "%Path:~-1%" == ";" ( 
		set "Path=%Path%;"
	)
	SET UserPath=%Path%%MERGEBATHY_DLLS%;%MERGEBATHY_EXE%
	SET Path=%UserPath%
	goto:endRun
	
rem Get directly from Windows registry the user PATH variable value.
:GetUserPath
	set "UserPath="
	setlocal enabledelayedexpansion
	for /F "skip=2 tokens=1,2*" %%N in ('%SystemRoot%\System32\reg.exe query "HKCU\Environment" /v "Path" 2^>nul') do (
		if /I "%%N" == "Path" (
			set "UserPath=%%P"
			echo Last char on Path = "!UserPath:~-1!"
			echo.
			if not "!UserPath:~-1!" == ";" ( 
				set "UserPath=!UserPath!;"
			)
			echo.
			echo The user PATH is = !UserPath!
			echo.
			goto:run
		)
	)
	endlocal
	echo.
	echo There is no user PATH defined.
	set "UserPath="
	echo.
	goto:run 

:run
	:: Assign usr environment variables
	SETX BAG_HOME "%BAG_HOME%"
	SETX MERGEBATHY_DLLS "%MERGEBATHY_DLLS%"
	SETX MERGEBATHY_EXE "%MERGEBATHY_EXE%"
	:: Set usr Path variable
	SET UserPath=%UserPath%%%MERGEBATHY_DLLS%%;%%MERGEBATHY_EXE%%
	SETX Path "%UserPath%"
	goto:endRun
	
:endRun
	::pause
	echo.
	echo BAG_HOME = %BAG_HOME%
	echo.
	echo MERGEBATHY_DLLS = %MERGEBATHY_DLLS%
	echo.
	echo MERGEBATHY_EXE = %MERGEBATHY_EXE%
	echo.
	echo UserPath = %UserPath%
	echo.
	echo Verify by typing (when setting globally, open new terminal):
	echo 		mergeBathy
	echo.
	cmd /k
	EXIT /B
	