**********************************************************************
CC0 License
**********************************************************************
MergeBathy - Tool to combine one or more bathymetric data files onto a single input grid.
Written in 2015 by Samantha J. Zambo (samantha.zambo@gmail.com) while employed by the U.S. Naval Research Laboratory.
Written in 2015 by Todd Holland while employed by the U.S. Naval Research Laboratory.
Written in 2015 by Nathaniel Plant while employed by the U.S. Naval Research Laboratory.
Written in 2015 by Kevin Duvieilh while employed by the U.S. Naval Research Laboratory.
Written in 2015 by Paul Elmore while employed by the U.S. Naval Research Laboratory.
Written in 2015 by Will Avera while employed by the U.S. Naval Research Laboratory.
Written in 2015 by Brian Bourgeois while employed by the U.S. Naval Research Laboratory.
Written in 2015 by A. Louise Perkins while employed by the U.S. Naval Research Laboratory.
Written in 2015 by David Lalejini while employed by the U.S. Naval Research Laboratory.
To the extent possible under law, the author(s) and the U.S. Naval Research Laboratory have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide. This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

**********************************************************************
This ReadMe.txt is for the MergeBathy Software Suite.
**********************************************************************
Developed in C++, MergeBathy (2015) is cross-platform and multi-threaded software suite for constructing digital bathymetric models.  
It provides the user with a set of modelling tools to construct custom bathymetric surfaces, including splines-in-tension routines for interpolation output or as an intermediate resampling step when merging multiple bathymetry data sets.  
Notable to MergeBathy is its user-friendly and flexible processing options made possible from its integrated bathymetric process framework.  
See OSP Publication MergeBathy (2015).

**********************************************************************
This file contains Quick Install instructions for MergeBathy_CPP.  
	The user should refer to mergeBathy_install.rtf in MergeBathy_Repos-mergeBathy_DOCS for detailed instructions.

**********************************************************************
Download & Install:
	Please download the following repositories in the same directory as /mergeBathy_CPP
	Documentation: https://github.com/Sammie-Jo/MergeBathy_Repos-mergeBathy_DOCS
	Examples: https://github.com/Sammie-Jo/MergeBathy_Repos-TEST_CENTER
	Data for Examples: https://github.com/Sammie-Jo/MergeBathy_Repos-DATA_CENTER

	We recommend downloading to a parent directory /MergeBathy_Repos/
		The structure will be:
		/MergeBathy_Repos/mergeBathy_CPP
		/MergeBathy_Repos/TEST_CENTER
		/MergeBathy_Repos/DATA_CENTER
		/MergeBathy_Repos/mergeBathy_DOCS
**********************************************************************
Set Environment Variables:
	Windows: double-click
		setenv.bat for
	Unix:
		>. ./setenv.sh
	
	This sets the usr environment variables.
	To set environment variables for a single terminal/cmd window, run locally when prompted [default].
	To set environment variables permanently, run globally when prompted.
	You must run setenv before making command-line calls to MergeBathy, unless the environment variables are set globally.
	Alternatively, scripts can be used as in the examples, where environment variables are set at the beginning of the script.
**********************************************************************
Verify by typing into terminal:
	mergeBathy
	
	If environment variables were set globally in setenv, you must restart terminal/cmd window.
	If set locally [default], then environment vars are only good for the current terminal/cmd window open.
**********************************************************************
Run Examples:
	To run examples:
	Windows: double-click
		./TEST_CENTER/Test_Set_10_OSP_Paper/Active_Testing_Site/bat_files/runInterpTests_e.bat
		./TEST_CENTER/Test_Set_10_OSP_Paper/Active_Testing_Site/bat_files/runInterpTests_e_meters.bat
	Windows: cmd
		>cd ../TEST_CENTER/Test_Set_10_OSP_Paper/Active_Testing_Site/bat_files
		>runInterpTests_e.bat
		>runInterpTests_e_meters.bat

	Unix: terminal
		>cd ../TEST_CENTER/Test_Set_10_OSP_Paper/Active_Testing_Site/bash_files
		>sh ./runInterpTests_e.sh
		>sh ./runInterpTests_e_meters.sh

	Outputs to:
		../TEST_CENTER/Test_Set_10_OSP_Paper/Active_Testing_Site/output_files
		../TEST_CENTER/Test_Set_10_OSP_Paper/Active_Testing_Site/command_line_execution

=============================================================================================
Possible Errors when compiling from source:

(1) error LNK1104: cannot open file '..\mergeBathy_CPP\x64\Release\mergeBathy.exe' 

	Build > Clean mergeBathy



(2) error LNK1123: failure during conversion to COFF: file invalid or corrupt

	The solution was found via googling at:
	http://stackoverflow.com/questions/10888391/error-link-fatal-error-lnk1123-failure-during-conversion-to-coff-file-inval

	In summary, you can either:
	(A) Reintstall VS 2010 SP1 and/or VS 2010 Compiler Pack

	(B) or the simplest, make sure your using cvtres.exe in the first location below and not the second.  This requires either adjusting your path or renaming the file in the second location. 
	"You may have two versions of the cvtres.exe utility in your path. 
	One at C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\BIN\cvtres.exe 
	and one at C:\Windows\Microsoft.NET\Framework\v4.0.30319\cvtres.exe"



