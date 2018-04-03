#!/bin/bash
# Author: Samantha J Zambo
# Date: 3/22/18
# For OSP Publication 'MergeBathy 2015'
#
# To Run:
# 	1) Open terminal
# 	2) Change directory to folder containing runALL.sh
# 	3) Type ". ./setenv.sh"
# 	4) if permission is denied type "chmod +rxw setenv.sh"
#
##
shopt -s expand_aliases
# Must set this option, else script will not expand aliases.

alias cls='printf "\033c"'
echo $(cls)

echo ""
echo "***********************************************"
echo "* setenv.bat"
echo "* This script will set the usr environment variables to run MergeBathy."
echo "* If set globally, environment variables are set permanently."
echo "* If set locally, environment variables are set only for the current session. Restart required to reset."
echo "* "
echo "* To work in both x86 and x64, and both Debug and Release,"
echo "* set environment variables locally for a single terminal/cmd instance."
echo "* If set locally, you must run this script before calling MergeBathy for every new terminal/cmd window."
echo "* Alternatively, make a script of MergeBathy calls and define the vars in the beginning of the script."
echo "***********************************************"
echo ""
#=====================================================================
# Find platform and configuration 
echo "Warning: GMT and MBZ pre-splining routines will crash in 32bit for large data sets."
read -p "Select (1) 32bit (2) 64bit (3) both : "  arg1
read -p "Select (1) Debug (2) Release (3) both : " arg2
read -p "Select (1) Locally (current instance only) (2) Globally (permanently in usr envars) : " arg3

if [ $arg1 ] ;then bitModeID=$arg1 ;fi
if [ $arg2 ] ;then configModeID=$arg2 ;fi
instanceModeID=1
if [ $arg3 ] ;then instanceModeID=$arg3 ;fi

if [ $bitModeID -eq 1 ] ;then bitFolder=x86
elif [ $bitModeID -eq 2 ] ;then bitFolder=x64
else echo "Invalid bit selection. Exit." ;exit $?
fi

#:config
if [ $configModeID -eq 1 ] ;then configFolder=Debug
elif [ $configModeID -eq 2 ] ;then configFolder=Release
else echo "Invalid configuration selection. Exit." ;exit $?
fi

#:setTemps
# Gets bash equivalent of %~dp0
cd `dirname $0`
homed=`pwd`
cd -
echo ""
BAG_HOME=$homed/configdata
MERGEBATHY_DLLS=.:/usr/local/lib:$homed/extlibs/libs_unix/bag/$bitFolder/$configFolder/lib:$homed/extlibs/libs_unix/libxml2/$bitFolder/$configFolder/lib:$homed/extlibs/libs_unix/hdf5/$bitFolder/$configFolder/lib:$homed/extlibs/libs_unix/beecrypt/$bitFolder/$configFolder/lib
MERGEBATHY_EXE=$homed/$bitFolder/$configFolder

#:instance
if [ $instanceModeID -eq 1 ] ;then GetUserPATH=0
elif [ $instanceModeID -eq 2 ] ;then GetUserPATH=1
else echo "Invalid instance selection. Local instance selected by default." ;GetUserPATH=0
fi

if [ $GetUserPATH -eq 0 ] ;then
	#:runTemps
	# Assign usr environment variables
	UserPATH=$PATH:$MERGEBATHY_DLLS:$MERGEBATHY_EXE
	export PATH=$UserPATH
	export BAG_HOME=$BAG_HOME
	export MERGEBATHY_DLLS=$MERGEBATHY_DLLS
	export MERGEBATHY_EXE=$MERGEBATHY_EXE
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MERGEBATHY_DLLS
else
	if [[ ! -z $PATH ]] ;then
		UserPATH=$PATH
		echo "The user PATH is = " ;echo $UserPATH ;echo
	else
		echo "There is no user PATH defined." ;echo
		UserPATH=
	fi

	#:run
	# delete vars (delete line if found)
	sed -i '/export BAG_HOME=/c\' ~/.bashrc
	sed -i '/export MERGEBATHY_DLLS=/c\' ~/.bashrc
	sed -i '/export MERGEBATHY_EXE=/c\' ~/.bashrc
	
	# delete vars appended to PATH if found
	sed -ie '/export PATH/s/"$MERGEBATHY_DLLS":"$MERGEBATHY_EXE"//g' ~/.bashrc
	sed -ie '/export LD_LIBRARY_PATH/s/:"$MERGEBATHY_DLLS"//g' ~/.bashrc	
	PATH=`echo $PATH | sed -e 's/:"$MERGEBATHY_DLLS":"$MERGEBATHY_EXE"$//g'`
	#echo PATH=$PATH

	num=$(grep -n 'export PATH' ~/.bashrc | awk -F: '{print $1}')
	# Assign usr environment variables
	if [[ ! -z "$num" ]] ;then
		sed -i "${num}i\export BAG_HOME=\"$BAG_HOME\"" ~/.bashrc 
		((num++))
		sed -i "${num}i\export MERGEBATHY_DLLS=\"$MERGEBATHY_DLLS\"" ~/.bashrc 
		((num++))
		sed -i "${num}i\export MERGEBATHY_EXE=\"$MERGEBATHY_EXE\"" ~/.bashrc 
	else
		echo 'export BAG_HOME="'"$BAG_HOME"'"' >> ~/.bashrc 
		echo 'export MERGEBATHY_DLLS="'"$MERGEBATHY_DLLS"'"' >> ~/.bashrc 
		echo 'export MERGEBATHY_EXE="'"$MERGEBATHY_EXE"'"' >> ~/.bashrc 
	 fi
	# Set usr PATH variable
	UserPATH=$PATH:$MERGEBATHY_DLLS:$MERGEBATHY_EXE
	
	# Append to PATH if exists else write PATH to .bashrc
	grep -q 'export PATH' ~/.bashrc && sed -i '/export PATH=/ s/$/:"$MERGEBATHY_DLLS":"$MERGEBATHY_EXE"/' ~/.bashrc  || echo 'export PATH=$PATH:"$MERGEBATHY_DLLS":"$MERGEBATHY_EXE"' >> ~/.bashrc 
	grep -q 'export LD_LIBRARY_PATH' ~/.bashrc && sed -i '/export LD_LIBRARY_PATH=/ s/$/:"$MERGEBATHY_DLLS"/' ~/.bashrc  || echo 'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":"$MERGEBATHY_DLLS"' >> ~/.bashrc 
	sed -i 's/::/:/g' ~/.bashrc
	sed -i 's/::/:/g' ~/.bashrc
fi

echo
echo BAG_HOME = ;echo $BAG_HOME ;echo
echo
echo MERGEBATHY_DLLS = ;echo $MERGEBATHY_DLLS ;echo
echo
echo MERGEBATHY_EXE = ;echo $MERGEBATHY_EXE ;echo
echo
echo UserPATH = ;echo $PATH ;echo
echo
echo "Verify by typing (when setting globally, open new terminal):"
echo "		mergeBathy"
echo
echo "if permission denied try: chmod 777 -R ./x64/* ./x86/*"
echo
echo

wait


