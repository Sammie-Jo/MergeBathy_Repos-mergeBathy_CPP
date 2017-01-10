##
# MergeBathy 4.0 C++ Makefile
#
# Created by Louise Perkins and Samantha Zambo
# 1/1/2015
#
# To run:
#
# >make clean
# >make BITFLAG=-m32 CONFIGFLAG=Debug
# >make BITFLAG=-m32 CONFIGFLAG=Release (this is the same as >make)
# >make BITFLAG=-m64 CONFIGFLAG=Debug
# >make BITFLAG=-m64 CONFIGFLAG=Release
##

#===============================================================================
# Define Compilers and Compiler Options.
#===============================================================================
# GNU C/C++ Compilers.
GCC = gcc -w  -fpermissive 
CPP = g++ -w  -fpermissive

# Bit and Configuration Flags.
BITFLAG = -m32
CONFIGFLAG = Debug
ifeq (${BITFLAG},$(filter $(BITFLAG),-m64))
	BITFOLDER=x64
else
	BITFOLDER=x86
endif

# Additional Flags.
GMT_FLAGS = -DGMT_SHARE_PATH=\".\" #-DNVLinux 

#Check Operating System and Architecture
ifeq ($(OS),Windows_NT) #Windows
    $(info Windows OS)
	LDLIBS = -lws2_32 -lopengl32 
	
	LDFLAGS = \
	-L./extlibs/libs_win/$(BITFOLDER) \
	-L./extlibs/libs_win/$(BITFOLDER)/$(CONFIGFLAG)

	CFLAGS = \
	-I./extlibs/libs_win/includes/beecrypt \
	-I./extlibs/libs_win/includes/hdf5 \
	-I./extlibs/libs_win/includes \
	-I./api \
	-I./mergeBathy/includes \
	-I./mergeBathy/GSF \
	-I./mergeBathy/MB_ZGrid \
	-I./mergeBathy/GMT_Surface \
	-I./mergeBathy/GMT_Surface/include \
	-I./mergeBathy/Error_Estimator 

    #PROCESSOR_ARCHITEW6432 
    CFLAGS += -D WIN32 -D _WIN32 -D NVWIN3X -D _CONSOLE #-D _MBCS -D BOOST_ALL_NO_LIB 
#    ifeq ($(PROCESSOR_ARCHITECTURE),AMD64) #AMD64 64-bit Architecture (x86_64, x64)
#        CFLAGS += -D AMD64 -D _WIN64
#        $(info AMD64 Processor)
#    endif
#    ifeq ($(PROCESSOR_ARCHITECTURE),x86) #IA32 Intel Architecture 32-bit (i386)
#        CFLAGS += -D IA32
#        $(info IA32 Processor)
#    endif
else #anything else

	LDLIBS = -pthread 
	GMTFLAGS += -DNVLinux 
#    UNAME_S := $(shell uname -s)
#    ifeq ($(UNAME_S),Linux) #Linux
#        CCFLAGS += -D LINUX
#        $(info Linux OS)
#    endif
#    ifeq ($(UNAME_S),Darwin) #Darwin
#        CCFLAGS += -D OSX
#        $(info Darwin OS)
#    endif
#    UNAME_P := $(shell uname -p)
#    ifeq ($(UNAME_P),x86_64) #AMD64 64-bit Architecture (x86_64, x64)
#        CCFLAGS += -D AMD64
#        $(info AMD64 Processor)
#    endif
#    ifneq ($(filter %86,$(UNAME_P)),) #IA32 Intel Architecture 32-bit (i386)
#        CCFLAGS += -D IA32
#        $(info IA32 Processor)
#    endif
#    ifneq ($(filter arm%,$(UNAME_P)),) #ARM Architecture
#        CCFLAGS += -D ARM
#        $(info ARM Architecture)
#    endif

    # Library Flags.
    LDLIBS += -lxml2 -lhdf5 -lbeecrypt #-lpthread

    # Non-Library Linker Flags. (-L = lib locations)
    LDFLAGS = \
    -L./extlibs/libs_unix/bag/$(BITFOLDER) \
    -L./extlibs/libs_unix/libxml2/$(BITFOLDER)/$(CONFIGFLAG)/lib \
    -L./extlibs/libs_unix/hdf5/$(BITFOLDER)/$(CONFIGFLAG)/lib \
    -L./extlibs/libs_unix/beecrypt/$(BITFOLDER)/$(CONFIGFLAG)/lib 

    # Compiler Flags. (-I = includes locations)
    CFLAGS = \
    -I./extlibs/libs_unix/beecrypt/$(BITFOLDER)/$(CONFIGFLAG)/include \
    -I./extlibs/libs_unix/hdf5/$(BITFOLDER)/$(CONFIGFLAG)/include \
    -I./extlibs/libs_unix/libxml2/$(BITFOLDER)/$(CONFIGFLAG)/include/libxml2 \
    -I./api \
    -I./mergeBathy/includes \
    -I./mergeBathy/GSF \
    -I./mergeBathy/MB_ZGrid \
    -I./mergeBathy/GMT_Surface \
    -I./mergeBathy/GMT_Surface/include \
    -I./mergeBathy/Error_Estimator  
endif

# Optimization and Debug Flags.
ifeq ($(findstring, Debug, $(CONFIGFLAG)),Debug) #Debug
	CFLAGS += -D_DEBUG" -Og -Z7 -ggdb -DBOOST_UBLAS_NDEBUG -D_DISABLE_3RDPARTY_WARNINGS=1 -oo
	LDFLAGS += -DEBUG
	baglib=-lbagd
	ifeq ($(OS),Windows_NT)
		CFLAGS += -MDd  
	endif
else #Release
	CFLAGS += -DNDEBUG -O2 
	baglib=-lbag
	ifeq ($(OS),Windows_NT)
		CFLAGS += -MD 
	endif
endif

CFLAGS += ${LDLIBS} ${LDFLAGS} 

#===============================================================================
# Define Output Directories.
#===============================================================================
# Directories
OUTPUT_FILE = "./$(BITFOLDER)/${CONFIGFLAG}/mergeBathy"
INTERMEDIATE_DIR = ./mergeBathy/$(BITFOLDER)/${CONFIGFLAG}
BINDIR = ./$(BITFOLDER)/${CONFIGFLAG}

# Print info.
$(info $$var is [${OUTPUT_FILE}])
$(info $$var is [${INTERMEDIATE_DIR}])
$(info $$var is [${BINDIR}])
$(info $$var is [${BITFLAG}])
#===============================================================================
# Define Object Files.
#===============================================================================
# OBJS
OBJS = \
mergeBathy.o \
fileReader.o \
fileWriter.o \
fileBagWriter.o \
xmlWriter.o \
grid.o \
LatLong-UTMconversion.o \
rng.o \
externalInterpolators.o \
mergeBathyOld.o \
computeOffset.o \
standardOperations.o \
bathyTool.o \
subSampleData.o \
consistentWeights.o \
regr_xzw.o \
scalecInterpTile.o \
scalecInterp.o \
scalecInterpPerturbations.o \
kriging.o \
MB_Threads.o

# GSF_OBJS
GSF_OBJS = \
geod.o \
geodesic.o \
gsf.o \
gsf_dec.o \
gsf_enc.o \
gsf_indx.o \
gsfReader.o \
gsfSensorName.o

# MB_ZGRID_OBJS
MB_ZGRID_OBJS = mb_zgrid.o

# GMT_SURF_OBJS
GMT_SURF_OBJS = \
gmt_bcr.o \
gmt_calclock.o \
gmt_cdf.o \
gmt_customio.o \
gmt_grdio.o \
gmt_init.o \
gmt_io.o \
gmt_map.o \
gmt_nc.o \
gmt_proj.o \
gmt_stat.o \
gmt_support.o \
gmt_vector.o \
surf.o \
surface.o \
processSurface.o

# ERR_EST_OBJS
ERR_EST_OBJS = \
geom.o \
mesh.o \
pingList.o \
pointList.o \
sHullDelaunay.o \
vincenty.o \
Bathy_Grid.o \
GradientGrid.o

# ALG_OBJS
ALG_OBJS = \
alglibmisc.o \
alglibinternal.o \
ap.o \
linalg.o \
dataanalysis.o \
diffequations.o \
fasttransforms.o \
integration.o \
interpolation.o \
optimization.o \
solvers.o \
specialfunctions.o \
statistics.o

# Append directories to objects.
OUT_OBJS=$(addprefix ${INTERMEDIATE_DIR}/,${OBJS})
OUT_GSF_OBJS=$(addprefix ${INTERMEDIATE_DIR}/,${GSF_OBJS})
OUT_MB_ZGRID_OBJS=$(addprefix ${INTERMEDIATE_DIR}/,${MB_ZGRID_OBJS})
OUT_SURF_OBJS=$(addprefix ${INTERMEDIATE_DIR}/,${GMT_SURF_OBJS})
OUT_ERR_EST_OBJS=$(addprefix ${INTERMEDIATE_DIR}/,${ERR_EST_OBJS})
OUT_ALG_OBJS=$(addprefix ${INTERMEDIATE_DIR}/,${ALG_OBJS})
#===============================================================================
# Define Rules.
#===============================================================================
# Build MergeBathy Target.
mergeBathy : ${BINDIR} ${INTERMEDIATE_DIR} \
	${OUT_GSF_OBJS} \
	${OUT_MB_ZGRID_OBJS} \
	${OUT_SURF_OBJS} \
	${OUT_ERR_EST_OBJS} \
	${OUT_OBJS} \
	${OUT_ALG_OBJS}
	
	${CPP} ${BITFLAG} ${CFLAGS} ${GMT_FLAGS} -lm ${OUT_GSF_OBJS} ${OUT_MB_ZGRID_OBJS} ${OUT_SURF_OBJS} ${OUT_ERR_EST_OBJS} ${OUT_OBJS} ${OUT_ALG_OBJS} $(baglib) -o ${OUTPUT_FILE}

# Clean 32bit object files.
clean-x86 :
	rm -f ./mergeBathy/x86/Debug/*.o
	rm -f ./mergeBathy/x86/Release/*.o

# Clean 64bit object files.
clean-x64 : 
	rm -f ./mergeBathy/x64/Debug/*.o
	rm -f ./mergeBathy/x64/Release/*.o
	
# Clean 32 and 64 bit object files.	
clean : 
	rm -f ./mergeBathy/x86/Debug/*.o
	rm -f ./mergeBathy/x86/Release/*.o
	rm -f ./mergeBathy/x64/Debug/*.o
	rm -f ./mergeBathy/x64/Release/*.o
	
# Makes the output directory.
${BINDIR} :
	if [ ! -d "$(BINDIR)" ];then     \
                mkdir -p $(BINDIR);           \
        fi

# Makes the intermediate directory.
${INTERMEDIATE_DIR} :
	if [ ! -d "$(INTERMEDIATE_DIR)" ];then     \
                mkdir -p $(INTERMEDIATE_DIR);           \
        fi

#ALG_OBJS
${INTERMEDIATE_DIR}/%.o : ./mergeBathy/ALG/%.cpp
	${CPP} ${BITFLAG} ${CFLAGS} -c $< -o $@

#GSF_OBJS
${INTERMEDIATE_DIR}/%.o : ./mergeBathy/GSF/%.c
	${GCC} ${BITFLAG} ${CFLAGS} -c $< -o $@

#MB_ZGRID_OBJS
${INTERMEDIATE_DIR}/%.o : ./mergeBathy/MB_ZGrid/%.c
	${GCC} ${BITFLAG} ${CFLAGS} -c $< -o $@

#GMT_SURF_OBJS
${INTERMEDIATE_DIR}/%.o : ./mergeBathy/GMT_Surface/%.c
	${GCC} ${BITFLAG} ${GMT_FLAGS} ${CFLAGS} -c $< -o $@

#ERR_EST_OBJS
${INTERMEDIATE_DIR}/%.o : ./mergeBathy/Error_Estimator/%.cpp
	${CPP} ${BITFLAG} ${CFLAGS} -c $< -o $@

#OBJS
${INTERMEDIATE_DIR}/%.o : ./mergeBathy/%.cpp
	${CPP} ${BITFLAG} ${CFLAGS} -c $< -o $@

