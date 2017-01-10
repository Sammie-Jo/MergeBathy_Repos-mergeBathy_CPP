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