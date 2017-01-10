//Standard Template Library
#include <iostream>
#include <ostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <Windows.h>
#include <vector>

////CGAL Library(Computational Geometry Algorithms Library)
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Delaunay_triangulation_2.h>
//#include <CGAL/Delaunay_triangulation_3.h>
//#include <CGAL/Interpolation_traits_2.h>
//#include <CGAL/Interpolation_gradient_fitting_traits_2.h>
//#include <CGAL/natural_neighbor_coordinates_2.h>
//#include <CGAL/interpolation_functions.h>
//#include <CGAL/sibson_gradient_fitting.h>
//
//#include "computational_geometry_algo_interp.h"
////#include "numerical_recipes_interp.h"
////
//////Shared
//typedef CGAL::Exact_predicates_inexact_constructions_kernel		K;
//typedef K::FT                                                   Coord_type;			//This is the Kernel FieldNumberType
//typedef K::Point_2												Point;
//typedef CGAL::Delaunay_triangulation_2<K>						Delaunay_triangulation;
//
////Interpolation specific
//typedef CGAL::Interpolation_traits_2<K>                         Traits;
//typedef CGAL::Interpolation_gradient_fitting_traits_2<K>        Traits_grad;
//
//
////Globals
//Delaunay_triangulation											T;
//std::map<Point,Coord_type, K::Less_xy_2>						function_values;
//std::map<Point, K::Vector_2, K::Less_xy_2>                      function_gradients;
//
////Input grid
//extern std::vector<double> grid_xvi;
//extern std::vector<double> grid_yvi;
//extern std::vector<double> grid_zvi;
//
////Output grid
//extern std::vector<double> grid_xv;
//extern std::vector<double> grid_yv;
//extern std::vector<double> grid_zv;
//
////Execute a maximum of THREAD_POOL_SIZE threads on a maximum
////of CPU_POOL_SIZE central processing units.
////#define THREAD_POOL_SIZE    8
//
//
//struct COLUMN_RNG_STRUCT
//{
//	int begin;
//	int end;
//};

//
///*
//Implemented from Computational Geometry Algorithms Library(CGAL) Reference Manual; http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Interpolation_ref/Chapter_intro.html
//Citations for algorithms
//[Far90]  G. Farin. Surfaces over Dirichlet tesselations. Comput. Aided Geom. Design, 7:281-292, 1990 
//[Sib80]  R. Sibson. A vector identity for the Dirichlet tesselation. Math. Proc. Camb. Phil. Soc., 87:151-155, 1980. 
//[Sib81]  R. Sibson. A brief description of natural neighbour interpolation. In Vic Barnet, editor, Interpreting Multivariate Data, pages 21-36. John Wiley & Sons, Chichester, 1981.  
//*/
//
///*
//The Kernel choice determines if the computation is geometrically exact or 
//an approximation. Tradeoff is speed versus exactness; 
//see section "11.2.6 Choosing a Kernel and Predefined Kernels"
//of user's manual. These will need to be options on the command line once migrated
//to merge bathy system.
//*/
//
//double vector_min(std::vector<double> vector_in)
//{
//	double  dmin = 1e10;
//	
//	for (int i=0; i < vector_in.size(); i++)
//	{
//		//Extract minimum values
//		if (vector_in.at(i) < dmin)
//		{
//			dmin = vector_in.at(i);
//		}
//	}
//	
//	return dmin;
//}
//
//double vector_max(std::vector<double> vector_in)
//{
//	double  dmax = -1e10;
//	
//	for (int i=0; i < vector_in.size(); i++)
//	{
//		//Extract maximum values
//		if (vector_in.at(i) > dmax)
//		{
//			dmax = vector_in.at(i);
//		}
//	}
//	
//	return dmax;
//}
//

//double estimate_cols_per_tile(double total_cols,double prcnt,double min_colseg)
//{
//	double value;
//
//    if (prcnt>100)
//		prcnt = 100; 
//
//    if (prcnt<0)
//		prcnt = 0;
//    
//    if (total_cols <= min_colseg)
//	{
//        value = total_cols;
//	}
//	else
//	{
//		value = total_cols * prcnt/100;
//    
//		if  (value <= min_colseg)
//			value = min_colseg;
//        
//		value = floor(value + .5);
//	}
//
//    return value;
//}
//
//void build_coltile_set(double ncols_total,double ncols_per_tile, std::vector<double>& col_strip_start, std::vector<double>& col_strip_end)
//{
//    double tile_count		= 0;
//    double nRemainderCols	= fmod(ncols_total,ncols_per_tile);
//	col_strip_start.empty();
//	col_strip_end.empty();
//	
//    for (double nColIndex = 1; nColIndex <= (ncols_total - nRemainderCols); nColIndex+=ncols_per_tile)
//	{
//        double col_start = (double)nColIndex;
//        double col_end   = (double)nColIndex+ncols_per_tile-1;
//        tile_count = tile_count + 1;
//		col_strip_start.push_back(col_start);
//        col_strip_end.push_back(col_end);
//	}
//
//    if (nRemainderCols > 0)
//	{
//        double col_start  = ncols_total-nRemainderCols+1;
//        double col_end    = ncols_total;
//        tile_count = tile_count + 1;
//        col_strip_start.push_back(col_start);
//        col_strip_end.push_back(col_end);
//	}
//    return;
//}
//
//
//void grid_segmentor(double num_columns,double num_threads, vector<double>& col_strip_start, vector<double>& col_strip_end)
//{
//	double prcnt_colseg = floor((num_columns/num_threads) + .5) / num_columns * (100.0);
//	std::cout << "prcnt_colseg = " << prcnt_colseg << "\n";
//    double min_colseg          = 10;
//    double ncols_per_tile      = estimate_cols_per_tile(num_columns,prcnt_colseg,min_colseg);
//	cout << "ncols_per_tile = " << ncols_per_tile << "\n";
//
//    build_coltile_set(num_columns,ncols_per_tile,col_strip_start,col_strip_end);
//	
//	double ncolTiles  = (double)col_strip_start.size();
//    
//    //If last segment has less than min_colseg then delete it and add the
//    //cols from this last segment to the previous segment
//	double last_col_cnt = col_strip_end.at(((int)(ncolTiles-1))) - col_strip_start.at((int)(ncolTiles-1))+1; 
//    if (last_col_cnt < min_colseg)
//	{
//        cout << "WARNING: Col segment count fault; merging columns with previous segment.\n";
//        
//        col_strip_end.at(((int)ncolTiles-2)) = col_strip_end.at(((int)ncolTiles-2)) + last_col_cnt;
//        
//		col_strip_start.pop_back();
//		col_strip_end.pop_back();
//	}   
//}

//int build_delaunay_triangulation(std::vector<double>& xvect, std::vector<double>& yvect, std::vector<double>& zvect)
//{
//	int status = 0;
//
//	for (int ndx = 0; ndx < xvect.size(); ndx++)
//	{
//		K::Point_2 p(xvect[ndx],yvect[ndx]);
//		T.insert(p);
//	}
//	
//	return status;
//}
//
//int build_function_values(std::vector<double>& xvect, std::vector<double>& yvect, std::vector<double>& zvect)
//{
//	int status = 0;
//
//	for (int ndx = 0; ndx < xvect.size(); ndx++)
//	{
//		K::Point_2 p(xvect[ndx],yvect[ndx]);
//		function_values.insert(std::make_pair(p,zvect[ndx]));
//	}
//	
//	return status;
//}
//
//DWORD WINAPI compute_farin_ntrlnbr_cont_inter_col_thread(LPVOID param)
//{
//	//Pointer to structure representing the begin and end of
//	//the columns to process with this thread
//	COLUMN_RNG_STRUCT* pRangeStruct = (COLUMN_RNG_STRUCT*)param;
//	
//	int xndx_start = pRangeStruct->begin;
//	int xndx_end   = pRangeStruct->end;
//
//	int status     = 0;
//	int xcnt = (int)grid_xv.size();
//	int ycnt = (int)grid_yv.size();
//
//	//Estimate the gradient on the 2D natural neighbor space
//	sibson_gradient_fitting_nn_2(T, std::inserter(function_gradients,function_gradients.begin()),CGAL::Data_access<std::map<Point,Coord_type,K::Less_xy_2>>(function_values),Traits_grad());
//	std::vector<std::pair<Point, Coord_type>> coords;
//		
//	for (int xndx = xndx_start; xndx <= xndx_end; xndx++)
//	{
//		for (int yndx = 0; yndx < ycnt; yndx++)
//		{
//			K::Point_2 p(grid_xv[xndx],grid_yv[yndx]);
//			std::vector<std::pair<Point, Coord_type>> coords;
//			Coord_type norm = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords)).second;
//
//			try
//			{
//				std::pair<Coord_type, bool> result_farin_c1  = CGAL::farin_c1_interpolation(coords.begin(), coords.end(), norm, p, CGAL::Data_access<std::map<Point, Coord_type,  K::Less_xy_2>>(function_values),CGAL::Data_access<std::map<Point, K::Vector_2, K::Less_xy_2>>(function_gradients), Traits_grad());
//
//				if (result_farin_c1.second == true)
//				{
//					grid_zv[yndx + (ycnt*xndx)] = result_farin_c1.first;
//				}
//				else
//				{
//					Coord_type result_farin = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, CGAL::Data_access<std::map<Point, Coord_type, K::Less_xy_2>>(function_values));
//					grid_zv[yndx + (ycnt*xndx)] = result_farin;
//				}
//			}
//			catch(...)
//			{
//				grid_zv[yndx + (ycnt*xndx)] = NaN; //NR implements IEEE Floating Point NaN.
//			}
//		}
//	}
//
//	status = 1;
//
//	return (DWORD)status;
//}
//
//DWORD WINAPI compute_sibson_ntrlnbr_cont_inter_col_thread(LPVOID param)
//{
//	COLUMN_RNG_STRUCT* pRangeStruct = (COLUMN_RNG_STRUCT*)param;
//
//	int xndx_start = pRangeStruct->begin;
//	int xndx_end   = pRangeStruct->end;
//
//	//Represents a strip from the main grid matrix
//	//xndx_start = start column for this thread
//	//xndx_end   = end column for this thread
//	int status = 0;
//	int xcnt   = (int)grid_xv.size();
//	int ycnt   = (int)grid_yv.size();
//
//	//Estimate the gradient on the 2D natural neighbor space
//	sibson_gradient_fitting_nn_2(T, std::inserter(function_gradients,function_gradients.begin()),CGAL::Data_access<std::map<Point,Coord_type,K::Less_xy_2>>(function_values),Traits_grad());
//	std::vector<std::pair<Point, Coord_type>> coords;
//
//	//create a ndx boundary test method; future
//	for (int xndx = xndx_start; xndx <= xndx_end; xndx++)
//	{
//		for (int yndx = 0; yndx < ycnt; yndx++)
//		{
//			K::Point_2 p(grid_xv[xndx],grid_yv[yndx]);
//			std::vector<std::pair<Point, Coord_type>> coords;
//			Coord_type norm = CGAL::natural_neighbor_coordinates_2(T, p, std::back_inserter(coords)).second;
//
//			std::pair<Coord_type, bool> result_sibson_c1 = CGAL::sibson_c1_interpolation_square(coords.begin(), coords.end(), norm, p,			\
//										CGAL::Data_access<std::map<Point, Coord_type,  K::Less_xy_2>>(function_values),                         \
//										CGAL::Data_access<std::map<Point, K::Vector_2, K::Less_xy_2>>(function_gradients),Traits_grad());
//
//			if (result_sibson_c1.second == true)
//			{
//				grid_zv[yndx + (ycnt*xndx)] = result_sibson_c1.first;
//			}
//			else
//			{
//				Coord_type result_sibson  = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, CGAL::Data_access<std::map<Point, Coord_type, K::Less_xy_2>>(function_values));
//				grid_zv[yndx + (ycnt*xndx)] = result_sibson;
//			}
//		}
//	}
//	
//	status = 1;
//
//	return (DWORD)status;
//}

//DWORD WINAPI compute_sibson_ntrlnbr_interp_col_thread(LPVOID param)
//{
//	COLUMN_RNG_STRUCT* pRangeStruct = (COLUMN_RNG_STRUCT*)param;
//
//	int xndx_start = pRangeStruct->begin;
//	int xndx_end   = pRangeStruct->end;
//
//	//Represents a strip from the main grid matrix
//	//xndx_start = start column for this thread
//	//xndx_end   = end column for this thread
//	int status = 0;
//	int xcnt   = (int)grid_xv.size();
//	int ycnt   = (int)grid_yv.size();
//
//	//create a ndx boundary test method; future
//	for (int xndx = xndx_start; xndx <= xndx_end; xndx++)
//	{
//		for (int yndx = 0; yndx < ycnt; yndx++)
//		{
//			K::Point_2 p(grid_xv[xndx],grid_yv[yndx]);
//			vector<pair<Point, Coord_type>> coords;
//			Coord_type norm = CGAL::natural_neighbor_coordinates_2(T, p, back_inserter(coords)).second;
//			Coord_type res  = CGAL::linear_interpolation(coords.begin(), coords.end(),
//				norm, CGAL::Data_access<map<Point, Coord_type, K::Less_xy_2>>(function_values));
//			grid_zv[yndx + (ycnt*xndx)] = res;
//		}
//	}
//	
//	status = 1;
//
//	return (DWORD)status;
//}

//int sibson_natural_neighbor_interp( )
//{
//	std::cout << "Executing Sibson's Natural Neighbor Interpolation...\n";
//
//	HANDLE      hThreads[THREAD_POOL_SIZE];
//	int         slot = 0;
//	DWORD       threadID;
//	DWORD       rc;
//		
//	int status			= 0;
//	int xcnt			= (int)grid_xv.size();
//	int ycnt			= (int)grid_yv.size();
//	int zcnt			= xcnt*ycnt;
//	double ncols_total	= xcnt;
//	std::vector<double> grid_zv_sub = std::vector<double>(ycnt);
//	grid_zv				= std::vector<double>(zcnt);
//	std::vector<double> col_strip_start;
//	std::vector<double> col_strip_end;
//
//	grid_segmentor(ncols_total,THREAD_POOL_SIZE,col_strip_start,col_strip_end);
//	
//	int    nTasks      = (int)col_strip_start.size();
//	COLUMN_RNG_STRUCT  param_rng_struct;
//
//	std::cout << "Number of threads available in pool = " << THREAD_POOL_SIZE << "\n";
//	std::cout << "Number of strips to process         = " << col_strip_start.size() << "\n";
//
//	for (int task_ndx = 1; task_ndx <= nTasks; task_ndx++)
//	{
//		double ncols_in_strip = col_strip_end.at(task_ndx-1) - col_strip_start.at(task_ndx-1) + 1;
//		std::cout << "Column Strip(TASK) # = " << task_ndx << "\t" << col_strip_start.at(task_ndx-1)-1 << "\t" << col_strip_end.at(task_ndx-1)-1 << ";(" << ncols_in_strip << ")\n";
//		
//		//There are no more threads available in the pool; wait for one to become available before creating another.
//		if (task_ndx > THREAD_POOL_SIZE)
//		{
//			std::cout << "\tWaiting for thread in pool to become available...\n";
//			rc	 = WaitForMultipleObjects(THREAD_POOL_SIZE,hThreads,FALSE,INFINITE);
//			slot = rc - WAIT_OBJECT_0;
//		}
//
//		param_rng_struct.begin = (int)col_strip_start.at(task_ndx-1)-1;
//		param_rng_struct.end   = (int)col_strip_end.at(task_ndx-1)-1;
//
//		hThreads[slot++] = CreateThread(NULL,0,compute_sibson_ntrlnbr_interp_col_thread,(LPVOID) &param_rng_struct,0,&threadID);
//		std::cout << "\tThread in pool available; slot " << slot << "\n";
//	}
//
//	//Ensure all threads in pool are finished
//	std::cout << "Waiting for all threads in pool to finish...\n";
//	rc = WaitForMultipleObjects(THREAD_POOL_SIZE,hThreads,TRUE,INFINITE);
//
//
//	//Terminate all threads in pool 
//	for(slot=0;slot<THREAD_POOL_SIZE;slot++)
//	{
//		CloseHandle(hThreads[slot]);
//		std::cout << "Terminatig thread; slot " << slot << "\n";
//	}
//
//
//	///////////////////
//	//Oringinal(Verified); single threaded version
//	//for (int xndx = 0; xndx < xcnt; xndx++)
//	//{
//	//	std::cout << "Processing column :" << xndx+1 << " of " << xcnt << "\n";
//
//		//Execute with a boost::thread
//		//boost::thread thrd1(boost::bind(&compute_sibson_ntrlnbr_interp_col,xndx));
//		//thrd1.join();
//		
//	    //Execute as a single main thread
//		//status = compute_sibson_ntrlnbr_interp_col(xndx);
//	//}
//	return status;
//}
//
//int sibson_natural_neighbor_cont_interp()
//{
//	std::cout << "Executing Sibson Continuous; Natural Neighbor Interpolation w/ Gradient Est...\n";
//
//	HANDLE hThreads[THREAD_POOL_SIZE];
//	int    slot = 0;
//	DWORD  threadID;
//	DWORD  rc;
//
//	int status = 0;
//	int xcnt = (int)grid_xv.size();
//	int ycnt = (int)grid_yv.size();
//	int zcnt = xcnt*ycnt;
//	double ncols_total = xcnt;
//	std::vector<double> grid_zv_sub = std::vector<double>(ycnt);
//	grid_zv = std::vector<double>(zcnt);
//	std::vector<double> col_strip_start;
//	std::vector<double> col_strip_end;
//
//	grid_segmentor(ncols_total,THREAD_POOL_SIZE,col_strip_start,col_strip_end);
//	int nTasks = (int)col_strip_start.size();
//	COLUMN_RNG_STRUCT	param_rng_struct;
//	
//	for (int task_ndx = 1; task_ndx <= nTasks; task_ndx++)
//	{
//		double ncols_in_strip = col_strip_end.at(task_ndx-1) - col_strip_start.at(task_ndx-1) + 1;
//		std::cout << "Column Strip(TASK) # = " << task_ndx << "\t" << col_strip_start.at(task_ndx-1)-1 << "\t" << col_strip_end.at(task_ndx-1)-1 << ";(" << ncols_in_strip << ")\n";
//		if (task_ndx > THREAD_POOL_SIZE)
//		{
//			std::cout << "\tWaiting for thread in pool to become available...\n";
//			rc = WaitForMultipleObjects(THREAD_POOL_SIZE,hThreads,FALSE,INFINITE);
//			slot = rc - WAIT_OBJECT_0;
//		}
//
//		param_rng_struct.begin = (int)col_strip_start.at(task_ndx-1)-1;
//		param_rng_struct.end   = (int)col_strip_end.at(task_ndx-1)-1;
//
//		LPSECURITY_ATTRIBUTES   lpThreadAttributes = NULL;
//		SIZE_T                  dwStackSize        = 0;
//		LPTHREAD_START_ROUTINE  lpStartAddress     = compute_sibson_ntrlnbr_cont_inter_col_thread;
//		LPVOID                  lpParameter        = (LPVOID)&param_rng_struct;
//		DWORD                   dwCreationFlags    = 0;
//		LPDWORD                 lpThreadId         = &threadID;
//		hThreads[slot++] = CreateThread(lpThreadAttributes,dwStackSize,lpStartAddress,lpParameter,dwCreationFlags,lpThreadId);
//		std::cout << "\tThread in pool available; slot " << slot << "\n";
//	}
//
//	rc = WaitForMultipleObjects(THREAD_POOL_SIZE,hThreads,TRUE,INFINITE);
//
//	//Terminate all threads 
//	for(slot = 0; slot<THREAD_POOL_SIZE;slot++)
//	{
//		CloseHandle(hThreads[slot]);
//		std::cout << "Terminating thread; slot " << slot << "\n";
//	}
//
//	std::cout << "Executing Sibson Continuous; Natural Neighbor Interpolation w/ Gradient Est...(DONE)\n";
//
//	return status;
//}
//
//int farin_natural_neighbor_cont_interp()
//{
//	std::cout << "Executing Farin's; Natural Neighbor Interpolation...\n";
//
//	HANDLE hThreads[THREAD_POOL_SIZE];
//	int    slot = 0;
//	DWORD  threadID;
//	DWORD  rc;
//
//	int status = 0;
//	int xcnt = (int)grid_xv.size();
//	int ycnt = (int)grid_yv.size();
//	int zcnt = xcnt*ycnt;
//	double ncols_total = xcnt;
//	std::vector<double> grid_zv_sub = std::vector<double>(ycnt);
//	grid_zv = std::vector<double>(zcnt);
//	std::vector<double> col_strip_start;
//	std::vector<double> col_strip_end;
//
//	grid_segmentor(ncols_total,THREAD_POOL_SIZE,col_strip_start,col_strip_end);
//	int nTasks = (int)col_strip_start.size();
//	COLUMN_RNG_STRUCT	param_rng_struct;
//	
//	for (int task_ndx = 1; task_ndx <= nTasks; task_ndx++)
//	{
//		double ncols_in_strip = col_strip_end.at(task_ndx-1) - col_strip_start.at(task_ndx-1) + 1;
//		std::cout << "Column Strip(TASK) # = " << task_ndx << "\t" << col_strip_start.at(task_ndx-1)-1 << "\t" << col_strip_end.at(task_ndx-1)-1 << ";(" << ncols_in_strip << ")\n";
//		if (task_ndx > THREAD_POOL_SIZE)
//		{
//			std::cout << "\tWaiting for thread in pool to become available...\n";
//			rc = WaitForMultipleObjects(THREAD_POOL_SIZE,hThreads,FALSE,INFINITE);
//			slot = rc - WAIT_OBJECT_0;
//		}
//
//		param_rng_struct.begin = (int)col_strip_start.at(task_ndx-1)-1;
//		param_rng_struct.end   = (int)col_strip_end.at(task_ndx-1)-1;
//
//		LPSECURITY_ATTRIBUTES   lpThreadAttributes = NULL;
//		SIZE_T                  dwStackSize        = 0;
//		LPTHREAD_START_ROUTINE  lpStartAddress     = compute_farin_ntrlnbr_cont_inter_col_thread;
//		LPVOID                  lpParameter        = (LPVOID)&param_rng_struct;
//		DWORD                   dwCreationFlags    = 0;
//		LPDWORD                 lpThreadId         = &threadID;
//		hThreads[slot++] = CreateThread(lpThreadAttributes,dwStackSize,lpStartAddress,lpParameter,dwCreationFlags,lpThreadId);
//		std::cout << "\tThread in pool available; slot " << slot << "\n";
//	}
//
//	rc = WaitForMultipleObjects(THREAD_POOL_SIZE,hThreads,TRUE,INFINITE);
//
//	//Terminate all threads 
//	for(slot = 0; slot<THREAD_POOL_SIZE;slot++)
//	{
//		CloseHandle(hThreads[slot]);
//		std::cout << "Terminating thread; slot " << slot << "\n";
//	}
//
//	std::cout << "Executing Farin's; Natural Neighbor Interpolation...(DONE)\n";
//
//	return status;
//}
//
//int build_input_grid_locs(const std::vector<double>& xvect,    \
//	                      const std::vector<double>& yvect,    \
//						  const std::vector<double>& zvect,    \
//						  std::vector<double>& grid_xvi,       \
//						  std::vector<double>& grid_yvi,       \
//						  std::vector<double>& grid_zvi,       \
//						  double dx, double dy)
//{
//	//User has input a gridded file.
//	int status = 0;
//
//	
//	try
//	{
//		int		zcnt = (int)zvect.size();
//		std::cout << "Building input grid locations...";
//		double dxmin_in = vector_min(xvect);
//		double dxmax_in = vector_max(xvect);
//		double dymin_in = vector_min(yvect);
//		double dymax_in = vector_max(yvect);
//
//		double  tempvalx= dxmin_in;
//		int     xcnt    = 0;
//		while (tempvalx <= dxmax_in)
//		{
//			grid_xvi.push_back(tempvalx);
//			tempvalx	= tempvalx + dx;
//			xcnt		= xcnt + 1;
//		}
//	
//		double  tempvaly= dymin_in;
//		int     ycnt	= 0;
//		while (tempvaly <= dymax_in)
//		{
//			grid_yvi.push_back(tempvaly);
//			tempvaly	= tempvaly + dy;
//			ycnt		= ycnt + 1;
//		}
//
//		int  xycnt = xcnt*ycnt;
//		grid_zvi		= std::vector<double>(xycnt);
//
//		std::cout << "(Done)\n";
//		std::cout << "\t" << xcnt  << "  X computed grid locations.\n";
//		std::cout << "\t" << ycnt  << "  Y computed grid locations.\n";
//		std::cout << "\t" << xycnt << " XY Product; Grid Sample Count.\n";
//		std::cout << "\t" << zcnt  << "  Z input.\n";
//		std::cout << "\t X min = " << dxmin_in << "\n";
//		std::cout << "\t X max = " << dxmax_in << "\n";
//		std::cout << "\t Y min = " << dymin_in << "\n";
//		std::cout << "\t Y max = " << dymax_in << "\n";
//		if (zcnt != xycnt)
//		{
//			std::cout << "Dataset input not a grid.\n";
//			status = -1;
//		}
//		else
//		{
//			std::cout << "Dataset input is a grid.\n";
//			status = 1;
//
//			for (int yndx = 0; yndx < ycnt; yndx++)
//			{
//				for (int xndx = 0; xndx < xcnt; xndx++)
//				{
//					grid_zvi[xndx + (xcnt*yndx)] = zvect[xndx + (xcnt*yndx)];
//				}
//			}
//		}
//	}
//	catch(...)
//	{
//		status = -1;
//		std::cout << "Exception in build_input_grid_locs()\n";
//	}
//
//	return status;
//
//}
//
//
//int build_output_grid_locs(const std::vector<double>& xvect,   \
//						   const std::vector<double>& yvect,   \
//						   double dx, double dy) 
//{
//	int status = 0;
//
//	try
//	{
//		std::cout << "Building output grid locations...";
//		double dxmin_in = vector_min(xvect);
//		double dxmax_in = vector_max(xvect);
//		double dymin_in = vector_min(yvect);
//		double dymax_in = vector_max(yvect);
//
//		double  tempvalx= dxmin_in;
//		int     xcnt    = 0;
//		while (tempvalx <= dxmax_in)
//		{
//			grid_xv.push_back(tempvalx);
//			tempvalx	= tempvalx + dx;
//			xcnt		= xcnt + 1;
//		}
//	
//		double  tempvaly= dymin_in;
//		int     ycnt	= 0;
//		while (tempvaly <= dymax_in)
//		{
//			grid_yv.push_back(tempvaly);
//			tempvaly	= tempvaly + dy;
//			ycnt		= ycnt + 1;
//		}
//
//		std::cout << "(Done)\n";
//		std::cout << "\t" << xcnt << " X grid locations.\n";
//		std::cout << "\t" << ycnt << " Y grid locations.\n";
//		std::cout << "\t X min = " << dxmin_in << "\n";
//		std::cout << "\t X max = " << dxmax_in << "\n";
//		std::cout << "\t Y min = " << dymin_in << "\n";
//		std::cout << "\t Y max = " << dymax_in << "\n";
//	}
//	catch(...)
//	{
//		status = -1;
//		std::cout << "Exception in build_output_grid_locs()\n";
//	}
//
//	status = 1;
//
//	return status;
//}