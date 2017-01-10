//computational_geometry_algo_interp.h

#ifndef COMPUTATIONAL_GEOMETRY_ALGO_INTERP_H
#define COMPUTATIONAL_GEOMETRY_ALGO_INTERP_H

//Standard Template Library
#include <iostream>
#include <ostream>
#include <fstream>
#include <time.h>
#include <math.h>
#include <windows.h>


//Irregular Grid
int build_delaunay_triangulation(std::vector<double>&, std::vector<double>&, std::vector<double>&);
int build_function_values(std::vector<double>&, std::vector<double>&, std::vector<double>&);
int build_input_grid_locs(const std::vector<double>&,const std::vector<double>&,const std::vector<double>&,std::vector<double>&,std::vector<double>&,std::vector<double>&,double, double);
int build_output_grid_locs(const std::vector<double>&,const std::vector<double>&,double,double);

double vector_max(std::vector<double>);
double vector_min(std::vector<double>);

DWORD WINAPI compute_sibson_ntrlnbr_interp_col(LPVOID);
DWORD WINAPI compute_sibson_ntrlnbr_cont_inter_col_thread(LPVOID);
DWORD WINAPI compute_farin_ntrlnbr_cont_inter_col_thread(LPVOID);

int sibson_natural_neighbor_interp();
int sibson_natural_neighbor_cont_interp();
int farin_natural_neighbor_cont_interp();

#endif