#include "Bathy_Grid.h"
#include "GradientGrid.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <math.h>
#include <list>
#include <vector>
#include "../grid.h"
#include "pointList.h"
#include "sHullDelaunay.h"

void Bathy_Grid::clear()
{
	if(ptl != NULL)
	{
		delete ptl;
		ptl = NULL;
	}

	if(tin != NULL)
	{
		delete tin;
		tin = NULL;
	}

	size_t sz = interpGrids.size();
    for(size_t i = 0; i < sz; ++i){
		delete interpGrids[i];
		interpGrids[i] = NULL;
	}
	interpGrids.clear();
	ensemble_X.clear();//sam added 12/28/15 don't know why it wasnt already here
	ensemble_Y.clear();//sam added 12/28/15 don't know why it wasnt already here
	ensemble_Z.clear();//sam added 12/28/15 don't know why it wasnt already here
	ensemble_H.clear();//sam added 12/28/15 don't know why it wasnt already here
	ensemble_V.clear();//sam added 12/28/15 don't know why it wasnt already here
	ensemble_U.clear();//sam added 12/28/15 don't know why it wasnt already here
	ensemble_U2.clear();
	ensemble_U3.clear();
	ensemble_U4.clear();
	ensemble_U5.clear();
}

void InterpGrid::clear()
{
	if(ptlSurface != NULL)
	{
		delete ptlSurface;
		ptlSurface = NULL;
	}

	if(grads != NULL)
	{
		delete grads;
		grads = NULL;
	}

	triangles.clear();
	e.clear();
	e2.clear();
	e3.clear();
	e4.clear();
	e5.clear();
	z.clear();
	s.clear();

	//For Raster/Bag Interpolation SJZ
	z0.clear();
	e0.clear();
	zK.clear();
	eK.clear();
	nei.clear();
	rei.clear();
}

void Bathy_Grid::Construct_TinRaster(vector<double> *xInit, vector<double> *yInit, vector<double> *zInit, vector<double> *eInit)
{
	//Each BathyGrid object gets its own copy.
	ptl = new PointList();
	ptl->setFromVectorsRaster(*xInit, *yInit, *zInit, *eInit/*, *z0Init, *e0Init, *zKInit, *eKInit*/);

	tin = new SHullDelaunay();
	std::cout << "Initializing triangle for Raster Conversion." << std::endl;
	tin->insert(*ptl);

	std::cout << "Done Converting to Raster." << std::endl;
}
void Bathy_Grid::Construct_TinRaster(vector<double> *xInit, vector<double> *yInit, vector<double> *zInit, vector<double> *eInit, vector<double> *neiInit, vector<double> *reiInit, vector<double> *z0Init, vector<double> *e0Init, vector<double> *zKInit, vector<double> *eKInit)
{
	//Each BathyGrid object gets its own copy.
	ptl = new PointList();
	ptl->setFromVectorsRaster(*xInit, *yInit, *zInit, *eInit, *neiInit, *reiInit, *z0Init, *e0Init, *zKInit, *eKInit);

	tin = new SHullDelaunay();
	std::cout << "Initializing triangle for Raster Conversion." << std::endl;
	tin->insert(*ptl);

	std::cout << "Done Converting to Raster." << std::endl;
}
//void InterpGrid::estimateRaster(vector<double> *xSurf, vector<double> *ySurf, vector<double> *zSurf, const double& sH, const double& alpha, const double& deltaMin, SHullDelaunay* tin, string interpMethod)
//{
//	ptlSurface = new PointList();
//	ptlSurface->setFromVectors(*xSurf, *ySurf, *zSurf, tin->getOffsetX(), tin->getOffsetY());
//
//	grads = new GradientGrid();
//
//	//triangles = vector<Triangle>();
//	triangles.reserve(ptlSurface->size());
//
//	//Compute Gradients if depths are known.
//	if(ptlSurface->getAvgZ() != 0.00)
//		(*grads).calc_GradientGrid(*xSurf, *ySurf, *zSurf);
//	//Compute Depths and Gradients
//	else
//	{
//		this->scaleFactor = sH;
//		this->alpha = alpha;
//		this->deltaMin = deltaMin;
//		vector<Point> pVec(tin->determinePingLocationInTriangle(*ptlSurface, sH, alpha, deltaMin,*grads, triangles, interpMethod));
//
//		e.reserve(pVec.size());
//		z.reserve(pVec.size());
//		//s.reserve(pVec.size());
//
//		for (int i = 0; i < (int) pVec.size(); i++)
//		{
//			e.push_back(pVec[i].u);
//			z.push_back(pVec[i].z);
//			//s.push_back(pVec[i].s);
//		}
//	}
//}

void Bathy_Grid::Construct_Tin(vector<double> *xInit, vector<double> *yInit, vector<double> *zInit, vector<double> *hInit, vector<double> *vInit)
{
	//Each BathyGrid object gets its own copy.
	ptl = new PointList();
	ptl->setFromVectors(*xInit, *yInit, *zInit, *hInit, *vInit);

	tin = new SHullDelaunay();
	//std::cout << "Initializing triangle for error computation" << std::endl;
	tin->insert(*ptl);

	//std::cout << "Done initializing" << std::endl;
}

void InterpGrid::estimate(vector<double> *xSurf, vector<double> *ySurf, vector<double> *zSurf, const double& sH, const double& alpha, const double& deltaMin, SHullDelaunay* tin, string depthInterpMethod, string errorInterpMethod, string extrapMethod)
{
	ptlSurface = new PointList();
	ptlSurface->setFromVectors(*xSurf, *ySurf, *zSurf, tin->getOffsetX(), tin->getOffsetY());

	grads = new GradientGrid();

	//triangles = vector<Triangle>();
	triangles.reserve(ptlSurface->size());

	//Compute Gradients before-hand if depths are known.
	//This catches if we pre-splined.
	if(ptlSurface->getAvgZ() != 0.00)
		(*grads).calc_GradientGrid(*xSurf, *ySurf, *zSurf);
	
	//Compute Depths, Gradients, and uncertainties
	this->scaleFactor = sH;
	this->alpha = alpha;
	this->deltaMin = deltaMin;
	vector<Point> pVec(tin->determinePingLocationInTriangle(*ptlSurface, sH, alpha, deltaMin,*grads, triangles, depthInterpMethod, errorInterpMethod, extrapMethod));

	e.reserve(pVec.size());
	e2.reserve(pVec.size());
	e3.reserve(pVec.size());
	e4.reserve(pVec.size());
	e5.reserve(pVec.size());
	z.reserve(pVec.size());
	s.reserve(pVec.size());

	nei.reserve(pVec.size());
	rei.reserve(pVec.size());
	z0.reserve(pVec.size());
	e0.reserve(pVec.size());
	eK.reserve(pVec.size());
	eK.reserve(pVec.size());

	for (int i = 0; i < (int) pVec.size(); i++)
	{
		e.push_back(pVec[i].u);
		e2.push_back(pVec[i].u2);
		e3.push_back(pVec[i].u3);
		e4.push_back(pVec[i].u4);
		e5.push_back(pVec[i].u5);
		z.push_back(pVec[i].z);
		s.push_back(pVec[i].s);

	//	e.push_back(pVec[i].e);
		nei.push_back(pVec[i].nei);
		rei.push_back(pVec[i].rei);
		z0.push_back(pVec[i].z0);
		e0.push_back(pVec[i].e0);
		zK.push_back(pVec[i].zK);
		eK.push_back(pVec[i].eK);
		
	}
}
//void InterpGrid::estimate(vector<double> *xSurf, vector<double> *ySurf, vector<double> *zSurf,vector<double> *eSurf,vector<double> *hSurf,vector<double> *vSurf, vector<double> *nmseiSurf,vector<double> *reiSurf,const double& sH, const double& alpha, const double& deltaMin, SHullDelaunay* tin, string interpMethod)
//{
//	ptlSurface = new PointList();
//	ptlSurface->setFromVectors(*xSurf, *ySurf, *zSurf, tin->getOffsetX(), tin->getOffsetY());
//
//	grads = new GradientGrid();
//
//	//triangles = vector<Triangle>();
//	triangles.reserve(ptlSurface->size());
//
//	//Compute Gradients if depths are known.
//	if(ptlSurface->getAvgZ() != 0.00)
//		(*grads).calc_GradientGrid(*xSurf, *ySurf, *zSurf);
//	//Compute Depths and Gradients
//	else
//	{
//		this->scaleFactor = sH;
//		this->alpha = alpha;
//		this->deltaMin = deltaMin;
//		vector<Point> pVec(tin->determinePingLocationInTriangle(*ptlSurface, sH, alpha, deltaMin,*grads, triangles, interpMethod));
//
//		e.reserve(pVec.size());
//		e2.reserve(pVec.size());
//		e3.reserve(pVec.size());
//		e4.reserve(pVec.size());
//		e5.reserve(pVec.size());
//		z.reserve(pVec.size());
//		s.reserve(pVec.size());
//
//		for (int i = 0; i < (int) pVec.size(); i++)
//		{
//			e.push_back(pVec[i].u);
//			e2.push_back(pVec[i].u2);
//			e3.push_back(pVec[i].u3);
//			e4.push_back(pVec[i].u4);
//			e5.push_back(pVec[i].u5);
//			z.push_back(pVec[i].z);
//			s.push_back(pVec[i].s);
//		}
//	}
//}
bool Bathy_Grid::ensemble()
{
	//This assumes that the correct grid is MBZ
	//and the extra leading rows in GMT need to 
	//be ignored for calculations.  These rows
	//are in the beginning and provide extra lats for 1 extra lon.
	//**This is believed to be fixed;verify before removing! SJZ 12/28/15.
	int i, k, j, diff = 0;
	double temp = 0;//,n,h;
	InterpGrid *kInd;
	int cnt = (const int)interpGrids.size();//3;// {MBZ,GMT,ALG} //
	double sigma = 0, sigma2 = 0, sigma3 = 0, sigma4 = 0, sigma5 = 0;
	vector<int> row, js, jsInit;
	vector<int> xcols, yrows;
	row.reserve(cnt);
	xcols.reserve(cnt);	//# rows (x) in each grid
	yrows.reserve(cnt);	//#cols (y) in each grid
	js.reserve(cnt);
	jsInit.reserve(cnt);	//last index pos at in the set of uniq. new starting index
	vector<InterpGrid*>::iterator it = interpGrids.begin();
	vector<InterpGrid*>::iterator kIt;
	k = (const int)(*it)->e.size();

	vector<double> sortGrid1;
	vector<double> sortGrid2;
	vector<double> sortGrid3;

	//find the max grid in order to know the number of extra rows to skip in GMT beginning
	for(it = interpGrids.begin()++; it != interpGrids.end(); it++)
	{
		j = (const int)(*it)->e.size();
		if(j < k)
		{
			k = j;
			diff = k - j;
			kInd = *&*it;
			kIt = it;
		}
		//xcols and yrows are indices starting at 0.  Add +1 to get dimensions.
		xcols.push_back((const int)(((*it)->ptlSurface->getPositions().back().x
			- (*it)->ptlSurface->getPositions().front().x)/(*it)->deltaMin)+1);
		yrows.push_back((const int)(abs(((*it)->ptlSurface->getPositions().back().y
			- (*it)->ptlSurface->getPositions().front().y)/(*it)->deltaMin))+1);
		js.push_back((const int)(*it)->e.size()-1);
		jsInit.push_back((const int)(*it)->e.size()-1);//beginning row index
	}
	vector<Point> mbzItXY;
	vector<Point> algItXY;
	vector<Point> gmtItXY;
	int mbzflag = 0;
	int masterGrid = 0;;
	for(it = interpGrids.begin(); it != interpGrids.end(); it++)
	{
		switch((*it)->gtype)
		{
			case ALGSpline:
				if(!mbzflag)
					masterGrid = 1; //because there is only GMT then ALG
				algItXY = (*it)->ptlSurface->getPositions();
				break;
			case GMT:
				gmtItXY = (*it)->ptlSurface->getPositions();
				break;//sorted by y min to max
			case MBZ:
				mbzflag = 1;
				masterGrid = MBZ; //should be 0
				mbzItXY = (*it)->ptlSurface->getPositions();
				break;
		}
	}

	it = interpGrids.begin();
	ensemble_H.reserve(k);
	ensemble_V.reserve(k);
	ensemble_U.reserve(k);
	ensemble_U2.reserve(k);
	ensemble_U3.reserve(k);
	ensemble_U4.reserve(k);
	ensemble_U5.reserve(k);
	vector<int> js2;
	int j2;
	js2.reserve(cnt);
	int counter = 0;
	int counter1 = 0, counter2 = 0;
	int counter3 = 0, counter4 = 0; 
	js[0] = 0; //last visited index
	jsInit[0] = 0; //index at beginning of current row visited
	js2.assign(js.begin(),js.end());
	//i = MBZ, k = GMT, j = same x,y loc in both
	for(i =0; i < k; i++)
	{
		for(it = interpGrids.begin(); it != interpGrids.end(); it++)
		{
			switch((*it)->gtype)
			{
				case ALGSpline:
					//ALGSpline j should be the same x,y as MBZ
					//counter = 2;		
					j = xcols[counter]*(counter3)+counter4;
					js[counter] = j;
					counter3++; 
					if(counter3==abs(yrows[masterGrid]))
					{
						counter3=0;
						counter4++;
					}
					if(counter4==xcols[masterGrid])
						counter4=0;

					break; //sorted by y max to min
				case GMT:
					//GMT j should be the same x,y as MBZ
					//counter = 1;				
					j = xcols[counter]*(abs(yrows[counter])-counter1)-(xcols[counter]-counter2);
					js[counter] = j;
					counter1++; 
					if(counter1==abs(yrows[masterGrid]))
					{
						counter1=0;
						counter2++;
					}
					if(counter2==xcols[masterGrid])
						counter2=0;

					break;//sorted by y min to max
				case MBZ:
					//j index is always i since MBZ is assumed to be correct.
					//counter = 0;
					j = i;
					js[counter] = j;
					break;//sorted by x min to max
			}
			if(mbzflag)
			{
				if(counter==ALGSpline){//2
					if(! (mbzItXY[js[masterGrid]].x == algItXY[js[ALGSpline]].x && mbzItXY[js[masterGrid]].y == algItXY[js[ALGSpline]].y))
					{
						std::cout<< "Error: ALG and MBZ X, Y locations do not align!" << std::endl;
						return false;
					}
				}
				if(counter==GMT){//1
					if(! (mbzItXY[js[masterGrid]].x == gmtItXY[js[GMT]].x && mbzItXY[js[masterGrid]].y == gmtItXY[js[GMT]].y))
					{
						std::cout<< "Error: GMT and MBZ X, Y locations do not align!" << std::endl;
						return false;
					}
				}
			}
			else
			{
				if(counter==masterGrid){
					if(! (algItXY[js[masterGrid]].x == gmtItXY[js[counter-1]].x && algItXY[js[masterGrid]].y == gmtItXY[js[counter-1]].y))
					{
						std::cout<< "Error: ALG and GMT X, Y locations do not align!" << std::endl;
						return false;
					}
				}
			}
			sigma += (*it)->e[j];
			sigma2 += (*it)->e2[j];
			sigma3 += (*it)->e3[j];
			sigma4 += (*it)->e4[j];
			sigma5 += (*it)->e5[j];
			counter++;
		}
		if(mbzflag)
		{
			ensemble_X.push_back(mbzItXY[js[masterGrid]].x);
			ensemble_Y.push_back(mbzItXY[js[masterGrid]].y);
			ensemble_Z.push_back(mbzItXY[js[masterGrid]].z);
		}
		else
		{
			ensemble_X.push_back(algItXY[js[masterGrid]].x);
			ensemble_Y.push_back(algItXY[js[masterGrid]].y);
			ensemble_Z.push_back(algItXY[js[masterGrid]].z);
		}
		ensemble_U.push_back(sigma/cnt);
		ensemble_U2.push_back(sigma2/cnt);
		ensemble_U3.push_back(sigma3/cnt);
		ensemble_U4.push_back(sigma4/cnt);
		ensemble_U5.push_back(sigma5/cnt);
		
		//SJZ
		temp = sqrt(pow(ensemble_U[i],2)/2);
		ensemble_H.push_back(temp);
		ensemble_V.push_back(temp);
		
		sigma = 0, sigma2 = 0, sigma3 = 0, sigma4 = 0, sigma5 = 0; counter = 0;
	}
	//printensemble("../Output_Files/ensemble_Test.txt");
	return true;
}

void Bathy_Grid::printensemble(string z_OutputFileName)
{
	std::ofstream outFile;
	//vector<InterpGrid*>::iterator it = interpGrids.begin();
	//vector<Point> ps;
	outFile.open(z_OutputFileName.c_str());
	outFile.precision(6);
	outFile.setf(std::ios::fixed, std::ios::floatfield);
	std::cout << "Writing file to " << z_OutputFileName.c_str() << std::endl;
	if (outFile.is_open())
	{
		//it = interpGrids.begin();
		//ps=(*it)->ptlSurface->getPositions();
		for(int i = 0; i < (int) ensemble_X.size(); i++){
			outFile << ensemble_X[i] << "\t" << ensemble_Y[i] << "\t" << ensemble_Z[i]
					<< "\t" << ensemble_U[i]  << "\t" << ensemble_U2[i]
					<< "\t" << ensemble_U3[i]<< "\t" << ensemble_U4[i]
					<< "\t" << ensemble_U5[i] << std::endl;
		}
	}
	outFile.close();
}

void Bathy_Grid::addToList( InterpGrid* g)
{
	interpGrids.push_back(&*g);
}

//Get the first interpGrid of type g
InterpGrid* Bathy_Grid::getGrid(GridType g)
{
	vector<InterpGrid*>::iterator it;
	for(it = interpGrids.begin(); it != interpGrids.end(); it++)
		if((*it)->gtype == g) return *it;
}

//Get all interpGrids
vector<InterpGrid*> Bathy_Grid::getGrids()
{
	return interpGrids;
}