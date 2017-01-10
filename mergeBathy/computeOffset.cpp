#include <exception>
#include "computeOffset.h"
#include "bathyTool.h"
/*computeOffset - Requires more than 1 data set in order to find the offsets between them. Constructs a grid x,y and computes the grid values for each data set by calling bathyTool.  BathyTool executes the heart of mergeBathy computations for generating a bathymetric surface. The n resulting gridded data sets from the n input data sets are compared to find the differences in their estimations. The computed offsets are removed from the input data before returning to the calling function.  This is a correction to the input data.
*/

void computeOffset(vector<TRUE_DATA> *inputData, double x0, double y0, double x1, double y1, vector<double> *xInterpVector, vector<double> *yInterpVector, double gridSpacingX, double gridSpacingY, vector<double> *zDataIn, vector<double> *eDataIn, vector<double> *hEDataIn, vector<double> *vEDataIn, string &kernel,  map<string, int> addOpts)
{
	//************************************************************************************
	//0. Declare and initialize local variables and objects for use with computeOffset.
	//************************************************************************************
	#pragma region	--Declare and Initialize local vars
	int inDataSize = (const int)(*inputData).size();

	vector<double> cOffset		= vector<double>(inDataSize);
	vector<double> cOffsetError = vector<double>(inDataSize);

	vector<double> xt = vector<double>();
	vector<double> yt = vector<double>();

	double minXi = (*xInterpVector)[0];
	double maxXi = (*xInterpVector)[0];
	double minYi = (*yInterpVector)[0];
	double maxYi = (*yInterpVector)[0];

	double localSpacingX = gridSpacingX / 2.00;
	double localSpacingY = gridSpacingY / 2.00;
	double val;

	#pragma region --Find Min and Max of Xi and Yi
	for (int i = 0; i < (const int)(*xInterpVector).size(); i++)
	{
		if ((*xInterpVector)[i] < minXi)
			minXi = (*xInterpVector)[i];
		else if ((*xInterpVector)[i] > maxXi)
			maxXi = (*xInterpVector)[i];

		if ((*yInterpVector)[i] < minYi)
			minYi = (*yInterpVector)[i];
		else if ((*yInterpVector)[i] > maxYi)
			maxYi = (*yInterpVector)[i];
	}

	double nx = ceil((maxXi-minXi) / localSpacingX);
	double ny = ceil((maxYi-minYi) / localSpacingY);
	#pragma endregion

	//A. Recompute the spacing
	localSpacingX = (maxXi-minXi) / nx;
	localSpacingY = (maxYi-minYi) / ny;
	#pragma endregion

	//************************************************************************************
	//I. Start computing grid locations and spacings for computing offset.
	//************************************************************************************
	#pragma region --Compute Grid locations and spacings
	double meanXt = 0.00;
	double meanYt = 0.00;
	//A. Calculate X Grid Interpolation Locations
	double locationValue = minXi;
	int xtSize = 0;
	while (locationValue <= maxXi) 
	{
		xt.push_back(locationValue);
		meanXt += locationValue;
		locationValue += localSpacingX;
		xtSize += 1;
	}
	//Get the mean for later use
	meanXt = meanXt / (double)xt.size();

	//B. Calculate Y Grid Interpolation Locations
	locationValue = minYi;
	int ytSize = 0;
	while (locationValue <= maxYi)
	{
		yt.push_back(locationValue);
		meanYt += locationValue;
		locationValue += localSpacingY;
		ytSize += 1;
	}
	//Get the mean for later use
	meanYt = meanYt / (double)yt.size();

	//C. set up the vectors and grids
	int xt_yt_Size = xtSize*ytSize;
	vector<double> xMeshVector2 = vector<double>(xt_yt_Size);
	vector<double> yMeshVector2 = vector<double>(xt_yt_Size);
	dgrid xMeshGrid2 = dgrid(ytSize, xtSize);
	dgrid yMeshGrid2 = dgrid(ytSize, xtSize);

	//D. Reform the data to a grid that we can use for calculations
	int clx = 0;
	int cly = 0;
	int currentLoc = 0;
	for (int i = 0; i < ytSize; i++)
	{
		currentLoc = i;
		for (int j = 0; j < xtSize; j++)
		{
			xMeshGrid2(i,j) = xt[clx] - meanXt;
			yMeshGrid2(i,j) = yt[cly] - meanYt;
			xMeshVector2[currentLoc] = xt[clx];
			yMeshVector2[currentLoc] = yt[cly];
			currentLoc = currentLoc + ytSize;
			clx += 1;
		}
		clx = 0;
		cly += 1;
	}

	OUTPUT_DATA xyzOut;
	boost::numeric::ublas::mapped_matrix<double> zData (inDataSize, xMeshVector2.size());
	boost::numeric::ublas::mapped_matrix<double> wData (inDataSize, xMeshVector2.size());

	dgrid xMeshTemp;
	dgrid yMeshTemp;
	vector<double> xtTemp;
	vector<double> ytTemp;
	double minCurXLoc;
	double minCurYLoc;

	//F. Make the local data structure for overlapping data
	vector<TRUE_DATA> localInputData = vector<TRUE_DATA>(inDataSize);
	for (int i = 0; i < inDataSize; i++)
	{
		for (int j = 0; j < (const int)(*inputData)[i].x.size(); j++)
		{
			if( ((*inputData)[i].x[j] > (minXi - gridSpacingX))
				&& ((*inputData)[i].x[j] < (maxXi + gridSpacingX)) )
			{
				if( ((*inputData)[i].y[j] > minYi - gridSpacingY)
					&& ((*inputData)[i].y[j] < (maxYi + gridSpacingY)) )
				{
					localInputData[i].x.push_back((*inputData)[i].x[j]);
					localInputData[i].y.push_back((*inputData)[i].y[j]);
					localInputData[i].depth.push_back((*inputData)[i].depth[j]);
					localInputData[i].error.push_back((*inputData)[i].error[j]);
					localInputData[i].v_Error.push_back((*inputData)[i].v_Error[j]);
					localInputData[i].h_Error.push_back((*inputData)[i].h_Error[j]);
				}
			}
		}
	}
	#pragma endregion

	//************************************************************************************
	//II. Run each input data set through bathy tool to compute the offset at that data set
	//************************************************************************************
	#pragma region --Compute Offset, Run bathyTool for each dataset
	int ii = 0; // SJZ to catch weird behavior
	for (int i = 0; i < inDataSize; i++)
	{
		ii = i+1;
		cout << "\tComputing Data Set Number: " << ii << " ...";
		minCurXLoc = localInputData[i].x[0];
		minCurYLoc = localInputData[i].y[0];
		//A. Get the min
		for (int j = 0; j < (const int)localInputData[i].x.size(); j++)
		{
			if(localInputData[i].x[j] < minCurXLoc)
				minCurXLoc = localInputData[i].x[j];
			if(localInputData[i].y[j] < minCurYLoc)
				minCurYLoc = localInputData[i].y[j];
		}
		xMeshTemp = dgrid(xMeshGrid2);
		yMeshTemp = dgrid(yMeshGrid2);
		xtTemp = vector<double>(xt);
		ytTemp = vector<double>(yt);

		//B. Call bathyTool
		bathyTool(&localInputData[i].x,&localInputData[i].y, &localInputData[i].depth, &localInputData[i].error, &localInputData[i].h_Error, &localInputData[i].v_Error, &xMeshTemp, &yMeshTemp, &xtTemp, &ytTemp, gridSpacingX, gridSpacingY, minCurXLoc, minCurYLoc, meanXt, meanYt, kernel, addOpts, 0.5, false, false, NEITOL_COMPUTE_OFFSET, &xyzOut);

		//Added 8/7/14 sqrt results to get errors
		for (int j = 0; j < (const int)xyzOut.nEi.size(); j++)
		{
			val = 1.00 - sqrt(xyzOut.nEi[j]);
			if (val > eps || val < -eps)
			{
				if(xyzOut.depth[j] > eps || xyzOut.depth[j] < -eps)
					zData(i,j) = xyzOut.depth[j];
				wData(i,j) = val;
			}
		}

		//C. Clear the variables
		xyzOut.depth.clear();
		xyzOut.error.clear();
		xyzOut.nEi.clear();
		xyzOut.rEi.clear();

		if(abs(addOpts.find("-propUncert")->second)==1)
		{
			xyzOut.depth0.clear();
			xyzOut.error0.clear();
		}
		if(abs(addOpts.find("-kalman")->second)==1)
		{
			xyzOut.depthK.clear();
			xyzOut.errorK.clear();
		}
		xMeshTemp.clear();
		yMeshTemp.clear();
		xtTemp.clear();
		ytTemp.clear();
		localInputData[i].x.clear();
		localInputData[i].y.clear();
		localInputData[i].depth.clear();
		localInputData[i].error.clear();
		localInputData[i].h_Error.clear();
		localInputData[i].v_Error.clear();
		cout << "... Done!" << endl;
	}

	localInputData.clear();
	#pragma endregion

	//************************************************************************************
	//III. Do some offsetting calculations based on the values returned by bathyTool.
	// Form the matrix equation: DZ(k,k') = [Delta(k) - Delta(k')]dz(k)=>D = [R]*dz
	// These are observed differences between data sets.
	//************************************************************************************
	#pragma region Offsetting Calculations

	int NR = (int)floor(xtSize * ytSize * (pow(inDataSize, 2.00)) / 2.00);
	vector<int> idVector;

	int curLoc = 0;
	double sumW = 0.00;

	//Create empty sparse matrices big enough to hold ...
	//all pairs
	boost::numeric::ublas::mapped_matrix<double> dzDData(1, NR, 0);
	//the correlation matrix
	boost::numeric::ublas::mapped_matrix<double> dzRData(inDataSize, NR, 0);
	//the joint weight, sqrt(wt(k)*wt(k'))
	boost::numeric::ublas::mapped_matrix<double> dzWData(1, NR, 0);

	//A. Form the matrix equation: DZ(k,k') = [Delta(k) - Delta(k')]dz(k)=>D = [R]*dz
	//Here are observed differences between data sets
	if(!idVector.empty())
		idVector.clear();

	//loop through all offset pairs in order to find the points that contribute to the offset.
	for (int i = 0; i < inDataSize; i++){//k
		for (int j = i + 1; j < inDataSize; j++){//kprime
			//get overlap
			for (int k = 0; k < (const int)wData.size2(); k++){//id
				//these are points that contribute to offset
				if(((wData(i,k) > eps || wData(i,k) < -eps) && (wData(j,k) > eps || wData(j,k) < -eps ))){
					dzDData(0, curLoc) = zData(i,k)-zData(j,k);
					dzRData(i, curLoc) = 1.0000;
					dzRData(j, curLoc) = -1.00 * dzRData(i,curLoc);
					dzWData(0, curLoc) = sqrt(wData(i,k)*wData(j,k));
					curLoc = curLoc + 1;
				}
			}
		}
		printf("Finished weighting %d of %d sets\n", i, inDataSize);
	}

	//PAE - 12/31/2009. If there is no overlap, which is possible, leave the subroutine.
	//MATLAB and CPP versions may not agree!!!
	//This check occurs later in MATLAB version, but we need to do this here to avoid run-time crash.
	//MATLAB version uses NaN's if no overlap, then exits later without making corrections.
	if(curLoc == 0)
		return ;

	//Create empty sparse matrices
	boost::numeric::ublas::mapped_matrix<double>  dzDData2 (1, curLoc, 0);
	boost::numeric::ublas::mapped_matrix<double>  dzRData2 (inDataSize, curLoc, 0);
	boost::numeric::ublas::mapped_matrix<double>  dzWData2 (1, curLoc, 0);

	//convert to weighted space (priestly p.315)
	sumW = 0;
	for (int i = 0; i < curLoc; i++)
	{
		dzDData2(0,i) = dzDData(0,i)*dzWData(0,i);
		dzWData2(0,i) = dzWData(0,i);
		for (int j = 0; j < inDataSize; j++)
		{
			dzRData2(j,i) = dzRData(j,i) * (dzWData2(0,i) + (eps));
		}
		sumW += dzWData2(0,i);
	}

	//Save as compressed sparse matrix for efficient computations
	boost::numeric::ublas::compressed_matrix<double> dzRData2_C (dzRData2);
	boost::numeric::ublas::compressed_matrix<double> dzDData2_C (dzDData2);
	boost::numeric::ublas::compressed_matrix<double> dzWData2_C (dzWData2);

	//model-model correlation (mult. by weights here)
	boost::numeric::ublas::compressed_matrix<double> rtRGrid = prod(dzRData2_C, trans(dzRData2_C))/sumW;
	//model-data correlation
	boost::numeric::ublas::compressed_matrix<double> rtDGrid = prod(dzDData2_C, trans(dzRData2_C))/sumW;

	//B. Get the values that are not zero
	//It is possible that there is no overlap for some data sets, remove
	for (int i = 0; i < (const int)rtRGrid.size2(); i++)
	{
		if (rtRGrid(i,i) > 0+eps)
			idVector.push_back(i);
	}

	if(idVector.size() == 0)
	{
		cerr << "No overlap found." << endl;
		return ;
	}

	int newFileSize = (const int)idVector.size();
	boost::numeric::ublas::mapped_matrix<double> temp  (rtRGrid);
	boost::numeric::ublas::mapped_matrix<double> temp2 (rtDGrid);

	rtRGrid.clear();
	rtDGrid.clear();
	//rtRGrid.resize (idVector.size() + 1, idVector.size() + 1); // produced !preserve error in UNIX // SJZ 12/10/14
	//rtDGrid.resize (1, idVector.size() + 1);
	rtRGrid = boost::numeric::ublas::compressed_matrix<double>(idVector.size() + 1, idVector.size() + 1);
	rtDGrid = boost::numeric::ublas::compressed_matrix<double>(1, idVector.size() + 1);

	//C. Do some matrix setting based on the values that were found in the idVector
	//Augment with Lagrange Mult., B
	int B = 0; // sum of offsets equals this
	for (int i = 0; i < (const int)idVector.size(); i++)
	{
		for (int j = 0; j < (const int)idVector.size(); j++)
		{
			rtRGrid(i,j) = temp(idVector[i], idVector[j]);
			if( i == idVector.size()-1)
				rtRGrid(i+1,j) = 1; //append row of ones
		}
		rtRGrid(i,idVector.size()) = 1; //append column of ones
		rtDGrid(0,i) = temp2(0, idVector[i]); //Keep id datasets
	}
	rtRGrid.erase_element(rtRGrid.size1()-1, rtRGrid.size2()-1);//= 0;
	rtDGrid(0,idVector.size()) = B; // append B

	int count = 0;
	boost::numeric::ublas::compressed_matrix<double> rtRInv;
	boost::numeric::ublas::compressed_matrix<double> dzSolution;
	boost::numeric::ublas::compressed_matrix<double> rtRTemp;

	//D. Invert the matrix
	while((dzSolution.nnz() == 0) && (count <= 10))
	{
		double p = 0;
		if (count > 0)
		{
			p = eps*pow(2.00,count);
		}
		rtRTemp = rtRGrid;
		for (int i = 0; i < newFileSize; i++)
			for (int j = 0; j < newFileSize; j++)
				rtRTemp(i,j) = rtRTemp(i,j) + p;

		rtRInv = boost::numeric::ublas::compressed_matrix<double> (rtRTemp.size1(), rtRTemp.size2());		
		//rtRInv.resize(rtRTemp.size1(), rtRTemp.size2()); // produced !preserve error in UNIX // SJZ 12/10/14
		InvertMatrix(rtRTemp, rtRInv);

		dzSolution = prod(rtRInv, trans(rtDGrid));
		count += 1;

		rtRInv.clear();
	}

	vector<int>::iterator it = idVector.begin();
	//E. Compute the offset
	double meanOffset = 0;
	for (int i = 0; i < newFileSize; i++)
	{
		if(abs(dzSolution(i,0)) > eps)
		{
			cOffset[*it] = dzSolution(i,0);
			meanOffset += cOffset[*it];
		}
		it++;
	}
	meanOffset /= (double)cOffset.size();

	boost::numeric::ublas::compressed_matrix<double> msz = (prod(dzDData2_C, trans(dzDData2_C))/sumW);
	boost::numeric::ublas::compressed_matrix<double> msrTemp = prod(trans(dzSolution), rtRGrid);
	boost::numeric::ublas::compressed_matrix<double> msrTemp2 = prod(msrTemp, dzSolution);
	boost::numeric::ublas::compressed_matrix<double> msr = msz - msrTemp2;

	//F. Compute the offset error
	double offErrTmp;
	double meanOffsetError = 0;
	it = idVector.begin();
	int m=0;
	for (int i = 0; i < newFileSize; i++){
		//default error
		while(m<*it){
			cOffsetError[m] = cOffsetError[m]+ sqrt(msz(0,0));
			meanOffsetError += cOffsetError[m];
			m++;
		}
		offErrTmp = (rtRInv(i,i)*msr(0,0));
		if ((sumW-2)>0){
			offErrTmp = offErrTmp/(sumW-2.00);
			if (offErrTmp >= 0)
				cOffsetError[*it] = sqrt(offErrTmp);
		}
		meanOffsetError += cOffsetError[*it];//this is not correct
		it++;m++;
	}
	meanOffsetError /= (double)cOffsetError.size();
	#pragma endregion

	//************************************************************************************
	//IV. Apply the offsets and clear the local variables
	//************************************************************************************
	#pragma region Apply Offset
	double pinDepthValue = 0;
	curLoc = 0;
	for (int i = 0; i < inDataSize; i++)
	{
		cOffset[i] = cOffset[i] - meanOffset;
//		cOffsetError[i] = cOffsetError[i] - meanOffsetError;

		if (abs(cOffset[i]) > (*inputData)[i].maximumDataOffset)
		{
			curLoc += 1;
			pinDepthValue += cOffset[i];
		}
	}
	if (curLoc != 0)
		pinDepthValue /= (double)curLoc;

	curLoc = 0;
	for (int i = 0; i < inDataSize; i++)
	{
		cOffset[i] -= pinDepthValue;
		cout << "\tApplying a " << cOffset[i] << " meter offset with " << cOffsetError[i] << " meters of combined uncertainty to data set " << i+1 << " ...";

		for (int j = 0; j < (const int)(*inputData)[i].depth.size(); j++)
		{
			(*zDataIn)[curLoc] = (*zDataIn)[curLoc] - cOffset[i];
			//(*eDataIn)[curLoc] = (*eDataIn)[curLoc] - cOffsetError[i];
			//(*hEDataIn)[curLoc] = (*hEDataIn)[curLoc] - cOffsetError[i];
			//(*vEDataIn)[curLoc] = (*vEDataIn)[curLoc] - cOffsetError[i];
			curLoc += 1;
		}
		cout << "... Done!" << endl;
	}
	#pragma endregion

	//msz.clear();
	//msr.clear();
	//temp.clear();
	//temp2.clear();
	//idVector.clear();
	////rtRInv.clear();
	//dzSolution.clear();
	//rtRTemp.clear();
	//rtRGrid.clear();
	//rtDGrid.clear();
	//dzDData.clear();
	//dzWData.clear();
	//dzRData.clear();
	//dzDData2.clear();
	//dzWData2.clear();
	//dzRData2.clear();
	//cOffset.clear();
	//cOffsetError.clear();
	//xt.clear();
	//yt.clear();
	//xMeshGrid2.clear();
	//yMeshGrid2.clear();
}

/* Matrix inversion routine.
Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
 template<class T>
bool InvertMatrix (const boost::numeric::ublas::compressed_matrix<T>& input, boost::numeric::ublas::compressed_matrix<T>& inverse) {
	using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;
	// create a working copy of the input
	matrix<T> A(input);
	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = (const int)lu_factorize(A,pm);
	if( res != 0 ) return false;

	if(determinant(pm)==0) // SJZ
	{
		cout << "ERROR: Determinant is 0!" << endl;
		cerr << "ERROR: Stop." << endl;
	}

	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<T>(A.size1()));

	// back substitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}

//This can be improved upon...There is a better one on the internets 
//but I didn't have time to implement it so this one was good for now. // SJZ
//This had not been checked! I only plugged the function in. 
 template <typename size_type, typename A>
 int determinant(const boost::numeric::ublas::permutation_matrix<size_type,A>& pm)
 {
     int pm_sign=1;
     size_type size=pm.size();
     for(size_type i = 0; i < size; ++i)
         if(i != pm(i))
             pm_sign = pm_sign*-1; // swap_rows would swap a pair of rows here, so we change sign
     return pm_sign;
 }


