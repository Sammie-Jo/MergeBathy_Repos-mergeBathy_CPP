#include "subSampleData.h"
#include "consistentWeights.h"
#include <math.h>

//************************************************************************************
//SUBROUTINE I: Partition2: This subroutine partitions the grid data for
//sorting and maintains each element's initial index when it is moved.  Used to sort the subsampled data.
//************************************************************************************
int Partition2( vector<double> *sX, vector<double> *sY, vector<double> *sZ, vector<double> *sE, vector<double> *sE2, int left, int right)
{
	double val1 = (*sX)[left];
	double val2 = (*sY)[left];

	double tempX;
	double tempY;
	double tempZ;
	double tempE;
	double tempE2;

	int lm = left - 1;
	int rm = right + 1;

	for(;;) {
		do{
			rm--;
		}while (((*sX)[rm] > val1));// && ((*sY)[rm] > val2));

		do {
			lm++;
		}while( ((*sX)[lm] < val1));// && ((*sY)[lm] < val2));

		if(lm < rm) {
			tempX = (*sX)[rm];
			tempY = (*sY)[rm];
			tempZ = (*sZ)[rm];
			tempE = (*sE)[rm];
			tempE2 = (*sE2)[rm];

			(*sX)[rm] = (*sX)[lm];
			(*sY)[rm] = (*sY)[lm];
			(*sZ)[rm] = (*sZ)[lm];
			(*sE)[rm] = (*sE)[lm];
			(*sE2)[rm] = (*sE2)[lm];

			(*sX)[lm] = tempX;
			(*sY)[lm] = tempY;
			(*sZ)[lm] = tempZ;
			(*sE)[lm] = tempE;
			(*sE2)[lm] = tempE2;
		}else
			return rm;
	}
}

//************************************************************************************
//SUBROUTINE II: SelectionSort2: This subroutine is a simple selection sort that
//is used when the grid data gets small enough.  It maintains each element's
//initial index when it is moved. Used to sort the subsampled data.
//************************************************************************************
void SelectionSort2(vector<double> *sX, vector<double> *sY, vector<double> *sZ, vector<double> *sE, vector<double> *sE2, int left, int right)
{
	double tempX;
	double tempY;
	double tempZ;
	double tempE;
	double tempE2;
	for(int i = left; i < right; i++) {
		int min = i;
		for(int j=i+1; j <= right; j++){
			if (((*sX)[j] < (*sX)[min]))// && ((*sY)[j] < (*sY)[min]))
			{
				min = j;
			}
		}
		tempX = (*sX)[min];
		tempY = (*sY)[min];
		tempZ = (*sZ)[min];
		tempE = (*sE)[min];
		tempE2 = (*sE2)[min];

		(*sX)[min] = (*sX)[i];
		(*sY)[min] = (*sY)[i];
		(*sZ)[min] = (*sZ)[i];
		(*sE)[min] = (*sE)[i];
		(*sE2)[min] = (*sE2)[i];

		(*sX)[i] = tempX;
		(*sY)[i] = tempY;
		(*sZ)[i] = tempZ;
		(*sE)[i] = tempE;
		(*sE2)[i] = tempE2;
	}
}

//************************************************************************************
//SUBROUTINE III: Quicksort2: This function is a simple quick sort function
//to sort data so it can be parsed. Used to sort the subsampled data.
//************************************************************************************
void Quicksort2( vector<double> *sX, vector<double> *sY, vector<double> *sZ, vector<double> *sE, vector<double> *sE2, int left, int right)
{
	if(left < (right-7))
	{
		int split_pt = Partition2(sX, sY, sZ, sE, sE2, left, right);
		Quicksort2(sX, sY, sZ, sE, sE2, left, split_pt);
		Quicksort2(sX, sY, sZ, sE, sE2, split_pt+1, right);
	}
	else SelectionSort2(sX, sY, sZ, sE, sE2, left, right);
}

//************************************************************************************
//SUBROUTINE IV: subsampleData:
//This subsamples and spaces the data points.
//************************************************************************************
int subsampleData(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, vector<double> *inputDataHErr, vector<double> *inputDataVErr, double gridSpacingX, double gridSpacingY, double &inX0, double &inY0, double meanX, double meanY, bool dispIntermResults, vector< vector<double> > *subsampledData)
{
	//************************************************************************************
	//0. Declare local variables and objects.
	//************************************************************************************
	//A. Initialize local variables
	int returnValue = 0;

	double wtol = 0.10000;
	//The array of number of cells in each dimension
	double jMaxX, jMaxY;
	//Scale parameters, indicating the step size in each dimension
	/*Matlab denotes the scale parameters as being an MxN2 array:
	N2 is either 1 (indicating constant smoothing length scales,
	or M, indicating variable length scales for each data point*/
	double DX = gridSpacingX;
	double DY = gridSpacingY;
	//The location of the first grid point (i.e., location at Ji=1)
	//Make nice integer values
	double x0 = (floor(inX0/DX))*(DX);
	double y0 = (floor(inY0/DY))*(DY);

	double s3Value;
	int jProd = 1;
	unsigned int i, j;
	unsigned int count;
	uint inputDataXSize = (uint)(*inputDataX).size();

	//The array of indices into each cell
	dgrid Ji(inputDataXSize,2,1);
	dgrid J(inputDataXSize,2,0);

	vector<int>::iterator iter;
	vector<double> weights = vector<double>(inputDataXSize, 2.00);

	if (dispIntermResults)
		printf("Subsampling Data: ");

	/*
	SAM notes:
	We need to find the ranges (or distances from 0)
	by subtracting min x and min y from x and y.

	X0 is the floor of min x and the floor of min y, in order to get an integer.
	We then subtract the mean from our values (X-repmat(X0,N,1)), scale by DX, and round.
	By scaling and rounding we will assign the same value in J(:,1) for all x
	that fall within the same range of x values. Likewise, we will assign the
	same value in J(:,2) for all y falling within the same range of y values.

	We now want to map these to index into a multi-D array of unique indices.
	Ji is our new multi-D array of unique indices we want to index into by mapping J.
	Multi-D array J is stored in an M*N x 1 vector Ji via Jprod to ensure uniqueness.
	We want a unique value representing all x AND y falling within the
	same x AND y range (a tile).  We need to make sure these values are
	unique so that tiles with inverse values for x range AND y range, such as
	(1,2) and (2,1), do not map to the same value.  This is done by multiplying 
	y range values(J(:,2)-1) by the max x range value Jmax and adding it to J(:,1).

	As DX grows, so does the number of observations that will be merged in to the same index
	i.e. the larger the DX,the more values that will be assigned the same index
	and values with the same index in Ji will be merged by calculating a single
	value for each index such that indices are unique
	*/

	//B. Initialize the Ji, and J dgrids with coordinate data
	//map data to scaled points
	//J = 1,1...,1 is location X0(1,1,...,1)
	//subtract the min to get the range and scale by Dx, then add 1 and round
	J(0,0) = round(1.00+(((*inputDataX)[0] - x0)/DX));
	J(0,1) = round(1.00+(((*inputDataY)[0] - y0)/DY));
	jMaxX = J(0,0);
	jMaxY = J(0,1);

	//C. Now loop to set the values and square the errors
	for (i = 0; i < inputDataXSize; i++){
		(*inputDataE)[i] = pow((*inputDataE)[i], 2);
		Ji(i,0) = 1;
		Ji(i,1) = i;
		J(i,0) = round(1.00+(((*inputDataX)[i] - x0)/DX));
		J(i,1) = round(1.00+(((*inputDataY)[i] - y0)/DY));

		if (J(i,0) > jMaxX)
			jMaxX = J(i,0);
		if (J(i,1) > jMaxY)
			jMaxY = J(i,1);
	}
	if (dispIntermResults)
		printf(".");

	//************************************************************************************
	//I. Begin subsampling the data
	//************************************************************************************
	//A. Compute the weights of each data point
	consistentWeights(inputDataZ, inputDataE, &wtol, &weights, &s3Value);

	//************************************************************************************
	//II. Fill Ji and get jProd
	//************************************************************************************
	//A. Map J to index of a multi-D array of unique indices. Jprod will store
	//rescaled versions of multi-D array J into a M*N x 1 vector Ji. We repeat
	//the multiplication process to ensure that different N-blocks of indices
	//have different numbers (i.e. no way to repeat an x and y index due to
	//resetting JProd on first go and then multiplication of y-indices by JProd
	//on second go. The whole point of this is to count the number of unique
	//DX cell distances.
	jProd = 1;

	for (i = 0; i < 2; i++){
		for (j = 0; j < (uint)inputDataXSize; j++){
			Ji(j,0) = Ji(j,0) + ((J(j,i) - 1) * (double)jProd);
		}
		if (i == 0)
			jProd = jProd * (int)jMaxX;
		else
			jProd = jProd * (int)jMaxY;
	}
	if (dispIntermResults)
		printf(".");

	//B. Sort data in ascending order and keep track of the indexes
	Quicksort(&Ji, 0, inputDataXSize - 1);
	
	if (dispIntermResults)
		printf("...");

	//************************************************************************************
	//III. Remove repeated points
	//************************************************************************************
	//A. Resize the vector and remove unique elements
	vector<int> JiVector = vector<int>(inputDataXSize,0);
	for (i = 0; i < inputDataXSize; i++){
		JiVector[i] = (int)Ji(i,0);
	}
	iter = unique(JiVector.begin(), JiVector.end());
	JiVector.resize(iter-JiVector.begin());

	if (dispIntermResults)
		printf(".");

	//B. Set up the computation vectors
	size_t	JiVectorSize = JiVector.size();
	zvector zi(JiVectorSize,0.0);
	zvector ni(JiVectorSize,0.0);		//The number of observations going into each cell
	zvector wi(JiVectorSize,0.0);
	zvector w2zi(JiVectorSize,0.0);		// weights against data
	zvector w2ei(JiVectorSize,0.0);		// weights against a priori errors
	zvector w2i(JiVectorSize,0.0);		// sum of squared weights
	zvector si(JiVectorSize,0.0);		// holds weighted sum of squares, initially
	dvector JiNew(JiVectorSize,0.0);	//Reset vector Ji for new use; The array of indices into each cell

	JiNew[0] = Ji(0,0);

	//C. Set up the output vectors
	//The mean position of the data in each cell
	(*subsampledData)[0] = vector<double>(JiVectorSize, 0.0); //x
	(*subsampledData)[1] = vector<double>(JiVectorSize, 0.0); //y
	//The mean value at each interp. cell
	(*subsampledData)[2] = vector<double>(JiVectorSize, 0.0); //z
	//The standard error (=std. dev./sqrt(n))
	//(or, if ni<3, insert average value of all other si values)
	//ONLY COMPUTED FOR z(:,1), others are weighted identically
	(*subsampledData)[3] = vector<double>(JiVectorSize, 0.0); //e
	(*subsampledData)[4] = vector<double>(JiVectorSize, 0.0); //w
	(*subsampledData)[5] = vector<double>(JiVectorSize, 0.0); //hu
	(*subsampledData)[6] = vector<double>(JiVectorSize, 0.0); //vu
	count = 0;
	int cur_position=0;
	double hval=0,vval=0;
	//************************************************************************************
	//IV. Subsample and apply weights
	//************************************************************************************
	//A. Loop through and set values
	for (i = 0; i < inputDataXSize; i++){
		if(Ji(i,0) > JiNew[count]){ //Jisort(i) > Ji(cnt)
			count = count + 1;
			JiNew[count] = Ji(i,0); //Ji(cnt) = Jisort(i)
		}

		cur_position = (int)Ji(i,1);
		vval=(*inputDataVErr)[(int)Ji(i,1)];
		hval=(*inputDataVErr)[(int)Ji(i,1)];

		// Ji(:,0) is Jisort, Ji(:,1) is sortid(i)
		ni[count] = ni[count].real() + 1;
		wi[count] =  wi[count].real() + weights[(int)Ji(i,1)];
		w2i[count] =  w2i[count].real() + pow(weights[(int)Ji(i,1)],2);
		zi[count] =  zi[count].real() + (*inputDataZ)[(int)Ji(i,1)] * weights[(int)Ji(i,1)];
		w2zi[count] =  w2zi[count].real() + (*inputDataZ)[(int)Ji(i,1)] * pow(weights[(int)Ji(i,1)],2);
		si[count] =  si[count].real() + pow((*inputDataZ)[(int)Ji(i,1)] * weights[(int)Ji(i,1)], 2);
		w2ei[count] = w2ei[count] + (*inputDataE)[(int)Ji(i,1)] * pow(weights[(int)Ji(i,1)],2);
		(*subsampledData)[0][count] = (*subsampledData)[0][count] + (*inputDataX)[(int)Ji(i,1)] * weights[(int)Ji(i,1)];
		(*subsampledData)[1][count] = (*subsampledData)[1][count] + (*inputDataY)[(int)Ji(i,1)] * weights[(int)Ji(i,1)];
		(*subsampledData)[5][count] = (*subsampledData)[5][count] + pow((*inputDataHErr)[(int)Ji(i,1)] * weights[(int)Ji(i,1)],2);
		(*subsampledData)[6][count] = (*subsampledData)[6][count] + pow((*inputDataVErr)[(int)Ji(i,1)] * weights[(int)Ji(i,1)],2);
	}
	if (dispIntermResults)
		printf("...");

	//B. Set the subdata values with the real data from the zgrid complex data
	for (i = 0; i < (uint)JiVectorSize; i++){
		// Handle division by 0.
		if(wi[i].real()==0){
			zi[i] = 0;
			(*subsampledData)[0][i] = 0;
			(*subsampledData)[1][i] = 0;
			(*subsampledData)[5][i] = 0;
			(*subsampledData)[6][i] = sqrt(DEFAULT_E);
		}
		else{
			zi[i] = zi[i].real() / wi[i].real();
			(*subsampledData)[0][i] = (*subsampledData)[0][i] / wi[i].real();
			(*subsampledData)[1][i] = (*subsampledData)[1][i] / wi[i].real();
			(*subsampledData)[5][i] = sqrt((*subsampledData)[5][i] / w2i[i].real());
			(*subsampledData)[6][i] = sqrt((*subsampledData)[6][i] / w2i[i].real());
		}

		si[i] = (si[i].real() - (2.00 * zi[i].real() * w2zi[i].real()) + (w2i[i].real() * pow(zi[i].real(),2)))/(ni[i].real());
		si[i] = (si[i].real())/(1.00 + ni[i].real());

		//Handle division by 0.
		if(w2i[i].real()==0)
			si[i] = ((ni[i].real() - 1.00) * si[i])/(ni[i].real());
		else
			si[i] = (((ni[i].real() - 1.00) * si[i]) + (w2ei[i].real()/w2i[i].real()))/(ni[i].real());

		si[i] = (std::complex<double>)sqrt(si[i]);

		(*subsampledData)[2][i] = (double)zi[i].real();
		(*subsampledData)[3][i] = (double)si[i].real();
		(*subsampledData)[4][i] = pow((double)si[i].real(),2);
	}

	//C. Do a second quicksort to get all of the X values in a row and maintain the rest of the data in the proper indexes
	//Quicksort2(&(*subsampledData)[0], &(*subsampledData)[1], &(*subsampledData)[2], &(*subsampledData)[3], &(*subsampledData)[4], 0, JiVectorSize);
	if (dispIntermResults)
		printf("..");

	inX0 = x0;
	inY0 = y0;
	
	//************************************************************************************
	//V. Clean up variables and return.
	//************************************************************************************
	Ji.clear();
	J.clear();
	weights.clear();
	JiVector.clear();
	zi.clear();
	ni.clear();
	wi.clear();
	w2zi.clear();
	w2ei.clear();
	w2i.clear();
	si.clear();
	JiNew.clear();

	if (dispIntermResults)
		printf(" Done Subsampling Data\n");

	return returnValue;
}

//************************************************************************************
//SUBROUTINE V: round: This function is used to properly round double variables
//************************************************************************************
double round(double x){
	if((ceil(x)-x) < .5){
		return ceil(x);
	}else{
		return floor(x);
	}
}

//************************************************************************************
//SUBROUTINE VI: Partition: This subroutine partitions the grid data for
//sorting and maintains each element's initial index when it is moved
//************************************************************************************
int Partition( dgrid *d, int left, int right)
{
	double val = (*d)(left,0);
	int lm = left - 1;
	int rm = right + 1;
	for(;;)
	{
		do{
			rm--;
		}while ((*d)(rm,0) > val);

		do {
			lm++;
		}while( (*d)(lm,0) < val);

		if(lm < rm)
		{
			double r = (*d)(rm,1);
			double tempr = (*d)(rm,0);
			(*d)(rm,0) = (*d)(lm,0);
			(*d)(lm,0) = tempr;

			(*d)(rm,1) = (*d)(lm,1);
			(*d)(lm,1) = r;
		}
		else
			return rm;
	}
}

//************************************************************************************
//SUBROUTINE VII: SelectionSort: This subroutine is a simple selection sort that
//is used when the grid data gets small enough.  It maintains each element's
//initial index when it is moved
//************************************************************************************
void SelectionSort(dgrid *data, int left, int right)
{
	double r;
	for(int i = left; i < right; i++) {
		int min = i;
		for(int j=i+1; j <= right; j++){
			if((*data)(j,0) < (*data)(min,0)) {
				min = j;
			}
		}

		double temp = (*data)(min,0);
		if(temp==144)
			int somethingw=0;

		double min0=(*data)(i,0);

		(*data)(min,0) = (*data)(i,0);
		(*data)(i,0) = temp;

		r = (*data)(min,1);

		double min1=(*data)(i,1);
		(*data)(min,1) = (*data)(i,1);
		(*data)(i,1) = r;
	}
}

//************************************************************************************
//SUBROUTINE VIII: Quicksort: This function is a simple quick sort function
//to sort data so it can be parsed.
//************************************************************************************
void Quicksort( dgrid *d, int left, int right)
{
	if(left < (right-7)) {
		int split_pt = Partition(d,left, right);
		Quicksort(d, left, split_pt);
		Quicksort(d, split_pt+1, right);
	}
	else SelectionSort(d, left, right);
}

// left is the index of the leftmost element of the subarray; right is one
// past the index of the rightmost element

//************************************************************************************
//SUBROUTINE IX: merge_helper: This function is a helper function for MergeSort.
//************************************************************************************
void merge_helper(dgrid *input, int left, int right, dgrid *scratch)
{
	/* base case: one element */
	if(right == left + 1)
		return;
	else
	{
		int i = 0;
		int length = right - left;
		int midpoint_distance = length/2;
		/* l and r are to the positions in the left and right subarrays */
		int l = left, r = left + midpoint_distance;

		/* sort each subarray */
		merge_helper(input, left, left + midpoint_distance, scratch);
		merge_helper(input, left + midpoint_distance, right, scratch);

		/* merge the arrays together using scratch for temporary storage */
		for(i = 0; i < length; i++)
		{
			/* Check to see if any elements remain in the left array; if so,
			 * we check if there are any elements left in the right array; if
			 * so, we compare them.  Otherwise, we know that the merge must
			 * use take the element from the left array */
			if(l < left + midpoint_distance &&
					(r == right || std::min((*input)(l,0), (*input)(r,0)) == (*input)(l,0)))
			{
				(*scratch)(i,0) = (*input)(l,0);
				(*scratch)(i,1) = (*input)(l,1);
				l++;
			}
			else
			{
				(*scratch)(i,0) = (*input)(r,0);
				(*scratch)(i,1) = (*input)(r,1);
				r++;
			}
		}
		/* Copy the sorted subarray back to the input */
		for(i = left; i < right; i++)
		{
			(*input)(i,0) = (*scratch)(i - left,0);
			(*input)(i,1) = (*scratch)(i - left,1);
		}
	}
}

//************************************************************************************
//SUBROUTINE X: MergeSort: This function is a simple merge sort function
//to sort data so it can be parsed. Elements are sorted in from least to greatest.
// To sort elements in reverse order from greatest to least
// replace mine with max.
//************************************************************************************
int MergeSort(dgrid *input, int size)
{
	// dgrid *scratch = (dgrid *)malloc(2*size * sizeof(double));
	dgrid *scratch = new dgrid(size, 2, 0);
	if(scratch != NULL)
	{
		merge_helper(input, 0, size, scratch);
		delete scratch;
		return 1;
	}
	else
		return 0;
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ORIGINAL FUNCTIONS FOR MONTE CARLO RUNS!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//************************************************************************************
//SUBROUTINE IV: subsampleData:
//This subsamples and spaces the data points.
//************************************************************************************
//int subsampleData_ORIGINAL(vector<double> *inputDataX, vector<double> *inputDataY, vector<double> *inputDataZ, vector<double> *inputDataE, double gridSpacingX, double gridSpacingY, double &inX0, double &inY0, double meanX, double meanY, bool dispIntermResults, vector< vector<double> > *subsampledData)
//{
//	//************************************************************************************
//	//0. Declare local variables and objects.
//	//************************************************************************************
//	//A. Initialize local variables
//	int returnValue = 0;
//
//	double wtol = 0.10000;
//	double jMaxX, jMaxY;
//	double DX = gridSpacingX;
//	double DY = gridSpacingY;
//	double x0 = (floor(inX0/DX))*(DX);
//	double y0 = (floor(inY0/DY))*(DY);
//
//	double s3Value;
//	int jProd = 1;
//	unsigned int i, j;
//	unsigned int count;
//	uint inputDataXSize = (uint)(*inputDataX).size();
//
//	dgrid Ji(inputDataXSize,2,1);
//	dgrid J(inputDataXSize,2,0);
//
//	vector<int>::iterator iter;
//	vector<double> weights = vector<double>(inputDataXSize, 2.00);
//
//	if (dispIntermResults)
//		printf("Subsampling Data: ");
//
//	//B. Initialize the Ji, and J dgrids with coordinate data
//	J(0,0) = round(1.00+(((*inputDataX)[0] - x0)/DX));
//	J(0,1) = round(1.00+(((*inputDataY)[0] - y0)/DY));
//	jMaxX = J(0,0);
//	jMaxY = J(0,1);
//
//	//C. Now loop to set the values and square the errors
//	for (i = 0; i < inputDataXSize; i++){
//		(*inputDataE)[i] = pow((*inputDataE)[i], 2);
//		Ji(i,0) = 1;
//		Ji(i,1) = i;
//		J(i,0) = round(1.00+(((*inputDataX)[i] - x0)/DX));
//		J(i,1) = round(1.00+(((*inputDataY)[i] - y0)/DY));
//
//		if (J(i,0) > jMaxX)
//			jMaxX = J(i,0);
//		if (J(i,1) > jMaxY)
//			jMaxY = J(i,1);
//	}
//	if (dispIntermResults)
//		printf(".");
//
//	//************************************************************************************
//	//I. Begin subsampling the data
//	//************************************************************************************
//	//A. Compute the weights of each data point
//	consistentWeights(inputDataZ, inputDataE, &wtol, &weights, &s3Value);
//
//	//************************************************************************************
//	//II. Fill Ji and get jProd
//	//************************************************************************************
//	//A. Map J to index of a multi-D array of unique indices. Jprod will store
//	//rescaled versions of multi-D array J into a M*N x 1 vector Ji. We repeat
//	//the multiplication process to ensure that different N-blocks of indicies
//	//have different numbers (i.e. no way to repeat an x and y index due to
//	//resetting JProd on first go and then multiplication of y-indices by JProd
//	//on second go. The whole point of this is to cound the number of unique
//	//DX cell distances.
//	jProd = 1;
//
//	for (i = 0; i < 2; i++){
//		for (j = 0; j < inputDataXSize; j++){
//			Ji(j,0) = Ji(j,0) + ((J(j,i) - 1) * (double)jProd);
//		}
//		if (i==0)
//			jProd = jProd * (int)jMaxX;
//		else
//			jProd = jProd * (int)jMaxY;
//	}
//	if (dispIntermResults)
//		printf(".");
//
//	//B. Sort data in ascending order and keep track of the indexes
//	Quicksort(&Ji, 0, inputDataXSize - 1);
//
//	if (dispIntermResults)
//		printf("...");
//
//	//************************************************************************************
//	//III. Remove repeated points
//	//************************************************************************************
//	//A. Resize the vector and remove unique elements
//	vector<int> JiVector = vector<int>(inputDataXSize,0);
//	for (i = 0; i < inputDataXSize; i++){
//		JiVector[i] = (const int)Ji(i,0);
//	}
//	iter = unique(JiVector.begin(), JiVector.end());
//	JiVector.resize(iter-JiVector.begin());
//
//	if (dispIntermResults)
//		printf(".");
//
//	//B. Set up the computation vectors
//	uint JiVectorSize = (uint)JiVector.size();
//	zvector zi(JiVectorSize,0.0);
//	zvector ni(JiVectorSize,0.0);
//	zvector wi(JiVectorSize,0.0);
//	zvector w2zi(JiVectorSize,0.0); // weights against data
//	zvector w2ei(JiVectorSize,0.0); // weights against a priori errors
//	zvector w2i(JiVectorSize,0.0); // sum of squared weights
//	zvector si(JiVectorSize,0.0); // holds weighted sum of squares, initially
//	dvector JiNew(JiVectorSize,0.0); //Reset vector Ji for new use.
//
//	JiNew[0] = Ji(0,0);
//
//	//C. Set up the output vectors
//	(*subsampledData)[0] = vector<double>(JiVectorSize, 0.0);
//	(*subsampledData)[1] = vector<double>(JiVectorSize, 0.0);
//	(*subsampledData)[2] = vector<double>(JiVectorSize, 0.0);
//	(*subsampledData)[3] = vector<double>(JiVectorSize, 0.0);
//	(*subsampledData)[4] = vector<double>(JiVectorSize, 0.0);
//	count = 0;
//
//	//************************************************************************************
//	//IV. Subsample and apply weights
//	//************************************************************************************
//	//A. Loop through and set values
//	for (i = 0; i < inputDataXSize; i++){
//		if(Ji(i,0) > JiNew[count]){
//			count = count + 1;
//			JiNew[count] = Ji(i,0);
//		}
//
//		ni[count] = ni[count].real() + 1;
//		wi[count] =  wi[count].real() + weights[(int)Ji(i,1)];
//		w2i[count] =  w2i[count].real() + pow(weights[(int)Ji(i,1)],2);
//		zi[count] =  zi[count].real() + (*inputDataZ)[(int)Ji(i,1)] * weights[(int)Ji(i,1)];
//		w2zi[count] =  w2zi[count].real() + (*inputDataZ)[(int)Ji(i,1)] * pow(weights[(int)Ji(i,1)],2);
//		si[count] =  si[count].real() + pow((*inputDataZ)[(int)Ji(i,1)] * weights[(int)Ji(i,1)], 2);
//		w2ei[count] = w2ei[count] + (*inputDataE)[(int)Ji(i,1)] * pow(weights[(int)Ji(i,1)],2);
//		(*subsampledData)[0][count] = (*subsampledData)[0][count] + (*inputDataX)[(int)Ji(i,1)] * weights[(int)Ji(i,1)];
//		(*subsampledData)[1][count] = (*subsampledData)[1][count] + (*inputDataY)[(int)Ji(i,1)] * weights[(int)Ji(i,1)];
//	}
//	if (dispIntermResults)
//		printf("...");
//
//	//B. Set the subdata values with the real data from the zgrid complex data
//	for (i = 0; i < JiVectorSize; i++){
//		zi[i] = zi[i].real() / wi[i].real();
//		(*subsampledData)[0][i] = (*subsampledData)[0][i] / wi[i].real();
//		(*subsampledData)[1][i] = (*subsampledData)[1][i] / wi[i].real();
//		(*subsampledData)[0][i] -= meanX;
//		(*subsampledData)[1][i] -= meanY;
//
//		si[i] = (si[i].real() - (2.00 * zi[i].real() * w2zi[i].real()) + (w2i[i].real() * pow(zi[i].real(),2)))/(ni[i].real());
//
//		si[i] = (si[i].real())/(1.00 + ni[i].real());
//
//		si[i] = (((ni[i].real() - 1.00) * si[i]) + (w2ei[i].real()/w2i[i].real()))/(ni[i].real());//sqrt( (((ni.get(i,0).real() - 1) * si.get(i,0).real()) + (w2ei.get(i,0).real()/ w2i.get(i,0).real()))/(ni.get(i,0).real()) ));
//		si[i] = (std::complex<double>)sqrt(si[i]);
//
//		(*subsampledData)[2][i] = (double)zi[i].real();
//		(*subsampledData)[3][i] = (double)si[i].real();
//		(*subsampledData)[4][i] = pow((double)si[i].real(),2);
//	}
//
//	//C. Do a second quicksort to get all of the X values in a row and maintain the rest of the data in the proper indexes
//	//Quicksort2(&(*subsampledData)[0], &(*subsampledData)[1], &(*subsampledData)[2], &(*subsampledData)[3], &(*subsampledData)[4], 0, JiVectorSize);
//	if (dispIntermResults)
//		printf("..");
//
//	inX0 = x0;
//	inY0 = y0;
//
//	//************************************************************************************
//	//V. Clean up variables and return.
//	//************************************************************************************
//	Ji.clear();
//	J.clear();
//	weights.clear();
//	JiVector.clear();
//	zi.clear();
//	ni.clear();
//	wi.clear();
//	w2zi.clear();
//	w2ei.clear();
//	w2i.clear();
//	si.clear();
//	JiNew.clear();
//
//	if (dispIntermResults)
//		printf(" Done Subsampling Data\n");
//
//	return returnValue;
//}