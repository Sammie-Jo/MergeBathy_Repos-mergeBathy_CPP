#include "kriging.h"
#include "outFileStructs.h"

bool krigingDebugValue = false;

int roundKrig(double x){
	if((ceil(x)-x) < .5){
		return (int)ceil(x);
	}else{
		return (int)floor(x);
	}
}

void ordinaryKrigingOfResiduals_PreCompute(vector<double> *subX, vector<double> *subY, vector<double> *residualObservations, vector<double> *twoGammaHatVector, vector<double> *distanceVectorBinCenters, vector<double> *aVectorFine, dgrid *invGammaDArray, dgrid *A)
{
	//Declare some variables
	dgrid twoGammaHat;
	dgrid twoGammaHatFine;

	dgrid distanceArray((const uint)(*residualObservations).size(),(const uint)(*residualObservations).size(), 0.00);
	dgrid angleArray(distanceArray);
	dgrid deltaRSquaredArray(distanceArray);

	dgrid distanceArray2;
	dgrid angleArray2;
	dgrid deltaRSquaredArray2;
	dgrid indexMax;
	dgrid workingGrid;

	vector<double> angleVector = vector<double>(8);
	vector<double> angleVectorPerpendicular = vector<double>(4);
	vector<double> distanceVector;
	vector<double> distanceVector2;
	vector<int> r95Indices;
	vector<double> twoGammaBar = vector<double>(8,0);
	vector<double> twoGammaHatPerpendicular;
	vector<double> bVector;

	double min1, max1, min2, max2;
	double tileSize, distanceBinSize;
	double tempSum, tempMax;//, tempMin;
	double tempTermDFT_Re, tempTermDFT_Im;
	double standardDeviationRi, meanRi;
	double angleMultiplier = 180.00 / PI;
	double alphaAlt = 1.00;
	double workingDouble1, workingDouble2, cValue, gamma;
	double gammaComputedX, gammaComputedY;
	double phi_2;
	double thetaM;

	int i, j, k;
	int riSize = (const int)(*residualObservations).size();
	int indexLoc;
	int gammaDZeroIndexLoc;

	angleVector[0] = -180.00;
	angleVector[1] = -135.00;
	angleVector[2] = -90.00;
	angleVector[3] = -45.00;
	angleVector[4] = 0.00;
	angleVector[5] = 45.00;
	angleVector[6] = 90.00;
	angleVector[7] = 135.00;

	min1 = max1 = (*subX)[0];
	min2 = max2 = (*subY)[0];

	standardDeviationRi = standardDeviation(residualObservations);
	meanRi = mean(residualObservations);

	//Get the minimum and maximum values of the sub data in the x and y directions
	//Also get the indices of the residuals that lie between the bounds
	for (i = 0; i < (const int)(*subX).size(); i++){
		if ((*subX)[i] > max1)
		{
			max1 = (*subX)[i];
		}else if ((*subX)[i] < min1)
		{
			min1 = (*subX)[i];
		}
		if ((*subY)[i] > max2)
		{
			max2 = (*subY)[i];
		}else if ((*subY)[i] < min2)
		{
			min2 = (*subY)[i];
		}
		if ( ( (*residualObservations)[i] > (-2.00 * standardDeviationRi + meanRi) ) && ( (*residualObservations)[i] < (2.00*standardDeviationRi+meanRi) ) ){
			r95Indices.push_back(i);
		}
	}

	//Establish tile size
	if ((max1-min1) > (max2-min2))
	{
		tileSize = (max1-min1);
	}else
	{
		tileSize = (max2-min2);
	}

	distanceBinSize = tileSize / 100.00;

	//Get the distances that we are computing
	j = roundKrig(tileSize/distanceBinSize);
	for (i = 0; i <= j; i++)
	{
		distanceVector.push_back((double)i * distanceBinSize);
	}

	//Calculate the array information
	//First get distances, angles and squared difference in residuals. Place
	//into n-by-n square-matrix arrays.
	for (i = 0; i < riSize; i++){
		for (j = 0; j < riSize; j++){
			distanceArray(i,j) = sqrt( pow((*subX)[i] - (*subX)[j], 2) + pow((*subY)[i] - (*subY)[j], 2) );
			angleArray(i,j) = angleMultiplier*atan2( (*subY)[i] - (*subY)[j],(*subX)[i] - (*subX)[j] );
			deltaRSquaredArray(i,j) = pow((*residualObservations)[i] - (*residualObservations)[j], 2);
		}
	}

	//Use Methods of Moments Estimator as coded in Appendix A below to get
	//empirical variogram, two_gamma_hat. Calder's Eqns. (15-17). Also in Cressie, 1993.
	methodsOfMomentsEstimator(&distanceVector, &angleVector, &distanceArray, &angleArray, &deltaRSquaredArray, &twoGammaHat);

	//Compute phi_2 and thetaM
	//Need to fit two_gamma_bar to second Fourier eigen function. Take FFT of two_gamma_bar.
	//Let Two_Gamma_Bar = fft(two_gamma_bar). Then, coefficients g0 and g2 in Eqn (18) of
	// Calder's paper are g0 = 2*abs(H_bar(1)) and g2 = abs(H_bar(3)). The phase
	// angle, phi_2 = atan2(imag(H_bar(3)), real(H_bar(3))) - pi/2.
	// The real-world angle of the axis of minimum variation is theta_m = (-phi_2)/2 + pi/2.
	k = 2;
	tempTermDFT_Re = 0;
	tempTermDFT_Im = 0;
	phi_2 = 0.00;
	thetaM = 0.00;
	for (i = 0; i < 8; i++){
		tempSum = 0.00;
		for (j = 0; j < (const int)twoGammaHat.cols(); j++){
			tempSum += twoGammaHat(i,j);
		}
		twoGammaBar[i] = (double)(1.00 / (double)distanceVector.size()) * tempSum;

		tempTermDFT_Re = tempTermDFT_Re + twoGammaBar[i]*sin(2.00*PI*((double)k)*((double)(i+1))/8.00);
		tempTermDFT_Im = tempTermDFT_Im + twoGammaBar[i]*cos(2.00*PI*((double)k)*((double)(i+1))/8.00);
		//k+=1;
	}

	phi_2 = atan2((tempTermDFT_Im), tempTermDFT_Re)/2.00;

	thetaM = ((-1.00*phi_2) + PI/2.00)*180.00/PI;

	// III. Fit the variograms along the two axes found in Section II to the
	// spherical model for variograms. Use Levenberg-Marquardt algorithm for
	// fitting parameters in the spherical model.
	// 1. Compute fine-scale empirical variograms in the orthogonal directions
	// theta_m and theta_m_perpendicular using method of moments estimator again
	// and angle bins pi/2.

	angleVectorPerpendicular[0] = -90.00-thetaM;
	angleVectorPerpendicular[1] = 0.00-thetaM;
	angleVectorPerpendicular[2] = 90.00-thetaM;
	angleVectorPerpendicular[3] = 180.00-thetaM;

	//a. Use distance_array2 = 2.5 times more resolution.
	j = roundKrig(tileSize/(distanceBinSize/2.5));
	for (i = 0; i <= j; i++)
	{
		distanceVector2.push_back((double)i * (distanceBinSize/2.5));
	}

	distanceArray2 = dgrid((const uint)r95Indices.size(),(const uint)r95Indices.size());
	angleArray2 = dgrid((const uint)r95Indices.size(),(const uint)r95Indices.size());
	deltaRSquaredArray2 = dgrid((const uint)r95Indices.size(),(const uint)r95Indices.size());
	for (i = 0; i < (const int)r95Indices.size(); i ++){
		for (j = 0; j < (const int)r95Indices.size(); j++){
			distanceArray2(i,j) = distanceArray.get(r95Indices.at(i), r95Indices.at(j));
			angleArray2(i,j) = angleArray.get(r95Indices.at(i), r95Indices.at(j));
			deltaRSquaredArray2(i,j) = deltaRSquaredArray.get(r95Indices.at(i), r95Indices.at(j));
		}
	}

	methodsOfMomentsEstimator(&distanceVector2, &angleVectorPerpendicular, &distanceArray2, &angleArray2, &deltaRSquaredArray2, &twoGammaHatFine);

	(*twoGammaHatVector) = vector<double> (twoGammaHatFine.cols());
	twoGammaHatPerpendicular = vector<double> (twoGammaHatFine.cols());
	if ( ((twoGammaHatFine.get(twoGammaHatFine.rows()-1,twoGammaHatFine.cols()-1) > twoGammaHatFine.get(twoGammaHatFine.rows()-2,twoGammaHatFine.cols()-1)) || (twoGammaHatFine.get(twoGammaHatFine.rows()-1,twoGammaHatFine.cols()-1) == 0)) && (twoGammaHatFine.get(twoGammaHatFine.rows()-2,twoGammaHatFine.cols()-1) > 0) ){
		max2 = twoGammaHatFine(twoGammaHatFine.rows()-2,0);
		for (i = 0; i < (const int)(twoGammaHatFine.cols()); i++){
			(*twoGammaHatVector)[i] = twoGammaHatFine(twoGammaHatFine.rows()-2,i);
			twoGammaHatPerpendicular[i] = twoGammaHatFine(twoGammaHatFine.rows()-1,i);
			if ((*twoGammaHatVector)[i] > max2)
				max2 = (*twoGammaHatVector)[i];
		}
	}else{
		max2 = twoGammaHatFine(twoGammaHatFine.rows()-1,0);
		for (i = 0; i < (const int)(twoGammaHatFine.cols()); i++){
			(*twoGammaHatVector)[i] = twoGammaHatFine(twoGammaHatFine.rows()-1,i);
			twoGammaHatPerpendicular[i] = twoGammaHatFine(twoGammaHatFine.rows()-2,i);
			if ((*twoGammaHatVector)[i] > max2)
				max2 = (*twoGammaHatVector)[i];
		}
	}

	(*distanceVectorBinCenters) = vector<double>(distanceVector2.size()-1, 0.00);
	for (i = 0; i < (const int)(*distanceVectorBinCenters).size(); i++){
		(*distanceVectorBinCenters)[i] = distanceVector2[i] + ((distanceVector2[0] + distanceVector2[1])/2.00);
	}

	dgrid indexMinMaxLimits((const uint)angleVectorPerpendicular.size(), 1, 0.00);
	for (i = 0; i < (const int)angleVectorPerpendicular.size(); i++){
		tempMax = twoGammaHatFine(i,0);
		for (j = 0; j < (const int)twoGammaHatFine.cols(); j++)
		{
			if(twoGammaHatFine(i,j) >  tempMax)
			{
				indexMinMaxLimits(i,0) = j;
				tempMax = twoGammaHatFine(i,j);
			}
		}
	}

	max1 = indexMinMaxLimits.max();
	min1 = indexMinMaxLimits.min();
	indexMinMaxLimits.clear();

	//IV. Transform variogram from Section III to final semivariogram in accord
	// with Calder's recipe below. Calder's Eqns. (20-22)
	if(distanceVector2[(int)min1] == 0.00){
		alphaAlt = 1.00;
	}else
	{
		alphaAlt = distanceVector2[(int)max1] / distanceVector2[(int)min1]; //distanceVector2.at(InDem);
	}
	indexLoc = (int)floor(distanceVector2.size()*max1 / twoGammaHatFine.cols());

	(*aVectorFine) = vector<double>(3);
	(*aVectorFine)[0] = 0.00;
	(*aVectorFine)[1] = 0.50*(max2 + (*twoGammaHatVector)[(*twoGammaHatVector).size()-1]);
	(*aVectorFine)[2] = distanceVector2[indexLoc];

	workingGrid = dgrid(2,2,0.00);
	workingGrid(0,0) = 1.00;
	workingGrid(1,1) = alphaAlt;

	//A. Compute A matrix, Calder's Eqn. (22), to account for anisotropy.
	(*A) = dgrid(trans(matmult(matmult(rotMtx(-thetaM),workingGrid),rotMtx(thetaM))));
	workingGrid.clear();

	//Initialize storage vector for Calder's change to distance to account for anisotropy. Calder's Eq. (21)
	dgrid distancePrimeArray = dgrid(riSize,riSize, 0);
	dgrid gammaDArray = dgrid(riSize+1, riSize+1, 1);

	gammaDZeroIndexLoc = 0;
	//cout << aVectorFine << endl;
	for (i = 0; i < riSize; i++)
	{
		for (j = 0; j < riSize; j++)
		{
			//Apply Calder's change to distance to account for anisotropy, Calder's Eqn. (21)
			workingDouble1 = (*subX)[i]-(*subX)[j];
			workingDouble2 = (*subY)[i]-(*subY)[j];

			gammaComputedX = workingDouble1*(*A)(0,0) + workingDouble2*(*A)(1,0);
			gammaComputedY = workingDouble1*(*A)(0,1) + workingDouble2*(*A)(1,1);

			tempSum = pow(gammaComputedX,2) + pow(gammaComputedY,2);
			if (tempSum <=0)
			{
				distancePrimeArray(i,j) = 0;
			}else
			{
				distancePrimeArray(i,j) = sqrt(tempSum);
			}

			//Compute semivariogram from data point to interpolation point.
			if (1 >= (distancePrimeArray(i,j) / (*aVectorFine)[2])){
				cValue = distancePrimeArray(i,j) / (*aVectorFine)[2];
				gamma =  (*aVectorFine)[0] +  (*aVectorFine)[1] * ( 1.5000 * cValue  - 0.5000 * pow(cValue,3) );
			}else
			{
				gamma =  (*aVectorFine)[0] +  (*aVectorFine)[1];
			}
			gammaDArray(i,j) = 0.5*gamma;
			//Since fit of modeled variogram may give negative numbers, fallback to emperical variogram if needed.
			if ((0.5*gamma) <= 0)
			{
				min1 = pow(((*distanceVectorBinCenters)[0] - distancePrimeArray(i,j)), 2);
				indexLoc = 0;
				for (k = 0; k < (const int)(*distanceVectorBinCenters).size(); k++){
					if (pow(((*distanceVectorBinCenters)[k] - distancePrimeArray(i,j)), 2) < min1){
						min1 = pow(((*distanceVectorBinCenters)[k] - distancePrimeArray(i,j)), 2);
						indexLoc = k;
					}
				}
				gammaDArray(i,j) = (*twoGammaHatVector)[indexLoc];
			}
		}
	}
	//Compute solution for the interpolated residual surface using ordinary
	// kriging as provided in Davis JC, Statistics and Data Analysis in Geology,
	// 3rd edition, pp. 420-1, New York: Wiley, 2002.
	gammaDArray(riSize,riSize) = 0.0;
	bVector = vector<double>(riSize+1, 1);

	// Get inverse of W for later use
	//try {
		(*invGammaDArray) = dgrid(inv(gammaDArray, 0, 'U', true, 'N'));
	/*}
	catch (exception &e) {
		cout << "An exception occurred. Exception Thrown: " << e.what() << '\n';
	}*/

	//Clean up variable space
	distancePrimeArray.clear();
	gammaDArray.clear();

	twoGammaHat.clear();
	twoGammaHatFine.clear();

	distanceArray.clear();
	angleArray.clear();
	deltaRSquaredArray.clear();
	distanceArray2.clear();
	angleArray2.clear();
	deltaRSquaredArray2.clear();
	indexMax.clear();
	workingGrid.clear();
	gammaDArray.clear();
	distancePrimeArray.clear();

	angleVector.clear();
	angleVectorPerpendicular.clear();
	distanceVector.clear();
	distanceVector2.clear();
	twoGammaBar.clear();
	r95Indices.clear();
	bVector.clear();
}

void ordinaryKrigingOfResiduals_PreComputeTile(vector<double> *subX, vector<double> *subY, vector<double> *residualObservations, vector<double> *interpX, vector<double> *interpY, double spacingX, double spacingY, vector<double> *outDepthKrig, vector<double> *outErrorKrig)
{
	//Initialize some variables
	double numberOfComputedTiles;
	double computedTilesPower;
	int i;//, j, k;
	int outerLoop, innerLoop;
	int kx, ky;
	double numStepsX, numStepsY;
	double numStepsXInterp, numStepsYInterp;
	double LMAX_x, LMAX_y;
	double minSubX, minSubY, maxSubX, maxSubY;
	double minInterpX, maxInterpX, minInterpY, maxInterpY;
	double zKriged, varZKriged;
	double newSpacingX, newSpacingY;

	vector<int> idx;
	vector<int> idxInterp;
	vector<double> subX_idx;
	vector<double> subY_idx;
	vector<double> ri_idx;

	double xmin, xmax, ymin, ymax;
	double xmin2, xmax2, ymin2, ymax2;
	vector<double> twoGammaHatVector;
	vector<double> distanceVectorBinCenters;
	vector<double> aVectorFine;
	dgrid invGammaDArray;
	dgrid AGrid;

	vector<double> localSubX = vector<double>((*subX));
	vector<double> localSubY = vector<double>((*subY));
	vector<double> localInterpX = vector<double>((*interpX));
	vector<double> localInterpY = vector<double>((*interpY));
	double meanSubX = 0.0000;
	double meanSubY = 0.0000;
	double meanX_std = 0.0000;
	double meanY_std = 0.0000;
	int subDataXLength = (const int)(*subX).size();
	int outputDepthSize = (const int)(*interpX).size();

	//int tempIn;
	//cout << "\n\nCHECK: ";
	//cin >> tempIn;

	//Compute the number of subtiles that we need to break this stuff down into
	i = 6;
	do
	{
		computedTilesPower = pow(2.0,i);
		i++;
	}while ((const int)(*residualObservations).size() >= computedTilesPower);

	//numberOfComputedTiles = (computedTilesPower / KRIGED_SIZE_THRESHOLD);

	for (i = 0; i < subDataXLength; i++)
	{
		meanSubX += localSubX[i];
		meanSubY += localSubY[i];
	}
	meanSubX /= (double)subDataXLength;
	meanSubY /= (double)subDataXLength;

	for (i = 0; i < outputDepthSize; i++)
	{
		localInterpX[i] -= meanSubX;
		localInterpY[i] -= meanSubY;
	}

	for (i = 0; i < subDataXLength; i++)
	{
		localSubX[i] -= meanSubX;
		localSubY[i] -= meanSubY;
		meanX_std = meanX_std + localSubX[i];
		meanY_std = meanY_std + localSubY[i];
	}
	meanX_std /= (double)subDataXLength;
	meanY_std /= (double)subDataXLength;

	//A. Compute the variance estimate
	double std_x = 0;
	double std_y = 0;
	for (i = 0; i < subDataXLength; i++){
		std_x = std_x + pow((localSubX[i] - meanX_std),2);
		std_y = std_y + pow((localSubY[i]- meanY_std),2);
	}
	std_x = std_x / (double)(subDataXLength - 1.00);
	std_y = std_y / (double)(subDataXLength - 1.00);
	std_x = sqrt(std_x);
	std_y = sqrt(std_y);

	//B. Scale the data and grid.
	localSubX /= std_x;
	localSubY /= std_y;
	localInterpX /= std_x;
	localInterpY /= std_y;

	newSpacingX = spacingX / std_x;
	newSpacingY = spacingY / std_y;

	minSubX = localSubX[0];
	minSubY = localSubY[0];
	maxSubX = localSubX[0];
	maxSubY = localSubY[0];
	for (i = 0; i < (const int)localSubX.size(); i++)
	{
		if (localSubX[i] < minSubX)
			minSubX = localSubX[i];
		else if (localSubX[i] > maxSubX)
			maxSubX = localSubX[i];
		if (localSubY[i] < minSubY)
			minSubY = localSubY[i];
		else if (localSubY[i] > maxSubY)
			maxSubY = localSubY[i];
	}

	minInterpX = localInterpX[0];
	minInterpY = localInterpY[0];
	maxInterpX = localInterpX[0];
	maxInterpY = localInterpY[0];
	for (i = 0; i < (const int)localInterpX.size(); i++)
	{
		if (localInterpX[i] < minInterpX)
			minInterpX = localInterpX[i];
		else if (localInterpX[i] > maxInterpX)
			maxInterpX = localInterpX[i];
		if (localInterpY[i] < minInterpY)
			minInterpY = localInterpY[i];
		else if (localInterpY[i] > maxInterpY)
			maxInterpY = localInterpY[i];
	}

	maxInterpX = maxInterpX*1.01;
	maxInterpY = maxInterpY*1.01;

	//THIS IS THE TRICKY STUFF
	numberOfComputedTiles = sqrt( (computedTilesPower) * (1.00 + abs(newSpacingX / (maxSubX - minSubX))) );

	if (numberOfComputedTiles < sqrt( (computedTilesPower) * (1.00 + abs(newSpacingY / (maxSubY - minSubY))) ))
		numberOfComputedTiles = sqrt( (computedTilesPower) * (1.00 + abs(newSpacingY / (maxSubY - minSubY))) );

	kx = (int) floor( sqrt( numberOfComputedTiles ));//was ceil
	ky = (int) floor( sqrt( numberOfComputedTiles ));
	numStepsX = (maxSubX - minSubX) / (double)kx;
	numStepsY = (maxSubY - minSubY) / (double)ky;
	numStepsXInterp = (maxInterpX - minInterpX) / (double)kx;
	numStepsYInterp = (maxInterpY - minInterpY) / (double)ky;
	LMAX_x = ((numStepsX)/double(kx)); //% specify max overlap between tiles, //SJZ because it showed less visual artifacts 2/17/16
	LMAX_y = ((numStepsY)/double(ky)); //% often = 10*length scale   //numberOfComputedTiles
	//LMAX_x = 5.00*(((maxSubX - minSubX)/(*subX).size())/double(numberOfComputedTiles)); //% specify max overlap between tiles,
	//LMAX_y = 5.00*(((maxSubY - minSubY)/(*subY).size())/double(numberOfComputedTiles)); //% often = 10*length scale

	//Initialize the tiling routine here
	for (outerLoop = 0; outerLoop < kx; outerLoop++)
	{
		//1. Compute appropriate overlap between tiles horizontally and get the
		//useful data. Recall that we scaled xi-array by 1./std(x). Then,
		//get indices of x-coordinates in x-array for the data to be
		//interpolated. Indices to be filtered further for useful y-coordinates.
		xmin = (minSubX + (outerLoop*numStepsX)) - LMAX_x;//(idxi(1))-LMAX(1); % find tile limits
		xmax = (minSubX + ((outerLoop+1.00)*numStepsX)) + LMAX_x;//(idxi(1))-LMAX(1); % find tile limits
		xmin2 = (minInterpX + (outerLoop*numStepsXInterp));//(idxi(1))-LMAX(1); % find tile limits
		xmax2 = (minInterpX + ((outerLoop+1.00)*numStepsXInterp));//(idxi(1))-LMAX(1); % find tile limits

		//DO XMIN2 XMAX2

		for (i = 0; i < (const int)localSubX.size(); i++){
			if ((localSubX[i] < xmax) && (localSubX[i] > xmin)){
				idx.push_back(i);
			}
		}
		//cout << "XMM " << xmin << " " << xmax << endl;
		//cout << "IDX: " << idx.size() << endl;

		for (i = 0; i < (const int)localInterpX.size(); i++){
			if ((localInterpX[i] < xmax2) && (localInterpX[i] >= xmin2))
			{
				idxInterp.push_back(i);
			}
		}

		if (idx.size() >= 4)
		{
			for (innerLoop = 0; innerLoop < ky; innerLoop++)
			{
				//2. Compute appropriate overlap between tiles horizontally and get the
				//useful data. Recall that we scaled xi-array by 1./std(x). Then,
				//get indices of x-coordinates in x-array for the data to be
				//interpolated. Indices to be filtered further for useful y-coordinates.
				ymin = (minSubY + (innerLoop*numStepsY)) - LMAX_y;//(idxi(1))-LMAX(1); % find tile limits
				ymax = (minSubY + ((innerLoop+1.00)*numStepsY)) + LMAX_y;//(idxi(1))-LMAX(1); % find tile limits
				ymin2 = (minInterpY + (innerLoop*numStepsYInterp));//(idxi(1))-LMAX(1); % find tile limits
				ymax2 = (minInterpY + ((innerLoop+1.00)*numStepsYInterp));//(idxi(1))-LMAX(1); % find tile limits

				for (i = 0; i < (const int)idx.size(); i++){
					if ((localSubY[idx[i]] < ymax) && (localSubY[idx[i]] > ymin)){
						subX_idx.push_back(localSubX[idx[i]]);
						subY_idx.push_back(localSubY[idx[i]]);
						ri_idx.push_back((*residualObservations)[idx[i]]);
					}
				}

				if (subX_idx.size() >= 4)
				{
					ordinaryKrigingOfResiduals_PreCompute(&subX_idx, &subY_idx, &ri_idx, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);

					//3. Krig the data
					//cout << "PROCESS CALL" << endl;
					for (i = 0; i < (const int)idxInterp.size(); i++)
					{
						if ( ((*interpX)[i] >= xmin && (*interpX)[i] <= xmax) && ((*interpY)[i] >= ymin && (*interpY)[i] <= ymax) )
						//if ( (localInterpY[idxInterp[i]] >= ymin2 && localInterpY[idxInterp[i]] < ymax2) )
						{
							zKriged = 0.00;
							varZKriged = 0.00;
							//Compute the residual value
							ordinaryKrigingOfResiduals_PostCompute(&subX_idx, &subY_idx, &ri_idx, &localInterpX[idxInterp[i]], &localInterpY[idxInterp[i]], &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);

							(*outDepthKrig)[idxInterp[i]] = zKriged;
							(*outErrorKrig)[idxInterp[i]] = varZKriged;
						}
					}

					twoGammaHatVector.clear();
					distanceVectorBinCenters.clear();
					aVectorFine.clear();
					invGammaDArray.clear();
					AGrid.clear();
				}
				subX_idx.clear();
				subY_idx.clear();
				ri_idx.clear();
			}
		}

		idx.clear();
		idxInterp.clear();
	}

	localSubX.clear();
	localSubY.clear();
	localInterpX.clear();
	localInterpY.clear();
}

/*void ordinaryKrigingOfResiduals_PreComputeTile_ORIG_NOT_ROBUST(vector<double> *subX, vector<double> *subY, vector<double> *residualObservations, vector<double> *interpX, vector<double> *interpY, double spacingX, double spacingY, vector<double> *outDepthKrig, vector<double> *outErrorKrig, double recursionDepth)
{
	//Initialize some variables
	int numberOfComputedTiles;
	int computedTilesPower;
	int i, j, k;
	int outerLoop, innerLoop;
	int kx, ky;
	double numStepsX, numStepsY;
	double LMAX_x, LMAX_y;
	double minSubX, minSubY, maxSubX, maxSubY;
	double zKriged, varZKriged;
	double newSpacingX, newSpacingY;

	vector<int> idx;
	vector<double> subX_idx;
	vector<double> subY_idx;
	vector<double> ri_idx;

	double xmin, xmax, ymin, ymax;
	vector<double> twoGammaHatVector;
	vector<double> distanceVectorBinCenters;
	vector<double> aVectorFine;
	dgrid invGammaDArray;
	dgrid AGrid;

	//int tempIn;
	//cout << "\n\nCHECK: ";
	//cin >> tempIn;

	//Comptue the number of subtiles that we need to break this stuff down into
	i = 6;
	do
	{
		computedTilesPower = pow(2.0,i);
		i++;
	}while ((*residualObservations).size() >= computedTilesPower);

	//numberOfComputedTiles = (computedTilesPower / KRIGED_SIZE_THRESHOLD);

	minSubX = (*subX)[0];
	minSubY = (*subY)[0];
	maxSubX = (*subX)[0];
	maxSubY = (*subY)[0];
	for (i = 0; i < (*subX).size(); i++)
	{
		if ((*subX)[i] < minSubX)
			minSubX = (*subX)[i];
		else if ((*subX)[i] > maxSubX)
			maxSubX = (*subX)[i];
		if ((*subY)[i] < minSubY)
			minSubY = (*subY)[i];
		else if ((*subY)[i] > maxSubY)
			maxSubY = (*subY)[i];
	}

	//THIS IS THE TRICKY STUFF
	numberOfComputedTiles = sqrt( ((double)computedTilesPower) * (1.00 + abs(spacingX / (maxSubX - minSubX))) );

	if (numberOfComputedTiles < sqrt( ((double)computedTilesPower) * (1.00 + abs(spacingY / (maxSubY - minSubY))) ))
		numberOfComputedTiles = sqrt( ((double)computedTilesPower) * (1.00 + abs(spacingY / (maxSubY - minSubY))) );

	kx = (int) ceil( sqrt( (double) numberOfComputedTiles ));
	ky = (int) ceil( sqrt( (double) numberOfComputedTiles ));
	numStepsX = (maxSubX - minSubX) / (double)kx;
	numStepsY = (maxSubY - minSubY) / (double)ky;
	LMAX_x = ((numStepsX)/double(kx)); //% specify max overlap between tiles,
	LMAX_y = ((numStepsY)/double(ky)); //% often = 10*length scale   //numberOfComputedTiles
	//LMAX_x = 5.00*(((maxSubX - minSubX)/(*subX).size())/double(numberOfComputedTiles)); //% specify max overlap between tiles,
	//LMAX_y = 5.00*(((maxSubY - minSubY)/(*subY).size())/double(numberOfComputedTiles)); //% often = 10*length scale

	//Initialize the tiling routine here
	for (outerLoop = 0; outerLoop < kx; outerLoop++)
	{
		//1. Compute appropriate overlap between tiles horizontally and get the
		//useful data. Recall that we scaled xi-array by 1./std(x). Then,
		//get indices of x-coordinates in x-array for the data to be
		//interpolated. Indices to be filtered further for useful y-coordinates.
		if (!idx.empty())
		{
			idx.clear();
			xmax = (minSubX + ((outerLoop+1)*numStepsX)) + LMAX_x;//(idxi(1))-LMAX(1); % find tile limits
		}else
		{
			xmin = (minSubX + (outerLoop*numStepsX)) - LMAX_x;//(idxi(1))-LMAX(1); % find tile limits
			xmax = (minSubX + ((outerLoop+1.00)*numStepsX)) + LMAX_x;//(idxi(1))-LMAX(1); % find tile limits
		}

		for (i = 0; i < (*subX).size(); i++){
			if (((*subX)[i] < xmax) && ((*subX)[i] > xmin)){
				idx.push_back(i);
			}
		}
		//cout << "XMM " << xmin << " " << xmax << endl;
		//cout << "IDX: " << idx.size() << endl;

		if (idx.size() >= 4*ky)
		{
			for (innerLoop = 0; innerLoop < ky; innerLoop++)
			{
				//2. Compute appropriate overlap between tiles horizontally and get the
				//useful data. Recall that we scaled xi-array by 1./std(x). Then,
				//get indices of x-coordinates in x-array for the data to be
				//interpolated. Indices to be filtered further for useful y-coordinates.
				if (!subX_idx.empty())
				{
					subX_idx.clear();
					subY_idx.clear();
					ri_idx.clear();
					ymax = (minSubY + ((innerLoop+1)*numStepsY)) + LMAX_y;//(idxi(1))-LMAX(1); % find tile limits
				}else
				{
					ymin = (minSubY + (innerLoop*numStepsY)) - LMAX_y;//(idxi(1))-LMAX(1); % find tile limits
					ymax = (minSubY + ((innerLoop+1.00)*numStepsY)) + LMAX_y;//(idxi(1))-LMAX(1); % find tile limits
				}

				for (i = 0; i < idx.size(); i++){
					if (((*subY)[idx[i]] < ymax) && ((*subY)[idx[i]] > ymin)){
						subX_idx.push_back((*subX)[idx[i]]);
						subY_idx.push_back((*subY)[idx[i]]);
						ri_idx.push_back((*residualObservations)[idx[i]]);
					}
				}

				if (subX_idx.size() >= 4)
				{
					//cout << "SX_IDX: " << subX_idx.size() << endl;
					if (subX_idx.size() > KRIGED_SIZE_THRESHOLD && recursionDepth <= 8)//*2
					{
						//cout << "RECURSIVE CALL" << endl;
						newSpacingX = spacingX * (recursionDepth + 1.00);
						newSpacingY = spacingY * (recursionDepth + 1.00);
						ordinaryKrigingOfResiduals_PreComputeTile(&subX_idx, &subY_idx, &ri_idx, interpX, interpY, newSpacingX, newSpacingY, outDepthKrig, outErrorKrig, recursionDepth+1.00);
					}else
					{
						//cout << "PRE-PROCESS CALL" << endl;
						ordinaryKrigingOfResiduals_PreCompute(&subX_idx, &subY_idx, &ri_idx, &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid);

						//3. Krig the data
						//cout << "PROCESS CALL" << endl;
						for (i = 0; i < (*interpX).size(); i++)
						{
							if ( ((*interpX)[i] >= xmin && (*interpX)[i] <= xmax) && ((*interpY)[i] >= ymin && (*interpY)[i] <= ymax) )
							{
								zKriged = 0.00;
								varZKriged = 0.00;
								//Compute the residual value
								ordinaryKrigingOfResiduals_PostCompute(&subX_idx, &subY_idx, &ri_idx, &(*interpX)[i], &(*interpY)[i], &twoGammaHatVector, &distanceVectorBinCenters, &aVectorFine, &invGammaDArray, &AGrid, &zKriged, &varZKriged);

								if ((*outDepthKrig)[i] == 0)
								{
									(*outDepthKrig)[i] = zKriged;
									(*outErrorKrig)[i] = varZKriged;
								}else if ((*outErrorKrig)[i] > varZKriged)
								{
										(*outDepthKrig)[i] = zKriged;
										(*outErrorKrig)[i] = varZKriged;
								}
							}
						}

						twoGammaHatVector.clear();
						distanceVectorBinCenters.clear();
						aVectorFine.clear();
						invGammaDArray.clear();
						AGrid.clear();
					}
				}else
				{
					continue;
				}
				subX_idx.clear();
				subY_idx.clear();
				ri_idx.clear();
			}
			subX_idx.clear();
			subY_idx.clear();
			ri_idx.clear();
		}else
		{
			continue;
		}

		idx.clear();
	}
	//cout << "DONE LOOP" << endl;
	idx.clear();
}*/

void ordinaryKrigingOfResiduals_PostCompute(vector<double> *subX, vector<double> *subY, vector<double> *residualObservations, double *xGridValue, double *yGridValue, vector<double> *twoGammaHatVector, vector<double> *distanceVectorBinCenters, vector<double> *aVectorFine, dgrid *invGammaDArray, dgrid *A, double *zKrigValue, double *varZKrigValue)
{
	//Set up the variables
	vector<double> B = vector<double>((*residualObservations).size()+1, 1.0);
	vector<double> Y = vector<double>((*residualObservations).size()+1, 0.0);
	double distancePrimeInterp;
	double xValueTemp, yValueTemp;
	double gammaComputedX, gammaComputedY;
	double tempSum;
	double cValue, gamma, tempGammaMin;
	double zComputed, varZComputed;
	int i, j, k;

	//Find ordinary kringing solution for residuals at the interpolated surfaces
	//Build B Vector, Eqn. (5.102)
	for (i = 0; i < (const int)(*residualObservations).size(); i++)
	{
		//Apply Calder's change to distance to account for anisotropy
		xValueTemp = (*subX)[i]-(*xGridValue);
		yValueTemp = (*subY)[i]-(*yGridValue);

		gammaComputedX = xValueTemp*(*A)(0,0) + yValueTemp*(*A)(1,0);
		gammaComputedY = xValueTemp*(*A)(0,1) + yValueTemp*(*A)(1,1);

		tempSum = pow(gammaComputedX,2) + pow(gammaComputedY,2);
		if (tempSum <=0)
		{
			distancePrimeInterp = 0;
		}else
		{
			distancePrimeInterp = sqrt(tempSum);
		}

		if (1 >= (distancePrimeInterp / (*aVectorFine)[2])){
			cValue = distancePrimeInterp / (*aVectorFine)[2];
			gamma =  (*aVectorFine)[0] +  (*aVectorFine)[1] * ( 1.5000 * cValue  - 0.5000 * pow(cValue,3) );
		}else
		{
			gamma =  (*aVectorFine)[0] +  (*aVectorFine)[1];
		}

		//Compute semivariogram from data point to interpolation point.
		B[i] = 0.50*gamma;

		//Since fit of modeled variogram may give negative numbers, fallback to emperical variogram if needed.
		if (B[i] <= 0)
		{
			k = 0;
			for (j = 0; j < (const int)(*distanceVectorBinCenters).size(); j++){
				if (j == 0){
					tempGammaMin = pow(((*distanceVectorBinCenters).at(j) - distancePrimeInterp),2);
					k = j;
				}else{
					if (pow(((*distanceVectorBinCenters)[j] - distancePrimeInterp),2) < tempGammaMin){
						tempGammaMin = pow(((*distanceVectorBinCenters)[j] - distancePrimeInterp),2);
						k = j;
					}
				}
			}
			B[i] = (*twoGammaHatVector)[k];
		}
		//Build Y Vector in Davis, Eqn, (5.103).
		Y[i] = (*residualObservations)[i];
	}

	//Davis's solution for ordinary kriging equations, Eqn (5.105).
	// Davis's solution for the variance of the interpolated point when semmivariogram is used, Eqn. (5.106).
	zComputed = 0;
	varZComputed = 0;
	tempSum = 0;
	for (i = 0; i < (const int)B.size(); i++)
	{
		tempSum = 0;
		for (j = 0; j < (const int)(*invGammaDArray).cols(); j++)
		{
			tempSum = tempSum + (*invGammaDArray)(i,j)*B[j];
		}
		//Davis's solution for ordinary kriging equations, Eqn (5.105).
		zComputed = zComputed + (Y[i]*tempSum);

		// Davis's solution for the variance of the interpolated point when semmivariogram is used, Eqn. (5.106).
		varZComputed = varZComputed + (B[i]*tempSum);
	}

	(*zKrigValue) = zComputed;
	(*varZKrigValue) = varZComputed;

	B.clear();
	Y.clear();
}

void methodsOfMomentsEstimator(vector<double> *distanceVector, vector<double> *angleVector, dgrid *distanceArray, dgrid *angleArray, dgrid *deltaRSquaredArray, dgrid *twoGammaHat)
{
	(*twoGammaHat) = dgrid((const uint)(*angleVector).size(), (const uint)(*distanceVector).size(), 0.00);

	int currentTempVecLoc = 0;
	int currentAngleVecLoc = 0;

	vector<double> tempVec;
	vector<double> filteredDistanceVector;
	vector<double> filteredDeltaRSquaredVector;
	vector<int> hatAngle;

	double angleBinSize = 360.00/((double)(*angleVector).size());
	double deltaRSquaredArraySum = 0.00;
	double nDeltaR = 0.00;
	double angleMultiplier = 180.00 / PI;
//	double workingDouble;
	double maxVal;

	int i,j,k;

	//Bin distances
	for (i = 0; i < (const int)(*angleVector).size(); i++){
		currentAngleVecLoc = 0;

		//Find and sum points in the bins
		if (i == 0){
			for (j = 0; j < (const int)(*angleArray).cols(); j++){
				for (k = 0; k < (const int)(*angleArray).rows(); k++){
					if ( ( ( (*angleVector)[ (*angleVector).size()-1 ] + (angleBinSize / 2.00) ) <= (*angleArray)(k,j) ) || ( (*angleArray)(k,j) < ( (*angleVector)[0] + (angleBinSize / 2.00) ) ) ){
						filteredDistanceVector.push_back((*distanceArray)(k,j));
						filteredDeltaRSquaredVector.push_back((*deltaRSquaredArray)(k, j));

						currentAngleVecLoc += 1;
					}
				}
			}
		}else{
			for (j = 0; j < (const int)(*angleArray).cols(); j++){
				for (k = 0; k < (const int)(*angleArray).rows(); k++){
					if ( ( ( (*angleVector)[i] - (angleBinSize / 2.00) ) <= (*angleArray)(k,j) ) && ( (*angleArray)(k,j) < ( (*angleVector)[i] + (angleBinSize / 2.00) ) ) ){
						filteredDistanceVector.push_back((*distanceArray)(k,j));
						filteredDeltaRSquaredVector.push_back((*deltaRSquaredArray)(k, j));

						currentAngleVecLoc += 1;
					}
				}
			}
		}

		deltaRSquaredArraySum = 0;
		nDeltaR = 0;

		//Calculate h_hat for the jj-th angle and the kk-th distance bin
		tempVec = vector<double>(filteredDistanceVector.size(),0.00);
		for (j = 0; j < ((const int)(*distanceVector).size() - 1); j++){
			currentTempVecLoc = 0;
			for (k = 0; k < (const int)filteredDistanceVector.size(); k++){
				if ( ((*distanceVector)[j] <= filteredDistanceVector[k]) && (filteredDistanceVector[k] < (*distanceVector)[j+1])){
					tempVec[currentTempVecLoc] = filteredDeltaRSquaredVector[k];
					deltaRSquaredArraySum += filteredDeltaRSquaredVector[k];
					currentTempVecLoc += 1;
				}
			}
			nDeltaR = nDeltaR + currentTempVecLoc;
			if ( (currentTempVecLoc > 0) && (deltaRSquaredArraySum > 0) ){
				(*twoGammaHat)(i,j+1) = (1.00/nDeltaR)*deltaRSquaredArraySum;
			}else if (j > 0){
				(*twoGammaHat)(i,j+1) = (*twoGammaHat)(i,j);
			}else{
				(*twoGammaHat)(i,j+1) = 0.00;
			}
		}

		for (j = 0; j < (const int)(*twoGammaHat).cols(); j++){
			if ((*twoGammaHat)(i,j) == 0 ){
				hatAngle.push_back(j);
			}
		}

		//E. Perfrom linear interpolation to get rid of zeros after first distance bin.
		if ( (hatAngle.size() > 1) && (hatAngle[hatAngle.size()-1] < ((const int)(*twoGammaHat).cols() - 1)) ){
			maxVal = (*twoGammaHat)(i,hatAngle[0]);

			for (j = 0; j < (const int)hatAngle.size(); j++){
				if((*twoGammaHat)(i,hatAngle[j]) > maxVal){
					maxVal = (*twoGammaHat)(i,hatAngle[j]);
				}
			}
			for (j = hatAngle[0]; j <= hatAngle[hatAngle.size() - 1]; j++){
				(*twoGammaHat)(i, j) = maxVal;
			}
		}
		filteredDistanceVector.clear();
		filteredDeltaRSquaredVector.clear();
		tempVec.clear();
		hatAngle.clear();
	}

	//************************************************************************************
	// 2. Clear all allocated memory from the function
	//************************************************************************************
	hatAngle.clear();
	tempVec.clear();
	filteredDistanceVector.clear();
	filteredDeltaRSquaredVector.clear();
}

dgrid rotMtx(double x)
{
	double y = PI * x / 180.00;
	dgrid R(2,2);
	R(0,0) = cos(y);
	R(0,1) = sin(y);
	R(1,0) = -sin(y);
	R(1,1) = cos(y);

	return R;
}