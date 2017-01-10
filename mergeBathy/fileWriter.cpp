#include "fileWriter.h"

//************************************************************************************
// SUBROUTINE I: Function call for writing flat files
//************************************************************************************
int writeFile(string &fName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, map<string, int> *additionalOptions)
{
	//************************************************************************************
	//0. Declare and initialize local variables and objects.  Open the file.
	//************************************************************************************
	ofstream fileOutM;
	string fileName = fName;
	fileName.append((".txt"));
	cout << "Writing Output Data to: " << fileName << endl;

	//1. Open the file
	fileOutM.open(fileName.c_str(), ios::out);
	fileOutM.precision(6);
	fileOutM.setf(std::ios::fixed, std::ios::floatfield);

	//add header
	fileOutM << "%lon\tlat\tdepth";
	if ((*additionalOptions).find("-noerr")->second == 0)
	{
		fileOutM << "\te";
	}
	if ((*additionalOptions).find("-nmsei")->second == 1)
	{
		fileOutM << "\tnEi";
	}
	if ((*additionalOptions).find("-msri")->second == 1)
	{
		fileOutM << "\trEi";
	}
	if ((*additionalOptions).find("-propUncert")->second == 1)
	{
		fileOutM << "\tz0\te0";
	}
	if ((*additionalOptions).find("-kalman")->second == 1)
	{
		fileOutM << "\tkz\tkvar";
	}
	fileOutM << "\n";

	//************************************************************************************
	//I. Loop through the data and output to file. Check additional arguments for special output
	//************************************************************************************
	for (int i = 0; i < (const int)(*x).size(); i++){
		fileOutM << (*x)[i] << "\t" << (*y)[i] << "\t" << (*z)[i];
		if ((*additionalOptions).find("-noerr")->second == 0)
		{
			fileOutM << "\t" << (*e)[i];
		}
		if ((*additionalOptions).find("-nmsei")->second == 1)
		{
			fileOutM << "\t" << (*nei)[i];
		}
		if ((*additionalOptions).find("-msri")->second == 1)
		{
			fileOutM << "\t" << (*rei)[i];
		}
		fileOutM << "\n";
	}

	//************************************************************************************
	//II. Close the output file and return.
	//************************************************************************************
	fileOutM.close();

	return SUCCESS;
}

//Print out all individual files.
int writeFiles(string &fName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, map<string, int> *additionalOptions, vector<double> *z0, vector<double> *e0, vector<double>* kz, vector<double>* kvar)
{
	//************************************************************************************
	//0. Declare and initialize local variables and objects.  Open the file.
	//************************************************************************************
	
	int cnt = 0;
	vector<double> *zTemp;
	vector<double> *eTemp;
	string fNameTemp;
	string header;
	while(cnt < 3)
	{
		//Add header
		header = "%lon\tlat";
		fNameTemp = fName;
		switch(cnt)
		{
			case 0: //MSE
				cnt++;
				if((*additionalOptions).find("-mse")->second == 1)
				{
					header.append("\tdepth");
					zTemp = z;
					if ((*additionalOptions).find("-noerr")->second == 0)
					{
						header.append("\te");
						eTemp = e;
					}
				} else continue;
				break;
			case 1: //Propagated Uncertainty
				cnt++;
				if((*additionalOptions).find("-propUncert")->second == 1)
				{
					fNameTemp.append("P");
					header.append("\tdepthP");
					zTemp = z0;
					if ((*additionalOptions).find("-noerr")->second == 0)
					{
						header.append("\teP");
						eTemp = e0;
					}
				}else continue;
				break;
			case 2: //Kalman
				cnt++;
				if((*additionalOptions).find("-kalman")->second == 1)
				{
					fNameTemp.append("K");
					header.append("\tdepthK");
					zTemp = kz;
					if ((*additionalOptions).find("-noerr")->second == 0)
					{
						header.append("\teK");
						eTemp = kvar;
					}
				}else continue;
				break;
		}
		if ((*additionalOptions).find("-nmsei")->second == 1)
		{
			header.append("\tnEi");
		}
		if ((*additionalOptions).find("-msri")->second == 1)
		{
			header.append("\trEi");
		}
		header.append("\n");

		ofstream fileOutM;
		string fileName = fNameTemp;
		fileName.append((".txt"));
		cout << "Writing Output Data to: " << fileName << endl;

		//1. Open the file
		fileOutM.open(fileName.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		fileOutM << header;
		//************************************************************************************
		//I. Loop through the data and output to file. Check additional arguments for special output
		//************************************************************************************
		for (int i = 0; i < (const int)(*x).size(); i++)
		{
			fileOutM << (*x)[i] << "\t" << (*y)[i] << "\t" << (*zTemp)[i];
			if ((*additionalOptions).find("-noerr")->second == 0)
			{
				fileOutM << "\t" << (*eTemp)[i];
			}
			if ((*additionalOptions).find("-nmsei")->second == 1)
			{
				fileOutM << "\t" << (*nei)[i];
			}
			if ((*additionalOptions).find("-msri")->second == 1)
			{
				fileOutM << "\t" << (*rei)[i];
			}
			fileOutM << "\n";
		}

		//************************************************************************************
		//II. Close the output file and return.
		//************************************************************************************
		fileOutM.close();
	}
	return SUCCESS;
}

//************************************************************************************
// SUBROUTINE II: Function call for writing flat file MSE
// followed by Kalman
//************************************************************************************
int writeFile(string &fName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, map<string, int> *additionalOptions, vector<double>* kz, vector<double>* kvar)
{
	//************************************************************************************
	//0. Declare and initialize local variables and objects.  Open the file.
	//************************************************************************************
	ofstream fileOutM;
	string fileName = fName;
	fileName.append((".txt"));
	cout << "Writing Output Data to: " << fileName << endl;

	//1. Open the file
	fileOutM.open(fileName.c_str(), ios::out);
	fileOutM.precision(6);
	fileOutM.setf(std::ios::fixed, std::ios::floatfield);

	//************************************************************************************
	//I. Loop through the data and output to file. Check additional arguments for special output
	//************************************************************************************
	for (int i = 0; i < (const int)(*x).size(); i++){
		fileOutM << (*x)[i] << "\t" << (*y)[i] << "\t" << (*z)[i];
		if ((*additionalOptions).find("-noerr")->second == 0)
		{
			fileOutM << "\t" << (*e)[i];
		}
		if ((*additionalOptions).find("-nmsei")->second == 1)
		{
			fileOutM << "\t" << (*nei)[i];
		}
		if ((*additionalOptions).find("-msri")->second == 1)
		{
			fileOutM << "\t" << (*rei)[i];
		}
		fileOutM << "\t" << (*kz)[i];
		if ((*additionalOptions).find("-noerr")->second == 0)
		{
			fileOutM << "\t" << (*kvar)[i];
		}
		fileOutM << "\n";
	}

	//************************************************************************************
	//II. Close the output file and return.
	//************************************************************************************
	fileOutM.close();

	return SUCCESS;
}

//************************************************************************************
// SUBROUTINE III: Function call for writing flat files that
// match matlab output format for debugging!!! 
//************************************************************************************
int writeFileMatch(string &fName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, map<string, int> *additionalOptions, vector<double>* kz, vector<double>* kvar)
{
	string header = "%latitude,longitude,depth_m,depth_ft,X,Y,Nerr,Rerr,Error,KZ,KVAR";
	if ((*additionalOptions).find("-noerr")->second == 1)
	{
		cout << "-noerr option found; print out will be the same as MSE" << endl;
	}

	//************************************************************************************
	//0. Declare and initialize local variables and objects.  Open the file.
	//************************************************************************************
	ofstream fileOutM;
	string fileName = fName;
	fileName.append((".txt"));
	cout << "Writing Output Data to: " << fileName << endl;

	//1. Open the file
	fileOutM.open(fileName.c_str(), ios::out);
	fileOutM.setf(std::ios::fixed, std::ios::floatfield);
	fileOutM << header << "\n";
	vector<double> zft=vector<double>((*z).size(),0.00);
	vector<double> xx=vector<double>((*z).size(),0.00);
	vector<double> yy=vector<double>((*z).size(),0.00);
	//************************************************************************************
	//I. Loop through the data and output to file. Check additional arguments for special output
	//************************************************************************************
	for (int i = 0; i < (const int)(*x).size(); i++){
		fileOutM << setprecision(8) << (*y)[i] <<  "\t" << setprecision(8) << (*x)[i] << "\t" << setprecision(4) << (*z)[i] << "\t" << setprecision(4) << (zft)[i] << "\t" << setprecision(4) << (xx)[i] << "\t" << setprecision(4) << (yy)[i]; 

		if ((*additionalOptions).find("-nmsei")->second == 1)
		{
			fileOutM << "\t" << setprecision(4) << (*nei)[i];
		}
		if ((*additionalOptions).find("-msri")->second == 1)
		{
			fileOutM << "\t" << setprecision(4) << (*rei)[i];
		}
		if ((*additionalOptions).find("-noerr")->second == 0)
		{
			fileOutM << "\t" << setprecision(4) << (*e)[i];
			fileOutM << "\t" << setprecision(4) << (*kz)[i];
			fileOutM << "\t" << setprecision(4) << (*kvar)[i];
		}
		fileOutM << "\n";
	}

	//************************************************************************************
	//II. Close the output file and return.
	//************************************************************************************
	fileOutM.close();

	return SUCCESS;
}

//************************************************************************************
// SUBROUTINE IV: Function call for writing ARC ASCII Raster
// files
//************************************************************************************
int writeRasterFile(string &fileName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, double spacingX, double spacingY, double xSizeTemp, double ySizeTemp, map<string, int> *additionalOptions, char ZoneRef[4], vector<double> *z0, vector<double> *e0, vector<double> *zK, vector<double> *eK)
{
	uint xSize = (uint)xSizeTemp;
	uint ySize = (uint)ySizeTemp;

	//************************************************************************************
	//0. Declare and initialize local variables and objects.  Open the file and initialize the headers.
	//************************************************************************************
	ofstream fileOutM;
	string fileNameTemp;
	uint curLoc = 0;

	//1. Open the file
	fileNameTemp = fileName;
	fileNameTemp.append("_depth");

	//2. Create projection file
	rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

	//3. Finish creating the main output file
	fileNameTemp.append((".asc"));
	cout << "Writing Output Data to: " << fileNameTemp << endl;

	fileOutM.open(fileNameTemp.c_str(), ios::out);
	fileOutM.precision(6);
	fileOutM.setf(std::ios::fixed, std::ios::floatfield);

	//4. Build the headers
	fileOutM << "ncols " << ySize << endl;  
	fileOutM << "nrows " << xSize << endl;
	fileOutM << "xllcenter " << (*x)[0] << endl;
	fileOutM << "yllcenter " << (*y)[0] << endl;
	fileOutM << "cellsize " << spacingX << endl;
	fileOutM << "nodata_value 9999999" << endl;

	dgrid flippedGrid = dgrid(ySize, xSize);
	for (uint i = 0; i < xSize; i++)
	{
		curLoc = i;
		for (uint j = 0; j < ySize; j++)
		{
			flippedGrid(j, (xSize - 1)-i) = (*z)[curLoc];
			curLoc += xSize;
		}
	}

	//************************************************************************************
	//I. Loop through the data and output to file.
	//************************************************************************************
	//A. Output Depth
	curLoc = 0;
	for (uint i = 0; i < xSize; i++)
	{
		for (uint j = 0; j < ySize; j++)
		{
			fileOutM << flippedGrid(j,i) << " ";
		}
		fileOutM << endl;
	}
	fileOutM.close();

	//B. Check if EI is output
	if ((*additionalOptions).find("-noerr")->second == 0)
	{
		fileNameTemp = fileName;
		fileNameTemp.append("_ei");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;
		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*e)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	}

	//C. Check if NEI is output
	if ((*additionalOptions).find("-nmsei")->second == 1)
	{
		fileNameTemp = fileName;
		fileNameTemp.append("_nei");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;
		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*nei)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	}
	
	//D. Check if REI is output
	if ((*additionalOptions).find("-msri")->second == 1)
	{
		fileNameTemp = fileName;
		fileNameTemp.append("_rei");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;

		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*rei)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	}

	//E. Check if propagated uncertainty is output
	if ((*additionalOptions).find("-propUncert")->second == 1)
	{
		//1. Propagated Uncertainty Depth
		fileNameTemp = fileName;
		fileNameTemp.append("_depthP");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;
		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*z0)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();

		//1. Propagated Uncertainty Error
		fileNameTemp = fileName;
		fileNameTemp.append("_eiP");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;
		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*e0)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	}

	//F. Check if Kalman is output
	if ((*additionalOptions).find("-kalman")->second == 1)
	{
		//1. Kalman Depth
		fileNameTemp = fileName;
		fileNameTemp.append("_depthK");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;
		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*zK)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	
		//1. Kalman Error
		fileNameTemp = fileName;
		fileNameTemp.append("_eiK");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;
		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*eK)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	}
	flippedGrid.clear();
	return SUCCESS;
}

int writeRasterFile(string &fileName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, double spacingX, double spacingY, double xSizeTemp, double ySizeTemp, map<string, int> *additionalOptions, char ZoneRef[4])
{
	uint xSize = (uint)xSizeTemp;
	uint ySize = (uint)ySizeTemp;

	//************************************************************************************
	//0. Declare and initialize local variables and objects.  Open the file and initialize the headers.
	//************************************************************************************
	ofstream fileOutM;
	string fileNameTemp;
	uint curLoc = 0;

	//1. Open the file
	fileNameTemp = fileName;
	fileNameTemp.append("_depth");

	//2. Create projection file
	rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

	//3. Finish creating the main output file
	fileNameTemp.append((".asc"));
	cout << "Writing Output Data to: " << fileNameTemp << endl;

	fileOutM.open(fileNameTemp.c_str(), ios::out);
	fileOutM.precision(6);
	fileOutM.setf(std::ios::fixed, std::ios::floatfield);

	//4. Build the headers
	fileOutM << "ncols " << ySize << endl;  
	fileOutM << "nrows " << xSize << endl;
	fileOutM << "xllcenter " << (*x)[0] << endl;
	fileOutM << "yllcenter " << (*y)[0] << endl;
	fileOutM << "cellsize " << spacingX << endl;
	fileOutM << "nodata_value 9999999" << endl;

	dgrid flippedGrid = dgrid(ySize, xSize);
	for (uint i = 0; i < xSize; i++)
	{
		curLoc = i;
		for (uint j = 0; j < ySize; j++)
		{
			flippedGrid(j, (xSize - 1)-i) = (*z)[curLoc];
			curLoc += xSize;
		}
	}

	//************************************************************************************
	//I. Loop through the data and output to file.
	//************************************************************************************
	//A. Output Depth
	curLoc = 0;
	for (uint i = 0; i < xSize; i++)
	{
		for (uint j = 0; j < ySize; j++)
		{
			fileOutM << flippedGrid(j,i) << " ";
		}
		fileOutM << endl;
	}
	fileOutM.close();

	//B. Check if EI is output
	if ((*additionalOptions).find("-noerr")->second == 0)
	{
		fileNameTemp = fileName;
		fileNameTemp.append("_ei");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;
		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*e)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	}

	//C. Check if NEI is output
	if ((*additionalOptions).find("-nmsei")->second == 1)
	{
		fileNameTemp = fileName;
		fileNameTemp.append("_nei");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;
		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*nei)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	}
	//D. Check if REI is output
	if ((*additionalOptions).find("-msri")->second == 1)
	{
		fileNameTemp = fileName;
		fileNameTemp.append("_rei");
		//2. Create projection file
		rasterWGS84_SpheroidCreator(fileNameTemp, ZoneRef);

		//3. Finish creating the main output file
		fileNameTemp.append((".asc"));
		cout << "Writing Output Data to: " << fileNameTemp << endl;

		fileOutM.open(fileNameTemp.c_str(), ios::out);
		fileOutM.precision(6);
		fileOutM.setf(std::ios::fixed, std::ios::floatfield);

		//4. Build the headers
		fileOutM << "ncols " << ySize << endl;
		fileOutM << "nrows " << xSize << endl;
		fileOutM << "xllcenter " << (*x)[0] << endl;
		fileOutM << "yllcenter " << (*y)[0] << endl;
		fileOutM << "cellsize " << spacingX << endl;
		fileOutM << "nodata_value 9999999" << endl;

		curLoc = 0;
		for (uint i = 0; i < xSize; i++)
		{
			curLoc = i;
			for (uint j = 0; j < ySize; j++)
			{
				flippedGrid(j,(xSize - 1)-i) = (*rei)[curLoc];
				curLoc += xSize;
			}
		}
		for (uint i = 0; i < xSize; i++)
		{
			for (uint j = 0; j < ySize; j++)
			{
				fileOutM << flippedGrid(j,i) << " ";
			}
			fileOutM << endl;
		}
		fileOutM.close();
	}

	flippedGrid.clear();
	return SUCCESS;
}

int rasterWGS84_SpheroidCreator(string &fileName, char ZoneRef[4])
{
	//************************************************************************************
	//0. Declare and initialize local variables and objects.  Open the file and initialize the headers.
	//************************************************************************************
	ofstream fileOutM;
	string fileNameTemp;
	string projectionData;

	//1. Open the file
	fileNameTemp = fileName;
	fileNameTemp.append((".prj"));
	cout << "Writing UTM Projection Output Data to: " << fileNameTemp << endl;
	fileOutM.open(fileNameTemp.c_str(), ios::out);
	fileOutM.precision(6);
	fileOutM.setf(std::ios::fixed, std::ios::floatfield);

	string projectionFileNorth[60];
	string projectionFileSouth[60];

	//2. Make all of the strings
	projectionFileNorth[0] = "PROJCS[\"WGS_1984_UTM_Zone_1N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-177.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32601]]";
	projectionFileNorth[1] = "PROJCS[\"WGS_1984_UTM_Zone_2N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-171.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32602]]";
	projectionFileNorth[2] = "PROJCS[\"WGS_1984_UTM_Zone_3N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-165.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32603]]";
	projectionFileNorth[3] = "PROJCS[\"WGS_1984_UTM_Zone_4N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-159.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32604]]";
	projectionFileNorth[4] = "PROJCS[\"WGS_1984_UTM_Zone_5N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-153.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32605]]";
	projectionFileNorth[5] = "PROJCS[\"WGS_1984_UTM_Zone_6N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-147.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32606]]";
	projectionFileNorth[6] = "PROJCS[\"WGS_1984_UTM_Zone_7N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-141.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32607]]";
	projectionFileNorth[7] = "PROJCS[\"WGS_1984_UTM_Zone_8N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-135.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32608]]";
	projectionFileNorth[8] = "PROJCS[\"WGS_1984_UTM_Zone_9N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-129.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32609]]";
	projectionFileNorth[9] = "PROJCS[\"WGS_1984_UTM_Zone_10N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-123.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32610]]";
	projectionFileNorth[10] = "PROJCS[\"WGS_1984_UTM_Zone_11N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-117.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32611]]";
	projectionFileNorth[11] = "PROJCS[\"WGS_1984_UTM_Zone_12N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-111.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32612]]";
	projectionFileNorth[12] = "PROJCS[\"WGS_1984_UTM_Zone_13N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-105.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32613]]";
	projectionFileNorth[13] = "PROJCS[\"WGS_1984_UTM_Zone_14N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-99.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32614]]";
	projectionFileNorth[14] = "PROJCS[\"WGS_1984_UTM_Zone_15N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-93.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32615]]";
	projectionFileNorth[15] = "PROJCS[\"WGS_1984_UTM_Zone_16N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-87.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32616]]";
	projectionFileNorth[16] = "PROJCS[\"WGS_1984_UTM_Zone_17N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-81.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32617]]";
	projectionFileNorth[17] = "PROJCS[\"WGS_1984_UTM_Zone_18N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-75.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32618]]";
	projectionFileNorth[18] = "PROJCS[\"WGS_1984_UTM_Zone_19N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-69.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32619]]";
	projectionFileNorth[19] = "PROJCS[\"WGS_1984_UTM_Zone_20N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-63.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32620]]";
	projectionFileNorth[20] = "PROJCS[\"WGS_1984_UTM_Zone_21N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-57.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32621]]";
	projectionFileNorth[21] = "PROJCS[\"WGS_1984_UTM_Zone_22N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-51.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32622]]";
	projectionFileNorth[22] = "PROJCS[\"WGS_1984_UTM_Zone_23N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-45.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32623]]";
	projectionFileNorth[23] = "PROJCS[\"WGS_1984_UTM_Zone_24N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-39.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32624]]";
	projectionFileNorth[24] = "PROJCS[\"WGS_1984_UTM_Zone_25N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-33.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32625]]";
	projectionFileNorth[25] = "PROJCS[\"WGS_1984_UTM_Zone_26N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-27.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32626]]";
	projectionFileNorth[26] = "PROJCS[\"WGS_1984_UTM_Zone_27N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-21.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32627]]";
	projectionFileNorth[27] = "PROJCS[\"WGS_1984_UTM_Zone_28N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-15.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32628]]";
	projectionFileNorth[28] = "PROJCS[\"WGS_1984_UTM_Zone_29N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-9.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32629]]";
	projectionFileNorth[29] = "PROJCS[\"WGS_1984_UTM_Zone_30N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",-3.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32630]]";
	projectionFileNorth[30] = "PROJCS[\"WGS_1984_UTM_Zone_31N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",3.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32631]]";
	projectionFileNorth[31] = "PROJCS[\"WGS_1984_UTM_Zone_32N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",9.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32632]]";
	projectionFileNorth[32] = "PROJCS[\"WGS_1984_UTM_Zone_33N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",15.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32633]]";
	projectionFileNorth[33] = "PROJCS[\"WGS_1984_UTM_Zone_34N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",21.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32634]]";
	projectionFileNorth[34] = "PROJCS[\"WGS_1984_UTM_Zone_35N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",27.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32635]]";
	projectionFileNorth[35] = "PROJCS[\"WGS_1984_UTM_Zone_36N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",33.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32636]]";
	projectionFileNorth[36] = "PROJCS[\"WGS_1984_UTM_Zone_37N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",39.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32637]]";
	projectionFileNorth[37] = "PROJCS[\"WGS_1984_UTM_Zone_38N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",45.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32638]]";
	projectionFileNorth[38] = "PROJCS[\"WGS_1984_UTM_Zone_39N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",51.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32639]]";
	projectionFileNorth[39] = "PROJCS[\"WGS_1984_UTM_Zone_40N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",57.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32640]]";
	projectionFileNorth[40] = "PROJCS[\"WGS_1984_UTM_Zone_41N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",63.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32641]]";
	projectionFileNorth[41] = "PROJCS[\"WGS_1984_UTM_Zone_42N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",69.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32642]]";
	projectionFileNorth[42] = "PROJCS[\"WGS_1984_UTM_Zone_43N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",75.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32643]]";
	projectionFileNorth[43] = "PROJCS[\"WGS_1984_UTM_Zone_44N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",81.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32644]]";
	projectionFileNorth[44] = "PROJCS[\"WGS_1984_UTM_Zone_45N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",87.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32645]]";
	projectionFileNorth[45] = "PROJCS[\"WGS_1984_UTM_Zone_46N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",93.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32646]]";
	projectionFileNorth[46] = "PROJCS[\"WGS_1984_UTM_Zone_47N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",99.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32647]]";
	projectionFileNorth[47] = "PROJCS[\"WGS_1984_UTM_Zone_48N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",105.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32648]]";
	projectionFileNorth[48] = "PROJCS[\"WGS_1984_UTM_Zone_49N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",111.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32649]]";
	projectionFileNorth[49] = "PROJCS[\"WGS_1984_UTM_Zone_50N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",117.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32650]]";
	projectionFileNorth[50] = "PROJCS[\"WGS_1984_UTM_Zone_51N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",123.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32651]]";
	projectionFileNorth[51] = "PROJCS[\"WGS_1984_UTM_Zone_52N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",129.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32652]]";
	projectionFileNorth[52] = "PROJCS[\"WGS_1984_UTM_Zone_53N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",135.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32653]]";
	projectionFileNorth[53] = "PROJCS[\"WGS_1984_UTM_Zone_54N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",141.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32654]]";
	projectionFileNorth[54] = "PROJCS[\"WGS_1984_UTM_Zone_55N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",147.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32655]]";
	projectionFileNorth[55] = "PROJCS[\"WGS_1984_UTM_Zone_56N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",153.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32656]]";
	projectionFileNorth[56] = "PROJCS[\"WGS_1984_UTM_Zone_57N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",159.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32657]]";
	projectionFileNorth[57] = "PROJCS[\"WGS_1984_UTM_Zone_58N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",165.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32658]]";
	projectionFileNorth[58] = "PROJCS[\"WGS_1984_UTM_Zone_59N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",171.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32659]]";
	projectionFileNorth[59] = "PROJCS[\"WGS_1984_UTM_Zone_60N\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",0.0],PARAMETER[\"Central_Meridian\",177.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32660]]";

	projectionFileSouth[0] = "PROJCS[\"WGS_1984_UTM_Zone_1S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-177.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32701]]";
	projectionFileSouth[1] = "PROJCS[\"WGS_1984_UTM_Zone_2S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-171.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32702]]";
	projectionFileSouth[2] = "PROJCS[\"WGS_1984_UTM_Zone_3S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-165.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32703]]";
	projectionFileSouth[3] = "PROJCS[\"WGS_1984_UTM_Zone_4S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-159.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32704]]";
	projectionFileSouth[4] = "PROJCS[\"WGS_1984_UTM_Zone_5S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-153.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32705]]";
	projectionFileSouth[5] = "PROJCS[\"WGS_1984_UTM_Zone_6S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-147.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32706]]";
	projectionFileSouth[6] = "PROJCS[\"WGS_1984_UTM_Zone_7S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-141.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32707]]";
	projectionFileSouth[7] = "PROJCS[\"WGS_1984_UTM_Zone_8S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-135.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32708]]";
	projectionFileSouth[8] = "PROJCS[\"WGS_1984_UTM_Zone_9S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-129.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32709]]";
	projectionFileSouth[9] = "PROJCS[\"WGS_1984_UTM_Zone_10S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-123.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32710]]";
	projectionFileSouth[10] = "PROJCS[\"WGS_1984_UTM_Zone_11S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-117.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32711]]";
	projectionFileSouth[11] = "PROJCS[\"WGS_1984_UTM_Zone_12S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-111.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32712]]";
	projectionFileSouth[12] = "PROJCS[\"WGS_1984_UTM_Zone_13S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-105.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32713]]";
	projectionFileSouth[13] = "PROJCS[\"WGS_1984_UTM_Zone_14S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-99.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32714]]";
	projectionFileSouth[14] = "PROJCS[\"WGS_1984_UTM_Zone_15S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-93.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32715]]";
	projectionFileSouth[15] = "PROJCS[\"WGS_1984_UTM_Zone_16S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-87.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32716]]";
	projectionFileSouth[16] = "PROJCS[\"WGS_1984_UTM_Zone_17S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-81.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32717]]";
	projectionFileSouth[17] = "PROJCS[\"WGS_1984_UTM_Zone_18S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-75.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32718]]";
	projectionFileSouth[18] = "PROJCS[\"WGS_1984_UTM_Zone_19S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-69.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32719]]";
	projectionFileSouth[19] = "PROJCS[\"WGS_1984_UTM_Zone_20S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-63.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32720]]";
	projectionFileSouth[20] = "PROJCS[\"WGS_1984_UTM_Zone_21S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-57.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32721]]";
	projectionFileSouth[21] = "PROJCS[\"WGS_1984_UTM_Zone_22S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-51.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32722]]";
	projectionFileSouth[22] = "PROJCS[\"WGS_1984_UTM_Zone_23S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-45.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32723]]";
	projectionFileSouth[23] = "PROJCS[\"WGS_1984_UTM_Zone_24S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-39.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32724]]";
	projectionFileSouth[24] = "PROJCS[\"WGS_1984_UTM_Zone_25S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-33.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32725]]";
	projectionFileSouth[25] = "PROJCS[\"WGS_1984_UTM_Zone_26S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-27.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32726]]";
	projectionFileSouth[26] = "PROJCS[\"WGS_1984_UTM_Zone_27S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-21.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32727]]";
	projectionFileSouth[27] = "PROJCS[\"WGS_1984_UTM_Zone_28S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-15.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32728]]";
	projectionFileSouth[28] = "PROJCS[\"WGS_1984_UTM_Zone_29S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-9.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32729]]";
	projectionFileSouth[29] = "PROJCS[\"WGS_1984_UTM_Zone_30S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",-3.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32730]]";
	projectionFileSouth[30] = "PROJCS[\"WGS_1984_UTM_Zone_31S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",3.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32731]]";
	projectionFileSouth[31] = "PROJCS[\"WGS_1984_UTM_Zone_32S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",9.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32732]]";
	projectionFileSouth[32] = "PROJCS[\"WGS_1984_UTM_Zone_33S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",15.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32733]]";
	projectionFileSouth[33] = "PROJCS[\"WGS_1984_UTM_Zone_34S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",21.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32734]]";
	projectionFileSouth[34] = "PROJCS[\"WGS_1984_UTM_Zone_35S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",27.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32735]]";
	projectionFileSouth[35] = "PROJCS[\"WGS_1984_UTM_Zone_36S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",33.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32736]]";
	projectionFileSouth[36] = "PROJCS[\"WGS_1984_UTM_Zone_37S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",39.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32737]]";
	projectionFileSouth[37] = "PROJCS[\"WGS_1984_UTM_Zone_38S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",45.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32738]]";
	projectionFileSouth[38] = "PROJCS[\"WGS_1984_UTM_Zone_39S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",51.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32739]]";
	projectionFileSouth[39] = "PROJCS[\"WGS_1984_UTM_Zone_40S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",57.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32740]]";
	projectionFileSouth[40] = "PROJCS[\"WGS_1984_UTM_Zone_41S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",63.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32741]]";
	projectionFileSouth[41] = "PROJCS[\"WGS_1984_UTM_Zone_42S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",69.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32742]]";
	projectionFileSouth[42] = "PROJCS[\"WGS_1984_UTM_Zone_43S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",75.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32743]]";
	projectionFileSouth[43] = "PROJCS[\"WGS_1984_UTM_Zone_44S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",81.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32744]]";
	projectionFileSouth[44] = "PROJCS[\"WGS_1984_UTM_Zone_45S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",87.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32745]]";
	projectionFileSouth[45] = "PROJCS[\"WGS_1984_UTM_Zone_46S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",93.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32746]]";
	projectionFileSouth[46] = "PROJCS[\"WGS_1984_UTM_Zone_47S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",99.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32747]]";
	projectionFileSouth[47] = "PROJCS[\"WGS_1984_UTM_Zone_48S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",105.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32748]]";
	projectionFileSouth[48] = "PROJCS[\"WGS_1984_UTM_Zone_49S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",111.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32749]]";
	projectionFileSouth[49] = "PROJCS[\"WGS_1984_UTM_Zone_50S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",117.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32750]]";
	projectionFileSouth[50] = "PROJCS[\"WGS_1984_UTM_Zone_51S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",123.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32751]]";
	projectionFileSouth[51] = "PROJCS[\"WGS_1984_UTM_Zone_52S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",129.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32752]]";
	projectionFileSouth[52] = "PROJCS[\"WGS_1984_UTM_Zone_53S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",135.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32753]]";
	projectionFileSouth[53] = "PROJCS[\"WGS_1984_UTM_Zone_54S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",141.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32754]]";
	projectionFileSouth[54] = "PROJCS[\"WGS_1984_UTM_Zone_55S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",147.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32755]]";
	projectionFileSouth[55] = "PROJCS[\"WGS_1984_UTM_Zone_56S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",153.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32756]]";
	projectionFileSouth[56] = "PROJCS[\"WGS_1984_UTM_Zone_57S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",159.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32757]]";
	projectionFileSouth[57] = "PROJCS[\"WGS_1984_UTM_Zone_58S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",165.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32758]]";
	projectionFileSouth[58] = "PROJCS[\"WGS_1984_UTM_Zone_59S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",171.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32759]]";
	projectionFileSouth[59] = "PROJCS[\"WGS_1984_UTM_Zone_60S\",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137.0,298.257223563]],PRIMEM[\"Greenwich\",0.0],UNIT[\"Degree\",0.0174532925199433]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"False_Easting\",500000.0],PARAMETER[\"False_Northing\",10000000.0],PARAMETER[\"Central_Meridian\",177.0],PARAMETER[\"Scale_Factor\",0.9996],PARAMETER[\"Latitude_Of_Origin\",0.0],UNIT[\"Meter\",1.0],AUTHORITY[\"EPSG\",32760]]";

	string hemisphere;
	string location;
	string temp;

	if (isalpha(ZoneRef[2]))
	{
		location = ZoneRef[0];
		temp = ZoneRef[1];
		location.append(temp);
		hemisphere = ZoneRef[2];
	}else
	{
		location = ZoneRef[0];
		hemisphere = ZoneRef[1];
	}

	int index = atoi(location.c_str());

	if (hemisphere >= "N")
	{
		fileOutM << projectionFileNorth[index - 1] << endl;
	}else
	{
		fileOutM << projectionFileSouth[index - 1] << endl;
	}

	fileOutM.close();

	return SUCCESS;
}