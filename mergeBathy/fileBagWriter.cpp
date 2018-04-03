/**********************************************************************
* CC0 License
**********************************************************************
* MergeBathy - Tool to combine one or more bathymetric data files onto a single input grid.
* Written in 2015 by Samantha J.Zambo(samantha.zambo@gmail.com) while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Todd Holland while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Nathaniel Plant while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Kevin Duvieilh while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Paul Elmore while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Will Avera while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by Brian Bourgeois while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by A.Louise Perkins while employed by the U.S.Naval Research Laboratory.
* Written in 2015 by David Lalejini while employed by the U.S.Naval Research Laboratory.
* To the extent possible under law, the author(s) and the U.S.Naval Research Laboratory have dedicated all copyright and related and neighboring rights to this software to the public domain worldwide.This software is distributed without any warranty.
* You should have received a copy of the CC0 Public Domain Dedication along with this software.If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
**********************************************************************/
#include "fileBagWriter.h"
#include "xmlWriter.h"

#include <stdio.h>
#include <stdlib.h>
#include "bag.h"
#include "hdf5.h"
//#include "bag_opt_surfaces.h"

#include <fstream>
#include <algorithm>
#include "grid.h"

#define DPRINT fprintf(stderr, "%s %s %d\n",__FILE__, __FUNCTION__,__LINE__); fflush(stderr);

#ifdef _MSC_VER
//Disable warnings since this is a Third-party file. -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( push )		//Save current warning state
	//#include "../WarningStates.h"	//Disable all Warnings!!!
	#pragma warning ( disable : 4996 )	//Deprecated call strncpy, sprintf, fopen
#endif
#endif

//#define XMLFile sample
#ifndef XMLFILE
#ifdef _MSC_VER
#define XMLFILE "sample.xml"
#else
#define XMLFILE "./sample.xml"
#endif
#endif

void CHECK_ERROR(bagError err);

void CHECK_ERROR(bagError err)
{
	if( err != BAG_SUCCESS )
	{
		char *errstr;
		if( bagGetErrorString( err, (u8**)&errstr ) == BAG_SUCCESS )
		{
			fprintf( stderr, "\nError create Bag: {%s}\n\n", errstr );
		}
	}
}

//************************************************************************************
// SUBROUTINE I.A: Function call for writing BAG files
//************************************************************************************
// added by L. Perkins, March 10, 2014 to write a bag file
// (HDF-5) as per JCGM-100 guidance
int writeBagFile(string &fName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei, int gridSpacingX, int gridSpacingY, int xSize, int ySize, map<string, int> *additionalOptions, char ZoneRef[4], vector<double> *z0, vector<double> *e0, vector<double> *zK, vector<double> *eK)
{
	//************************************************************************************
	//0. Prepare datasets for bags.	Flip the grids.
	//************************************************************************************
	dgrid flippedGrid_z = dgrid(ySize, xSize);
	dgrid flippedGrid_e = dgrid(ySize, xSize);
	dgrid flippedGrid_nei = dgrid(ySize, xSize);
	dgrid flippedGrid_rei = dgrid(ySize, xSize);
	dgrid flippedGrid_z0 = dgrid(ySize, xSize);
	dgrid flippedGrid_e0 = dgrid(ySize, xSize);
	dgrid flippedGrid_zK = dgrid(ySize, xSize);
	dgrid flippedGrid_eK = dgrid(ySize, xSize);

	int curLoc = 0;
	for (int i = 0; i < xSize; i++)
	{
		curLoc = i;
		for (int j = 0; j < ySize; j++)
		{

			flippedGrid_z(j, (xSize - 1)-i)   = ((*z)[curLoc]==NaN)?NULL_GENERIC:(*z)[curLoc]*-1;
			flippedGrid_e(j, (xSize - 1)-i)   = ((*e)[curLoc]==NaN)?NULL_GENERIC:(*e)[curLoc];
			flippedGrid_nei(j, (xSize - 1)-i) = ((*nei)[curLoc]==NaN)?NULL_GENERIC:(*nei)[curLoc];
			flippedGrid_rei(j, (xSize - 1)-i) = ((*rei)[curLoc]==NaN)?NULL_GENERIC:(*rei)[curLoc];
			flippedGrid_z0(j, (xSize - 1)-i)  = ((*z0)[curLoc]==NaN)?NULL_GENERIC:(*z0)[curLoc]*-1;
			flippedGrid_e0(j, (xSize - 1)-i)  = ((*e0)[curLoc]==NaN)?NULL_GENERIC:(*e0)[curLoc];
			flippedGrid_zK(j, (xSize - 1)-i)  = ((*zK)[curLoc]==NaN)?NULL_GENERIC:(*zK)[curLoc]*-1;
			flippedGrid_eK(j, (xSize - 1)-i)  = ((*eK)[curLoc]==NaN)?NULL_GENERIC:(*eK)[curLoc];
			curLoc += xSize;
		}
	}
	
	//************************************************************************************
	//0. Declare and initialize local variables and objects.	Open the file.
	//************************************************************************************
	const int SEP_SIZE = 3;
	int	nrows = xSize, ncols = ySize;

	float	*surf			 = new float[ncols];
	float	*uncert			 = new float[ncols];
	float	*nominal_depth	 = new float[ncols];
	bagVerticalCorrector *sep_depth = new bagVerticalCorrector[SEP_SIZE*SEP_SIZE];
	//bagVerticalCorrector *sep_depth_row = new bagVerticalCorrector[SEP_SIZE];

	float	surfRange[2];
	float	uncertRange[2];
	float	nominal_depthRange[2];
	float	sep_depthRange[2];

	/* Supporting variables */
	int	i;
	int	j;
	bagError err;
	bagData	data;
	bagHandle bagHandle;		 /* Primary Pointer to BAG object */
	char outFileName[256];
	char xmlFileName[256];	 

	memset (&data, 0, sizeof(data));

	int bagerr = 0;
	//bagerr = baginit(XMLFile, fName, (*z), (*e), (*x), (*y));					// invoke modified example from NAVO

	//************************************************************************************
	//I. Loop through the data and output to file. Check additional arguments for special output
	//************************************************************************************
	int DEFAULT_MSE=0;
	if (((*additionalOptions).find("-mse")->second == -1 && (*additionalOptions).find("-propUncert")->second == -1) && (*additionalOptions).find("-kalman")->second == -1)
		DEFAULT_MSE=1;
		
	string fNameTemp;
	int count=0;
	while(count!=3)
	{
		//A. Get z and e dataset and file name for bag
		cout<<fName<<endl;
		if(count==0)
		{
			if ((*additionalOptions).find("-mse")->second == 1 || DEFAULT_MSE)
			{
				*z=flippedGrid_z.vec();
				*e=flippedGrid_e.vec();
				fNameTemp = fName;
			}else count++;
		}
		if(count==1)
		{
			if ((*additionalOptions).find("-propUncert")->second == 1)
			{
				*z=flippedGrid_z0.vec();
				*e=flippedGrid_e0.vec();
				fNameTemp = fName.append("_P");
			}else count++;
		}
		if(count==2)
		{
			if ((*additionalOptions).find("-kalman")->second == 1)
			{
				*z=flippedGrid_zK.vec();
				*e=flippedGrid_eK.vec();
				fNameTemp = fName.append("_K");
			}else break;
		}
		count++;

		//B. write xml file for BAG
		string fNameXML = fNameTemp;
		int status = writeXMLFile(fNameXML, x, y, z, e, nei, rei, gridSpacingX, gridSpacingY, nrows, ncols, additionalOptions, ZoneRef, z0, e0, zK, eK);//xSize, ySize,

		if(status)
		{
			fprintf(stderr,"Error: writeXMLFile() unsuccessful.\n\n");
			return(-1);
		}

		strncpy( xmlFileName, fNameXML.c_str(), 255 );						 /* Store the XML fileName */
		strncpy( outFileName, (fNameTemp.append(".bag")).c_str(), 255 );	 /* Store the BAG fileName to write */
		remove ( outFileName );
		
		//use debug DPRINT from Stacy Johnson
		DPRINT fNameXML;

		//C. Find mins and maxes 
		surfRange[0] = (f32)*std::min_element(z->begin(), z->end());
		surfRange[1] = (f32)*std::max_element(z->begin(), z->end());
		uncertRange[0] = (f32)*std::min_element(e->begin(), e->end());
		uncertRange[1] = (f32)*std::max_element(e->begin(), e->end());
		nominal_depthRange[0] = 20L;
		nominal_depthRange[1] = (float)(20.0 + (float)((xSize-1)*(xSize-1)+xSize)/20.0);
		sep_depthRange[0] = (float)0.3333;
		sep_depthRange[1] = (float)103.333;

		printf( "Attempting to initialize a BAG!\n" );

		data.min_elevation = (f32) surfRange[0];
		data.max_elevation = (f32) surfRange[1];
		data.min_uncertainty = (f32) uncertRange[0];
		data.max_uncertainty = (f32) uncertRange[1];

		printf( "Creating the BAG, " );	

		//D. Initialize bag file.
		err = bagInitDefinitionFromFile(&data, xmlFileName);
		CHECK_ERROR(err);
		printf( "	ErrorCode for bagInitDefinition = %d\n", err );
		
		data.compressionLevel = 1;
		
		//E. Create bag file.
		err = bagFileCreate((u8*)outFileName, &data, &bagHandle);
		CHECK_ERROR(err);

		printf( "	finished initial creation of the bag, errCode = %d\n", err );
		printf( "Dims from XML r,c = [%d, %d]\n", 
				bagGetDataPointer(bagHandle)->def.nrows,
				bagGetDataPointer(bagHandle)->def.ncols );
	
		//************************************************************************************
		//II. Mandatory Data Set: Elevation.
		//************************************************************************************

		int rows = bagGetDataPointer(bagHandle)->def.nrows;
		int cols = bagGetDataPointer(bagHandle)->def.ncols;
		int cnt = 0;		

		// Load a row into surf and write
		for( i = 0; i < rows; i++ )
		{
			for (j = 0; j < cols; j++)
			{
				surf[j]	= (f32)(*z)[cnt]; 
				cnt++;
			}
			err = bagWriteRow( bagHandle, i, 0, cols-1, Elevation, (void *)surf );
		}
		CHECK_ERROR(err);

		//Update elevation surface
		err = bagUpdateSurface( bagHandle, Elevation );
		CHECK_ERROR(err);

		//************************************************************************************
		//II. Mandatory Data Set: Uncertainty.
		//************************************************************************************

		cnt = 0;
		// Load a row into uncert and write
		for( i = 0; i < rows; i++ )
		{
			for (j = 0; j < cols; j++)
			{
				uncert[j]	= (f32)((*e)[cnt]);
				cnt++;
			}
			err = bagWriteRow( bagHandle, i, 0, cols-1, Uncertainty, (void *)uncert );
		}
		CHECK_ERROR(err);

		//Update uncertainty surface
		err = bagUpdateSurface( bagHandle, Uncertainty );
		CHECK_ERROR(err);

		//************************************************************************************
		//II. Optional Data Set: Nominal_Depth
		//************************************************************************************

		/* adding optional nominal elevation dataset */
		bagGetDataPointer(bagHandle)->opt[Nominal_Elevation].datatype = H5T_NATIVE_FLOAT; 
		bagGetDataPointer(bagHandle)->opt[Nominal_Elevation].nrows = rows;
		bagGetDataPointer(bagHandle)->opt[Nominal_Elevation].ncols = cols;
			
		err = bagCreateOptionalDataset (bagHandle, bagGetDataPointer(bagHandle), Nominal_Elevation);
		CHECK_ERROR(err);

		bagAllocArray (bagHandle, 0, 0, rows-1, cols-1, Nominal_Elevation);
		cnt = 0;
		for( i = 0; i < rows; i++ )
		{
			for (j = 0; j < cols; j++)
			{
				nominal_depth[j] = (f32)((*z)[cnt]);//i*nrows+j]);
				cnt++;
			}
			err = bagWriteRow( bagHandle, i, 0, cols-1, Nominal_Elevation, (void *)nominal_depth );
		}
		CHECK_ERROR(err);

		err = bagUpdateSurface (bagHandle, Nominal_Elevation);
		CHECK_ERROR(err);

		bagFreeArray (bagHandle, Nominal_Elevation);

		//************************************************************************************
		//II. Optional Data Set: Surface_Correction 
		//************************************************************************************

		/* adding optional sep elevation dataset */
		bagGetDataPointer(bagHandle)->opt[Surface_Correction].nrows = SEP_SIZE;
		bagGetDataPointer(bagHandle)->opt[Surface_Correction].ncols = SEP_SIZE;

		err = bagCreateCorrectorDataset (bagHandle, bagGetDataPointer(bagHandle), 2, BAG_SURFACE_IRREGULARLY_SPACED);
		CHECK_ERROR(err);

		err = bagWriteCorrectorVerticalDatum (bagHandle, 1, (u8 *)"Test");
		CHECK_ERROR(err);

		err = bagWriteCorrectorVerticalDatum (bagHandle, 2, (u8 *)"Unknown");
		CHECK_ERROR(err);

		for( i = 0; i < SEP_SIZE; i++ )
		{
			/*for (j = 0; j < SEP_SIZE; j++)
			{
				sep_depth_row[j] = (f32)((*sep_depth)[i*nrows+j]);
			}
			err = bagWriteRow( bagHandle, i, 0, SEP_SIZE-1, Surface_Correction, (void *)&sep_depth_row );*/
			err = bagWriteRow( bagHandle, i, 0, SEP_SIZE-1, Surface_Correction, (void *)&sep_depth[i*SEP_SIZE] );
		}
		CHECK_ERROR(err);

		err = bagUpdateSurface (bagHandle, Surface_Correction);
		CHECK_ERROR(err);

		bagFreeArray (bagHandle, Surface_Correction);

		//************************************************************************************
		//II. Close the output file and return. 
		//************************************************************************************
		err = bagFileClose( bagHandle );	
		free (data.metadata);
		CHECK_ERROR(err);
		printf(	"Excellent... our bag is cooked!, Final ErrorCode = %d\n", err );
	}

	delete [] surf;
	delete [] uncert;
	delete [] nominal_depth;
	delete [] sep_depth;
	//delete [] sep_depth_row;

	return EXIT_SUCCESS;
}

 



















int readBagFile(string &fName, vector<double> *x, vector<double> *y, vector<double> *z, vector<double> *e, vector<double> *nei, vector<double> *rei,int gridSpacingX, int gridSpacingY, int xSize, int ySize, map<string, int> *additionalOptions)
{
	bagHandle hnd;
	u32 i, j, k;
	bagError stat;
	f32 *data = NULL;
	bagVerticalCorrector *vdata = NULL;
	char gisExportFileName[512], xmlFileName[512], bagFileName[512];
	char datum[55];
	int summaryOnly;
	FILE *oFile;
	s32 num_opt_datasets; 
	int opt_dataset_entities[10];
	bagData	xml_data;
	int NominalFound =0, sepFound = 0;

	//gisExportFileName[0] = '\0';
	xmlFileName[0]		 = '\0';
	//bagFileName[0]		 = '\0';
	datum[0]		 = '\0';
	summaryOnly = 0;	/* By default verbosly output information to standard out */

	 /* if ( argc < 2 )
	{
		printf("usage:	%s [-e arcgisAsciiFileName ] [-xml xmlFileName ] <bagFilename>\n", argv[0]);
		return EXIT_FAILURE;
	}*/

	int DEFAULT_MSE=0;
	if (((*additionalOptions).find("-mse")->second == -1 && (*additionalOptions).find("-propUncert")->second == -1) && (*additionalOptions).find("-kalman")->second == -1)
		DEFAULT_MSE=1;

	string fNameTemp;
	int count=0;
	while(count!=3)
	{
		//A. Get z and e dataset and file name for bag
		cout<<fName<<endl;
		if(count==0)
		{
			if ((*additionalOptions).find("-mse")->second == 1 || DEFAULT_MSE)
			{
				fNameTemp = fName;
			}else count++;
		}
		if(count==1)
		{
			if ((*additionalOptions).find("-propUncert")->second == 1)
			{
				fNameTemp = fName.append("_P");
			}else count++;
		}
		if(count==2)
		{
			if ((*additionalOptions).find("-kalman")->second == 1)
			{
				fNameTemp = fName.append("_K");
			}else break;
		}
		count++;

	string fName2=fNameTemp;
	summaryOnly = 1;		
	strncpy( gisExportFileName, (fNameTemp.append("Out")).c_str(), 255 );	 /* Store the Raster export fileName */
	//strncpy( xmlFileName, "sampleOut.xml", 255 );	 /* Store the XML fileName */
	strncpy( bagFileName, (fName2.append(".bag")).c_str(), 255 );	 /* Store the BAG fileName to write */

	//ProcessCommandInput( argc, argv, gisExportFileName, xmlFileName, bagFileName, &summaryOnly );

	if( bagFileName[0] == '\0' )
	{
		printf( "Error: No bag input file specified. Exiting.\n" );
		exit(-1);
	}
	else
	{
		printf( "Input BAG file: %s\n", bagFileName );
	}

	fprintf( stdout, "trying to open: {%s}...\n", bagFileName );
	fflush( stdout );

	stat = bagFileOpen (&hnd, BAG_OPEN_READ_WRITE, (u8*)bagFileName);
	if (stat != BAG_SUCCESS)
	{
		fprintf(stderr, "bag file unavailable ! %d\n", stat);
		fflush(stderr);
		return EXIT_FAILURE;
	}
	printf("bagFileOpen status = %d\n", stat);

	if (xmlFileName[0] != '\0')
		stat = bagInitDefinitionFromFile(&xml_data, xmlFileName);

	switch(xml_data.def.depthCorrectionType)
	{
		case Unknown_Correction:
		default:
			printf("Depth Correction type is unknown.\n");
			break;
		case True_Depth:
			printf("Depths in elevation dataset are TRUE\n");
			break;
		case Nominal_Depth_Meters:
			printf("Depths in elevation dataset are Nominal at 1500m\\s.\n");
			break;
		case Nominal_Depth_Feet:
			printf("Depths in elevation dataset are Nominal at 4800ft\\s.\n");
			break;
		case Corrected_Carters:
			printf("Depths in elevation dataset were corrected via Carter's tables\n");
			break;
		case Corrected_Matthews:
			printf("Depths in elevation dataset were corrected via Matthew's tables\n");
			break;
	}
		

	printf("BAG row/column extents: %dx%d\n", bagGetDataPointer(hnd)->def.nrows, bagGetDataPointer(hnd)->def.ncols);
	printf("BAG South-West Corner : %lf, %lf\n", bagGetDataPointer(hnd)->def.swCornerX, bagGetDataPointer(hnd)->def.swCornerY );
	printf("BAG Node Spacing		: %lf, %lf\n", bagGetDataPointer(hnd)->def.nodeSpacingX, bagGetDataPointer(hnd)->def.nodeSpacingY );
	printf("BAG Horiz. CoordSys	 : %s\n", bagGetDataPointer(hnd)->def.referenceSystem.horizontalReference);
	printf("BAG Vert.	CoordSys	 : %s\n", bagGetDataPointer(hnd)->def.referenceSystem.verticalReference);
	
	
/*
	stat = bagReadXMLStream( hnd );
	printf("status for bagReadDataset(Metadata) = %d\n", stat);
*/
	/* read the BAG file to determine whether or not any optional datasets exist */
	stat = bagGetOptDatasets(&hnd, &num_opt_datasets, opt_dataset_entities);

	if(num_opt_datasets == 0)
	{
		printf("\n Only mandatory datasets were found in the BAG.\n");
		if(xml_data.def.depthCorrectionType == Nominal_Depth_Meters || 
			xml_data.def.depthCorrectionType == Nominal_Depth_Feet)
				printf("\nThe elevation dataset however includes only Nominal depth not True depth\n\n");
	}
	else
	{
		printf("\nOptional datasets have been found\n");
	}

	for(i = 0; i < (u32) num_opt_datasets; i++)
	{
		if(opt_dataset_entities[i] == Nominal_Elevation) //(i.e. Nomimal_Depth_Meters, Nominal_Depth_Feet)SJZ
		{
			printf("\n Nominal data has been found in the bag and is contained in the optional dataset\n\n");
			NominalFound = 1;
		}
		if(opt_dataset_entities[i] == Surface_Correction) //depthCorrectionType (i.e. Corrected_Carters, Corrected_Matthews, Unknown_Correction)SJZ
		{
			printf("\n Vertical Datum Corrections have been found in the bag and are contained in an optional dataset\n\n");
			sepFound = 1;
		}
	}

	if( !summaryOnly )	/* Don't display if just the summary is requested */
	{
		printf("metadata = {%s}\n\n", bagGetDataPointer(hnd)->metadata);
	}

	if( xmlFileName[0] != '\0' )	 /* Check request to export the XML Metadata to a file */
	{
		if( (oFile = fopen( xmlFileName, "wb" )) == NULL )
		{
			fprintf( stderr, "ERROR: Opening file %s for XML export\n", xmlFileName );
			exit(-1);
		}
		fprintf( oFile, "%s", bagGetDataPointer(hnd)->metadata );
		fclose( oFile );
	}

	/* Read data using the ReadRow method */
	if( !summaryOnly )		/* Display heights and uncertainties	*/
	{
		/* if optional datasets are present read them in */
		if(num_opt_datasets > 0)
		{
			if(NominalFound) //(i.e. Nomimal_Depth_Meters, Nominal_Depth_Feet)SJZ
			{
				fprintf(stdout, "Nominal Elevation:=	{\n\t");
				fflush(stdout);

				bagGetOptDatasetInfo(&hnd, Nominal_Elevation);

								
				data = (f32*)calloc (bagGetDataPointer(hnd)->opt[Nominal_Elevation].ncols, sizeof(f32));
				for (i=0; i < bagGetDataPointer(hnd)->opt[Nominal_Elevation].nrows; i++)
				{
					bagReadRow (hnd, i, 0, bagGetDataPointer(hnd)->opt[Nominal_Elevation].ncols-1, Nominal_Elevation, data);
					for (j=0; j < bagGetDataPointer(hnd)->def.ncols; j++)
					{
						fprintf(stdout, "%0.3f\t", data[j]);
					}
					fprintf(stdout, "\n\t");
					fflush(stdout); 
				}
				
				free(data);
				bagFreeInfoOpt (hnd);
				fprintf(stdout, "}\n\t");
				fflush(stdout);
			}

			if(sepFound) //depthCorrectionType (i.e. Corrected_Carters, Corrected_Matthews, Unknown_Correction)SJZ
			{
				fprintf(stdout, "Vertical Datum Correctors :=	{\n\t");
				fflush(stdout);

				bagGetOptDatasetInfo(&hnd, Surface_Correction);

								
				vdata = (bagVerticalCorrector*)calloc (bagGetDataPointer(hnd)->opt[Surface_Correction].ncols, sizeof(bagVerticalCorrector));
				for (i=0; i < bagGetDataPointer(hnd)->opt[Surface_Correction].nrows; i++)
				{
					bagReadRow (hnd, i, 0, bagGetDataPointer(hnd)->opt[Surface_Correction].ncols-1, Surface_Correction, vdata);
					for (j=0; j < bagGetDataPointer(hnd)->def.ncols; j++)
					{
						u32 limit;
						bagGetNumSurfaceCorrectors	(hnd, &limit);
						for (k=0; k < limit; k++)
						{
							fprintf(stdout, "Z%d=%0.3lf ", k, vdata[j].z[k]); 
						}
						fprintf(stdout, "X=%0.3lf Y=%0.3lf\t", vdata[j].x, vdata[j].y);
					}
					fprintf(stdout, "\n\t");
					fflush(stdout); 

					fprintf(stdout, "ROW %d\n", i);
				}
				
				free(vdata);
				bagFreeInfoOpt (hnd);
				fprintf(stdout, "\t}\n");
				fflush(stdout);
			}
		}
	} 
	
	if( gisExportFileName[0] != '\0' )
	{
		/* This is similar to the code above but I choose to implement the export
			 separately to keep the code easier to read which is suitable for a sample
			 demonstration application. 
		*** */
		char outName[512];
		/* Export the height data to <file>_heights.asc */
		sprintf( outName, "%s_heights.asc", gisExportFileName );
		if( (oFile = fopen( outName, "wb" )) == NULL )
		{
			fprintf( stderr, "ERROR: Opening file %s for height data export\n", outName );
			exit(-1);
		}

		fprintf( oFile, "ncols %d\n", bagGetDataPointer(hnd)->def.ncols );
		fprintf( oFile, "nrows %d\n", bagGetDataPointer(hnd)->def.nrows );
		fprintf( oFile, "xllcenter %lf\n", bagGetDataPointer(hnd)->def.swCornerX );
		fprintf( oFile, "yllcenter %lf\n", bagGetDataPointer(hnd)->def.swCornerY);
		fprintf( oFile, "cellsize %lf\n", bagGetDataPointer(hnd)->def.nodeSpacingX );
		fprintf( oFile, "nodata_value %f\n", NULL_GENERIC );
		data = (f32*)calloc (bagGetDataPointer(hnd)->def.ncols, sizeof(f32));
		for (i=0; i < bagGetDataPointer(hnd)->def.nrows; i++)
		{
			bagReadRow (hnd, i, 0, bagGetDataPointer(hnd)->def.ncols-1, Elevation, data);
			for (j=0; j < bagGetDataPointer(hnd)->def.ncols; j++)
			{
				fprintf(oFile, "%0.3f ", data[j]);
			}
			fprintf(oFile, "\n");
		}
		free(data);
		fclose( oFile );

		/* Export the uncertainty data to <file>_heights.asc */
		sprintf( outName, "%s_uncrt.asc", gisExportFileName );
		if( (oFile = fopen( outName, "wb" )) == NULL )
		{
			fprintf( stderr, "ERROR: Opening file %s for height data export\n", outName );
			exit(-1);
		}
		fprintf( oFile, "ncols %d\n", bagGetDataPointer(hnd)->def.ncols );
		fprintf( oFile, "nrows %d\n", bagGetDataPointer(hnd)->def.nrows );
		fprintf( oFile, "xllcenter %lf\n", bagGetDataPointer(hnd)->def.swCornerX );
		fprintf( oFile, "yllcenter %lf\n", bagGetDataPointer(hnd)->def.swCornerY);
		fprintf( oFile, "cellsize %lf\n", bagGetDataPointer(hnd)->def.nodeSpacingX );
		fprintf( oFile, "nodata_value %f\n", NULL_GENERIC );
		data = (f32*)calloc (bagGetDataPointer(hnd)->def.ncols, sizeof(f32));
		for (i=0; i < bagGetDataPointer(hnd)->def.nrows; i++)
		{
			bagReadRow (hnd, i, 0, bagGetDataPointer(hnd)->def.ncols-1, Uncertainty, data);
			for (j=0; j < bagGetDataPointer(hnd)->def.ncols; j++)
			{
				fprintf(oFile, "%0.3f ", data[j]);
			}
			fprintf(oFile, "\n");
		}
		free(data);
		fclose( oFile );

		/* Export the uncertainty data to <file>_nominal.asc */
		if(num_opt_datasets > 0)
		{
			if(NominalFound)
			{
				sprintf( outName, "%s_nominal.asc", gisExportFileName );
				if( (oFile = fopen( outName, "wb" )) == NULL )
				{
					fprintf( stderr, "ERROR: Opening file %s for nominal data export\n", outName );
					exit(-1);
				}
				fprintf( oFile, "ncols %d\n", bagGetDataPointer(hnd)->def.ncols );
				fprintf( oFile, "nrows %d\n", bagGetDataPointer(hnd)->def.nrows );
				fprintf( oFile, "xllcenter %lf\n", bagGetDataPointer(hnd)->def.swCornerX );
				fprintf( oFile, "yllcenter %lf\n", bagGetDataPointer(hnd)->def.swCornerY);
				fprintf( oFile, "cellsize %lf\n", bagGetDataPointer(hnd)->def.nodeSpacingX );
				fprintf( oFile, "nodata_value %f\n", NULL_GENERIC );
				bagGetOptDatasetInfo(&hnd, Nominal_Elevation);

				data = (f32*)calloc (bagGetDataPointer(hnd)->def.ncols, sizeof(f32));
				for (i=0; i < bagGetDataPointer(hnd)->def.nrows; i++)//breaking here
				{
					bagReadRow (hnd, i, 0, bagGetDataPointer(hnd)->def.ncols-1, Nominal_Elevation, data);
					for (j=0; j < bagGetDataPointer(hnd)->def.ncols; j++)
					{
						fprintf(oFile, "%0.3f ", data[j]);
					}
					fprintf(oFile, "\n");
				}
				free(data);
				fclose( oFile );
			}
		}
	}

	/* uses ReadDataset - reads the data arrays in one shot - an alternate method to read the data */
	bagUpdateSurface (hnd, Elevation);
	bagUpdateSurface (hnd, Uncertainty);


	fprintf(stdout, "min_elv	%0.3f, max_elv	%0.3f\n", 
			bagGetDataPointer(hnd)->min_elevation, 
			bagGetDataPointer(hnd)->max_elevation);
	fflush(stdout);
	fprintf(stdout, "min_uncert %0.3f, max_uncert %0.3f\n", 
			bagGetDataPointer(hnd)->min_uncertainty, 
			bagGetDataPointer(hnd)->max_uncertainty);
	fflush(stdout);

	stat =	bagFileClose( hnd );
	printf("stat for bagFileClose = %d\n", stat);
	}
	return EXIT_SUCCESS;

} /* main */

//int ProcessCommandInput( int argc, char **argv, char *gisFile, char *xmlFile, char *bagFile, int *summaryOnly )
//{
//	int i;
//
//	i = 1;
//	while( i < argc )	 /*	Process until done */
//	{
//		if( strcmp( argv[i], "-e") == 0 )
//		{
//			if( i < (argc-1) )
//			{
//				strncpy( gisFile, argv[i+1], 511 );
//			}
//			else
//			{
//				fprintf( stderr, "ERROR: Missing filename for the -e option\n" );
//				exit(-1);
//			}
//			i += 2;
//		}
//		else if( strcmp( argv[i], "-xml" ) == 0 )
//		{
//			if( i < (argc-1) )
//			{
//				strncpy( xmlFile, argv[i+1], 511 );
//			}
//			else
//			{
//				fprintf( stderr, "ERROR: Missing filename for the -xml option\n" );
//				exit(-1);
//			}
//			i += 2;			
//		}
//		else if( strcmp( argv[i], "-summary" ) == 0 )
//		{
//			*summaryOnly = 1;	/* Set flag to indicate that only a summary of output is desired */
//			i++;
//		}
//		else if( strcmp( argv[i], "-h" ) == 0 )
//		{
//			printf("usage: bagread [-summary ] [-e arcgisAsciiBaseName ] [-xml xmlFileName ] <bagFilename>\n" );
//			printf("	If -summary is given the XML and data arrays will not be displayed on the screen. The\n" );
//			printf("-e option writes out two ArcGIS compatible ascii files for the height and uncertainty\n" );
//			printf(" data. The files use the supplied filename and add _height.asc and _uncrt.asc for the\n" );
//			printf("two exported.grids. The -xml option writes the XML metadata to the specified filename.\n\n" );
//			exit( 0 );
//		}
//		else
//		{
//			if( argv[i][0] == '-' )
//			{
//				fprintf( stderr, "ERROR: Unexpected command line parameter: %s\n", argv[i] );
//				exit( -1 );
//			}
//			strncpy( bagFile, argv[i], 511 );
//			i++;
//		}
//	}
//	return EXIT_SUCCESS;
//}

#ifdef _MSC_VER
//Restore warning state -SJZ
#if _DISABLE_3RDPARTY_WARNINGS
	#pragma warning ( pop )
#endif 
#endif 
