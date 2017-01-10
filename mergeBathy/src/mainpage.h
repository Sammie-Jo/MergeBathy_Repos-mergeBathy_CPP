 /** @mainpage
 * @section version Version
 * 3.5.0
 * @section version_date Version Date
 * 04 August 2011
 * @section intro Introduction
 * Welcome to the mergeBathy documentation main page!
 * @n The mergeBathy utility is used to smooth bathymetric surfaces from mulitilple data sets according to a user defined smoothing scale.
 * @n The following documentation describes the calling conventions for mergeBathy and the arguments necessary for its operation.  
 * @n Each function of mergeBathy is also documented in the "Files" tab for reference by software developers.
 * @n
 * @section args Arguments
 * The main function takes 8 mandatory command line arguments and is capable of having additional arguments for separate processing.
 * @n Arguments encapsulated by <> are mandatory arguments. 
 * @n Arguments encapsulated by [] are optional arguments.
 * @n Some optional arguments have mandatory arguments that must be passed after the optional argument is invoked.  These mandatory arguments are indented below the optional argument that requries them and are encapsulated by <>.
 * 
 * @subsection mandatory_args Mandatory Arguments
 * <table border="1">
 * <tr><td><b><output_depth_file></b></td><td>A text file consisting of one line for each output point at which a valid depth was computed.  Each line consists of four columns, corresponding to Longitude, Latitude, Depth, and Uncertainty (the mean square interpolation error estimate) respectively.  Additional columns can be present depending on other additional command line options provided at runtime.</td></tr> 
 * <tr><td><b><grid_spacing></b></td><td>Computational grid spacing in the X direction defined in meters.</td></tr>
 * <tr><td align="right"><b>[grid_spacing_Y]</b></td><td align="left">Computational grid spacing in the Y direction defined in meters.  This argument is optional and if it is not provided the grid spacing in the Y direction will be assumed to be the same as the spacing in the X direction.</td></tr>
 * <tr><td><b><kernel_name></b></td><td>Interpolation window (i.e. filter window) to use, the choices are 'quadloess', 'linloess', 'hanning' or 'hann', and 'boxcar'. Hanning is the standard smoothing window.</td></tr>
 * <tr><td><b><input_file_list></b></td><td>Name of a textfile containing names of individual data files containing data to be merged.  Each line of the listfile should have a single filename.</td></tr>
 * <tr><td><b><ref_lon></b></td><td>Center longitude of the given computational area.</td></tr>
 * <tr><td><b><ref_lat></b></td><td>enter latitude of the given computational area.</td></tr>
 * <tr><td><b><rotation_angle></b></td><td>An angle [in degrees] of arbitrary rotation.  This allows the computation X and Y space to undergo a rotation, such that the along-shore or cross-shore direction can be oriented along a given axis.  This allows the user to more easily select unequal arbitrary smoothing length scales.</td></tr>
 * <tr><td><b><num_MC_runs(-1 if not MC)></b></td><td>Integer number of Monte Carlo iterations to perform. Use -1 if no Monte Carlo simulations are to be run.</td></tr>
 * </table>
 * @subsection optional_args Optional Arguments
 * <table border="1">
 * <tr><td><b>[-noerr]</b></td><td>Do not output the uncertainty estimate (the mean square interpolation error estimate) in the output file.</td></tr>
 * <tr><td><b>[-nmsei]</b></td><td>Write out the normalized mean square error estimate to the output file. This appears as an additional column in the ASCII file.</td></tr>
 * <tr><td><b>[-msri]</b></td><td>Write out the mean square of the residuals to the output file. This appears as an additional column in the ASCII file. It is given after the msei-column if both are given.</td></tr>
 * <tr><td><b>[-inputInMeters]</b></td><td>Specifies that all input files are using an X/Y coordinate system in meters instead of the typical Longitude/Latitude in degrees.</td></tr>
 * <tr><td><b>[-kriging]</b></td><td>Indicates that kriging will be used on to correct oversmoothing of the interpolation surface at the input data points. This step restores finer details lost from smoothing alone, but will add to the computation time required.</td></tr>
  * <tr><td><b>[-preInterpolatedLocations]</b></td><td>Allows a file to be provided to determine the exact Longitude and Latitude locations for the interpolation.  This allows irregularly gridded output at set points instead of points determined from a user defined grid spacing.  Using this option slows down computation speed significantly and therefore it is only recommended for small data sets; however, this option does function properly on large data sets.</td></tr>
 * <tr><td align="right"><b><interpolation_location_file_name></b></td><td align="left">The file name of a file containing Longitude and Latitude columns of data points that specify output interpolation grid locations.</td></tr>
 * <tr><td><b>[-computeOffset]</b></td><td>Computes the depth offset between input data sets.</td></tr>
 * <tr><td><b>[-outputRasterFile]</b></td><td>Produce an ARC ASCII Raster output file instead of the standard output file.</td></tr>
 * <tr><td><b>[-multiThread]</b></td><td>Multi-thread the mergeBathy application.</td></tr>
 * <tr><td align="right"><b><num_threads></b></td><td align="left">The number of threads to use in the processing routine of mergeBathy.</td></tr>
 * <tr><td><b>[-msmooth]</b></td><td>Specify a smoothing scale for the interpolation in meters.</td></tr>
 * <tr><td align="right"><b><smoothing_scale_x></b></td><td align="left">Smoothing scale in meters in the X direction.</td></tr>
 * <tr><td align="right"><b><smoothing_scale_y></b></td><td align="left">Smoothing scale in meters in the Y direction.</td></tr>
 * <tr><td><b>[-llsmooth]</b></td><td>Specify a smoothing scale for the interpolation in terms of longitued and latitude.</td></tr>
 * <tr><td align="right"><b><smoothing_scale_longitude (X)></b></td><td align="left">Smoothing scale in meters in the Longitude/X direction.</td></tr>
 * <tr><td align="right"><b><smoothing_scale_latitude (Y)></b></td><td align="left">Smoothing scale in meters in the Latitude/Y direction.</td></tr>
 * <tr><td><b>[-llgrid]</b></td><td>Specify that the grid spacing is defined in Longitude/Latitude instead of meters.</td></tr>
 * <tr><td><b>[-boundingBox]</b></td><td>Specify a bounding box that is used to cut data from the input data sets.</td></tr>
 * <tr><td align="right"><b><upper_bound></b></td><td align="left">The upper bound.</td></tr>
 * <tr><td align="right"><b><lower_bound></b></td><td align="left">The lower bound.</td></tr>
 * <tr><td align="right"><b><right_bound></b></td><td align="left">The right bound.</td></tr>
 * <tr><td align="right"><b><left_bound></b></td><td align="left">The left bound.</td></tr>
 * </table>
 * @subsubsection optional_args_ZG Optional Arguments For Using MB_ZGrid External Interpolator
 * <table border="1">
 * <tr><td><b>[-ZGrid]</b></td><td>Forces the included Z_Grid software package to run.  This function was NOT created by the design team of newMergeBathy but is included under the GNU Public License Agreement.</td></tr>
 * <tr><td align="right"><b><grid_spacing_X></b></td><td align="left">omputational grid spacing in the X direction defined in meters that will be used in MB_ZGrid computation.</td></tr>
 * <tr><td align="right"><b><grid_spacing_Y></b></td><td align="left">Computational grid spacing in the Y direction defined in meters that will be used in MB_ZGrid computation.</td></tr>
 * <tr><td align="right"><b><Z_Grid_Output_File_Name></b></td><td align="left">MB_ZGrid output file.  Each line consists of three columns, corresponding to Longitude, Latitude, and Depth respectively.</td></tr>
 * <tr><td align="right"><b><Tension_Factor></b></td><td align="left">Sets the tension of the interpolation.  A value of 0.0 yields a pure Laplace (minimum curvature) solution and a value of infinity yields a pure thin plate spline solution. A value of 1e10 value has commonly been used to yield spline solutions.</td></tr>  
 * <tr><td align="right"><b><Usage></b></td><td align="left">A value of 0 will perform the MB_ZGrid interpolation and write the results to the specified output file.  A value of 1 will perform the MB_ZGrid interpolation, write the results to the specified output file, and use the computed X,Y, and Z values as input for mergeBathy.  If used as input then the data will take the place of the data read from the input files.  This allows for a pre-smoothing effect before mergeBathy is run.</td></tr>
 * </table>
 * @subsubsection optional_args_GMT Optional Arguments For Using GMT Surfacel External Interpolator
 * <table border="1">
 * <tr><td><b>[-GMTSurface]</b></td><td>Forces the included GMT Surface software package to run.  This function was NOT created by the design team of newMergeBathy.</td></tr>
 * <tr><td align="right"><b><grid_spacing_X></b></td><td align="left">Computational grid spacing in the X direction defined in meters that will be used in GMT Surface computation.</td></tr>
 * <tr><td align="right"><b><grid_spacing_Y></b></td><td align="left">Computational grid spacing in the Y direction defined in meters that will be used in GMT Surface computation.</td></tr>
 * <tr><td align="right"><b><GMT_Surface_Output_File_Name></b></td><td align="left">MB_ZGrid output file.  Two versions of this file will be created.  The first creates a file where each line consists of three columns, corresponding to Longitude, Latitude, and Depth respectively.  The second creates a file where each line consists of four columns, corresponding to Longitude, Latitude, Depth, and Error respectively.  The second file has "_includeError.txt" appended to the end of the specified file name.</td></tr>
 * <tr><td align="right"><b><Tension_Factor></b></td><td align="left">The Laplacian tension operator between 0 and 1.  Typcially 0.1</td></tr>
 * <tr><td align="right"><b><scale_factor></b></td><td align="left">The multiplier value for a Confidence Interval to be used in error calculation. A value of 1.96 is typically used for a 95% Confidence Interval.</td></tr>
 * <tr><td align="right"><b><alpha></b></td><td align="left">This is the alpha value for error computation. Typically 2.0.</td></tr>
 * <tr><td align="right"><b><Usage></b></td><td align="left">A value of 0 will perform the GMT Surface interpolation and write the results to the specified output file.  A value of 1 will perform the GMT Surface interpolation, write the results to the specified output file, and use the computed X,Y,Z, and E values as input for mergeBathy.  If used as input then the data will take the place of the data read from the input files.  This allows for a pre-smoothing effect before mergeBathy is run.</td></tr>
 * </table>
 * @n
 * @section os_support Supported Operating Systems
 * Windows x86 (32 Bit)
 * @n Windows x64 (64 Bit)
 * @n Linux x86 (32 Bit)
 * @n Linux x64 (64 Bit)
 * @n
 * @section compiling Compiling Instructions
 * The tag ${HOME} references the location of the root directory of mergeBathy.
 * @subsection compiling_windows Compiling for the Windows Environment
 * -# The file "mergeBathy.sln" in the ${HOME}/mergeBathy directory is a Microsoft Visual Studio 2010 solution file that can be used to build mergeBathy.  Open this solution file to load all of the mergeBathy project files used in compilation.  
 * -# Using the solution platform dropdown menu, select either "Win32" for a 32 bit version of mergeBathy or "x64" for a 64 bit version of mergeBathy.  (x64 is the default).
 * -# Select "Build Solution" from the build menu to compile the sources into a windows executable.
 * -# The executable will be located in "${HOME}/mergeBathy/x86/Release" for the 32 bit version and "${HOME}/mergeBathy/x64/Release" for the 64 bit version.
 * @subsection compiling_linux Compiling for the Linux Environment
 * -# The make file in the ${HOME}/mergeBathy directory will compile all of the sources necessary to create the mergeBathy executable.
 * -# From the ${HOME}/mergeBathy directory type "make mergeBathy BITFLAG=-m32" to compile a 32 bit version of mergeBathy, or "make mergeBathy BITFLAG=-m64" to compile a 64 bit version of mergeBathy.  (BITFLAG=-m64 is the default).
 * -# The Linux executable will be located in "${HOME}/mergeBathy/x86/Release" for the 32 bit version and "${HOME}/mergeBathy/x64/Release" for the 64 bit version.
 * @n
 * @section about Contact information
 * Kevin Duvieilh, kevin.duvieilh@nrlssc.navy.mil
 * @n Paul Elmore, paul.elmore@nrlssc.navy.mil
 */