//============================================================================
// Name        : docu.h
// Author      : Jan 19, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

//This file exists purely to hold the text for the doxygen main page

/*! \mainpage jiba++
 *
 * \section Introduction
 *
 * jiba++ is a C++ framework for joint inversion of different types of geophysical data for sub-basalt, sub-salt and other
 * challenging imaging problems. The general idea of the joint inversion is similar to the jiba program by Heincke et al.,
 * however it is developed from scratch to provide maximal modularity and avoid lincence clashes between different parts of
 * the program.
 *
 * Note that this is research code with a focus on numerical correctness and maintainability. This has two implications:
 *   - User input is rarely checked for correctness. The code performs some checks for consistency, but non-sensical values
 *     and incoherent model descriptions in the input files can cause crashes or random results.
 *   - While I have followed basic performance guidelines and used state of the art methods wherever possible, I have
 *     not performed any detailed speed optimizations of the code. Also, if there is a choice between clarity and speed
 *     I usually opt for clarity. However, some basic profiling indicates that there are no major bottlenecks caused by this.
 *
   \section Requirements

 The code uses a number of libraries and tools all of which are freely available,
 but have to be installed before compilation. Please refer to the respective
 documentation on how to install these components.

   - Boost Jam
   - Boost
 These two are available from http://www.boost.org/users/download/

   - Boost numeric bindings from http://mathema.tician.de/dl/software/boost-numeric-bindings
   - netcdf from http://www.unidata.ucar.edu/software/netcdf/
   - cuda from http://www.nvidia.com/object/cuda_home.html#

These two libraries are optional and have to be activated in the code through appropriate pre-processor flags (see below).
   - atlas from http://math-atlas.sourceforge.net/
   - lapack from http://www.netlib.org/lapack/

 Alternatively you can install all but cuda through the package manager on ubuntu 08/10.
 Depending on the path you chose to install the libraries to, you might
 have to modify the Jamroot file. Refer to http://www.boost.org/doc/tools/build/index.html
 for a description of the syntax.

 \section Compilation

 Once you have installed everything and set the correct path you should only
 have to run "bjam" from the jiba++ directory. Depending on your compiler
 and whether you are doing a debug or release compilation, the programs and
 libraries will be put in different directories. These will be output during
 compilation.

 This code has been compiled with

 gcc 4.2 under SuSE Linux 10.1 and 11.0,
 gcc 4.3 under ubuntu 8/10,
 gcc 4.4 under ubuntu 10/04
 intel 11.0 under ubuntu 8/10

 Sun 5.10 is known to fail, this is an issue with partial template specialization in
 the boost libraries.

 Also, Visual Studio currently does not compile the code as we are using the c++ interface
 for the netcdf library that has not been ported to Windows, yet. An possible alternative
 is cygwin and the code should compile there. However, this has not been tested.

 *
 * \section pflags Preprocessor flags
 *
 * There are a number of flags for the preprocessor that determine with which features the code will be compiled. The following table summarizes
 * the flags and what happens if they are defined.
 *
 * <table>
 * <tr> <td>TESTR </td> <td>Run the unit tests for the R interface. Some versions of the boost unit-test system have  problems with the system call necessary for this test, so it is optional. </td></tr>
 * <tr> <td>HAVEUBCCODE </td> <td>Is the UBC gravity code present on the system. If defined we run a unit test that compares our results with that code.  </td></tr>
 * <tr> <td>HAVEATLAS </td> <td>This flag only affects FullSensitivityGravityCalculator.h This is to avoid dependency on blas etc. for the gravity code and makes using it on other machines easier. </td></tr>
 * </table>
 *
 *
 * \section History
 *
 * \version Jan 2008
 * This version is the first milestone and contains the deliverables for scalar gravity and FTG forward modeling and associated tools.
 *
 * \version Jan 2010
 * With the technical part of the extension to 3D complete, the forward codes for gravity and seismics and the joint inversion methods
 * have been delivered to the project sponsors.
 *
 *
 * \section License
 * The joint inversion and gravity forward modeling codes, as well as the associated tools are distributed to the sponsors of the JIBA project under the terms
 * and conditions defined in the project contract.
 *
 * The seismic modeling code written by Podvin and Lecomte is being distributed as part of the jiba++ code, but is licensed under the LGPL (see http://www.gnu.org/licenses/lgpl-2.1.html
 * for details). Note that this license is different from the public version of time3d and the result of a special agreement with the authors.
 */
#ifndef DOCU_H_
#define DOCU_H_

#endif /* DOCU_H_ */
