//============================================================================
// Name        : docu.h
// Author      : Jan 19, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

//This file exists purely to hold the text for the doxygen main page

/*! \mainpage jif3D
 *
 * \section Introduction
 *
 * jif3D is a C++ framework for joint inversion of different types of geophysical data for sub-basalt, sub-salt and other
 * challenging imaging problems. The general idea of the joint inversion is similar to the jiba program by Heincke et al.,
 * however it is developed from scratch to provide maximal modularity and avoid lincence clashes between different parts of
 * the program.
 *
 * Note that this is research code with a focus on numerical correctness and maintainability. This has two implications:
 *   - User input is rarely checked for correctness. The code performs some checks for consistency, but nonsensical values
 *     and incoherent model descriptions in the input files can cause crashes or random results.
 *   - While I have followed basic performance guidelines and used state of the art methods wherever possible, I have
 *     not performed any detailed speed optimizations of the code. Also, if there is a choice between clarity and speed
 *     I usually opt for clarity. However, some basic profiling indicates that there are no major bottlenecks caused by this.
 *   - Documentation is incomplete at the moment and a user manual does not exist, yet. We are working on
 *     this, but you currently have to look at the sources and this documentation to figure out what
 *     the code does and how it can be used.
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

 \section Getting the code

 You can check out the latest development version through svn via svn://svn.code.sf.net/p/jif3d/jif3dsvn/trunk/jif3d
 Alternatively you can download a release from http://sourceforge.net/projects/jif3d/

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
 gcc 4.6 under ubuntu 12/04
 intel 11.0 under ubuntu 8/10
 intel 12.0.2 under ubuntu 10/04 and ubuntu 12/04

 Sun 5.10 is known to fail, this is an issue with partial template specialization in
 the boost libraries.

 Sun 5.11 compiles most of the code, CUDA does not support the sun compilers though. However this
 can be quite easily fixed by removing the Gravity classes that use CUDA.

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
 * \version Mar 2011
 * This is the final version that is part of the first phase of the JIBA project. The joint inversion code should be relatively
 * stable and well documented. Also, we fixed a number of bugs and added new features for the last year.
 *
 * \version May 2013
 * The initial open-source beta release. This version contains many new and experimental features. The documentation
 * has been improved, but is still incomplete.
 *
 * \section License
 * All code is distributed under the GNU General Public License (GPL) v3. The full text of the license
 * can be found at http://opensource.org/licenses/GPL-3.0
 *
 * The seismic modeling code written by Podvin and Lecomte is being distributed as part of the jif3d code and is also licensed under the GPL.
 *
 * \section Papers
 * Two papers describe different aspects of jif3d
 *   - Moorkamp, M., Heincke, B., Jegen, M., Roberts, A.W., Hobbs, R.W. A framework for 3-D joint inversion of MT, gravity and seismic refraction data
        (2011) Geophysical Journal International, 184 (1), pp. 477-493.
      explains the overall concept of out framework and compares different coupling mechanisms.

    - Moorkamp, M., Jegen, M., Roberts, A., Hobbs, R. Massively parallel forward modeling of scalar
      and tensor gravimetry data (2010) Computers and Geosciences, 36 (5), pp. 680-686
      evaluates the performance of our GPGPU gravity forward modelling code.
   If you use this code in an academic context, we kindly ask you to cite the relevant publications.

   \section Acknowledgements

   The initial development of the framework was sponsored by Chevron, ExxonMobil, Nexen, RWE Dea,
   Shell, Statoil and Wintershall within the JIBA consortium. P. Podvin kindly made his eikonal
   solver publicly available. A. Avdeeva supported us with implementing the gradient calculation
   for the MT part of the joint inversion.
 */
#ifndef DOCU_H_
#define DOCU_H_

#endif /* DOCU_H_ */
