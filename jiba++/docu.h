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
 */
#ifndef DOCU_H_
#define DOCU_H_

#endif /* DOCU_H_ */
