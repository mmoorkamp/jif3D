//============================================================================
// Name        : NumUtil.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef NUMUTIL_H_
#define NUMUTIL_H_

/*! \file This file contains a collection of various simple numerical routines. 
 */

//! Determine the sign of a numeric type by comparing it with zero
template <class NumericType>
int sign (NumericType n)
{
   if (n >= 0) return 1;
   return -1;
}

#endif /*NUMUTIL_H_*/
