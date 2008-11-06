//============================================================================
// Name        : NumUtil.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef NUMUTIL_H_
#define NUMUTIL_H_

/*! \file NumUtil.h
 * This file contains a collection of various simple numerical routines.
 */
namespace jiba
  {
    /** \addtogroup util General utility routines */
    /* @{ */

    //! Determine the sign of a numeric type by comparing it with zero
    template<class NumericType>
    int sign(NumericType n)
      {
        if (n >= 0)
          return 1;
        return -1;
      }

    //! Generate a sequence of integers, for use with generate and generate_n algorithms, taken from Josuttis, p. 296
    /*! The constructor takes the starting value and each call to operator() returns the current value and increases
     * it by 1.
     */
    class IntSequence
      {
    private:
      int value;
    public:
      IntSequence(int start) :
        value(start)
        {
        }
      int operator()()
        {
          return value++;
        }
      };
  /* @} */
  }
#endif /*NUMUTIL_H_*/
