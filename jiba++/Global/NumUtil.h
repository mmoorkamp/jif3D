//============================================================================
// Name        : NumUtil.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef NUMUTIL_H_
#define NUMUTIL_H_
#include <cmath>

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
     * it by 1. This class is useful to generate values in a vector with std::generate.
     */
    class IntSequence
      {
    private:
      int value;
    public:
      //! The constructor takes the first value of the sequence
      IntSequence(int start) :
        value(start)
        {
        }
      int operator()()
        {
          return value++;
        }
      };

    //! Compare number1 and number2 and return true if the absolute value of number1 is smaller than the absolute value of number2, otherwise return false.
    /*! This template is useful in conjunction with the standard library. For example we can write
     * std::sort(Vector.begin(),Vector.end(),absLess<double>()) to sort a vector of numbers by their absolute values.
     *
     */
    template<class N1, class N2> struct absLess: public std::binary_function<
        N1, N2, bool>
      {
      //! Functors for the standard library have their functionality in operator()
      bool operator()(N1 number1, N2 number2) const
        {
          return std::abs(number1) < std::abs(number2);
        }
      };

  /* @} */
  }
#endif /*NUMUTIL_H_*/
