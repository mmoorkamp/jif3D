//============================================================================
// Name        : NumUtil.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef NUMUTIL_H_
#define NUMUTIL_H_

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>

#include "Jif3DGlobal.h"

/*! \file NumUtil.h
 * This file contains a collection of various simple numerical routines.
 */
namespace jif3D
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

    //! Convenience function for squaring numbers, more efficient than std::pow(x,2)
    template<class N1>
    N1 pow2(N1 x)
      {
        return x * x;
      }
    //! Convenience function for squaring numbers, more efficient than std::pow(x,2)
    template<class N1>
    N1 pow3(N1 x)
      {
        return x * x * x;
      }
    //! Compare number1 and number2 and return true if the absolute value of number1 is smaller than the absolute value of number2, otherwise return false.
    /*! This template is useful in conjunction with the standard library. For example we can write
     * std::sort(Vector.begin(),Vector.end(),absLess<double,double>()) to sort a vector of numbers by their absolute values.
     *
     */
    template<class N1, class N2> struct absLess: public std::binary_function<N1, N2, bool>
      {
      //! Functors for the standard library have their functionality in operator()
      bool operator()(N1 number1, N2 number2) const
        {
          return std::abs(number1) < std::abs(number2);
        }
      };
    //! Compare two numbers within a given tolerance
    /*! When comparing two floating point numbers, we typically
     * cannot simply compare them, due to rounding issues, so we
     * compare within the  tolerance specified by the numeric_limits::epsilon
     * specified for the type.
     */
    template<class N1, class N2> struct roughlyEqual: public std::binary_function<N1, N2,
        bool>
      {
      N1 epsilon;
      //! Functors for the standard library have their functionality in operator()
      bool operator()(N1 number1, N2 number2) const
        {
          return std::abs(number1 - number2) < epsilon;
        }
      roughlyEqual(N1 tol = std::numeric_limits < N1 > ::epsilon()) :
          epsilon(tol)
        {

        }
      };

    //! Checks whether the input is a power of two
    /*! Is there a positive m with 2^m=n~?
     * @param n A positive integer n
     * @return Is n a power of 2
     */
    inline bool IsPowerOfTwo(const size_t n)
      {
        size_t next = 1;
        while (next < n)
          next *= 2;
        return next == n;
      }
  /* @} */
  }
#endif /*NUMUTIL_H_*/
