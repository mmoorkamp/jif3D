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

    template<class N1, class N2> struct absLess: public std::binary_function<
        N1, N2, bool>
      {
      bool operator()(N1 number1, N2 number2) const
        {
          return abs(number1) < abs(number2);
        }
      };

    template<> struct absLess<double, double> : public std::binary_function<
        double, double, bool>
      {
      bool operator()(double number1, double number2) const
        {
          return fabs(number1) < fabs(number2);
        }
      };

    template<> struct absLess<float, float> : public std::binary_function<
        float, float, bool>
      {
      bool operator()(float number1, float number2) const
        {
          return fabs(number1) < fabs(number2);
        }
      };

    template<> struct absLess<long double, long double> : public std::binary_function<
        long double, long double, bool>
      {
      bool operator()(long double number1, long double number2) const
        {
          return fabs(number1) < fabs(number2);
        }
      };
  /* @} */
  }
#endif /*NUMUTIL_H_*/
