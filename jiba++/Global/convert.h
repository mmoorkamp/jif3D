//============================================================================
// Name        : convert.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef CONVERT_H_
#define CONVERT_H_

// File: convert.h
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <stdexcept>
#include <locale>

/*! \file convert.h
 * This file contains routines to convert streamable types to a string and versions of toupper and tolower that
 * are easier to handle in generic algorithms like transform etc.
 */

namespace jiba
  {
    /*! \addtogroup util General utility routines
     * These routines are of very general use and do not belong to a certain geophysical method. This includes
     * general mathematical methods such as wavelet transform or interpolation, but also small helper functions
     * and classes for conversion and comparison.
     */
    /* @{ */
    //! This exception signals that something went wrong during string conversion.
    class BadConversion: public std::runtime_error
      {
    public:
      //! The constructor takes the error description and simply passes it on to std::runtime_error.
      BadConversion(const std::string& s) :
        std::runtime_error(s)
        {
        }
      };
    //! Convert a streamable type into a string
    /*! This function can be used to make a string from any type
     * that has an overload for the streaming operator <<. This is
     * very useful for the conversion of numbers to strings to construct
     * filenames, for example.
     * @param x The input value that should be converted to a string
     * @return The string containing "x"
     */
    template<typename T> inline std::string stringify(const T& x)
      {
        std::ostringstream o;
        if (!(o << x))
          throw BadConversion(std::string("stringify(") + typeid(x).name() + ")");
        return o.str();
      }

    //! Convert from string to a streamable type
    /*! Convert a string to a type T that has the streaming operator >> overloaded.
     * @param s The string that should be converted
     * @param x The variable that holds the result of the conversion
     * @param failIfLeftoverChars Is it an error if there are more characters in the string than can be converted ? Default is true.
     */
    template<typename T> inline void convert(const std::string& s, T& x,
        bool failIfLeftoverChars = true)
      {
        std::istringstream i(s);
        char c;
        if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
          throw BadConversion(s);
      }
    //! Convert a character to upper case
    /*! The standard toupper function needs a locale
     * specified in the call that makes it awkward to
     * use with algorithms like transform. This class
     * simply takes the locale as an arguments to the
     * constructor and provides a unary overload of operator()
     * that is ready for use with standard algorithms.
     */
    class mytoupper
      {
    private:
      //! Storage for the locale, gets assigned in constructor to make use with algorithms easier
      const std::locale& loc;
    public:
      //! The constructor needs the locale to convert special characters properly
      mytoupper(const std::locale& l) :
        loc(l)
        {
        }
      //! Convert the character lowc to upper case given the locale of the constructor
      char operator()(const char lowc)
        {
          return std::toupper(lowc, loc);
        }
      };
    //! Convert a character to lower case
    /*! The standard tolower function needs a locale
     * specified in the call that makes it awkward to
     * use with algorithms like transform. This class
     * simply takes the locale as an arguments to the
     * constructor and provides a unary overload of operator()
     * that is ready for use with standard algorithms.
     */
    class mytolower
      {
    private:
      //! Storage for the locale, gets assigned in constructor to make use with algorithms easier
      const std::locale& loc;
    public:
      //! The constructor needs the locale to convert special characters properly
      mytolower(const std::locale& l) :
        loc(l)
        {
        }
      //! Convert the character lowc to lower case given the locale of the constructor
      char operator()(const char lowc)
        {
          return std::tolower(lowc, loc);
        }
      };
  /* @} */
  }
#endif /*CONVERT_H_*/
