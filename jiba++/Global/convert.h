#ifndef CONVERT_H_
#define CONVERT_H_

// File: convert.h
#include <iostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <stdexcept>
#include <locale>

class BadConversion : public std::runtime_error
  {
public:
    BadConversion(const std::string& s) :
      std::runtime_error(s)
      {
      }
  };
//! Convert a streamable type into a string
template<typename T> inline std::string stringify(const T& x)
  {
    std::ostringstream o;
    if (!(o << x))
      throw BadConversion(std::string("stringify(")
          + typeid(x).name() + ")");
    return o.str();
  }

//! Convert from string to a streamable type
template<typename T> inline void convert(const std::string& s, T& x,
    bool failIfLeftoverChars = true)
  {
    std::istringstream i(s);
    char c;
    if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
      throw BadConversion(s);
  }
//! Convert a character to upper case
class mytoupper
  {
private:
    const std::locale& loc;
public:
    mytoupper(const std::locale& l) :
      loc(l)
      {
      }
    char operator()(const char lowc)
      {
        return std::toupper(lowc, loc);
      }
  };
//! Convert a character to lower case
class mytolower
  {
private:
    const std::locale& loc;
public:
    mytolower(const std::locale& l) :
      loc(l)
      {
      }
    char operator()(const char lowc)
      {
        return std::tolower(lowc, loc);
      }
  };
#endif /*CONVERT_H_*/
