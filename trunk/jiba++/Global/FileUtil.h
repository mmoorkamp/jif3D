//============================================================================
// Name        : FileUtil.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef FILEUTIL_H_
#define FILEUTIL_H_

#include "FatalException.h"
#include <string>
#include <fstream>
#include <iostream>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>
#include <iostream>

/*! \file FileUtil.h
 * Utilities associated with the general handling of files and filenames.
 */

namespace jiba
  {
    /** \addtogroup util General utility routines */
    /* @{ */

    //! Given a filename, return the file extension including the dot, i.e. for test.t it will return .t
    /*! It is a standard task to extract the filename extension, i.e. the name without the ending, for example
     * to determine the expected type of a file. If several dots are present in a file the function returns
     * only the part including and after the last dot.
     * @param filename The complete filename to process.
     * @return The ending of the file
     */
    inline std::string GetFileExtension(const std::string &filename)
      {
        std::string ending;
        size_t dotpos = filename.find_last_of('.');
        if (dotpos != std::string::npos)
          ending = filename.substr(dotpos);
        return ending;
      }

    //! Given an open input filestream search for the string token in the file and return the line containing that token, throws if fails
    /*! This function performs a case-sensitive search for the content of the string token. If a line in the file contains
     * this string the function returns this line. If the token cannot be found it throws a FatalException.
     * @param filestream An open infile stream, will be pointing to the next line after finding the token
     * @param token The token to search for
     * @return The line in the file that contains the token
     */
    inline std::string FindToken(std::ifstream &filestream,
        const std::string &token)
      {
        bool found = false;
        bool end = false;
        std::string line;
        while (!found && !end)
          {
            if (std::getline(filestream, line))
              {
                found = boost::algorithm::contains(line, token);
              }
            else
              {
                end = true;
              }
          }
        if (found)
          {
            return line;
          }
        else
          {
            throw FatalException("Token " + token + " not found !");
          }
      }

    inline std::string AskFilename(const std::string &prompt)
      {
        std::cout << prompt;
        std::string filename;
        std::cin >> filename;
        if (!boost::filesystem::exists(filename))
          {
            throw jiba::FatalException("File " + filename + " does not exist");
          }
        return filename;
      }
  /* @} */
  }
#endif /* FILEUTIL_H_ */
