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
#include <boost/filesystem/convenience.hpp>
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

    //! Instead of always writing the same things in the main program, we can use this function to ask the use for a filename
    /*! This function provides a very simple way to get the name of a file from a use, by default it also checks
     * whether this file exists.
     * @param prompt The message to show to the user.
     * @param checkexists Shall the function check whether the file exists ? Default true
     * @return The name of the file as typed in by the user.
     */
    inline std::string AskFilename(const std::string &prompt, const bool checkexists = true)
      {
        std::cout << prompt;
        std::string filename;
        std::cin >> filename;
        if (checkexists && !boost::filesystem::exists(filename))
          {
            throw jiba::FatalException("File " + filename + " does not exist");
          }
        return filename;
      }

    //! Checks whether the file called abort exists in the current directory to signal the program that we want to stop
    inline bool WantAbort()
      {
        return boost::filesystem::exists("abort");
      }

    //! Remove the abort file in the current directory
    inline void RemoveAbort()
      {
        boost::filesystem::remove("abort");
      }
  /* @} */
  }
#endif /* FILEUTIL_H_ */
