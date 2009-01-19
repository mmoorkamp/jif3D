//============================================================================
// Name        : FileUtil.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef FILEUTIL_H_
#define FILEUTIL_H_

#include <string>
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
  /* @} */
  }
#endif /* FILEUTIL_H_ */
