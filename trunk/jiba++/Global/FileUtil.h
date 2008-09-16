//============================================================================
// Name        : FileUtil.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef FILEUTIL_H_
#define FILEUTIL_H_

#include <string>
namespace jiba
  {
    //! Given a filename, return the file extension including the dot, i.e. for test.t it will return .t
    inline std::string GetFileExtension(const std::string &filename)
      {
        std::string ending;
        size_t dotpos = filename.find_last_of('.');
        if (dotpos != std::string::npos)
          ending = filename.substr(dotpos);
        return ending;
      }
  }
#endif /* FILEUTIL_H_ */
