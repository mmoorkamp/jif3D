//============================================================================
// Name        : parrun.cpp
// Author      : Sep 21, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <boost/filesystem/operations.hpp>
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <typeinfo>
//#include <cstdlib>

namespace fs = boost::filesystem;

const std::string runext = "_run";
const std::string dirext = "_dir";

template<typename T> inline std::string stringify(const T& x)
  {
    std::ostringstream o;
    if (!(o << x))
      throw std::runtime_error(std::string("stringify(") + typeid(x).name()
          + ")");
    return o.str();
  }

void CleanFiles(const std::string &NameRoot)
  {
    fs::remove_all(NameRoot + dirext);
    fs::remove_all(NameRoot + runext);
  }

void MakeRunFile(const std::string &NameRoot)
  {
    std::string DirName = NameRoot + dirext;
    std::string RunFileName = NameRoot + runext;
    fs::create_directory(DirName);
    std::ofstream runfile;
    runfile.open(RunFileName.c_str());
    runfile << "#!/bin/bash\n";
    runfile << "cd " << DirName << "\n";
    runfile << "x3d\n";
    runfile << "cd .." << std::endl;
    runfile.close();
    int result = std::system((std::string("chmod u+x ./") + RunFileName).c_str());
    if (result)
      throw std::runtime_error("Cannot make script executable !");
  }

void Run(const std::string &NameRoot)
  {
    const std::string runname = "./" + NameRoot + runext;
    int result = std::system(runname.c_str());
    if (result)
      throw std::runtime_error("Cannot execute run script: " + runname
          + " Error Code: " + stringify(result));
  }

int main()
  {

    std::cout << "Starting " << std::endl;
    const int nruns = 50;
#pragma omp parallel for
    for (int i = 0; i < nruns; ++i)
      {
        std::string name = "test" + stringify(i);
#pragma omp critical
          {
            MakeRunFile(name);
          }
        Run(name);
        CleanFiles(name);

      }

  }

