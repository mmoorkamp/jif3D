//============================================================================
// Name        : ThreeDMTModel.cpp
// Author      : Apr 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "ThreeDMTModel.h"
#include "../Global/convert.h"
#include <fstream>
namespace jiba
  {

    ThreeDMTModel::ThreeDMTModel()
      {


      }

    ThreeDMTModel::~ThreeDMTModel()
      {

      }

    void ThreeDMTModel::WriteProjectFile()
      {
        const size_t nfreq = Frequencies.size();
        std::ofstream outfile("a.project");
        outfile << "  Version_of_X3D code (yyyy-mm-dd)\n";
        outfile << "  2006-06-06\n\n";
        outfile << "  Type_of_problem (0 - MT, 1 - CSMT, 2 - EDIP, 3 - MDIP)\n";
        outfile << "  1\n";
        outfile << "  Frequency (Hz)    File_with_results     File with 3D formation      File with 1st source      File with 2nd source\n";
        for (size_t i = 0; i < nfreq; ++i)
          {
            outfile << "  " << Frequencies.at(i) << "   work/csmt01" << jiba::stringify(i)
            <<"                csmt01.model                csmt01_1.source                csmt01_2.source\n";
          }
        outfile << "$\n";
        outfile << "  Threshold_of relative residual norm  (default value = 0.003)\n";
        outfile << "  0.003\n\n";
        outfile << "  Maximum_number (Nmax) of iterates (default value = 500)\n";
        outfile << "  1000\n";
        outfile << "  Name_for_QAA (default value = value of Generic_name)\n";
        outfile << "  csmtQaa\n\n";
        outfile << "  Far_Zone (Y or N; default value = N; if Y it defines slow,but accurate calculus of EM fields in far zone from a dipole source)\n";
        outfile << "  Y\n\n";
        outfile << "  Format_of_Output (0 - mfo, 1 - ASCII, 2 - mfo+ASCII, 3 - mfo+ASCII+mfa; default value = 0)\n";
        outfile << "  3\n";
      }
  }
