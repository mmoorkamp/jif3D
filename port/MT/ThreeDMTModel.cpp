//============================================================================
// Name        : ThreeDMTModel.cpp
// Author      : Apr 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "ThreeDMTModel.h"
#include <fstream>


namespace jif3D
  {

    ThreeDMTModel::ThreeDMTModel() :
        Frequencies()
      {

      }

    ThreeDMTModel::~ThreeDMTModel()
      {

      }

    ThreeDMTModel::ThreeDMTModel(const ThreeDMTModel &source) :
        ThreeDModelBase(source), Frequencies(source.Frequencies)
      {

      }

    ThreeDMTModel& ThreeDMTModel::operator=(const ThreeDMTModel& source)
      {
        if (&source != this)
          {
            //apart from copying the contents of the base class
            //we have to copy the vector of frequencies which
            //is the only additional data in this class
            ThreeDModelBase::operator=(source);
            DistortionParameters = source.DistortionParameters;
            Frequencies = source.Frequencies;
          }
        return *this;
      }

    ThreeDMTModel& ThreeDMTModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
      }

    void ThreeDMTModel::WriteModEM(const std::string filename)
      {
        std::ofstream outfile(filename.c_str());
        outfile << "# 3D MT model written by jif3D in WS format\n";
        const size_t nx = GetXCellSizes().size();
        const size_t ny = GetYCellSizes().size();
        const size_t nz = GetZCellSizes().size();
        outfile << nx << " " << ny << " " << nz << " 0 LOGE\n";
        for (double cellx : GetXCellSizes())
          {
            outfile << cellx << " ";
          }
        outfile << "\n";
        for (double celly : GetYCellSizes())
          {
            outfile << celly << " ";
          }
        outfile << "\n";
        for (double cellz : GetZCellSizes())
          {
            outfile << cellz << " ";
          }
        outfile << "\n";
        for (size_t i = 0; i < nz; ++i)
          {
            for (size_t j = 0; j < ny; ++j)
              {
                for (int k = nx - 1; k >= 0; --k)
                  {

                    outfile << std::log(1.0 / GetData()[k][j][i]) << " ";
                  }
                outfile << "\n";
              }
            outfile << "\n\n";
          }
        outfile << "\n  0.0 0.0 0.0\n 0.0 \n";
      }

    void ThreeDMTModel::ReadModEM(const std::string filename)
      {
        std::ifstream infile(filename.c_str());
        //swallow the first description line
        char dummy[1024];
        infile.getline(dummy, 1024);
        int nx, ny, nz;
        infile >> nx >> ny >> nz;
        infile.getline(dummy, 1024);
        SetXCellSizes().resize(nx);
        SetYCellSizes().resize(ny);
        SetZCellSizes().resize(nz);
        SetData().resize(boost::extents[nx][ny][nz]);
        for (int i = 0; i < nx; ++i)
          infile >> SetXCellSizes()[i];
        for (int i = 0; i < ny; ++i)
          infile >> SetYCellSizes()[i];
        for (int i = 0; i < nz; ++i)
          infile >> SetZCellSizes()[i];
        double value;
        for (int i = 0; i < nz; ++i)
          for (int j = 0; j < ny; ++j)
            for (int k = nx - 1; k >= 0; --k)
              {
                infile >> value;
                SetData()[k][j][i] = std::exp(-value);
              }
      }
  }
