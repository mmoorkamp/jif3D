//============================================================================
// Name        : ThreeDMTModel.cpp
// Author      : Apr 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <fstream>
#include "ThreeDMTModel.h"

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

    void ThreeDMTModel::ReadModEM(const std::string filename)
      {
        std::ifstream infile(filename.c_str());
        //swallow the first description line
        char dummy[1024];
        infile.getline(dummy,1024);
        int nx, ny, nz;
        infile >> nx >> ny >> nz;
        infile.getline(dummy,1024);
        this->SetXCellSizes().resize(nx);
        this->SetYCellSizes().resize(ny);
        this->SetZCellSizes().resize(nz);
        this->SetData().resize(boost::extents[nx][ny][nz]);
        for (int i = 0; i < nx; ++i)
          infile >> this->SetXCellSizes()[i];
        for (int i = 0; i < ny; ++i)
          infile >> this->SetYCellSizes()[i];
        for (int i = 0; i < nz; ++i)
          infile >> this->SetZCellSizes()[i];
        for (int i = 0; i < nx; ++i)
          for (int j = 0; j < ny; ++j)
            for (int k = 0; k < nz; ++k)
              infile >> this->SetData()[i][j][k];
      }
  }
