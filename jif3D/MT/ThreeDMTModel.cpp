//============================================================================
// Name        : ThreeDMTModel.cpp
// Author      : Apr 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "ThreeDMTModel.h"
#include <fstream>
#include "MTData.h"
#include "TipperData.h"

namespace jif3D
  {

    void MTDistortionSetter::operator()(ThreeDMTModel &Model, MTData &Data,
        const std::vector<double> &Dist)
      {
        Data.SetDistortion(Dist);
        // Model.SetDistortionParameters(Dist);
      }

    void MTDistortionSetter::operator()(ThreeDMTModel &Model, TipperData &Data,
        const std::vector<double> &Dist)
      {
        throw jif3D::FatalException("Cannot have extra inversion parameters with Tipper",
            __FILE__, __LINE__);
      }

    ThreeDMTModel::ThreeDMTModel()
      {

      }

    ThreeDMTModel::~ThreeDMTModel()
      {

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
        outfile << "\n  " << GetXCoordinates()[0] << " " << GetYCoordinates()[0] << " "
            << GetZCoordinates()[0] << "\n 0.0 \n";
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
        ThreeDModelBase::t3DModelDim XCS, YCS, ZCS;
        XCS.resize(nx);
        YCS.resize(ny);
        ZCS.resize(nz);
        SetData().resize(boost::extents[nx][ny][nz]);
        for (int i = 0; i < nx; ++i)
          infile >> XCS[i];
        for (int i = 0; i < ny; ++i)
          infile >> YCS[i];
        for (int i = 0; i < nz; ++i)
          infile >> ZCS[i];
        double value;
        for (int i = 0; i < nz; ++i)
          for (int j = 0; j < ny; ++j)
            for (int k = nx - 1; k >= 0; --k)
              {
                infile >> value;
                SetData()[k][j][i] = std::exp(-value);
              }
        SetXCellSizes(XCS);
        SetYCellSizes(YCS);
        SetZCellSizes(ZCS);
        double XOrigin, YOrigin, ZOrigin;
        infile >> XOrigin >> YOrigin >> ZOrigin;
        SetOrigin(XOrigin, YOrigin, ZOrigin);
      }
  }
