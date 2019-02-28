//============================================================================
// Name        : vtk2netcdf.cpp
// Author      : 15 Jun 2013
// Version     :
// Copyright   : 2013, mm489
//============================================================================

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "../MT/X3DModel.h"
#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"

int main()
  {
    jif3D::X3DModel Model;
    std::string Infilename = jif3D::AskFilename("Filename: ");
    jif3D::ThreeDModelBase::t3DModelDim XCellSizes, YCellSizes, ZCellSizes;
    jif3D::Read3DModelFromVTK(Infilename, XCellSizes, YCellSizes, ZCellSizes,
        Model.SetConductivities());
    Model.SetHorizontalCellSize(XCellSizes[1], YCellSizes[1], XCellSizes.size(),
        YCellSizes.size());
    Model.SetZCellSizes(ZCellSizes);
    const size_t nz = Model.GetZCellSizes().size();
    std::vector<double> bg_cond(nz), bg_thick(nz);
    //std::copy(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
    //		bg_thick.end());

    for (size_t i = 0; i < nz; ++i)
      {
        bg_thick.at(i) = Model.GetZCellSizes()[i];
        bg_cond.at(i) = Model.GetConductivities()[0][0][i];
      }

    Model.SetBackgroundConductivities(bg_cond);
    Model.SetBackgroundThicknesses(bg_thick);
    Model.WriteNetCDF(Infilename + ".nc");

    std::ofstream ascfile((Infilename + ".asc").c_str());
    std::vector<double> XPos, YPos, ZPos;
    std::partial_sum(XCellSizes.begin(), XCellSizes.end(), std::back_inserter(XPos));
    std::partial_sum(YCellSizes.begin(), YCellSizes.end(), std::back_inserter(YPos));
    std::partial_sum(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
        std::back_inserter(ZPos));
    const size_t nx = XPos.size();
    const size_t ny = YPos.size();

    for (size_t i = 0; i < nx; ++i)
      {
        for (size_t j = 0; j < ny; ++j)
          {
            for (size_t k = 0; k < nz; ++k)
              {
                ascfile << std::fixed;
                ascfile << std::setw(15) << XPos[i] - XCellSizes[i] / 2.0;
                ascfile << std::setw(15) << YPos[j] - YCellSizes[j] / 2.0;
                ascfile << std::setw(15) << ZPos[k] - Model.GetZCellSizes()[k] / 2.0;
                ascfile << std::setw(15) << 1.0 / Model.GetConductivities()[i][j][k]
                    << "\n";
              }

          }
      }

  }

