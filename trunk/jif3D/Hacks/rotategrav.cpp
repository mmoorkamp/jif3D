//============================================================================
// Name        : rotatetomo.cpp
// Author      : Feb 21, 2014
// Version     :
// Copyright   : 2014, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../Global/VecMat.h"
#include "../Global/ReadWriteSparseMatrix.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ScalarGravityData.h"
#include "../ModelBase/VTKTools.h"

int main()
  {
    std::string ncfilename = jif3D::AskFilename("Name of netcdf file: ");

    jif3D::ScalarGravityData Data;
    Data.ReadNetCDF(ncfilename);

    std::cout << "Rotation angle [degree]: ";
    double dangle = 0.0;
    std::cin >> dangle;
    double rangle = dangle / 180.0 * boost::math::constants::pi<double>();
    const size_t ndata =Data.GetMeasPosX().size();
    std::vector<double> newx(ndata), newy(ndata);

    for (size_t i = 0; i < ndata; ++i)
      {
        newx.at(i) = Data.GetMeasPosX().at(i) * cos(rangle)
            + Data.GetMeasPosY().at(i) * sin(rangle);
        newy.at(i) = Data.GetMeasPosY().at(i) * cos(rangle)
            - Data.GetMeasPosX().at(i) * sin(rangle);
      }
    Data.SetMeasurementPoints(newx, newy, Data.GetMeasPosZ());
    Data.WriteNetCDF(ncfilename + ".rot.nc");
    Data.WriteVTK(ncfilename + "_rot.statpos.vtk");

  }
