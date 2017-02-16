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
#include "../ModelBase/VTKTools.h"

int main()
  {
    std::string ncfilename = jif3D::AskFilename("Name of netcdf file: ");


    jif3D::rvec Data, Error;
    jif3D::ThreeDModelBase::tMeasPosVec MeasX, MeasY, MeasZ;
    jif3D::ReadScalarGravityMeasurements(ncfilename,Data,MeasX,MeasY,MeasZ,Error);

    std::cout << "Rotation angle [degree]: ";
    double dangle = 0.0;
    std::cin >> dangle;
    double rangle = dangle / 180.0 * boost::math::constants::pi<double>();

    for (size_t i = 0; i < MeasX.size(); ++i)
      {
        double newx = MeasX.at(i)  * cos(rangle) + MeasY.at(i)  * sin(rangle);
        double newy = MeasY.at(i)  * cos(rangle) - MeasX.at(i)  * sin(rangle);
        MeasX.at(i) = newx;
        MeasY.at(i) = newy;

      }
    jif3D::SaveScalarGravityMeasurements(ncfilename+".rot.nc",Data,MeasX,MeasY,MeasZ,Error);


    jif3D::Write3DDataToVTK(ncfilename + "_rot.statpos.vtk", "Station",
        jif3D::rvec(MeasX.size(), 1.0), MeasX, MeasY, MeasZ);

  }
