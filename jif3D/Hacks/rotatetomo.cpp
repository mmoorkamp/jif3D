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
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyData.h"
#include "../ModelBase/VTKTools.h"

int main()
  {
    std::string ncfilename = jif3D::AskFilename("Name of netcdf file: ");

    jif3D::TomographyData Data;
    Data.ReadNetCDF(ncfilename);


    std::cout << "Rotation angle [degree]: ";
    double dangle = 0.0;
    std::cin >> dangle;
    double rangle = dangle / 180.0 * boost::math::constants::pi<double>();

    std::vector<double> SourceX(Data.GetSourcePosX());
    std::vector<double> SourceY(Data.GetSourcePosY());
    std::vector<double> SourceZ(Data.GetSourcePosZ());

    std::vector<double> RecX(Data.GetMeasPosX());
    std::vector<double> RecY(Data.GetMeasPosY());
    std::vector<double> RecZ(Data.GetMeasPosZ());


    for (size_t i = 0; i < SourceX.size(); ++i)
      {
        double newx = SourceX.at(i) * cos(rangle) + SourceY.at(i) * sin(rangle);
        double newy = SourceY.at(i) * cos(rangle) - SourceX.at(i) * sin(rangle);
        SourceX.at(i) = newx;
        SourceY.at(i) = newy;
      }

    for (size_t i = 0; i < RecX.size(); ++i)
      {
        double newx = RecX.at(i) * cos(rangle) + RecY.at(i) * sin(rangle);
        double newy = RecY.at(i) * cos(rangle) - RecX.at(i) * sin(rangle) ;
        RecX.at(i) = newx;
        RecY.at(i) = newy;
      }
    Data.SetMeasurementPoints(RecX,RecY,Data.GetMeasPosZ());
    Data.SetSourcePoints(SourceX,SourceY,Data.GetSourcePosZ());
    Data.WriteNetCDF(ncfilename + ".rot.nc");
    Data.WriteSourcePoints(ncfilename + "_rot.sor.vtk");
    Data.WriteMeasurementPoints(ncfilename + "_rot.rec.vtk");

  }
