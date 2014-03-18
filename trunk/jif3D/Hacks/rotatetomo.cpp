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
#include "../Tomo/ThreeDSeismicModel.h"
#include "../ModelBase/VTKTools.h"

int main()
  {
    std::string ncfilename = jif3D::AskFilename("Name of netcdf file: ");


    jif3D::rvec Data, Error;
    jif3D::ThreeDSeismicModel Model;
    jif3D::ReadTraveltimes(ncfilename,Data,Error,Model);

    std::cout << "Rotation angle [degree]: ";
    double dangle = 0.0;
    std::cin >> dangle;
    double rangle = dangle / 180.0 * boost::math::constants::pi<double>();

    jif3D::ThreeDModelBase::tMeasPosVec SourceX(Model.GetSourcePosX());
    jif3D::ThreeDModelBase::tMeasPosVec SourceY(Model.GetSourcePosY());
    jif3D::ThreeDModelBase::tMeasPosVec SourceZ(Model.GetSourcePosZ());

    jif3D::ThreeDModelBase::tMeasPosVec RecX(Model.GetMeasPosX());
    jif3D::ThreeDModelBase::tMeasPosVec RecY(Model.GetMeasPosY());
    jif3D::ThreeDModelBase::tMeasPosVec RecZ(Model.GetMeasPosZ());

    jif3D::ThreeDSeismicModel::tIndexVec SourceIndices(Model.GetSourceIndices());
    jif3D::ThreeDSeismicModel::tIndexVec ReceiverIndices(Model.GetReceiverIndices());

    Model.ClearMeasurementPoints();
    Model.ClearSourcePos();
    for (size_t i = 0; i < SourceX.size(); ++i)
      {
        double newx = SourceX.at(i)  * cos(rangle) - SourceY.at(i)  * sin(rangle);
        double newy = SourceX.at(i)  * sin(rangle) + SourceY.at(i)  * cos(rangle);
        SourceX.at(i) = newx;
        SourceY.at(i) = newy;
        Model.AddSource(newx,newy,SourceZ[i]);
      }

    for (size_t i = 0; i < RecX.size(); ++i)
      {
        double newx = RecX.at(i)  * cos(rangle) - RecY.at(i)  * sin(rangle);
        double newy = RecX.at(i)  * sin(rangle) + RecY.at(i)  * cos(rangle);
        Model.AddMeasurementPoint(newx,newy,RecZ[i]);
      }

    for (size_t i = 0; i < SourceIndices.size(); ++i)
    {
    	Model.AddMeasurementConfiguration(SourceIndices[i],ReceiverIndices[i]);
    }



    jif3D::Write3DDataToVTK(ncfilename + "_rot.statpos.vtk", "Station",
        jif3D::rvec(SourceX.size(), 1.0), SourceX, SourceY, SourceZ);
    jif3D::SaveTraveltimes(ncfilename+".rot.nc",Data,Error,Model);
  }
