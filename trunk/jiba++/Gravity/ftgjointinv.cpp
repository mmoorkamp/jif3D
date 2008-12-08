//============================================================================
// Name        : ftgjointinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include "../Global/convert.h"
#include "ThreeDGravityModel.h"
#include "../Inversion/MatrixTools.h"
#include "ReadWriteGravityData.h"
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include "../Inversion/LinearInversion.h"
#include "../ModelBase/VTKTools.h"

namespace atlas = boost::numeric::bindings::atlas;
namespace ublas = boost::numeric::ublas;

int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(true, true);

    jiba::ThreeDGravityModel::tTensorMeasVec FTGData;
    jiba::ThreeDGravityModel::tScalarMeasVec ScalarData;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string modelfilename, ftgdatafilename, scalardatafilename;
    std::cout << "Mesh Filename: ";
    std::cin >> modelfilename;
    //we read in a complete modelfile, but we only use the mesh information
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the ftg data and read it in
    std::cout << "FTG Data Filename: ";
    std::cin >> ftgdatafilename;
    jiba::ReadTensorGravityMeasurements(ftgdatafilename, FTGData, PosX, PosY,
        PosZ);

    //get the name of the file containing the scalar data and read it in
    //we assume the same measurement positions for ftg and scalar data
    std::cout << "Scalar Data Filename: ";
    std::cin >> scalardatafilename;
    jiba::ReadScalarGravityMeasurements(scalardatafilename, ScalarData, PosX,
        PosY, PosZ);
    //create a vector for scalar and ftg data
    const size_t nmeas = ScalarData.size();
    jiba::rvec DataVector(nmeas * 10), DataError(nmeas * 10);
    //copy the ftg data to the beginning
    for (size_t i = 0; i < nmeas; ++i)
      {
        copy(FTGData.at(i).data().begin(), FTGData.at(i).data().end(),
            DataVector.begin() + i * 9);
      }
    copy(ScalarData.begin(), ScalarData.end(), DataVector.begin() + 9 * nmeas);
    //std::fill_n(DataError.begin(),DataError.size(),1.0);
    //transform(DataVector.begin(),DataVector.end(),DataError.begin(),boost::bind(std::multiplies<double>(),_1,0.02));
    for (size_t i = 0; i < nmeas * 10; ++i)
      {
        DataError( i) = 0.02 * DataVector(i);
      }
    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }
    //calculate the response, we don't actually care about the densities
    //and the model response, but we need the sensitivity matrix
    Model.CalcTensorGravity();
    Model.CalcGravity();
    Model.SaveTensorMeasurements(modelfilename + ".new.nc");
    jiba::rmat FTGSensitivities(Model.GetTensorSensitivities());
    jiba::rmat ScalarSensitivities(Model.GetScalarSensitivities());
    const size_t nmod = FTGSensitivities.size2();

    jiba::rmat JointSensitivities(DataVector.size(), nmod);
    ublas::matrix_range<jiba::rmat> FTGPart(JointSensitivities, ublas::range(0,
        FTGSensitivities.size1()), ublas::range(0, nmod));
    ublas::matrix_range<jiba::rmat> ScalarPart(JointSensitivities,
        ublas::range(nmeas * 9, nmeas * 10), ublas::range(0, nmod));

    FTGPart = FTGSensitivities;
    ScalarPart = ScalarSensitivities;

    jiba::rvec ModelWeight(nmod);
    for (size_t i = 0; i < nmod; ++i)
      {
        ModelWeight( i) = 1.0;
      }

    jiba::rvec InvModel;

    jiba::DataSpaceInversion()(JointSensitivities, DataVector, ModelWeight,
        DataError, 1.0, InvModel);

    std::copy(InvModel.begin(), InvModel.end(), Model.SetDensities().origin());
    Model.CalcTensorGravity();
    Model.CalcGravity();
    Model.SaveTensorMeasurements(modelfilename + ".inv_ftg.nc");
    Model.SaveScalarMeasurements(modelfilename + ".inv_scal.nc");
    Model.WriteVTK(modelfilename + ".inv.vtk");
  }
