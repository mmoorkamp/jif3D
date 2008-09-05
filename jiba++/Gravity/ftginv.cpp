//============================================================================
// Name        : ftginv.cpp
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

void SumSensitivitiesAndPlot(const jiba::rmat &Sens, const size_t startindex,
    const std::string &filename, const jiba::ThreeDGravityModel &Model)
  {
    const size_t nmod = Sens.size2();
    const size_t ndata = Sens.size1();
    jiba::rvec SummedSens(nmod);
    std::fill_n(SummedSens.begin(),nmod,0.0);
    for (size_t i = startindex; i < ndata; i += 9)
      {
        SummedSens += boost::numeric::ublas::matrix_row<const jiba::rmat>(Sens,
            i);
      }
    //for output we copy the summed sensitivities into a model type structure
    jiba::ThreeDModelBase::t3DModelData
        SensModel(
            boost::extents[Model.GetDensities().shape()[0]][Model.GetDensities().shape()[1]][Model.GetDensities().shape()[2]]);
    std::copy(SummedSens.begin(), SummedSens.end(), SensModel.data());

    jiba::Write3DModelToVTK(filename, "summed_sens", Model.GetXCellSizes(),
        Model.GetYCellSizes(), Model.GetZCellSizes(), SensModel);
  }

int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(false, true);

    jiba::ThreeDGravityModel::tTensorMeasVec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string modelfilename, datafilename;
    std::cout << "Mesh Filename: ";
    std::cin >> modelfilename;
    //we read in a complete modelfile, but we only use the mesh information
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;
    jiba::ReadTensorGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
    jiba::rvec DataVector(Data.size() * 9), DataError(Data.size() * 9);
    for (size_t i = 0; i < Data.size(); ++i)
      {
        copy(Data.at(i).data().begin(), Data.at(i).data().end(),
            DataVector.begin() + i * 9);
      }
    for (size_t i = 0; i < DataVector.size(); ++i)
      {
        DataError(i) = 0.02 * DataVector(i);
      }
    const size_t nmeas = PosX.size();

    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }
    //calculate the response, we don't actually care about the densities
    //and the model response, but we need the sensitivity matrix
    Model.CalcTensorGravity();
    Model.SaveTensorMeasurements(modelfilename + ".new.nc");
    jiba::rmat Sensitivities(Model.GetTensorSensitivities());
    //write out sensitivities for the 9 tensor elements
    for (size_t i = 0; i < 9; ++i)
      {
        SumSensitivitiesAndPlot(Sensitivities, i, modelfilename + stringify(i)
            + ".vtk", Model);
      }
    //depth weighting
    double z0;
    std::cout << "Enter z0: ";
    std::cin >> z0;
    const size_t nmod = Sensitivities.size2();
    const size_t zsize = Model.GetDensities().shape()[2];
    jiba::rvec WeightVector(zsize), ModelWeight(nmod);
    jiba::ConstructDepthWeighting(Model.GetXCellSizes(), Model.GetYCellSizes(),
        Model.GetZCellSizes(), z0, WeightVector);
    std::ofstream weightfile("weights.out");
    std::copy(WeightVector.begin(), WeightVector.end(),
        std::ostream_iterator<double>(weightfile, "\n"));
    for (size_t i = 0; i < nmod; ++i)
      {
        ModelWeight(i) =WeightVector(i % zsize);
      }

    //we can play around with the threshold for the included eigenvalues
    double evalthresh;
    std::cout << "Eigenvalue threshold: ";
    std::cin >> evalthresh;
    jiba::rvec InvModel;

    jiba::DataSpaceInversion()(Sensitivities, DataVector, ModelWeight,DataError,
        evalthresh, 1.0, InvModel);

    std::copy(InvModel.begin(), InvModel.end(), Model.SetDensities().origin());
    Model.CalcTensorGravity();
    Model.SaveTensorMeasurements(modelfilename + ".inv_data.nc");
    Model.PlotTensorMeasurements(modelfilename + ".inv.plot");
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
  }
