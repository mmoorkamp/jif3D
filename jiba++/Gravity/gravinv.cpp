//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "../Inversion/LinearInversion.h"
#include "../ModelBase/VTKTools.h"
int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(true);

    jiba::ThreeDGravityModel::tScalarMeasVec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string modelfilename, datafilename;
    std::cout << "Mesh Filename: ";
    std::cin >> modelfilename;
    //we read in a complete modelfile, but we only use the mesh information
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;
    jiba::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
    const size_t nmeas = PosX.size();
    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }
    //calculate the response, we don't actually care about the densities
    //and the model response, but we need the sensitivity matrix
    Model.CalcGravity();
    Model.SaveScalarMeasurements(modelfilename + ".new.nc");
    jiba::rmat Sensitivities(Model.GetScalarSensitivities());

    double z0;
    std::cout << "Enter z0: ";
    std::cin >> z0;

    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    const size_t nmod = xsize * ysize * zsize;
    jiba::rvec WeightVector(zsize), ModelWeight(nmod), SummedSens(nmod);
    SummedSens *= 0.0;
    for (size_t i = 0; i < nmeas; ++i)
      {
        SummedSens += boost::numeric::ublas::matrix_row<jiba::rmat>(Sensitivities, i);
      }
    jiba::ThreeDModelBase::t3DModelData SensModel(
            boost::extents[xsize][ysize][zsize]);
    std::copy(SummedSens.begin(),SummedSens.end(),SensModel.data());

    jiba::Write3DModelToVTK(modelfilename + ".sens.vtk", "summed_sens",
                Model.GetXCellSizes(), Model.GetYCellSizes(),
                Model.GetZCellSizes(), SensModel);

    jiba::ConstructDepthWeighting(Model.GetXCellSizes(), Model.GetYCellSizes(),
        Model.GetZCellSizes(), z0, WeightVector);
    std::ofstream weightfile("weights.out");
    std::copy(WeightVector.begin(),WeightVector.end(),std::ostream_iterator<double>(weightfile,"\n"));
    for (size_t i = 0; i < nmod; ++i)
      {
        ModelWeight(i) =WeightVector(i % zsize);
      }

    //we can play around with the threshold for the included eigenvalues
    double evalthresh;

    std::cout << "Eigenvalue threshold: ";
    std::cin >> evalthresh;
    jiba::rvec InvModel;
    jiba::rvec DataVec(Data.size());
    std::copy(Data.begin(),Data.end(),DataVec.begin());
    jiba::DataSpaceInversion(Sensitivities,DataVec,ModelWeight,evalthresh,InvModel);

    std::copy(InvModel.begin(), InvModel.end(), Model.SetDensities().origin());
    Model.CalcGravity();
    Model.SaveScalarMeasurements(modelfilename + ".inv_data.nc");
    Model.PlotScalarMeasurements(modelfilename + ".inv.plot");
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
  }
