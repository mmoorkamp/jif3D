//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file Invert scalar gravity data. The program reads in a model file that specifies the inversion mesh and
 * a file with the input data. Two parameters have to be specified for the inversion: z0 determines the depth scaling
 * applied to counteract the decay of the sensisitivities. The eigenvalue threshold is the relative threshold with respect
 * to the largest eigenvalue of the sensitivity matrix below which the eigenvalues are considered 0.
 *
 * The program outputs the calculated data and inversion model in netcdf and vtk file format.
 */

#include <iostream>
#include <fstream>
#include <string>
#include "../Global/convert.h"
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "../Inversion/LinearInversion.h"
#include "../ModelBase/VTKTools.h"
#include "DepthWeighting.h"

namespace ublas = boost::numeric::ublas;

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

    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    const size_t nmod = xsize * ysize * zsize;
    jiba::rmat AllSens(Model.GetScalarSensitivities());
    jiba::rmat Sensitivities(ublas::matrix_range<jiba::rmat>(AllSens,
            ublas::range(0, nmeas), ublas::range(0, nmod)));


    jiba::rvec WeightVector(zsize), ModelWeight(AllSens.size2()), DataError(nmeas);
    jiba::rvec DataVec(nmeas);
    std::copy(Data.begin(), Data.end(), DataVec.begin());
    //here we calculate the sensitivity summed over all data
    //this is to find an appropriate depth scaling and might be removed in the future

    for (size_t i = 0; i < nmeas; ++i)
      {
        DataError(i) = 0.02 * DataVec(i);
      }
    jiba::rvec MiddleSens(boost::numeric::ublas::matrix_row<jiba::rmat>(
        Sensitivities, nmeas / 2));
    double z0 = FitZ0(MiddleSens, Model, jiba::WeightingTerm(-3));
    std::cout << "Estimated z0: " << z0 << std::endl;

    //calculate the depth scaling
    jiba::ConstructDepthWeighting(Model.GetXCellSizes(), Model.GetYCellSizes(),
        Model.GetZCellSizes(), z0, WeightVector, jiba::WeightingTerm(-3));
    //and output the scaling weights
    std::ofstream weightfile("weights.out");
    std::copy(WeightVector.begin(), WeightVector.end(),
        std::ostream_iterator<double>(weightfile, "\n"));
    //the WeightVector only has length zsize, one entry for each depth level
    //the inversion routine needs a vector with a weight for each model parameter
    // the weights only depend on the depth of the cell
    std::fill_n(ModelWeight.begin(),ModelWeight.size(),0.0);
    for (size_t i = 0; i < nmod; ++i)
      {
        ModelWeight(i) =WeightVector(i % zsize);
      }

    //we can play around with the threshold for the included eigenvalues
    double evalthresh;
    std::cout << "Eigenvalue threshold: ";
    std::cin >> evalthresh;

    //here comes the core inversion
    jiba::rvec InvModel;
    jiba::DataSpaceInversion Inversion;
    Inversion(AllSens, DataVec, ModelWeight, DataError, evalthresh, 1.0,
        InvModel);

    jiba::ThreeDModelBase::t3DModelData SensModel(
            boost::extents[xsize][ysize][zsize]);
    for (size_t i = 0; i < nmeas; ++i)
      {
        boost::numeric::ublas::matrix_column<const jiba::rmat> filcolumn(
            Inversion.GetFilteredSens(), i);
        boost::numeric::ublas::matrix_row<jiba::rmat> senscolumn(Sensitivities,
            i);
        std::copy(filcolumn.begin(), filcolumn.end(), SensModel.data());
        jiba::Write3DModelToVTK(modelfilename + ".sensfil_data"
            + jiba::stringify(i) + ".vtk", "filtered_sens",
            Model.GetXCellSizes(), Model.GetYCellSizes(),
            Model.GetZCellSizes(), SensModel);
        std::copy(senscolumn.begin(), senscolumn.end(), SensModel.data());
        jiba::Write3DModelToVTK(modelfilename + ".sens_data" + jiba::stringify(
            i) + ".vtk", "filtered_sens", Model.GetXCellSizes(),
            Model.GetYCellSizes(), Model.GetZCellSizes(), SensModel);
      }

    //copy the inversion model in the Model class for plotting
    std::copy(InvModel.begin(), InvModel.begin() + nmod, Model.SetDensities().origin());
    //calculate the predicted data
    jiba::ThreeDGravityModel::tScalarMeasVec InvData(Model.CalcGravity());
    //and write out the data and model
    jiba::Write3DDataToVTK(modelfilename + ".inv_data.vtk", "grav_accel", InvData,
        Model.GetMeasPosX(), Model.GetMeasPosY(),
        Model.GetMeasPosZ());
    Model.SaveScalarMeasurements(modelfilename + ".inv_data.nc");
    Model.PlotScalarMeasurements(modelfilename + ".inv.plot");
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
  }
