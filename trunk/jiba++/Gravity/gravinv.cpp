//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


/*! \file gravinv.cpp
 * Invert scalar gravity data. The program reads in a model file that specifies the inversion mesh and
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
#include "FullSensitivityGravityCalculator.h"
#include "../Inversion/LinearInversion.h"
#include "../ModelBase/VTKTools.h"
#include "DepthWeighting.h"
#include "../Global/FatalException.h"

namespace ublas = boost::numeric::ublas;

int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(true);

    jiba::rvec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

    std::string modelfilename, datafilename;
    std::cout << "Starting model Filename: ";
    std::cin >> modelfilename;
    //we read in the starting modelfile
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;
    jiba::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
    if (Data.empty())
      throw jiba::FatalException("No measurements defined");
    const size_t nmeas = PosX.size();
    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }
    //calculate the response of the starting model
    boost::shared_ptr<jiba::FullSensitivityGravityCalculator>
        ScalarCalculator(jiba::CreateGravityCalculator<
            jiba::FullSensitivityGravityCalculator>::MakeScalar());
    jiba::rvec StartingData(ScalarCalculator->Calculate(Model));

    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    const size_t nmod = xsize * ysize * zsize;

    jiba::rmat AllSens(ScalarCalculator->GetSensitivities());
    jiba::rmat Sensitivities(ublas::matrix_range<jiba::rmat>(AllSens,
        ublas::range(0, nmeas), ublas::range(0, nmod)));
    std::cout << "Gridded model size: " << nmod << " Complete size: "
        << AllSens.size2() << std::endl;
    jiba::rvec WeightVector(zsize), ModelWeight(AllSens.size2());
    jiba::rvec DataDiffVec(nmeas), DataError(nmeas);
    std::transform(Data.begin(), Data.end(), StartingData.begin(),
        DataDiffVec.begin(), std::minus<double>());
    const double errorlevel = 0.02;
    std::transform(Data.begin(), Data.end(), DataError.begin(), boost::bind(
        std::multiplies<double>(), _1, errorlevel));

    jiba::rvec SensProfile;
    jiba::ExtractMiddleSens(Model, Sensitivities, 1, SensProfile);
    double z0 = FitZ0(SensProfile, Model.GetZCellSizes(), jiba::WeightingTerm(
        -3));
    std::cout << "Estimated z0: " << z0 << std::endl;

    //calculate the depth scaling
    jiba::ConstructDepthWeighting(Model.GetZCellSizes(), z0, WeightVector,
        jiba::WeightingTerm(-3));
    //and output the scaling weights
    std::ofstream weightfile("weights.out");
    std::copy(WeightVector.begin(), WeightVector.end(), std::ostream_iterator<
        double>(weightfile, "\n"));
    //the WeightVector only has length zsize, one entry for each depth level
    //the inversion routine needs a vector with a weight for each model parameter
    // the weights only depend on the depth of the cell
    // we also have the background layers, which get a weight of 0, so they are not changed in the inversion
    std::fill_n(ModelWeight.begin(), ModelWeight.size(), 0);
    for (size_t i = 0; i < nmod; ++i)
      {
        ModelWeight( i) = WeightVector(i % zsize);
      }

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;

    //here comes the core inversion
    jiba::rvec InvModel(AllSens.size2());
    std::fill(InvModel.begin(), InvModel.end(), 0.0);
    jiba::DataSpaceInversion Inversion;
    Inversion(AllSens, DataDiffVec, ModelWeight, DataError, lambda,
        InvModel);

    //add the result of the inversion to the starting model
    //we only add the gridded part, the  background is always 0 due to the weighting
    std::transform(InvModel.begin(), InvModel.begin() + nmod,
        Model.SetDensities().origin(), Model.SetDensities().origin(),
        std::plus<double>());
    //calculate the predicted data
    jiba::rvec InvData(ScalarCalculator->Calculate(Model));
    //and write out the data and model
    jiba::Write3DDataToVTK(modelfilename + ".inv_data.vtk", "grav_accel",
        InvData, Model.GetMeasPosX(), Model.GetMeasPosY(), Model.GetMeasPosZ());
    jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_data.nc",
        InvData, Model.GetMeasPosX(), Model.GetMeasPosY(),
        Model.GetMeasPosZ());
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");

    std::cout << std::endl;
    jiba::ThreeDModelBase::t3DModelData SensModel(
        boost::extents[xsize][ysize][zsize]);
    for (size_t i = 0; i < nmeas; ++i)
      {
        //the filtered sensitivities include the background layers
        boost::numeric::ublas::matrix_row<jiba::rmat> filrow(AllSens, i);
        boost::numeric::ublas::matrix_row<jiba::rmat> sensrow(Sensitivities, i);
        //we are only interested in the sensitivities for the gridded part
        std::copy(filrow.begin(), filrow.begin() + nmod, SensModel.data());
        jiba::Write3DModelToVTK(modelfilename + ".sensfil_data"
            + jiba::stringify(i) + ".vtk", "filtered_sens",
            Model.GetXCellSizes(), Model.GetYCellSizes(),
            Model.GetZCellSizes(), SensModel);
        std::copy(sensrow.begin(), sensrow.begin() + nmod, SensModel.data());
        jiba::Write3DModelToVTK(modelfilename + ".sens_data" + jiba::stringify(
            i) + ".vtk", "raw_sens", Model.GetXCellSizes(),
            Model.GetYCellSizes(), Model.GetZCellSizes(), SensModel);
      }
  }
