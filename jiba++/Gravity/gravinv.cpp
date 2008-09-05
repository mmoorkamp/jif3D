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

double DepthTerm(const double z, const double z0)
  {
    return pow(z + z0, -3);
  }
double DerivTerm(const double z, const double z0)
  {
    return -3 * pow(z + z0, -4);
  }

double FitZ0(const jiba::rvec &SummedSens,
    const jiba::ThreeDGravityModel &Model)
  {
    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    //the index of the center horizontal cell
    size_t startindex = zsize * ysize * xsize / 2 + zsize * ysize / 2;
    jiba::rvec sensprofile(zsize), zvalues(zsize);
    partial_sum(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
        zvalues.begin());
    copy(SummedSens.begin() + startindex, SummedSens.begin() + startindex
        + zsize, sensprofile.begin());
    transform(sensprofile.begin(), sensprofile.end(), sensprofile.begin(),
        boost::bind(std::divides<double>(), _1, *std::max_element(
            sensprofile.begin(), sensprofile.end())));
    double startz = zvalues(zsize - 1) * 1.1;

    const size_t ndata = zvalues.size();
    jiba::rmat sens(zsize, 1);
    jiba::rvec weights(1);
    jiba::rvec InvModel(1), DeltaModel(1);
    weights(0) = 1.0;
    InvModel(0) = startz;
    const size_t iterations = 100;
    jiba::rvec error(ndata), delta(ndata), calculated(ndata);
    std::fill_n(error.begin(), zsize, 1.0);
    const double evalthresh = 1e-6;
    std::ofstream outfile("match.out");
    double stepsize = 1e6;
    size_t i = 0;
    while (stepsize > 0.1 && i < iterations)
      {
        for (size_t j = 0; j < zsize; ++j)
          {
            calculated(j) = DepthTerm(zvalues(j),InvModel(0));
            delta(j) = sensprofile(j) - calculated(j);
          }
        double maximum =
            *std::max_element(calculated.begin(), calculated.end());
        transform(calculated.begin(), calculated.end(), calculated.begin(),
            boost::bind(std::divides<double>(), _1, maximum));
        for (size_t j = 0; j < zsize; ++j)
          {
            sens(j, 0) = DerivTerm(zvalues(j), InvModel(0)) / maximum;
            delta(j) = sensprofile(j) - calculated(j);
          }
        for (size_t j = 0; j < zsize; ++j)
          {
            outfile << zvalues(j) << " " << sensprofile(j) << " "
                << calculated(j) << std::endl;
          }

        jiba::DataSpaceInversion()(sens, delta, weights, error, evalthresh,0.0,
            DeltaModel);

        stepsize = boost::numeric::ublas::norm_2(DeltaModel);
        InvModel -= DeltaModel;
        std::cout << InvModel << " " << DeltaModel << " " << stepsize
            << std::endl;

        outfile << std::endl << std::endl;
        ++i;
      }

    return InvModel(0);
  }

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
    //read in the depth parameter for the sensitivity depth scaling

    //std::cout << "Enter z0: ";
    //std::cin >> z0;

    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    const size_t nmod = xsize * ysize * zsize;
    jiba::rvec WeightVector(zsize), ModelWeight(nmod), SummedSens(nmod),
        DataError(nmeas);
    jiba::rvec DataVec(nmeas);
    std::copy(Data.begin(), Data.end(), DataVec.begin());
    //here we calculate the sensitivity summed over all data
    //this is to find an appropriate depth scaling and might be removed in the future
    std::fill_n(SummedSens.begin(),nmod,0.0);
    for (size_t i = 0; i < nmeas; ++i)
      {
        DataError(i) = 0.02 * DataVec(i);
        SummedSens += boost::numeric::ublas::matrix_row<jiba::rmat>(
            Sensitivities, i);
      }
    jiba::rvec MiddleSens(boost::numeric::ublas::matrix_row<jiba::rmat>(
        Sensitivities, nmeas/2));
    double z0 = FitZ0(MiddleSens, Model);
    std::cout << "Estimated z0: " << z0 << std::endl;
    //for output we copy the summed sensitivities into a model type structure
    jiba::ThreeDModelBase::t3DModelData SensModel(
        boost::extents[xsize][ysize][zsize]);
    std::copy(SummedSens.begin(), SummedSens.end(), SensModel.data());

    jiba::Write3DModelToVTK(modelfilename + ".sens.vtk", "summed_sens",
        Model.GetXCellSizes(), Model.GetYCellSizes(), Model.GetZCellSizes(),
        SensModel);
    //calculate the depth scaling
    jiba::ConstructDepthWeighting(Model.GetXCellSizes(), Model.GetYCellSizes(),
        Model.GetZCellSizes(), z0, WeightVector);
    //and output the scaling weights
    std::ofstream weightfile("weights.out");
    std::copy(WeightVector.begin(), WeightVector.end(),
        std::ostream_iterator<double>(weightfile, "\n"));
    //the WeightVector only has length zsize, one entry for each depth level
    //the inversion routine needs a vector with a weight for each model parameter
    // the weights only depend on the depth of the cell
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
    Inversion(Sensitivities, DataVec, ModelWeight, DataError, evalthresh,1.0,
        InvModel);

    //write out the scaled sensitivities
    std::fill_n(SummedSens.begin(),nmod,0.0);
    for (size_t i = 0; i < nmeas; ++i)
      {
        SummedSens += boost::numeric::ublas::matrix_column<const jiba::rmat>(
            Inversion.GetFilteredSens(), i);
      }
    std::copy(SummedSens.begin(), SummedSens.end(), SensModel.data());
    jiba::Write3DModelToVTK(modelfilename + ".sens_fil.vtk", "filtered_sens",
        Model.GetXCellSizes(), Model.GetYCellSizes(), Model.GetZCellSizes(),
        SensModel);
    for (size_t i = 0; i < nmeas; ++i)
      {
        boost::numeric::ublas::matrix_column<const jiba::rmat> filcolumn(
            Inversion.GetFilteredSens(), i);
        boost::numeric::ublas::matrix_row<jiba::rmat> senscolumn(Sensitivities,
            i);
        std::copy(filcolumn.begin(), filcolumn.end(), SensModel.data());
        jiba::Write3DModelToVTK(modelfilename + ".sensfil_data" + stringify(i)
            + ".vtk", "filtered_sens", Model.GetXCellSizes(),
            Model.GetYCellSizes(), Model.GetZCellSizes(), SensModel);
        std::copy(senscolumn.begin(), senscolumn.end(), SensModel.data());
        jiba::Write3DModelToVTK(modelfilename + ".sens_data" + stringify(i)
            + ".vtk", "filtered_sens", Model.GetXCellSizes(),
            Model.GetYCellSizes(), Model.GetZCellSizes(), SensModel);
      }

    //copy the inversion model in the Model class for plotting
    std::copy(InvModel.begin(), InvModel.end(), Model.SetDensities().origin());
    //calculate the predicted data
    Model.CalcGravity();
    //and write out the data and model
    Model.SaveScalarMeasurements(modelfilename + ".inv_data.nc");
    Model.PlotScalarMeasurements(modelfilename + ".inv.plot");
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
  }
