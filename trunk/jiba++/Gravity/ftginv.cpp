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
#include "DepthWeighting.h"

namespace atlas = boost::numeric::bindings::atlas;

void SumSensitivitiesAndPlot(const jiba::rmat &Sens, const size_t startindex,
    const std::string &filename, const jiba::ThreeDGravityModel &Model)
  {
    const size_t nmod = Sens.size2();
    const size_t ndata = Sens.size1();
    jiba::rvec SummedSens(nmod);
    std::fill_n(SummedSens.begin(), nmod, 0.0);
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
    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];

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
        SumSensitivitiesAndPlot(Sensitivities, i, modelfilename
            + jiba::stringify(i) + ".vtk", Model);
      }
    //depth weighting
    const double midx = Model.GetXCoordinates()[Model.GetXCoordinates().size()
        - 1] / 2.0;
    const double midy = Model.GetYCoordinates()[Model.GetYCoordinates().size()
        - 1] / 2.0;

    jiba::rvec distances(nmeas);
    for (size_t i = 0; i < nmeas; ++i)
      {
        distances(i) = sqrt(pow(Model.GetMeasPosX()[i] - midx,2)+pow(Model.GetMeasPosY()[i] - midy,2));
      }
    const size_t midindex = distance(distances.begin(), std::min_element(
        distances.begin(), distances.end()));

    boost::array<jiba::ThreeDModelBase::t3DModelData::index,3> modelindex(
        Model.FindAssociatedIndices(Model.GetMeasPosX()[midindex],
            Model.GetMeasPosY()[midindex], 0.0));

    jiba::rvec MiddleSens(boost::numeric::ublas::matrix_row<jiba::rmat>(
        Sensitivities, midindex / 2 * 9));

    jiba::rvec SensProfile(zsize);
    for (size_t i = 0; i < zsize; ++i)
      {
        SensProfile(i) = MiddleSens(zsize*modelindex[1] + (zsize*ysize)*modelindex[0] + i );
      }

    const double decayexponent = -4.0;
    double z0 = FitZ0(SensProfile, Model.GetZCellSizes(), jiba::WeightingTerm(decayexponent));
    std::cout << "Estimated z0: " << z0 << std::endl;
    const size_t nmod = Sensitivities.size2();

    jiba::rvec WeightVector(zsize), ModelWeight(nmod);
    jiba::ConstructDepthWeighting(Model.GetXCellSizes(), Model.GetYCellSizes(),
        Model.GetZCellSizes(), z0, WeightVector, jiba::WeightingTerm(decayexponent));
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

    jiba::DataSpaceInversion()(Sensitivities, DataVector, ModelWeight,
        DataError, evalthresh, 1.0, InvModel);

    std::copy(InvModel.begin(), InvModel.end(), Model.SetDensities().origin());
    //we want to save the resulting predicted data
    jiba::ThreeDGravityModel::tTensorMeasVec FTGData(Model.CalcTensorGravity());
    Model.SaveTensorMeasurements(modelfilename + ".inv_data.nc");
    Model.PlotTensorMeasurements(modelfilename + ".inv.plot");
    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");

    jiba::Write3DTensorDataToVTK(modelfilename + ".ftgdata.vtk", "U", FTGData,
        Model.GetMeasPosX(), Model.GetMeasPosY(),
        Model.GetMeasPosZ());
  }
