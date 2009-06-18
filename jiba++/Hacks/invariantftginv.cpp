//============================================================================
// Name        : ftginv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/bind.hpp>
#include "../Global/convert.h"
#include "../Global/Interpolate.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Inversion/MatrixTools.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Inversion/LinearInversion.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../Gravity/DepthWeighting.h"
#include "../Gravity/FullSensitivityGravityCalculator.h"

namespace atlas = boost::numeric::bindings::atlas;
namespace ublas = boost::numeric::ublas;

//write the sensitivities for each measurement to a file
void WriteSensitivities(const std::string &nameroot,
    const std::string &sensname, const jiba::rmat &Sens,
    const jiba::ThreeDGravityModel &Model)
  {
    //create a data structure that mimics the geometry of the models
    jiba::ThreeDModelBase::t3DModelData
        SensModel(
            boost::extents[Model.GetDensities().shape()[0]][Model.GetDensities().shape()[1]][Model.GetDensities().shape()[2]]);
    const size_t ngrid = Model.GetDensities().num_elements();
    //for each measurement
    for (size_t i = 0; i < Model.GetMeasPosX().size(); ++i)
      {
        //extract the corresponding row of the sensitivity matrix
        boost::numeric::ublas::matrix_row<const jiba::rmat> sensrow(Sens, i);
        //copy to the Model structure to map to its geometric position
        std::copy(sensrow.begin(), sensrow.begin() + ngrid, SensModel.data());
        //write out a .vtk and a netcdf file
        jiba::Write3DModelToVTK(nameroot + sensname + jiba::stringify(i)
            + ".vtk", sensname, Model.GetXCellSizes(), Model.GetYCellSizes(),
            Model.GetZCellSizes(), SensModel);
        jiba::Write3DModelToNetCDF(nameroot + sensname + jiba::stringify(i)
            + ".nc", sensname, " ", Model.GetXCellSizes(),
            Model.GetYCellSizes(), Model.GetZCellSizes(), SensModel);
      }

  }

double CalcInvariant(const jiba::rvec &Data, const size_t index)
  {
    return Data(index) * Data(index + 4) + Data(index + 4) * Data(index + 8)
        + Data(index) * Data(index + 8) - Data(index + 3) * Data(index + 1)
        - Data(index + 7) * Data(index + 5) - Data(index + 2) * Data(index + 6);
  }

double CalcMisfit(const jiba::rvec &MeasData, const jiba::rvec &CurrData,
    const jiba::rvec &Model, const double lambda)
  {
    jiba::rvec DiffVector(MeasData.size());
    std::transform(MeasData.begin(), MeasData.end(), CurrData.begin(),
        DiffVector.begin(), std::minus<double>());
    std::transform(DiffVector.begin(),DiffVector.end(),MeasData.begin(),DiffVector.begin(),std::divides<double>());
    return sqrt(boost::numeric::ublas::norm_2(DiffVector) + lambda
        * boost::numeric::ublas::norm_2(Model));
  }

void LineSearch(jiba::rvec &DeltaModel, jiba::ThreeDGravityModel &Model,
    jiba::rvec &DataVector, boost::shared_ptr<
        jiba::FullSensitivityGravityCalculator> GravityCalculator,
    const double lambda)
  {
    const size_t nmod = Model.GetDensities().size();
    jiba::rvec LastModel(nmod);
    std::copy(Model.GetDensities().origin(), Model.GetDensities().origin()
        + nmod, LastModel.begin());
    std::transform(DeltaModel.begin(), DeltaModel.begin() + nmod,
        Model.SetDensities().origin(), Model.SetDensities().origin(),
        std::plus<double>());
    jiba::rvec FTGData(GravityCalculator->Calculate(Model));

    jiba::rvec CurrModel(nmod);
    std::copy(Model.SetDensities().origin(), Model.SetDensities().origin()
        + nmod, CurrModel.begin());
    size_t ndata = DataVector.size();
    jiba::rvec CurrData(ndata), DiffVector(ndata);
    for (size_t i = 0; i < ndata; ++i)
      {
        CurrData( i) = CalcInvariant(FTGData, i * 9);
      }

    double misfit = CalcMisfit(DataVector, CurrData, CurrModel, lambda);

    double mu1 = 10;
    std::transform(LastModel.begin(), LastModel.begin() + nmod,
        DeltaModel.begin(), Model.SetDensities().origin(), boost::bind(
            std::plus<double>(), _1, boost::bind(std::multiplies<double>(), _2,
                mu1)));
    std::copy(Model.SetDensities().origin(), Model.SetDensities().origin()
        + nmod, CurrModel.begin());

    FTGData = GravityCalculator->Calculate(Model);
    for (size_t i = 0; i < ndata; ++i)
      {
        CurrData( i) = CalcInvariant(FTGData, i * 9);
      }
    double misfit2 = CalcMisfit(DataVector, CurrData, CurrModel, lambda);

    double mu2 = 0.0001;
    std::transform(LastModel.begin(), LastModel.begin() + nmod,
        DeltaModel.begin(), Model.SetDensities().origin(), boost::bind(
            std::plus<double>(), _1, boost::bind(std::multiplies<double>(), _2,
                mu2)));
    FTGData = GravityCalculator->Calculate(Model);
    for (size_t i = 0; i < ndata; ++i)
      {
        CurrData( i) = CalcInvariant(FTGData, i * 9);
      }
    std::copy(Model.SetDensities().origin(), Model.SetDensities().origin()
        + nmod, CurrModel.begin());

    double misfit3 = CalcMisfit(DataVector, CurrData, CurrModel, lambda);
    double optstep = jiba::QuadraticInterpolation(mu2, misfit3, 1.0, misfit,
        mu1, misfit2);
    std::cout << "In Line search: ";
    std::cout << mu2 << " " << 1.0 << " " << mu1 << std::endl;
    std::cout << misfit3 << " " << misfit << " " << misfit2 << std::endl;
    std::cout << "Optimum step length: " << optstep << std::endl;

    std::transform(DeltaModel.begin(), DeltaModel.begin() + nmod,
        DeltaModel.begin(), boost::bind(std::multiplies<double>(), _1, optstep));
  }

int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model;
    boost::shared_ptr<jiba::FullSensitivityGravityCalculator>
        GravityCalculator(jiba::CreateGravityCalculator<
            jiba::FullSensitivityGravityCalculator>::MakeTensor());
    jiba::rvec Data;
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
    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    const size_t ndata = Data.size();
    const size_t nmeas = Data.size() / 9;
    jiba::rvec DataVector(nmeas), DataError(nmeas);
    for (size_t i = 0; i < nmeas; ++i)
      {
        DataVector( i) = CalcInvariant(Data, i * 9);
      }

    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    const size_t ngrid = xsize * ysize * zsize;
    const size_t nmod = ngrid + Model.GetBackgroundDensities().size();

    jiba::rvec InvModel(nmod);
    std::copy(Model.GetDensities().origin(), Model.GetDensities().origin()
        + ngrid, InvModel.begin());
    std::copy(Model.GetBackgroundDensities().begin(),
        Model.GetBackgroundDensities().end(), InvModel.begin() + ngrid);

    jiba::rvec StartingModel(InvModel);
    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }

    const double errorlevel = 0.02;
    std::fill(DataError.begin(),DataError.end(),errorlevel);


    const size_t maxiter = 30;
    size_t iter = 0;
    double gradnorm = 1e10;
    const double mingradnorm = 0.01;
    double datamisfit = 1e10;
    const double minmisfit = sqrt(boost::numeric::ublas::norm_2(DataError));
    jiba::rmat InvarSens(nmeas, nmod);
    jiba::rvec StartingDataVector(nmeas), DataDiffVector(nmeas);

    std::ofstream misfitfile("misfit.out");

    while (iter < maxiter && gradnorm > mingradnorm && datamisfit > minmisfit)
      {

        std::cout << "\n\n";
        std::cout << " Iteration: " << iter << std::endl;
        //calculate the response of the starting model
        jiba::rvec StartingData(GravityCalculator->Calculate(Model));
        const jiba::rmat &AllSens = GravityCalculator->GetSensitivities();

        for (size_t i = 0; i < nmeas; ++i)
          {
            StartingDataVector( i) = CalcInvariant(StartingData, i * 9);

            for (size_t j = 0; j < nmod; ++j)
              {
                InvarSens(i, j) = AllSens(i * 9, j) * StartingData(i * 9 + 4)
                    + AllSens(i * 9 + 4, j) * StartingData(i * 9)

                + AllSens(i * 9 + 4, j) * StartingData(i * 9 + 8) + AllSens(i
                    * 9 + 8, j) * StartingData(i * 9 + 4)

                + AllSens(i * 9, j) * StartingData(i * 9 + 8) + AllSens(i * 9
                    + 8, j) * StartingData(i * 9)

                - 2.0 * AllSens(i * 9 + 3, j) * StartingData(i * 9 + 3) - 2.0
                    * AllSens(i * 9 + 7, j) * StartingData(i * 9 + 7) - 2.0
                    * AllSens(i * 9 + 2, j) * StartingData(i * 9 + 2);
              }
            boost::numeric::ublas::matrix_row<jiba::rmat> CurrentRow(
                                       InvarSens, i);
            CurrentRow /= DataVector(i);
          }
        WriteSensitivities(modelfilename, "raw_sens",
                            InvarSens, Model);
        std::transform(DataVector.begin(), DataVector.end(),
            StartingDataVector.begin(), DataDiffVector.begin(), std::minus<
                double>());
        std::transform(DataDiffVector.begin(),DataDiffVector.end(),DataVector.begin(),DataDiffVector.begin(),std::divides<double>());
        datamisfit = sqrt(boost::numeric::ublas::norm_2(DataDiffVector));
        std::cout << "Data Misfit: " << datamisfit << std::endl;
        //jiba::rvec CurrModel(ngrid);
        //std::copy(Model.SetDensities().origin(),Model.SetDensities().origin()+ngrid,CurrModel.begin());
        double totalmisfit = CalcMisfit(DataVector, StartingDataVector,
            InvModel, lambda);
        std::cout << "Total Misfit: " << totalmisfit << std::endl;
        misfitfile << iter << " " << datamisfit << " " << totalmisfit
            << std::endl;
        //depth weighting
        jiba::rvec SensProfile;
        jiba::ExtractMiddleSens(Model, InvarSens, 1,SensProfile);

        const double decayexponent = -3.0;
        double z0 = FitZ0(SensProfile, Model.GetZCellSizes(),
            jiba::WeightingTerm(decayexponent));
        std::cout << "Estimated z0: " << z0 << std::endl;

        jiba::rvec WeightVector(zsize), ModelWeight(InvarSens.size2());
        jiba::ConstructDepthWeighting(Model.GetZCellSizes(), z0, WeightVector,
            jiba::WeightingTerm(decayexponent));

        std::fill_n(ModelWeight.begin(), ModelWeight.size(), 0.0);
        for (size_t i = 0; i < ngrid; ++i)
          {
            ModelWeight( i) = WeightVector(i % zsize);
          }
        jiba::rvec LastModel(InvModel);
        InvModel -= StartingModel;
        jiba::DataSpaceInversion()(InvarSens, DataDiffVector, ModelWeight,
            DataError, lambda, InvModel);
        //LineSearch(InvModel, Model, DataVector, GravityCalculator, lambda);
        //add the result of the inversion to the starting model
        //we only add the gridded part, the  background is always 0 due to the weighting
        InvModel += StartingModel;
        gradnorm = sqrt(boost::numeric::ublas::norm_2(InvModel-LastModel));
        std::cout << "Delta Model norm: " << gradnorm << std::endl;
        std::cout << "Data Error: " << minmisfit << std::endl;
        std::copy(InvModel.begin(), InvModel.begin() + ngrid,
            Model.SetDensities().origin());
        WriteSensitivities(modelfilename, "fill_sens",
                            InvarSens, Model);
        //std::fill(InvModel.begin(), InvModel.end(), 0.0);

        iter++;
      }
    //we want to save the resulting predicted data
    jiba::rvec FTGData(GravityCalculator->Calculate(Model));
    jiba::rvec relerror(nmeas);
    for (size_t i = 0; i < nmeas; ++i)
      {
        StartingDataVector( i) = CalcInvariant(FTGData, i * 9);
      }
    std::transform(DataVector.begin(), DataVector.end(),
        StartingDataVector.begin(), DataDiffVector.begin(),
        std::minus<double>());
    std::transform(DataDiffVector.begin(), DataDiffVector.end(),
        DataVector.begin(), relerror.begin(), std::divides<double>());

    jiba::SaveScalarGravityMeasurements(modelfilename + ".meas_invariant.nc",
        DataVector, Model.GetMeasPosX(), Model.GetMeasPosY(),
        Model.GetMeasPosZ());
    jiba::SaveScalarGravityMeasurements(modelfilename + ".inv_invariant.nc",
        StartingDataVector, Model.GetMeasPosX(), Model.GetMeasPosY(),
        Model.GetMeasPosZ());
    jiba::Write3DDataToVTK(modelfilename + ".inv_ftgmisfit.vtk",
        "MisfitInvariant", relerror, Model.GetMeasPosX(), Model.GetMeasPosY(),
        Model.GetMeasPosZ());

    //and we save the models and the tensor elements
    jiba::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc", FTGData,
        Model.GetMeasPosX(), Model.GetMeasPosY(), Model.GetMeasPosZ());

    Model.WriteVTK(modelfilename + ".inv_ftg.vtk");

    jiba::Write3DTensorDataToVTK(modelfilename + ".inv_ftgdata.vtk", "U",
        FTGData, Model.GetMeasPosX(), Model.GetMeasPosY(), Model.GetMeasPosZ());
  }
