//============================================================================
// Name        : test_X3DObjective.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE X3DObjective test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/filesystem.hpp>
#include <boost/make_shared.hpp>

#include "../Inversion/ThreeDModelObjective.h"
#include "../Global/Noise.h"
#include "X3DModel.h"
#include "OneDMTCalculator.h"
#include "OneDMTObjective.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "ReadWriteImpedances.h"
#include "MTTransforms.h"
#include "MTData.h"

BOOST_AUTO_TEST_SUITE (X3DObjective_Suite)

void MakeMTModel(jif3D::X3DModel &Model, jif3D::MTData &Data)
  {
    const size_t xsize = 4;
    const size_t ysize = 3;
    const size_t zsize = 2;
    const size_t nbglayers = 2;
    const size_t nmod = xsize * ysize * zsize;

    Model.SetMeshSize(xsize, ysize, zsize);

    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

    const double deltax = 100.0;
    const double deltay = 100.0;
    const double deltaz = 100.0;
    Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize, deltaz);
    ZCS[1] = deltaz * 2.0;
    Model.SetZCellSizes(ZCS);
    std::fill_n(Model.SetConductivities().origin(), nmod, 0.02);

    Model.SetConductivities()[0][0][0] = 0.025;
    std::fill_n(bg_conductivities.begin(), nbglayers, 0.02);
    for (size_t i = 0; i < nbglayers; ++i)
      {
        bg_conductivities[i] *= 1.05 + i * 0.1;
      }
    bg_thicknesses = ZCS;

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Data.SetFrequencies(
      { 1.0, 2.0});
    for (size_t i = 1; i < xsize - 1; ++i)
      {
        for (size_t j = 1; j < ysize - 1; ++j)
          {
            double currx = Model.GetXCoordinates()[i] + deltax / 3.0;
            double curry = Model.GetYCoordinates()[j] + deltay / 4.0;
            double currz = (j - 1) * deltaz;
            Data.AddMeasurementPoint(currx, curry, currz);
          }
      }
    Data.CompleteObject();
  }

bool Between(const double limit1, const double limit2, const double value)
  {
    const double upper = std::max(limit1, limit2);
    const double lower = std::min(limit1, limit2);
    return (lower <= value) && (upper >= value);
  }

void MakeTitanModel(jif3D::X3DModel &Model, jif3D::MTData &Data)
  {
    const size_t xsize = 4;
    const size_t ysize = 3;
    const size_t zsize = 2;
    const size_t nbglayers = 2;
    const size_t nmod = xsize * ysize * zsize;

    Model.SetMeshSize(xsize, ysize, zsize);
    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);
    std::vector<int> ExIndices(2), EyIndices(2), HIndices(2);

    const double deltax = 50.0;
    const double deltay = 50.0;
    const double deltaz = 50.0;
    Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);

    jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize, deltaz);
    Model.SetZCellSizes(ZCS);
    std::fill_n(Model.SetConductivities().origin(), nmod, 0.02);

    Model.SetConductivities()[0][0][0] = 0.025;
    std::fill_n(bg_conductivities.begin(), nbglayers, 0.02);
    for (size_t i = 0; i < nbglayers; ++i)
      {
        bg_conductivities[i] *= 1.05 + i * 0.1;
      }
    std::fill_n(bg_thicknesses.begin(), nbglayers, deltaz);

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Data.SetFrequencies(
      { 10.0 });
    //Model.SetFrequencies().push_back(2.0);
    //Model.SetFrequencies().push_back(5.0);
    //Model.SetFrequencies().push_back(10.0);
    for (size_t i = 0; i < xsize - 1; ++i)
      {
        for (size_t j = 0; j < ysize - 1; ++j)
          {
            double currx = Model.GetXCoordinates()[i] + deltax / 3.0;
            double curry = Model.GetYCoordinates()[j] + deltay / 4.0;
            double currz = j * deltaz;
            Data.AddMeasurementPoint(currx, curry, currz);
          }
      }

    std::fill_n(HIndices.begin(), HIndices.size(), 0);
    std::fill_n(EyIndices.begin(), EyIndices.size(), ysize - 1);
    // HIndices[0]= (ysize - 1) - 1;  HIndices[1]= 3*(ysize - 1) - 1;
    // EyIndices[0]= (ysize - 1) - 1;  EyIndices[1]= 3*(ysize - 1) - 1;
    ExIndices[0] = (ysize - 1) - 1;
    ExIndices[1] = 3 * (ysize - 1) - 1;
    Data.SetFieldIndices(ExIndices, EyIndices, HIndices, HIndices);
    Data.CompleteObject();
  }

BOOST_AUTO_TEST_CASE (X3D_fail_test)
  {
    jif3D::X3DModel Model;
    jif3D::MTData Observed;
    jif3D::X3DMTCalculator Calculator;
    jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> Objective(Calculator);
    BOOST_CHECK_THROW(Objective.SetObservedData(Observed), jif3D::FatalException);
    BOOST_CHECK_THROW(Objective.SetCoarseModelGeometry(Model), jif3D::FatalException);

    MakeMTModel (Model, Observed);
    BOOST_CHECK_NO_THROW(Objective.SetObservedData(Observed));

    BOOST_CHECK_NO_THROW(Objective.SetCoarseModelGeometry(Model));
    Observed.ClearMeasurementPoints();
    Observed.AddMeasurementPoint(10.0, 12.0, 0.0);
    Observed.AddMeasurementPoint(13.0, 14.0, 30.0);
    Objective.SetCoarseModelGeometry(Model);

    BOOST_CHECK_THROW(
        Objective.CalcMisfit(jif3D::rvec(Model.GetConductivities().num_elements())),
        jif3D::FatalException);
  }

BOOST_AUTO_TEST_CASE (X3D_3D_deriv_test)
  {

    jif3D::X3DModel Model;
    jif3D::MTData Data;
    MakeMTModel (Model, Data);

    const size_t nmod = Model.GetConductivities().num_elements();

    jif3D::X3DModel TrueModel(Model);
    std::fill_n(TrueModel.SetConductivities().origin(), nmod, 0.01);

    //we want to test the distortion correction as well
    boost::filesystem::path TDir = boost::filesystem::current_path();
    jif3D::X3DMTCalculator Calculator(TDir, "x3d", true);
    jif3D::rvec Observed = Calculator.Calculate(TrueModel, Data);
    std::ofstream impfile("impedance.out");
    std::copy(Observed.begin(), Observed.end(),
        std::ostream_iterator<double>(impfile, "\n"));

    jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> Objective(Calculator);
    std::vector<double> Error(
        jif3D::ConstructMTError(std::vector<double>(Observed.begin(), Observed.end()),
            0.02));
    Data.SetDataAndErrors(std::vector<double>(Observed.begin(),Observed.end()),Error);
    Objective.SetObservedData(Data);
    Objective.SetCoarseModelGeometry(Model);

    Objective.SetDataError(Error);
    const size_t nstat = Data.GetMeasPosX().size();
    jif3D::rvec ModelVec(nmod + 4 * nstat);
    // jif3D::rvec ModelVec(nmod);
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin() + nmod, ModelVec.begin());
    for (size_t i = 0; i < nstat; ++i)
      {
        ModelVec(nmod + i * 4) = 1.2;
        ModelVec(nmod + i * 4 + 1) = 0.1;
        ModelVec(nmod + i * 4 + 2) = -0.2;
        ModelVec(nmod + i * 4 + 3) = 0.8;

      }
    std::vector<double> C(ModelVec.begin() + nmod, ModelVec.end());
    std::cout << "C: ";
    std::copy(C.begin(), C.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    Data.SetDistortion(C);
    double misfit = Objective.CalcMisfit(ModelVec);
    BOOST_CHECK(misfit > 0.0);
    jif3D::rvec Gradient = Objective.CalcGradient(ModelVec);

    jif3D::X3DModel GradientModel(Model);
    std::copy(Gradient.begin(), Gradient.begin() + nmod,
        GradientModel.SetConductivities().origin());
    GradientModel.WriteNetCDF("gradmod.nc");

    std::ofstream outfile("grad3d.comp");

    for (size_t index = 0; index < ModelVec.size(); ++index)
      {
        double delta = 0.005;
        jif3D::rvec Forward(ModelVec);
        jif3D::rvec Backward(ModelVec);
        Forward(index) += delta;
        Backward(index) -= delta;
        double ForFDGrad = (Objective.CalcMisfit(Forward) - misfit) / (delta);
        double BackFDGrad = (misfit - Objective.CalcMisfit(Backward)) / delta;
        double CentFDGrad = (ForFDGrad + BackFDGrad) / 2.0;
        bool OK = Between(ForFDGrad, BackFDGrad, Gradient(index))
        || fabs((BackFDGrad - Gradient(index)) / BackFDGrad) < 0.01
        || fabs((ForFDGrad - Gradient(index)) / ForFDGrad) < 0.01
        || fabs((CentFDGrad - Gradient(index)) / CentFDGrad) < 0.03
        || fabs(Gradient(index)) < 4000.0;
        BOOST_CHECK(OK);
        if (!OK)
          {
            std::cout << "Comparison Gradient-FD ";
            std::cout << "Component: " << index << " " << ForFDGrad << " "
            << BackFDGrad << " " << CentFDGrad << " " << Gradient(index) << "\n"
            << std::endl;
          }

        outfile << index << " " << ForFDGrad << " " << BackFDGrad << " "
        << (ForFDGrad + BackFDGrad) / 2.0 << " " << Gradient(index) << std::endl;
      }
    std::cout << "\n\n\n";
    Objective.CalcMisfit(ModelVec);
    jif3D::rvec Diff = Objective.GetDataDifference();

    BOOST_CHECK(misfit > 0.0);
    jif3D::rvec Gradient2 = Objective.CalcGradient(ModelVec);
    Calculator.Calculate(Model, Data);
    jif3D::rmat Sens = Calculator.SensitivityMatrix(Model, Data, Diff);
    jif3D::rvec SensGrad = 2.0 * ublas::prec_prod(ublas::trans(Sens), Diff);
    std::ofstream diffile("diff.out");
    std::copy(Diff.begin(), Diff.end(), std::ostream_iterator<double>(diffile, "\n"));

    jif3D::rvec SensDat = ublas::prec_prod(Sens, ModelVec);
    std::ofstream outfile2("sens.comp");
    std::ofstream sensfile("sens.out");
    for (size_t i = 0; i < Sens.size1(); ++i)
      {
        ublas::matrix_row<jif3D::rmat> Row(Sens, i);
        std::copy(Row.begin(), Row.end(),
            std::ostream_iterator<double>(sensfile, " "));
        sensfile << std::endl;
      }

    for (size_t index = 0; index < ModelVec.size(); ++index)
      {
        double diff = fabs((Gradient(index) - SensGrad(index)) / Gradient(index));
        bool close = diff < 0.02;
        BOOST_CHECK(close);

        if (!close)
          {
            std::cout << "Comparison Gradient-Sensitivity ";
            std::cout << "Component: " << index << " " << Gradient(index) << " "
            << SensGrad(index) << " " << diff << std::endl;
          }
        outfile2 << index << " " << Gradient(index) << " " << Gradient2(index) << " "
        << SensGrad(index) << " " << diff << std::endl;
      }
  }

BOOST_AUTO_TEST_CASE (X3D_Titan_deriv_test)
  {

    jif3D::X3DModel Model;
    jif3D::MTData Data;
    MakeTitanModel (Model, Data);
    const size_t nmod = Model.GetConductivities().num_elements();

    jif3D::X3DModel TrueModel(Model);
    std::fill_n(TrueModel.SetConductivities().origin(), nmod, 0.01);

    //we want to test the distortion correction as well
    boost::filesystem::path TDir = boost::filesystem::current_path();
    jif3D::X3DMTCalculator Calculator(TDir, "x3d", true);
    jif3D::rvec Observed = Calculator.Calculate(TrueModel, Data);
    std::ofstream impfile("titandata.out");
    std::copy(Observed.begin(), Observed.end(),
        std::ostream_iterator<double>(impfile, "\n"));
    impfile.close();
    std::cout << "ExIndices: ";
    std::copy(Data.GetExIndices().begin(),Data.GetExIndices().end(),
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "EyIndices: ";
    std::copy(Data.GetEyIndices().begin(), Data.GetEyIndices().end(),
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    std::cout << "HIndices: ";
    std::copy(Data.GetHxIndices().begin(), Data.GetHxIndices().end(),
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;

    jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> Objective(Calculator);

    Objective.SetCoarseModelGeometry(Model);
    std::vector<double> Error(
        jif3D::ConstructMTError(std::vector<double>(Observed.begin(), Observed.end()),
            0.02));
    Data.SetDataAndErrors(std::vector<double>(Observed.begin(),Observed.end()),Error);
    Data.WriteNetCDF("titanobs.nc");
    Objective.SetDataError(Error);
    Objective.SetObservedData(Data);
    //const size_t nstat = Model.GetExIndices().size() / Model.GetFrequencies().size();
    jif3D::rvec ModelVec(nmod);
    // jif3D::rvec ModelVec(nmod);
    /*std::copy(Model.GetConductivities().origin(),
     Model.GetConductivities().origin() + nmod, ModelVec.begin());
     for (size_t i = 0; i < nstat; ++i)
     {
     ModelVec(nmod + i * 4) = 1.2;
     ModelVec(nmod + i * 4 + 1) = 0.1;
     ModelVec(nmod + i * 4 + 2) = -0.2;
     ModelVec(nmod + i * 4 + 3) = 0.8;

     }
     std::vector<double> C(ModelVec.begin() + nmod, ModelVec.end());
     std::cout << "C: ";
     std::copy(C.begin(), C.end(), std::ostream_iterator<double>(std::cout, " "));
     std::cout << std::endl;
     Model.SetDistortionParameters(C);*/
    double misfit = Objective.CalcMisfit(ModelVec);
    BOOST_CHECK(misfit > 0.0);
    jif3D::rvec Gradient = Objective.CalcGradient(ModelVec);

    jif3D::X3DModel GradientModel(Model);
    std::copy(Gradient.begin(), Gradient.begin() + nmod,
        GradientModel.SetConductivities().origin());
    GradientModel.WriteNetCDF("gradtitanmod.nc");

    std::ofstream outfile("gradtitan.comp");

    for (size_t index = 0; index < ModelVec.size(); ++index)
      {
        double delta = 0.001;
        jif3D::rvec Forward(ModelVec);
        jif3D::rvec Backward(ModelVec);
        Forward(index) += delta;
        Backward(index) -= delta;
        double ForFDGrad = (Objective.CalcMisfit(Forward) - misfit) / (delta);
        double BackFDGrad = (misfit - Objective.CalcMisfit(Backward)) / delta;
        double CentFDGrad = (ForFDGrad + BackFDGrad) / 2.0;
        bool OK = Between(ForFDGrad, BackFDGrad, Gradient(index))
        || fabs((BackFDGrad - Gradient(index)) / BackFDGrad) < 0.01
        || fabs((ForFDGrad - Gradient(index)) / ForFDGrad) < 0.01
        || fabs((CentFDGrad - Gradient(index)) / CentFDGrad) < 0.03;
        BOOST_CHECK(OK);
        if (!OK)
          {
            std::cout << "Comparison Gradient-FD ";
            std::cout << "Component: " << index << " " << ForFDGrad << " "
            << BackFDGrad << " " << CentFDGrad << " " << Gradient(index) << "\n"
            << std::endl;
          }

        outfile << index << " " << ForFDGrad << " " << BackFDGrad << " "
        << (ForFDGrad + BackFDGrad) / 2.0 << " " << Gradient(index) << std::endl;
      }
  }

BOOST_AUTO_TEST_CASE (X3D_3D_deriv_trans_test)
  {

    jif3D::X3DModel Model;
    jif3D::MTData Data;
    MakeMTModel (Model, Data);
    const size_t nmod = Model.GetConductivities().num_elements();

    jif3D::X3DModel TrueModel(Model);
    std::fill_n(TrueModel.SetConductivities().origin(), nmod, 0.01);

    //we want to test the distortion correction as well
    boost::filesystem::path TDir = boost::filesystem::current_path();
    auto MTTrans = boost::make_shared<jif3D::ComplexLogTransform>();
    jif3D::X3DMTCalculator Calculator(TDir, "x3d", true);
    jif3D::rvec Observed = Calculator.Calculate(TrueModel, Data);

    std::vector<double> Error(
        jif3D::ConstructMTError(std::vector<double>(Observed.begin(), Observed.end()),
            0.02));
    Data.SetDataAndErrors(std::vector<double>(Observed.begin(),Observed.end()),Error);
    jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> Objective(Calculator);
    Objective.SetDataTransform(MTTrans);
    Objective.SetObservedData(Data);


    Objective.SetDataError(Error);
    Objective.SetCoarseModelGeometry(Model);
    jif3D::rvec ModelVec(nmod);
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin() + nmod, ModelVec.begin());
    double misfit = Objective.CalcMisfit(ModelVec);
    BOOST_CHECK(misfit > 0.0);
    jif3D::rvec Gradient = Objective.CalcGradient(ModelVec);

    jif3D::X3DModel GradientModel(Model);
    std::copy(Gradient.begin(), Gradient.begin() + nmod,
        GradientModel.SetConductivities().origin());
    GradientModel.WriteNetCDF("gradmod.nc");

    std::ofstream outfile("grad3d_rhophi.comp");

    for (size_t index = 0; index < ModelVec.size(); ++index)
      {
        double delta = 0.01;
        jif3D::rvec Forward(ModelVec);
        jif3D::rvec Backward(ModelVec);
        Forward(index) += delta;
        Backward(index) -= delta;
        double ForFDGrad = (Objective.CalcMisfit(Forward) - misfit) / (delta);
        double BackFDGrad = (misfit - Objective.CalcMisfit(Backward)) / delta;
        double CentFDGrad = (ForFDGrad + BackFDGrad) / 2.0;
        bool OK = Between(ForFDGrad, BackFDGrad, Gradient(index))
        || fabs((BackFDGrad - Gradient(index)) / BackFDGrad) < 0.01
        || fabs((ForFDGrad - Gradient(index)) / ForFDGrad) < 0.01
        || fabs((CentFDGrad - Gradient(index)) / CentFDGrad) < 0.03;
        BOOST_CHECK(OK);
        if (!OK)
          {
            std::cout << "Comparison Gradient-FD ";
            std::cout << "Component: " << index << " " << ForFDGrad << " "
            << BackFDGrad << " " << CentFDGrad << " " << Gradient(index) << "\n"
            << std::endl;
          }

        outfile << index << " " << ForFDGrad << " " << BackFDGrad << " "
        << (ForFDGrad + BackFDGrad) / 2.0 << " " << Gradient(index) << std::endl;
      }
    std::cout << "\n\n\n";
    Objective.CalcMisfit(ModelVec);
    jif3D::rvec Diff = Objective.GetDataDifference();

  }

BOOST_AUTO_TEST_SUITE_END()
