//============================================================================
// Name        : test_X3DObjective.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE X3DTipperObjective test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/filesystem.hpp>

#include "../Inversion/ThreeDModelObjective.h"
#include "../Global/Noise.h"
#include "X3DModel.h"
#include "OneDMTCalculator.h"
#include "OneDMTObjective.h"
#include "X3DTipperCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "ReadWriteImpedances.h"
#include "MTTransforms.h"
#include "TipperData.h"

BOOST_AUTO_TEST_SUITE (X3DTipperObjective_Suite)

void MakeTipperModel(jif3D::X3DModel &Model, jif3D::TipperData &Data)
  {
    const size_t xsize = 12;
    const size_t ysize = 9;
    const size_t zsize = 2;
    const size_t nbglayers = 2;
    const size_t nmod = xsize * ysize * zsize;

    Model.SetMeshSize(xsize, ysize, zsize);

    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

    const double deltax = 50.0;
    const double deltay = 50.0;
    const double deltaz = 50.0;
    Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize, deltaz);
    Model.SetZCellSizes(ZCS);
    std::fill_n(Model.SetConductivities().origin(), nmod, 0.01);

    Model.SetConductivities()[xsize / 2 + 1][ysize / 2][0] = 0.025;
    Model.SetConductivities()[xsize / 2][ysize / 2][0] = 0.025;
    Model.SetConductivities()[xsize / 2 - 1][ysize / 2][0] = 0.025;
    Model.SetConductivities()[xsize / 2 + 1][ysize / 2 + 1][0] = 0.025;
    Model.SetConductivities()[xsize / 2][ysize / 2 + 1][0] = 0.025;
    Model.SetConductivities()[xsize / 2 - 1][ysize / 2 + 1][0] = 0.025;
    Model.SetConductivities()[xsize / 2 + 1][ysize / 2 - 1][0] = 0.025;
    Model.SetConductivities()[xsize / 2][ysize / 2 - 1][0] = 0.025;
    Model.SetConductivities()[xsize / 2 - 1][ysize / 2 - 1][0] = 0.025;
    std::fill_n(bg_conductivities.begin(), nbglayers, 0.01);
    for (size_t i = 0; i < nbglayers; ++i)
      {
        bg_conductivities[i] *= 1.05 + i * 0.1;
      }
    std::fill_n(bg_thicknesses.begin(), nbglayers, deltaz);

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    //Model.SetFrequencies().push_back(0.1);
    //Model.SetFrequencies().push_back(2.0);
    //Model.SetFrequencies().push_back(5.0);
    Data.SetFrequencies(
      { 7.0 });
    const size_t nsites = 2;
    const double startx = 205.0;
    const double starty = 205.0;
    for (size_t i = 1; i < nsites; ++i)
      {

        Data.AddMeasurementPoint(startx + i * deltax, starty + i * deltay, 0.0);
      }
    std::vector<double> Dummy(Data.GetMeasPosX().size() * 4);
    Data.SetDataAndErrors(Dummy, Dummy);
    Data.CompleteObject();
  }

bool Between(const double limit1, const double limit2, const double value)
  {
    const double upper = std::max(limit1, limit2);
    const double lower = std::min(limit1, limit2);
    return (lower <= value) && (upper >= value);
  }

BOOST_AUTO_TEST_CASE (X3D_3D_deriv_test)
  {

    jif3D::X3DModel Model;
    jif3D::TipperData Data;
    MakeTipperModel(Model, Data);
    const size_t xsize = Model.GetData().shape()[0];
    const size_t ysize = Model.GetData().shape()[1];
    const size_t zsize = Model.GetData().shape()[2];
    const size_t nmod = xsize * ysize * zsize;

    jif3D::X3DModel TrueModel(Model);
    std::fill_n(TrueModel.SetConductivities().origin(), nmod, 0.01);
    TrueModel.SetConductivities()[xsize / 2 + 1][ysize / 2][0] = 0.0025;
    TrueModel.SetConductivities()[xsize / 2][ysize / 2][0] = 0.0025;
    TrueModel.SetConductivities()[xsize / 2 - 1][ysize / 2][0] = 0.0025;
    TrueModel.SetConductivities()[xsize / 2 + 1][ysize / 2 + 1][0] = 0.0025;
    TrueModel.SetConductivities()[xsize / 2][ysize / 2 + 1][0] = 0.0025;
    TrueModel.SetConductivities()[xsize / 2 - 1][ysize / 2 + 1][0] = 0.0025;
    TrueModel.SetConductivities()[xsize / 2 + 1][ysize / 2 - 1][0] = 0.0025;
    TrueModel.SetConductivities()[xsize / 2][ysize / 2 - 1][0] = 0.0025;
    TrueModel.SetConductivities()[xsize / 2 - 1][ysize / 2 - 1][0] = 0.0025;
    //we want to test the distortion correction as well
    boost::filesystem::path TDir = boost::filesystem::current_path();
    jif3D::X3DTipperCalculator Calculator(TDir, "x3d", true);
    jif3D::rvec Observed = Calculator.Calculate(TrueModel, Data);
    std::ofstream impfile("tipper.out");
    std::copy(Observed.begin(), Observed.end(),
        std::ostream_iterator<double>(impfile, "\n"));
    std::vector<double> Freq(Data.GetFrequencies());

    //jif3D::WriteImpedancesToNetCDF("gra1dimp.nc", Freq, TrueModel.GetMeasPosX(),
    //    TrueModel.GetMeasPosY(), TrueModel.GetMeasPosZ(), Observed);

    jif3D::X3DModel FineModel;
    FineModel.SetMeshSize(3 * xsize, 3 * ysize, 3 * zsize);
    FineModel.SetHorizontalCellSize(100.0, 100.0, 3 * xsize, 3 * ysize);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(3 * zsize, 100.0);

    FineModel.SetZCellSizes(ZCS);
    FineModel.SetBackgroundConductivities(Model.GetBackgroundConductivities());
    FineModel.SetBackgroundThicknesses(Model.GetBackgroundThicknesses());

    std::vector<double> Error(Observed.size(), 0.02);
    Data.SetDataAndErrors(std::vector<double>(Observed.begin(), Observed.end()),Error);

    jif3D::ThreeDModelObjective<jif3D::X3DTipperCalculator> Objective(Calculator);
    Objective.SetObservedData(Data);
    Objective.SetCoarseModelGeometry(Model);
//        Objective.SetFineModelGeometry(FineModel);

    Objective.SetDataError(Error);
    const size_t nstat = Data.GetMeasPosX().size();
    jif3D::rvec ModelVec(nmod);
    // jif3D::rvec ModelVec(nmod);
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin() + nmod, ModelVec.begin());

    double misfit = Objective.CalcMisfit(ModelVec);
    BOOST_CHECK(misfit > 0.0);
    jif3D::rvec Gradient = Objective.CalcGradient(ModelVec);

    jif3D::X3DModel GradientModel(Model);
    std::copy(Gradient.begin(), Gradient.begin() + nmod,
        GradientModel.SetConductivities().origin());
    GradientModel.WriteNetCDF("gradtip.nc");

    std::ofstream outfile("grad3dtip.comp");

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
  }

BOOST_AUTO_TEST_SUITE_END()
