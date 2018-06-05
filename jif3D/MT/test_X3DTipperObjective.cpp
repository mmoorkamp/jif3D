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
#include "../Titan24/ReadWriteTitanData.h"
#include "X3DModel.h"
#include "OneDMTCalculator.h"
#include "OneDMTObjective.h"
#include "X3DTipperCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"
#include "ReadWriteImpedances.h"
#include "MTTransforms.h"

BOOST_AUTO_TEST_SUITE( X3DTipperObjective_Suite )

    void MakeMTModel(jif3D::X3DModel &Model)
      {
        const size_t xsize = 4;
        const size_t ysize = 3;
        const size_t zsize = 2;
        const size_t nbglayers = 2;
        const size_t nmod = xsize * ysize * zsize;

        Model.SetMeshSize(xsize, ysize, zsize);
        Model.SetZCellSizes().resize(zsize);

        std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

        const double deltax = 100.0;
        const double deltay = 100.0;
        const double deltaz = 100.0;
        Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);

        std::fill_n(Model.SetZCellSizes().begin(), zsize, deltaz);
        Model.SetZCellSizes()[1] = deltaz * 2.0;
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
        Model.SetFrequencies().push_back(1.0);
        Model.SetFrequencies().push_back(2.0);
        Model.SetFrequencies().push_back(5.0);
        Model.SetFrequencies().push_back(10.0);
        for (size_t i = 1; i < xsize - 1; ++i)
          {
            for (size_t j = 1; j < ysize - 1; ++j)
              {
                double currx = Model.GetXCoordinates()[i] + deltax / 3.0;
                double curry = Model.GetYCoordinates()[j] + deltay / 4.0;
                double currz = (j - 1) * deltaz;
                Model.AddMeasurementPoint(currx, curry, currz);
              }
          }
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
        MakeMTModel(Model);
        const size_t xsize = Model.GetXCoordinates().size();
        const size_t ysize = Model.GetYCoordinates().size();
        const size_t zsize = Model.GetZCoordinates().size();
        const size_t nmod = xsize * ysize * zsize;

        jif3D::X3DModel TrueModel(Model);
        std::fill_n(TrueModel.SetConductivities().origin(), nmod, 0.01);

        //we want to test the distortion correction as well
        boost::filesystem::path TDir = boost::filesystem::current_path();
        jif3D::X3DTipperCalculator Calculator(TDir, "x3d", true);
        jif3D::rvec Observed = Calculator.Calculate(TrueModel);
        std::ofstream impfile("tipper.out");
        std::copy(Observed.begin(), Observed.end(),
            std::ostream_iterator<double>(impfile, "\n"));
        std::vector<double> Freq(TrueModel.GetFrequencies());

        //jif3D::WriteImpedancesToNetCDF("gra1dimp.nc", Freq, TrueModel.GetMeasPosX(),
        //    TrueModel.GetMeasPosY(), TrueModel.GetMeasPosZ(), Observed);

        jif3D::ThreeDModelObjective<jif3D::X3DTipperCalculator> Objective(Calculator);
        Objective.SetObservedData(Observed);
        Objective.SetCoarseModelGeometry(Model);
        jif3D::rvec Error(jif3D::ConstructMTError(Observed, 0.02));

        Objective.SetDataError(Error);
        const size_t nstat = Model.GetMeasPosX().size();
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