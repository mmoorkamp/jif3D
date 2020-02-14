/*
 * test_surfacewave_gradient.cpp
 *
 *  Created on: 20 Nov 2019
 *      Author: bweise
 */

#define BOOST_TEST_MODULE SurfaceWaveModel test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "../SurfaceWaves/SurfaceWaveCalculator.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/ThreeDModelObjective.h"
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

BOOST_AUTO_TEST_SUITE( SW_Gradient_Test_Suite )

    jif3D::rvec CheckGradient(jif3D::ObjectiveFunction &Objective,
        const jif3D::rvec &Model)
      {
        Objective.CalcMisfit(Model);
        jif3D::rvec Gradient = Objective.CalcGradient(Model);

        std::ofstream gradfile("SW_vs_grad.out");
        jif3D::rvec FDGrad(3 * Model.size(), 0.0);
        //double Misfit = Objective.CalcMisfit(Model);
        for (size_t i = 0; i < Model.size(); ++i)
          {
            double delta = Model(i) * 0.001;
            jif3D::rvec Forward(Model);
            jif3D::rvec Backward(Model);
            Forward(i) += delta;
            Backward(i) -= delta;
            FDGrad(i) = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))
                / (2.0 * delta);
            //we have a problem here with small gradients
            //have to investigate more, switched off test for these cases for now
            gradfile << i << " " << FDGrad(i) << " " << Gradient(i) << std::endl;
            BOOST_CHECK_CLOSE(FDGrad(i), Gradient(i), 0.001);
          }
        return FDGrad;
      }
    BOOST_AUTO_TEST_CASE (derivative_test)
      {
        jif3D::SurfaceWaveModel TomoModel;
        //const size_t xsize = 5;
        //const size_t ysize = 6;
        //const size_t zsize = 7;
        std::vector<double> xcoords_m =
          { 0.0, 10000.0, 20000.0, 30000.0, 40000.0, 50000.0 };
        std::vector<double> ycoords_m =
          { 500000.0, 500000.0 + 10000.0, 500000.0 + 20000.0, 500000.0 + 30000.0, 500000.0
              + 40000.0, 500000.0 + 50000.0, 500000.0 + 60000.0 };
        std::vector<double> zcoords_m =
          { 0.0, 500.0, 1050.0, 1655.0, 2320.5, 3052.55, 3857.805, 4743.5855 };
        const int zsize = zcoords_m.size();

        TomoModel.SetCellCoords(xcoords_m, ycoords_m, zcoords_m);
        const double firstdepth = TomoModel.GetZCoordinates()[0];
        const double bottomdepth = TomoModel.GetZCoordinates()[zsize - 1];
        const double topvel = 1000.0;
        const double bottomvel = 5000.0;
        const double topdens = 2700;
        const double bottomdens = 4000;
        for (size_t i = 0; i < TomoModel.GetData().num_elements(); ++i)
          {
            double Depth = TomoModel.GetZCoordinates()[i % (zsize - 1)];
            double Velocity = topvel
                + Depth * (bottomvel - topvel) / (bottomdepth - firstdepth);
            double Density = topdens
                + Depth * (bottomdens - topdens) / (bottomdepth - firstdepth);

            TomoModel.SetData().origin()[i] = Velocity;
            TomoModel.SetVp().origin()[i] = Velocity * sqrt(3);
            TomoModel.SetDens().origin()[i] = Density;
          }

        jif3D::rvec InvModel(TomoModel.GetData().num_elements());
        std::copy(TomoModel.GetData().origin(),
            TomoModel.GetData().origin() + TomoModel.GetData().num_elements(),
            InvModel.begin());
        TomoModel.WriteVTK("SWgradtest.vtk");

        jif3D::SurfaceWaveData SWData;
        std::vector<double> xcoords =
          { 25000.0, 25000.0 };
        std::vector<double> ycoords =
          { 500000.0 + 35000.0, 500000.0 + 45000.0 };
        std::vector<double> zcoords =
          { 0.0, 0.0 };
        SWData.SetMeasurementPoints(xcoords, ycoords, zcoords);
        std::vector<double> dtp(1, 10.0);
        std::vector<double> err(1, 1.0);
        SWData.SetDataAndErrors(dtp, err);

        std::vector<double> epx =
          { 25000.0 };
        std::vector<double> epy =
          { 500000.0 - 100000.0 };
        std::vector<double> epz =
          { 0.0 };
        SWData.SetEventPositions(epx, epy, epz);
        std::vector<double> T =
          { 5.0 };
        SWData.SetPeriods(T);
        SWData.SetDummy(-999.9);
        SWData.SetLonCentr(-123.0);
        std::vector<double> esc =
          { 0 };
        std::vector<double> sc =
          { 0, 1 };
        SWData.SetEventStatCmb(esc);
        SWData.SetStatComb(sc);

        jif3D::SurfaceWaveCalculator Calculator;
        Calculator.set_data_err(err);
        jif3D::rvec ObservedTimes(Calculator.Calculate(TomoModel, SWData));

        jif3D::ThreeDModelObjective<jif3D::SurfaceWaveCalculator> TomoObjective(
            Calculator);
        std::copy(ObservedTimes.begin(), ObservedTimes.end(), dtp.begin());
        SWData.SetDataAndErrors(dtp, err);
        TomoObjective.SetObservedData(SWData);
        TomoObjective.SetCoarseModelGeometry(TomoModel);

        std::vector<double> TomoCovar(ObservedTimes.size());
        //we assume a general error of 0.5 s for the seismic data
        std::fill(TomoCovar.begin(), TomoCovar.end(), 0.5);
        TomoObjective.SetDataError(TomoCovar);
        //TomoObjective->SetPrecondDiag(PreCond);
        double ZeroMisfit = TomoObjective.CalcMisfit(InvModel);
        //we used the same model to calculate the observed data so the misfit should be 0
        BOOST_CHECK(ZeroMisfit == 0.0);

        //for the same model the synthetic data should equal the observed data
        jif3D::rvec SynthData = TomoObjective.GetSyntheticData();
        BOOST_CHECK(ObservedTimes.size() == SynthData.size());
        BOOST_CHECK(
            std::equal(ObservedTimes.begin(), ObservedTimes.end(), SynthData.begin()));
        //check the gradient by perturbing the travel times
        ObservedTimes *= 1.1;
        std::copy(ObservedTimes.begin(), ObservedTimes.end(), dtp.begin());
        SWData.SetDataAndErrors(dtp, err);
        TomoObjective.SetObservedData(SWData);
        double Misfit = TomoObjective.CalcMisfit(InvModel);
        BOOST_CHECK(Misfit > 0.0);

        jif3D::rvec Gradient = TomoObjective.CalcGradient(InvModel);
        std::copy(Gradient.begin(), Gradient.begin() + TomoModel.GetNModelElements(),
            TomoModel.SetData().origin());
        TomoModel.WriteVTK("SW_vs_grad.vtk");

        jif3D::rvec FDGrad = CheckGradient(TomoObjective, InvModel);
        std::copy(FDGrad.begin(), FDGrad.begin() + TomoModel.GetNModelElements(),
            TomoModel.SetData().origin());
        TomoModel.WriteVTK("SWfd.vtk");

        std::vector<double> dcdvs = Calculator.GetDcdvs();
        std::copy(dcdvs.begin(), dcdvs.begin() + TomoModel.GetNModelElements(),
            TomoModel.SetData().origin());
        TomoModel.WriteVTK("dcdvs.vtk");

      }

    BOOST_AUTO_TEST_SUITE_END()