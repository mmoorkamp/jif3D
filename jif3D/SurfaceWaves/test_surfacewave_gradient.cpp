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

BOOST_AUTO_TEST_SUITE (SW_Gradient_Test_Suite)

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
            BOOST_CHECK_CLOSE(FDGrad(i), Gradient(i), 0.1);
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
          { 0.0, 10000.0, 20000.0, 30000.0, };
        std::vector<double> ycoords_m =
          { 500000.0, 500000.0 + 10000.0, 500000.0 + 20000.0, 500000.0 + 30000.0, 500000.0
              + 40000.0};
        std::vector<double> zcoords_m =
        //{ 0.0, 2.0, 4.3, 6.8, 9.6, 12.8, 20.0 };
            { 0.0, 5000.0, 10000.0, 15000.0};
              //{ 0.0, 5000.0 };
        const int zsize = zcoords_m.size();

        TomoModel.SetCellCoords(xcoords_m, ycoords_m, zcoords_m);
        /*const double firstdepth = TomoModel.GetZCoordinates()[0];
         const double bottomdepth = TomoModel.GetZCoordinates()[zsize - 1];
         const double topvel = 1000.0;
         const double bottomvel = 5000.0;
         const double topdens = 2700;
         const double bottomdens = 4000;*/
        std::vector<double> alpha =
         { 6000.0, 8000.0, 8000.0 };
         std::vector<double> beta =
         { 3500.0, 4500.0, 4500.0 };
         std::vector<double> rho =
         { 2700.0, 3300.0, 3300.0 };
        /*std::vector<double> alpha =
         { 650.0, 750.0, 1400.0, 1800.0, 2150.0, 2800.0 };
         std::vector<double> beta =
         { 194.0, 270.0, 367.0, 485.0, 603.0, 740.0 };
         std::vector<double> rho =
         { 1820.0, 1860.0, 1910.0, 1960.0, 2020.0, 2090.0 };*/
        /*std::vector<double> alpha =
          { 5000.0 };
        std::vector<double> beta =
          { 3000.0 };
        std::vector<double> rho =
          { 2700.0 };*/
        auto exten =
            boost::extents[xcoords_m.size() - 1][ycoords_m.size() - 1][zcoords_m.size()
                - 1];
        jif3D::SurfaceWaveModel::t3DModelData Velocity(exten), Vp(exten), Density(exten);
        for (size_t i = 0; i < TomoModel.GetData().num_elements(); ++i)
          {
            /*double Depth = TomoModel.GetZCoordinates()[i % (zsize - 1)];
             double Velocity = topvel
             + Depth * (bottomvel - topvel) / (bottomdepth - firstdepth);
             double Density = topdens
             + Depth * (bottomdens - topdens) / (bottomdepth - firstdepth);*/
            Velocity.origin()[i] = beta[i % (zsize - 1)];
            Vp.origin()[i] = alpha[i % (zsize - 1)];
            Density.origin()[i] = rho[i % (zsize - 1)];

            TomoModel.SetData() = Velocity;
            TomoModel.SetVp(Vp);
            TomoModel.SetDens(Density);
          }

        jif3D::rvec InvModel(TomoModel.GetData().num_elements());
        std::copy(TomoModel.GetData().origin(),
            TomoModel.GetData().origin() + TomoModel.GetData().num_elements(),
            InvModel.begin());
        TomoModel.WriteVTK("SWgradtest.vtk");

        std::vector<double> ycoords
          { 517000.0, 535000.0 };
        std::vector<double> zcoords
          { 0.0, 0.0 };
        std::vector<double> xcoords
          { 25000.0, 15000.0 };

        jif3D::SurfaceWaveData SWData;
        SWData.SetMeasurementPoints(xcoords, ycoords, zcoords);

        std::vector<double> eventlat =
          { 49.0, 49.0 };
        std::vector<double> eventlon =
          { 150.0, 330.0 };
        std::vector<double> epz =
          { 0.0, 0.0 };
        SWData.SetEventPositions(eventlat, eventlon, epz);
        std::vector<double> T =
          { 1, 10, 20, 30, 40, 50, 60, 70, 80 };

        SWData.SetPeriods(T);
        SWData.SetLonCentr(-123.0);
        std::vector<int> stat1 =
          { 0 };
        std::vector<int> stat2 =
          { 1 };
        SWData.SetStationPairs(stat1, stat2);
        std::vector<int> NDataPerT =
          { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        SWData.SetNDataPerT(NDataPerT);
        std::vector<int> PairIndex =
          { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        std::vector<int> EventIndex =
          { 1, 1, 1, 1, 1, 1, 1, 1, 1 };
        std::vector<int> PeriodIndex =
          { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
        std::vector<double> dtp(PairIndex.size(), 1.0);
        std::vector<double> err(dtp.size(), 1.0);
        SWData.SetIndexMap(PairIndex, EventIndex, PeriodIndex);
        SWData.SetDataAndErrors(dtp, err);

        jif3D::SurfaceWaveCalculator Calculator;
//Calculator.set_data_err(err);
        jif3D::rvec ObservedTimes(Calculator.Calculate(TomoModel, SWData));
        std::copy(ObservedTimes.begin(),ObservedTimes.end(),std::ostream_iterator<double>(std::cout," "));
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

        /*std::vector<double> dcdvs = Calculator.GetDcdvs();
         std::copy(dcdvs.begin(), dcdvs.begin() + TomoModel.GetNModelElements(),
         TomoModel.SetData().origin());
         TomoModel.WriteVTK("dcdvs.vtk");
         std::vector<double> dcdvp = Calculator.GetDcdvp();
         std::copy(dcdvp.begin(), dcdvp.begin() + TomoModel.GetNModelElements(),
         TomoModel.SetData().origin());
         TomoModel.WriteVTK("dcdvp.vtk");
         std::vector<double> dcdrho = Calculator.GetDcdrho();
         std::copy(dcdrho.begin(), dcdrho.begin() + TomoModel.GetNModelElements(),
         TomoModel.SetData().origin());
         TomoModel.WriteVTK("dcdrho.vtk");*/

      }

    BOOST_AUTO_TEST_SUITE_END()
