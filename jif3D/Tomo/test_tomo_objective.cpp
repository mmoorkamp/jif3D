
#define BOOST_TEST_MODULE SeismicModel test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "../Tomo/TomographyCalculator.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/ThreeDModelObjective.h"
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

BOOST_AUTO_TEST_SUITE( Tomo_Objective_Test_Suite )

jif3D  ::rvec CheckGradient(jif3D::ObjectiveFunction &Objective, const jif3D::rvec &Model)
    {
      Objective.CalcMisfit(Model);
      jif3D::rvec Gradient = Objective.CalcGradient(Model);

      std::ofstream gradfile("tomograd.out");
      jif3D::rvec FDGrad(Model.size(),0.0);
      double Misfit = Objective.CalcMisfit(Model);
      for (size_t i = 0; i < Gradient.size(); ++i)
        {
          double delta = Model(i) * 0.0001;
          jif3D::rvec Forward(Model);
          jif3D::rvec Backward(Model);
          Forward(i) += delta;
          Backward(i) -= delta;
          FDGrad(i) = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))/(2.0 * delta);
          //we have a problem here with small gradients
          //have to investigate more, switched off test for these cases for now
          gradfile << i << " " << FDGrad(i) << " " << Gradient(i) << std::endl;
          BOOST_CHECK_CLOSE(FDGrad(i),Gradient(i),0.001);
        }
      return FDGrad;
    }

  BOOST_AUTO_TEST_CASE (derivative_test)
    {
      jif3D::ThreeDSeismicModel TomoModel;
      const size_t xsize = 5;
      const size_t ysize = 6;
      const size_t zsize = 7;
      TomoModel.SetCellSize(100, xsize, ysize, zsize);
      const double firstdepth = TomoModel.GetZCoordinates()[0];
      const double bottomdepth = TomoModel.GetZCoordinates()[zsize - 1];
      const double topvel = 1000.0;
      const double bottomvel = 5000.0;
      for (size_t i = 0; i < TomoModel.GetSlownesses().num_elements(); ++i)
        {
          double Depth = TomoModel.GetZCoordinates()[i % zsize];
          double Velocity = topvel + (Depth - firstdepth)
          * (bottomvel - topvel) / (bottomdepth - firstdepth);
          TomoModel.SetSlownesses().origin()[i] = 1.0 / Velocity;
        }

      const double minx = 50;
      const double miny = 50;
      const double maxx = 450;
      const double maxy = 550;
      const double deltax = 100;
      const double deltay = 100;
      const double measz = 50.0;
      const double sourcez = 650;
      const size_t nmeasx = boost::numeric_cast<size_t>((maxx - minx) / deltax);
      const size_t nmeasy = boost::numeric_cast<size_t>((maxy - miny) / deltay);
      for (size_t i = 0; i <= nmeasx; ++i)
        {
          for (size_t j = 0; j <= nmeasy; ++j)
            {
              TomoModel.AddSource(minx + i * deltax, miny + j * deltay, sourcez);
              TomoModel.AddMeasurementPoint(minx + i * deltax, miny + j * deltay, measz);
            }
        }

      const size_t nsource = TomoModel.GetSourcePosX().size();
      const size_t nmeas = TomoModel.GetMeasPosX().size();
      for (size_t i = 0; i < nmeas; ++i)
        {
          for (size_t j = 0; j < nsource; ++j)
            {
              if (j != i)

                {
                  TomoModel.AddMeasurementConfiguration(j, i);
                }
            }
        }

      jif3D::rvec InvModel(TomoModel.GetSlownesses().num_elements());
      std::copy(TomoModel.GetSlownesses().origin(),
          TomoModel.GetSlownesses().origin()
          + TomoModel.GetSlownesses().num_elements(), InvModel.begin());

      jif3D::TomographyCalculator Calculator(true);
      jif3D::rvec ObservedTimes(Calculator.Calculate(TomoModel));

      jif3D::ThreeDModelObjective<jif3D::TomographyCalculator> TomoObjective(
          Calculator);
      TomoObjective.SetObservedData(ObservedTimes);
      TomoObjective.SetFineModelGeometry(TomoModel);
      TomoObjective.SetCoarseModelGeometry(TomoModel);
      TomoModel.WriteVTK("tomotest.vtk");
      jif3D::rvec TomoCovar(ObservedTimes.size());
      //we assume a general error of 50 ms for the seismic data
      std::fill(TomoCovar.begin(), TomoCovar.end(), 0.05);
      TomoObjective.SetDataError(TomoCovar);
      //TomoObjective->SetPrecondDiag(PreCond);
      double ZeroMisfit = TomoObjective.CalcMisfit(InvModel);
      //we used the same model to calculate the observed data so the misfit should be 0
      BOOST_CHECK(ZeroMisfit == 0.0);

      //for the same model the synthetic data should equal the observed data
      jif3D::rvec SynthData = TomoObjective.GetSyntheticData();
      BOOST_CHECK(ObservedTimes.size() == SynthData.size());
      BOOST_CHECK(std::equal(ObservedTimes.begin(),ObservedTimes.end(),SynthData.begin()));
      //check the gradient by perturbing the travel times
      ObservedTimes *= 1.1;
      TomoObjective.SetObservedData(ObservedTimes);
      double Misfit = TomoObjective.CalcMisfit(InvModel);
      BOOST_CHECK(Misfit > 0.0);

      jif3D::rvec Gradient = TomoObjective.CalcGradient(InvModel);
      std::copy(Gradient.begin(),Gradient.end(),TomoModel.SetSlownesses().origin());
      TomoModel.WriteVTK("tomograd.vtk");

      jif3D::rvec FDGrad = CheckGradient(TomoObjective, InvModel);
      std::copy(FDGrad.begin(),FDGrad.end(),TomoModel.SetSlownesses().origin());
      TomoModel.WriteVTK("tomofd.vtk");

    }

  BOOST_AUTO_TEST_SUITE_END()
