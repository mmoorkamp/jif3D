
#define BOOST_TEST_MODULE DCObjective test
#define BOOST_TEST_MAIN ...
#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#endif
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "DCResistivityCalculator.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/ThreeDModelObjective.h"
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include "DCResistivityData.h"



BOOST_AUTO_TEST_SUITE( DC_Objective_Test_Suite )

    jif3D::rvec CheckGradient(jif3D::ObjectiveFunction &Objective,
        const jif3D::rvec &Model)
      {
        Objective.CalcMisfit(Model);
        jif3D::rvec Gradient = Objective.CalcGradient(Model);

        std::ofstream gradfile("dcgrad.out");
        jif3D::rvec FDGrad(Model.size(), 0.0);
        Objective.CalcMisfit(Model);
        for (size_t i = 0; i < Gradient.size(); ++i)
          {
            double delta = Model(i) * 0.01;
            jif3D::rvec Forward(Model);
            jif3D::rvec Backward(Model);
            Forward(i) += delta;
            Backward(i) -= delta;
            FDGrad(i) = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))
                / (2.0 * delta);
            gradfile << i << " " << FDGrad(i) << " " << Gradient(i) << std::endl;
            BOOST_CHECK_CLOSE(FDGrad(i), Gradient(i), 0.01);
          }
        return FDGrad;
      }

    BOOST_AUTO_TEST_CASE (derivative_test)
      {
        jif3D::ThreeDDCResistivityModel DCModel;
        jif3D::DCResistivityData DCData;
        const size_t xsize = 7;
        const size_t ysize = 8;
        const size_t zsize = 7;

        jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize,1.0);

        DCModel.SetMeshSize(xsize, ysize, zsize);
        DCModel.SetHorizontalCellSize(1, 1, xsize, ysize);
        DCModel.SetZCellSizes(ZCS);

        const size_t ngrid = xsize * ysize * zsize;
        //generate a random mesh resistivity between 1 and 3000 Ohmm
        const double rho = exp((1.0 + drand48()) * 4);
        std::fill_n(DCModel.SetResistivities().origin(), ngrid, rho);

        const double minx = 2.5;
        const double miny = 2.5;
        const double maxx = 4.5;
        const double maxy = 5.5;
        const double deltax = 2.0;
        const double deltay = 2.0;
        const double sourcez = 0.0;
        const size_t nmeasx = boost::numeric_cast<size_t>((maxx - minx) / deltax);
        const size_t nmeasy = boost::numeric_cast<size_t>((maxy - miny) / deltay);

        for (size_t i = 0; i <= nmeasx; ++i)
          {
            for (size_t j = 0; j <= nmeasy; ++j)
              {
                double sourcex = minx + i * deltax;
                double sourcey = miny + j * deltay;
                int sourceindex = lrand48() % (nmeasx * nmeasy);
                DCData.AddSource(sourcex, sourcey, sourcez, sourcex + 0.1, sourcey + 0.1,
                    sourcez);
                DCData.AddAllMeasurementPoint(sourcex + 0.2, sourcey, sourcez, sourcex + 0.2,
                    sourcey + 0.2, sourcez, sourceindex);
              }
          }

        jif3D::rvec InvModel(DCModel.GetResistivities().num_elements());
        std::copy(DCModel.GetResistivities().origin(),
        		DCModel.GetResistivities().origin() + ngrid, InvModel.begin());

        jif3D::DCResistivityCalculator Calculator;
        jif3D::rvec AppRes(Calculator.Calculate(DCModel, DCData));

        jif3D::rvec DCCovar(AppRes.size());
        std::fill(DCCovar.begin(), DCCovar.end(), 0.5);

        jif3D::ThreeDModelObjective<jif3D::DCResistivityCalculator> DCObjective(
            Calculator);


        DCData.SetDataAndErrors(std::vector<double>(AppRes.begin(),AppRes.end()),std::vector<double>(DCCovar.begin(),DCCovar.end()));



        DCObjective.SetObservedData(DCData);
        DCObjective.SetCoarseModelGeometry(DCModel);
        DCModel.WriteVTK("DCtest.vtk");

        DCObjective.SetDataError(std::vector<double>(DCCovar.begin(),DCCovar.end()));
        //TomoObjective->SetPrecondDiag(PreCond);
        double ZeroMisfit = DCObjective.CalcMisfit(InvModel);
        //we used the same model to calculate the observed data so the misfit should be 0
        BOOST_CHECK(ZeroMisfit == 0.0);

        //for the same model the synthetic data should equal the observed data
        jif3D::rvec SynthData = DCObjective.GetSyntheticData();
        BOOST_CHECK(AppRes.size() == SynthData.size());
        BOOST_CHECK(std::equal(AppRes.begin(), AppRes.end(), SynthData.begin()));
        //check the gradient by perturbing the resistivities times
        AppRes *= 1.1;
        DCData.SetDataAndErrors(std::vector<double>(AppRes.begin(),AppRes.end()),std::vector<double>(DCCovar.begin(),DCCovar.end()));

        DCObjective.SetObservedData(DCData);
        double Misfit = DCObjective.CalcMisfit(InvModel);
        BOOST_CHECK(Misfit > 0.0);

        jif3D::rvec Gradient = DCObjective.CalcGradient(InvModel);
        std::copy(Gradient.begin(), Gradient.end(), DCModel.SetResistivities().origin());
        DCModel.WriteVTK("dcgrad.vtk");

        jif3D::rvec FDGrad = CheckGradient(DCObjective, InvModel);
        std::copy(FDGrad.begin(), FDGrad.end(), DCModel.SetResistivities().origin());
        DCModel.WriteVTK("dcfd.vtk");

      }

    BOOST_AUTO_TEST_SUITE_END()
