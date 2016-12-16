//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Regularization test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "../Global/Jif3DPlatformHelper.h"
#include "../Gravity/test_common.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "MinDiffRegularization.h"
#include "GradientRegularization.h"
#include "CurvatureRegularization.h"
#include "MinimumSupport.h"
#include "CrossGradient.h"
#include "DotStructureConstraint.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Regularization_Test_Suite )

    void CheckGradient(jif3D::ObjectiveFunction &Objective, const jif3D::rvec &Model)
      {
        Objective.CalcMisfit(Model);
        jif3D::rvec Gradient = Objective.CalcGradient(Model);
        for (size_t i = 0; i < Gradient.size(); ++i)
          {
            double delta = Model(i) * 0.001;
            jif3D::rvec Forward(Model);
            jif3D::rvec Backward(Model);
            Forward(i) += delta;
            Backward(i) -= delta;
            double FDGrad = (Objective.CalcMisfit(Forward)
                - Objective.CalcMisfit(Backward)) / (2 * delta);
            BOOST_CHECK_CLOSE(FDGrad, Gradient(i), 0.01);
          }
      }

    BOOST_AUTO_TEST_CASE (mindiff_test)
      {
        srand((unsigned int)time(nullptr));
        jif3D::ThreeDGravityModel GravModel;
        const size_t nx = 5;
        const size_t ny = 4;
        const size_t nz = 3;
        GravModel.SetMeshSize(nx,ny,nz);


        const size_t msize = nx * ny * nz;
        jif3D::rvec StartModel(msize), PertModel(msize);
        std::generate(StartModel.begin(), StartModel.end(), rand);
        std::generate(PertModel.begin(), PertModel.end(), rand);

        jif3D::MinDiffRegularization Regularization(GravModel);
        Regularization.SetReferenceModel(StartModel);
        jif3D::rvec Diff = StartModel - PertModel;
        double Misfit = Regularization.CalcMisfit(PertModel);
        BOOST_CHECK_CLOSE(Misfit, ublas::inner_prod(Diff, Diff), 0.001);
        CheckGradient(Regularization, PertModel);
      }

    BOOST_AUTO_TEST_CASE (gradreg_test)
      {
        jif3D::ThreeDGravityModel GravModel;
        const size_t nx = 5;
        const size_t ny = 4;
        const size_t nz = 3;
        GravModel.SetMeshSize(nx,ny,nz);

        const size_t msize = GravModel.GetDensities().num_elements();
        jif3D::rvec StartModel(msize), PertModel(msize);
        jif3D::rvec ConstMod(msize);
        std::fill(ConstMod.begin(), ConstMod.end(), 1.0);
        std::generate(StartModel.begin(), StartModel.end(), rand);
        std::generate(PertModel.begin(), PertModel.end(), rand);

        jif3D::GradientRegularization Regularization(GravModel, 0.0);
        Regularization.SetReferenceModel(StartModel);
        Regularization.SetDataError(StartModel);
        Regularization.SetXWeight(5.0);
        Regularization.SetYWeight(4.0);
        Regularization.SetZWeight(3.0);
        double zero = Regularization.CalcMisfit(StartModel + ConstMod);
        BOOST_CHECK_CLOSE(zero, 0.0, 0.0001);
        Regularization.CalcMisfit(PertModel);
        CheckGradient(Regularization, PertModel);
      }

    BOOST_AUTO_TEST_CASE (minsupp_test)
      {
        srand((unsigned int)time(nullptr));
        jif3D::ThreeDGravityModel GravModel;
        const size_t nx = 5;
        const size_t ny = 4;
        const size_t nz = 3;
        GravModel.SetMeshSize(nx,ny,nz);


        const size_t msize = nx * ny * nz;
        jif3D::rvec StartModel(msize), PertModel(msize);
        std::generate(StartModel.begin(), StartModel.end(), rand);
        std::generate(PertModel.begin(), PertModel.end(), rand);

        boost::shared_ptr<jif3D::MatOpRegularization> Regularization(
            new jif3D::MinDiffRegularization(GravModel));
        Regularization->SetReferenceModel(StartModel);
        double beta = std::accumulate(StartModel.begin(), StartModel.end(), 0.0)
            / StartModel.size();
        jif3D::MinimumSupport MinSupp(Regularization, beta);

        /*double Misfit = */ MinSupp.CalcMisfit(PertModel);

        CheckGradient(MinSupp, PertModel);
      }

    BOOST_AUTO_TEST_CASE (mingradsupp_test)
      {
        srand((unsigned int)time(nullptr));
        jif3D::ThreeDGravityModel GravModel;
        const size_t nx = 5;
        const size_t ny = 4;
        const size_t nz = 3;
        GravModel.SetMeshSize(nx,ny,nz);


        const size_t msize = nx * ny * nz;
        jif3D::rvec StartModel(msize), PertModel(msize);
        std::generate(StartModel.begin(), StartModel.end(), rand);
        std::generate(PertModel.begin(), PertModel.end(), rand);

        boost::shared_ptr<jif3D::MatOpRegularization> Regularization(
            new jif3D::GradientRegularization(GravModel));
        Regularization->SetReferenceModel(StartModel);
        double beta = std::accumulate(StartModel.begin(), StartModel.end(), 0.0)
            / StartModel.size();
        jif3D::MinimumSupport MinSupp(Regularization, beta);

        /*double Misfit = */ MinSupp.CalcMisfit(PertModel);

        CheckGradient(MinSupp, PertModel);
      }

//this needs to be extended and refined
    BOOST_AUTO_TEST_CASE (curvreg_test)
      {
        jif3D::ThreeDGravityModel GravModel;
        jif3D::ThreeDSeismicModel GradModel;
        const size_t nx = 5;
        const size_t ny = 6;
        const size_t nz = 7;
        const double cellsize = 100;
        GravModel.SetMeshSize(nx,ny,nz);
        GradModel.SetCellSize(cellsize, nx, ny, nz);

        const size_t msize = GravModel.GetDensities().num_elements();
        jif3D::rvec StartModel(msize), PertModel(msize), GradModelVec(msize);
        jif3D::rvec ConstMod(msize);
        std::fill(ConstMod.begin(), ConstMod.end(), 1.0);
        std::generate(StartModel.begin(), StartModel.end(), rand);
        std::generate(PertModel.begin(), PertModel.end(), rand);

        jif3D::CurvatureRegularization Regularization(GravModel, 0.0);
        Regularization.SetReferenceModel(StartModel);
        Regularization.SetDataError(StartModel);
        Regularization.SetXWeight(5.0);
        Regularization.SetYWeight(4.0);
        Regularization.SetZWeight(3.0);
        double zero = Regularization.CalcMisfit(StartModel + ConstMod);
        BOOST_CHECK_SMALL(zero, 1e-11);

        double topslow = 1.0 / 1000.0;
        double bottomslow = 1.0 / 5000.0;

        const double firstdepth = GradModel.GetZCoordinates()[0];
        const double bottomdepth = GradModel.GetZCoordinates()[nz - 1];
        for (size_t i = 0; i < GradModel.GetSlownesses().num_elements(); ++i)
          {
            double Depth = GradModel.GetZCoordinates()[i % nz];
            double Slowness = topslow
                + (Depth - firstdepth) * (bottomslow - topslow)
                    / (bottomdepth - firstdepth);
            GradModelVec(i) = Slowness;
          }
        zero = Regularization.CalcMisfit(StartModel + GradModelVec);
        BOOST_CHECK_SMALL(zero, 2e-11);

        Regularization.CalcMisfit(PertModel);
        CheckGradient(Regularization, PertModel);
      }

//this needs to be extended and refined
    BOOST_AUTO_TEST_CASE (curvreg_tear_test)
      {
        jif3D::ThreeDGravityModel GravModel;
        jif3D::ThreeDSeismicModel GradModel;
        jif3D::ThreeDSeismicModel TearX, TearY, TearZ;
        const size_t nx = 5;
        const size_t ny = 6;
        const size_t nz = 7;
        const double cellsize = 100;
        GravModel.SetMeshSize(nx,ny,nz);
        GradModel.SetCellSize(cellsize, nx, ny, nz);

        MakeTearModel(GravModel, TearX);
        MakeTearModel(GravModel, TearY);
        MakeTearModel(GravModel, TearZ);
        jif3D::platform::srand48((unsigned int)time(nullptr));
        const double fraction = 0.1;
        for (size_t i = 0; i < TearX.GetNModelElements(); ++i)
          {
            if (jif3D::platform::drand48() > fraction)
              TearX.SetSlownesses().data()[i] = 0.0;
            if (jif3D::platform::drand48() > fraction)
              TearY.SetSlownesses().data()[i] = 0.0;
            if (jif3D::platform::drand48() > fraction)
              TearZ.SetSlownesses().data()[i] = 0.0;
          }
        const size_t msize = GravModel.GetDensities().num_elements();
        jif3D::rvec StartModel(msize), PertModel(msize), GradModelVec(msize);
        jif3D::rvec ConstMod(msize);
        std::fill(ConstMod.begin(), ConstMod.end(), 1.0);
        std::generate(StartModel.begin(), StartModel.end(), rand);
        std::generate(PertModel.begin(), PertModel.end(), rand);

        jif3D::CurvatureRegularization Regularization(GravModel, TearX, TearY, TearZ);
        Regularization.SetReferenceModel(StartModel);
        Regularization.SetDataError(StartModel);
        Regularization.SetXWeight(5.0);
        Regularization.SetYWeight(4.0);
        Regularization.SetZWeight(3.0);
        double zero = Regularization.CalcMisfit(StartModel + ConstMod);
        BOOST_CHECK_SMALL(zero, 1e-11);

        double topslow = 1.0 / 1000.0;
        double bottomslow = 1.0 / 5000.0;

        const double firstdepth = GradModel.GetZCoordinates()[0];
        const double bottomdepth = GradModel.GetZCoordinates()[nz - 1];
        for (size_t i = 0; i < GradModel.GetSlownesses().num_elements(); ++i)
          {
            double Depth = GradModel.GetZCoordinates()[i % nz];
            double Slowness = topslow
                + (Depth - firstdepth) * (bottomslow - topslow)
                    / (bottomdepth - firstdepth);
            GradModelVec(i) = Slowness;
          }
        zero = Regularization.CalcMisfit(StartModel + GradModelVec);
        BOOST_CHECK_SMALL(zero, 1e-11);

        Regularization.CalcMisfit(PertModel);
        CheckGradient(Regularization, PertModel);
      }

    BOOST_AUTO_TEST_CASE (crossgrad_test)
      {
        jif3D::ThreeDGravityModel GravModel;
        GravModel.SetDensities().resize(boost::extents[3][3][3]);
        jif3D::platform::srand48((unsigned int)time(nullptr));
        const int msize = GravModel.GetDensities().num_elements();
        jif3D::rvec PertModel(msize * 2);
        for (int i = 0; i < msize; ++i)
          {
            PertModel(i) = i + 1;
            PertModel(i + msize) = 1.0 + double(i % 2 == 0) * (i + 1);
          }

        jif3D::CrossGradient Regularization(GravModel);
        //if the two models are scaled versions of each other
        //the cross-gradient should be zero
        jif3D::rvec ZeroModel(msize * 2);
        for (int i = 0; i < msize; ++i)
          {
            ZeroModel(i) = jif3D::platform::drand48();
            ZeroModel(i + msize) = 3.2 * ZeroModel(i);
          }
        double zero = Regularization.CalcMisfit(ZeroModel);
        //practically it is very small
        BOOST_CHECK(zero < 1e-10);
        Regularization.CalcMisfit(PertModel);
        CheckGradient(Regularization, PertModel);
      }
//
//    BOOST_AUTO_TEST_CASE (dotgrad_test)
//      {
//        jif3D::ThreeDGravityModel GravModel;
//        GravModel.SetDensities().resize(boost::extents[3][3][3]);
//        srand48(time(NULL));
//        const int msize = GravModel.GetDensities().num_elements();
//        jif3D::rvec PertModel(msize * 2);
//        for (int i = 0; i < msize; ++i)
//          {
//            PertModel(i) = drand48();
//            PertModel(i + msize) = drand48();
//          }
//
//        jif3D::CrossGradient CrossReg(GravModel);
//        jif3D::DotStructureConstraint DotReg(GravModel);
//        //if the two models are scaled versions of each other
//        //the cross-gradient should be zero
//        jif3D::rvec ZeroModel(msize * 2);
//        for (int i = 0; i < msize; ++i)
//          {
//            ZeroModel(i) = drand48();
//            ZeroModel(i + msize) = 3.2 * ZeroModel(i);
//          }
//        double zero = DotReg.CalcMisfit(ZeroModel);
//        //practically it is very small
//        BOOST_CHECK(zero < 1e-10);
//
//        double cross = CrossReg.CalcMisfit(PertModel);
//        double dot = DotReg.CalcMisfit(PertModel);
//        BOOST_CHECK_CLOSE(cross, dot, 0.0001);
//        DotReg.CalcMisfit(PertModel);
//        CheckGradient(DotReg, PertModel);
//      }

    BOOST_AUTO_TEST_CASE (gradjoint_test)
      {
        srand((unsigned int)time(nullptr));
        jif3D::ThreeDGravityModel GravModel;
        GravModel.SetDensities().resize(boost::extents[5][4][3]);

        const size_t msize = GravModel.GetDensities().num_elements();
        jif3D::rvec StartModel(msize), PertModel(msize);

        std::generate(StartModel.begin(), StartModel.end(), rand);
        std::generate(PertModel.begin(), PertModel.end(), rand);
        boost::shared_ptr<jif3D::GradientRegularization> GradReg(
            new jif3D::GradientRegularization(GravModel));

        GradReg->SetReferenceModel(StartModel);

        boost::shared_ptr<jif3D::MinDiffRegularization> DiffReg(
            new jif3D::MinDiffRegularization(GravModel));
        DiffReg->SetReferenceModel(StartModel);
        jif3D::JointObjective Objective;
        Objective.AddObjective(GradReg,
            boost::shared_ptr<jif3D::ModelCopyTransform>(new jif3D::ModelCopyTransform),
            0.05);
        Objective.AddObjective(DiffReg,
            boost::shared_ptr<jif3D::ModelCopyTransform>(new jif3D::ModelCopyTransform),
            1.23);
        CheckGradient(Objective, PertModel);
      }
    BOOST_AUTO_TEST_SUITE_END()
