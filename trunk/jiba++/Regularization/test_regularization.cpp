//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Inversion test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "../Gravity/test_common.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "MinDiffRegularization.h"
#include "GradientRegularization.h"
#include "CurvatureRegularization.h"
#include "CrossGradient.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Regularization_Test_Suite )

void  CheckGradient(jiba::ObjectiveFunction &Objective, const jiba::rvec &Model)
    {
      Objective.CalcMisfit(Model);
      jiba::rvec Gradient = Objective.CalcGradient(Model);
      for (size_t i = 0; i < Gradient.size(); ++i)
        {
          double delta = Model(i) * 0.00001;
          jiba::rvec Forward(Model);
          jiba::rvec Backward(Model);
          Forward(i) += delta;
          Backward(i) -= delta;
          double FDGrad = (Objective.CalcMisfit(Forward) - Objective.CalcMisfit(Backward))/(2*delta);
          BOOST_CHECK_CLOSE(FDGrad,Gradient(i),0.001);
        }
    }

  BOOST_AUTO_TEST_CASE (mindiff_test)
    {
      srand(time(NULL));
      const size_t msize = 18;
      jiba::rvec StartModel(msize), PertModel(msize);
      std::generate(StartModel.begin(),StartModel.end(),rand);
      std::generate(PertModel.begin(),PertModel.end(),rand);

      jiba::MinDiffRegularization Regularization(StartModel);
      jiba::rvec Diff = StartModel - PertModel;
      double Misfit = Regularization.CalcMisfit(PertModel);
      BOOST_CHECK_CLOSE(Misfit,ublas::inner_prod(Diff,Diff),0.001);
      CheckGradient(Regularization,PertModel);
    }

  BOOST_AUTO_TEST_CASE (gradreg_test)
    {
      jiba::ThreeDGravityModel GravModel;
      GravModel.SetDensities().resize(boost::extents[5][4][3]);

      const size_t msize = GravModel.GetDensities().num_elements();
      jiba::rvec StartModel(msize), PertModel(msize);
      jiba::rvec ConstMod(msize);
      std::fill(ConstMod.begin(),ConstMod.end(),1.0);
      std::generate(StartModel.begin(),StartModel.end(),rand);
      std::generate(PertModel.begin(),PertModel.end(),rand);

      jiba::GradientRegularization Regularization(GravModel,0.0);
      Regularization.SetReferenceModel(StartModel);
      Regularization.SetXWeight(5.0);
      Regularization.SetYWeight(4.0);
      Regularization.SetZWeight(3.0);
      double zero = Regularization.CalcMisfit(StartModel+ConstMod);
      BOOST_CHECK_CLOSE(zero,0.0,0.0001);
      Regularization.CalcMisfit(PertModel);
      CheckGradient(Regularization,PertModel);
    }

  //this needs to be extended and refined
  BOOST_AUTO_TEST_CASE (curvreg_test)
    {
      jiba::ThreeDGravityModel GravModel;
      jiba::ThreeDSeismicModel GradModel;
      const size_t nx = 5;
      const size_t ny = 6;
      const size_t nz = 7;
      const double cellsize = 100;
      GravModel.SetDensities().resize(boost::extents[nx][ny][nz]);
      GradModel.SetCellSize(cellsize, nx, ny, nz);


      const size_t msize = GravModel.GetDensities().num_elements();
      jiba::rvec StartModel(msize), PertModel(msize), GradModelVec(msize);
      jiba::rvec ConstMod(msize);
      std::fill(ConstMod.begin(),ConstMod.end(),1.0);
      std::generate(StartModel.begin(),StartModel.end(),rand);
      std::generate(PertModel.begin(),PertModel.end(),rand);

      jiba::CurvatureRegularization Regularization(GravModel,0.0);
      Regularization.SetReferenceModel(StartModel);
      Regularization.SetXWeight(5.0);
      Regularization.SetYWeight(4.0);
      Regularization.SetZWeight(3.0);
      double zero = Regularization.CalcMisfit(StartModel+ConstMod);
      BOOST_CHECK_CLOSE(zero,0.0,0.0001);

      double topslow = 1.0/1000.0;
      double bottomslow = 1.0/5000.0;

      const double firstdepth = GradModel.GetZCoordinates()[0];
      const double bottomdepth = GradModel.GetZCoordinates()[nz -1];
      for (size_t i = 0; i < GradModel.GetSlownesses().num_elements(); ++i)
        {
          double Depth = GradModel.GetZCoordinates()[i % nz];
          double Slowness = topslow + (Depth - firstdepth) * (bottomslow - topslow)/(bottomdepth - firstdepth);
          GradModelVec(i) = Slowness;
        }
      zero = Regularization.CalcMisfit(StartModel+GradModelVec);
      BOOST_CHECK_CLOSE(zero,0.0,0.0001);

      Regularization.CalcMisfit(PertModel);
      CheckGradient(Regularization,PertModel);
    }

  BOOST_AUTO_TEST_CASE (crossgrad_test)
    {
      jiba::ThreeDGravityModel GravModel;
      GravModel.SetDensities().resize(boost::extents[3][3][3]);
      srand48(time(NULL));
      const int msize = GravModel.GetDensities().num_elements();
      jiba::rvec PertModel(msize *2);
      for (int i = 0; i < msize; ++i)
        {
          PertModel(i) = i + 1;
          PertModel(i + msize) = 1.0 + double ( i %2 == 0) * (i+1);
        }

      jiba::CrossGradient Regularization(GravModel);
      //if the two models are scaled versions of each other
      //the cross-gradient should be zero
      jiba::rvec ZeroModel(msize*2);
      for (int i = 0; i < msize; ++i)
        {
          ZeroModel(i) = drand48();
          ZeroModel(i + msize) = 3.2 * ZeroModel(i);
        }
      double zero = Regularization.CalcMisfit(ZeroModel);
      //practically it is very small
      BOOST_CHECK(zero < 1e-10);
      Regularization.CalcMisfit(PertModel);
      CheckGradient(Regularization,PertModel);
    }

  BOOST_AUTO_TEST_CASE (gradjoint_test)
    {
      srand(time(NULL));
      jiba::ThreeDGravityModel GravModel;
      GravModel.SetDensities().resize(boost::extents[5][4][3]);

      const size_t msize = GravModel.GetDensities().num_elements();
      jiba::rvec StartModel(msize), PertModel(msize);

      std::generate(StartModel.begin(),StartModel.end(),rand);
      std::generate(PertModel.begin(),PertModel.end(),rand);
      boost::shared_ptr<jiba::GradientRegularization> GradReg(new jiba::GradientRegularization(GravModel));

      GradReg->SetReferenceModel(StartModel);

      boost::shared_ptr<jiba::MinDiffRegularization> DiffReg(new jiba::MinDiffRegularization(StartModel));
      jiba::JointObjective Objective;
      Objective.AddObjective(GradReg,boost::shared_ptr<jiba::ModelCopyTransform>(new jiba::ModelCopyTransform),0.05);
      Objective.AddObjective(DiffReg,boost::shared_ptr<jiba::ModelCopyTransform>(new jiba::ModelCopyTransform),1.23);
      CheckGradient(Objective,PertModel);
    }
  BOOST_AUTO_TEST_SUITE_END()
