//============================================================================
// Name        : SetupRegularization.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Regularization/GradientRegularization.h"
#include "../Regularization/CurvatureRegularization.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "SetupRegularization.h"

namespace jiba
  {

    void SetTearModel(const po::variables_map &vm,
        const std::string &OptionName, const ThreeDModelBase &StartModel,
        jiba::ThreeDSeismicModel &TearModel)
      {
        if (vm.count(OptionName))
          {
            TearModel.ReadNetCDF(vm[OptionName].as<std::string> ());
          }
        else
          {
            TearModel.SetCellSize(StartModel.GetXCellSizes()[0],
                StartModel.GetXCellSizes().size(),
                StartModel.GetYCellSizes().size(),
                StartModel.GetZCellSizes().size());
            std::fill_n(TearModel.SetSlownesses().origin(),
                TearModel.GetSlownesses().num_elements(), 1.0);
          }

      }

    SetupRegularization::SetupRegularization()
      {
      }

    SetupRegularization::~SetupRegularization()
      {
      }

    po::options_description SetupRegularization::SetupOptions()
      {
        //add all possible options for regularization to an
        //options_description object that can be used for
        //parsing
        po::options_description desc("Regularization options");

        desc.add_options()("xreg", po::value(&xweight)->default_value(1.0),
            "The weight for the regularization in x-direction")("yreg",
            po::value(&yweight)->default_value(1.0),
            "The weight for the regularization in y-direction")("zreg",
            po::value(&zweight)->default_value(1.0),
            "The weight for the regularization in z-direction")("curvreg",
            "Use model curvature for regularization. If not set use gradient.")("tearmodx", po::value<
            std::string>(),
            "Filename for a model containing information about tear zones in x-direction.")(
            "tearmody", po::value<std::string>(),
            "Filename for a model containing information about tear zones in y-direction.")(
            "tearmodz", po::value<std::string>(),
            "Filename for a model containing information about tear zones in z-direction.")(
            "beta", po::value(&beta)->default_value(0.0),
            "The weight for the model parameter minimization in the regularization")(
            "substart", po::value(&substart)->default_value(false),
            "Substract the starting model when calculating the roughness");

        return desc;
      }

    boost::shared_ptr<jiba::MatOpRegularization> SetupRegularization::SetupObjective(
        const po::variables_map &vm, const ThreeDModelBase &StartModel,
        boost::shared_ptr<jiba::GeneralModelTransform> Transform,
        const jiba::rvec &CovModVec)
      {
        //setup possible tearing for the regularization for the three directions
        jiba::ThreeDSeismicModel TearModX, TearModY, TearModZ;
        SetTearModel(vm,"tearmodx",StartModel,TearModX);
        SetTearModel(vm,"tearmody",StartModel,TearModY);
        SetTearModel(vm,"tearmodz",StartModel,TearModZ);

        assert(StartModel.GetNModelElements() == TearModX.GetNModelElements());
        assert(StartModel.GetNModelElements() == TearModY.GetNModelElements());
        assert(StartModel.GetNModelElements() == TearModZ.GetNModelElements());

        //decide whether we want to use gradient base regularization
        //or curvature based regularization, the default is gradient
        boost::shared_ptr<jiba::MatOpRegularization> Regularization;
        if (vm.count("curvreg"))
          {
            Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
                new jiba::CurvatureRegularization(StartModel, TearModX,
                    TearModY, TearModZ, beta));

          }
        else
          {
            Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
                new jiba::GradientRegularization(StartModel, TearModX,
                    TearModY, TearModZ, beta));
          }
        //We either pass an empty covariance vector then the regularization class
        //takes care of setting the covariance to 1
        //or we set a covariance vector that has one value for each model cell
        if (!CovModVec.empty())
          {
            //as we treat each direction separately, the covariance
            //vector has to have a length 3 times the length of the model vector
            //here we copy the covariances to the appropriate position
            const size_t ngrid = StartModel.GetNModelElements();
            assert(CovModVec.size() == ngrid);
            jiba::rvec Cov(ngrid * 3);
            ublas::subrange(Cov, 0, ngrid) = CovModVec;
            ublas::subrange(Cov, ngrid, 2 * ngrid) = CovModVec;
            ublas::subrange(Cov, 2 * ngrid, 3 * ngrid) = CovModVec;
            Regularization->SetDataCovar(Cov);
          }

        //we can directly use the values for the weights without checking
        //the options because we set the default value to 1
        Regularization->SetXWeight(xweight);
        Regularization->SetYWeight(yweight);
        Regularization->SetZWeight(zweight);

        return Regularization;
      }
  }
