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

    SetupRegularization::SetupRegularization()
      {
      }

    SetupRegularization::~SetupRegularization()
      {
      }

    po::options_description SetupRegularization::SetupOptions()
      {
        po::options_description desc("Regularization options");

        desc.add_options()("xreg", po::value<double>(),
            "The weight for the regularization in x-direction")("yreg",
            po::value<double>(),
            "The weight for the regularization in y-direction")("zreg",
            po::value<double>(),
            "The weight for the regularization in z-direction")("curvreg",
            "Use model curvature for regularization.")("tearmodx", po::value<
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
        jiba::ThreeDSeismicModel TearModX, TearModY, TearModZ;
        if (vm.count("tearmodx"))
          {
            TearModX.ReadNetCDF(vm["tearmodx"].as<std::string> ());
          }
        else
          {
            TearModX.SetCellSize(StartModel.GetXCellSizes()[0],
                StartModel.GetXCellSizes().size(),
                StartModel.GetYCellSizes().size(),
                StartModel.GetZCellSizes().size());
            std::fill_n(TearModX.SetSlownesses().origin(),
                TearModX.GetSlownesses().num_elements(), 1.0);
          }

        if (vm.count("tearmody"))
          {
            TearModY.ReadNetCDF(vm["tearmody"].as<std::string> ());
          }
        else
          {
            TearModY.SetCellSize(StartModel.GetXCellSizes()[0],
                StartModel.GetXCellSizes().size(),
                StartModel.GetYCellSizes().size(),
                StartModel.GetZCellSizes().size());
            std::fill_n(TearModY.SetSlownesses().origin(),
                TearModY.GetSlownesses().num_elements(), 1.0);
          }

        if (vm.count("tearmodz"))
          {

            TearModZ.ReadNetCDF(vm["tearmodz"].as<std::string> ());

          }
        else
          {
            TearModZ.SetCellSize(StartModel.GetXCellSizes()[0],
                StartModel.GetXCellSizes().size(),
                StartModel.GetYCellSizes().size(),
                StartModel.GetZCellSizes().size());
            std::fill_n(TearModZ.SetSlownesses().origin(),
                TearModZ.GetSlownesses().num_elements(), 1.0);
          }

        assert(StartModel.GetNModelElements() == TearModX.GetNModelElements());
        assert(StartModel.GetNModelElements() == TearModY.GetNModelElements());
        assert(StartModel.GetNModelElements() == TearModZ.GetNModelElements());

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

        if (!CovModVec.empty())
          {
            const size_t ngrid = StartModel.GetNModelElements();
            assert(CovModVec.size() == ngrid);
            jiba::rvec Cov(ngrid * 3);
            ublas::subrange(Cov, 0, ngrid) = CovModVec;
            ublas::subrange(Cov, ngrid, 2 * ngrid) = CovModVec;
            ublas::subrange(Cov, 2 * ngrid, 3 * ngrid) = CovModVec;
            Regularization->SetDataCovar(Cov);
          }
        if (vm.count("xreg"))
          {
            const double xreg = vm["xreg"].as<double> ();
            std::cout << "Setting xreg: " << xreg << std::endl;
            Regularization->SetXWeight(xreg);
          }
        if (vm.count("yreg"))
          {
            const double yreg = vm["yreg"].as<double> ();
            std::cout << "Setting yreg: " << yreg << std::endl;
            Regularization->SetYWeight(yreg);
          }
        if (vm.count("zreg"))
          {
            const double zreg = vm["zreg"].as<double> ();
            std::cout << "Setting zreg: " << zreg << std::endl;
            Regularization->SetZWeight(zreg);
          }

        return Regularization;
      }
  }
