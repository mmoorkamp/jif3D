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
            "Use model curvature for regularization.")("tearmod", po::value<
            std::string>(),
            "Filename for a model containing information about tear zones.");

        return desc;
      }

    boost::shared_ptr<jiba::MatOpRegularization> SetupRegularization::SetupObjective(
        const po::variables_map &vm, const ThreeDModelBase &StartModel,
        boost::shared_ptr<jiba::GeneralModelTransform> Transform,
        const jiba::rvec &CovModVec)
      {
        jiba::ThreeDSeismicModel TearMod, NoTearMod;
        if (vm.count("tearmod"))
          {
            NoTearMod.SetCellSize(StartModel.GetXCellSizes()[0],
                StartModel.GetXCellSizes().size(),
                StartModel.GetYCellSizes().size(),
                StartModel.GetZCellSizes().size());
            std::fill_n(TearMod.SetSlownesses().origin(),
                TearMod.GetSlownesses().num_elements(), 1.0);
            TearMod.ReadNetCDF(vm["tearmod"].as<std::string> ());
            assert(StartModel.GetNModelElements() == TearMod. GetNModelElements());
          }

        boost::shared_ptr<jiba::MatOpRegularization> Regularization;
        if (vm.count("curvreg"))
          {
            if (vm.count("tearmod"))
              {
                Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
                    new jiba::CurvatureRegularization(StartModel, NoTearMod,
                        NoTearMod, TearMod));
              }
            else
              {
                Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
                    new jiba::CurvatureRegularization(StartModel));
              }
          }
        else
          {
            if (vm.count("tearmod"))
              {
                Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
                    new jiba::GradientRegularization(StartModel, NoTearMod,
                        NoTearMod, TearMod));
              }
            else
              {
                Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
                    new jiba::GradientRegularization(StartModel));
              }
          }

        /*if (vm.count("refmod"))
         {
         Regularization->SetReferenceModel(InvModel);
         }
         else
         {
         jiba::rvec ZeroMod(InvModel.size());
         ZeroMod.clear();
         Regularization->SetReferenceModel(ZeroMod);
         }*/
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
