//============================================================================
// Name        : SetupRegularization.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Regularization/GradientRegularization.h"
#include "../Regularization/CurvatureRegularization.h"
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
            "Use model curvature for regularization.")("refmod",
            "Substract the starting model as a reference model for the regularization.");

        return desc;
      }

    boost::shared_ptr<jiba::MatOpRegularization> SetupRegularization::SetupObjective(const po::variables_map &vm,
        const ThreeDSeismicModel &StartModel, boost::shared_ptr<jiba::GeneralModelTransform> Transform)
      {

        boost::shared_ptr<jiba::MatOpRegularization> Regularization;
        if (vm.count("curvreg"))
          {
            Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
                new jiba::CurvatureRegularization(StartModel));
          }
        else
          {
            Regularization = boost::shared_ptr<jiba::MatOpRegularization>(
                new jiba::GradientRegularization(StartModel));
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
        //Regularization->SetDataCovar(ModCov);
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
