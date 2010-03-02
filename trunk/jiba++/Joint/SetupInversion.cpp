//============================================================================
// Name        : SetupInversion.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include "SetupInversion.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/NonLinearConjugateGradient.h"

namespace jiba
  {

    SetupInversion::SetupInversion()
      {
      }

    SetupInversion::~SetupInversion()
      {
      }

    po::options_description SetupInversion::SetupOptions()
      {
        po::options_description desc("Inversion options");
        desc.add_options()("corrpairs", po::value<int>(),
            "The number correction pairs for L-BFGS")("nlcg",
            "Use NLCG optimization");
        return desc;
      }

    void SetupInversion::ConfigureInversion(const po::variables_map &vm,boost::shared_ptr<
        jiba::GradientBasedOptimization> &Optimizer, boost::shared_ptr<
        jiba::ObjectiveFunction> ObjFunction)
      {
        int correctionpairs = 5;
        if (vm.count("corrpairs"))
          {
            correctionpairs = vm["corrpairs"].as<int> ();
          }

        if (vm.count("nlcg"))
          {
            Optimizer = boost::shared_ptr<jiba::GradientBasedOptimization>(
                new jiba::NonLinearConjugateGradient(ObjFunction));
          }
        else
          {
            Optimizer = boost::shared_ptr<jiba::GradientBasedOptimization>(
                new jiba::LimitedMemoryQuasiNewton(ObjFunction, correctionpairs));
          }
      }
  }
