//============================================================================
// Name        : SetupInversion.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include <algorithm>
#include <boost/bind.hpp>
#include "SetupInversion.h"
#include "../Global/FatalException.h"
#include "../Global/convert.h"
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
        //setup the desctription object for the inversioj options
        po::options_description desc("Inversion options");
        desc.add_options()("corrpairs",
            po::value(&corrpairs)->default_value(5),
            "The number correction pairs for L-BFGS")("nlcg",
            "Use NLCG optimization, otherwise use L-BFGS");
        return desc;
      }

    boost::shared_ptr<jiba::GradientBasedOptimization> SetupInversion::ConfigureInversion(
        const po::variables_map &vm,
        boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction,
        const jiba::rvec &InvModel, const jiba::rvec &CovModVec)
      {
        //we can either use nlcg or L-BFGS for the optimizer
        boost::shared_ptr<jiba::GradientBasedOptimization> Optimizer;
        if (vm.count("nlcg"))
          {
            Optimizer = boost::shared_ptr<jiba::GradientBasedOptimization>(
                new jiba::NonLinearConjugateGradient(ObjFunction));
          }
        else
          {
            //for L-BFGS we check whether the number of correction pairs is positive
            if (corrpairs < 0)
              throw jiba::FatalException(
                  "Negative number of correction pairs specified !");
            Optimizer
                = boost::shared_ptr<jiba::GradientBasedOptimization>(
                    new jiba::LimitedMemoryQuasiNewton(ObjFunction,
                        corrpairs));
          }
        //if the model covariance is empty we let the optimizer object
        //take care of setting the values to 1
        //otherwise we perform some checks here
        if (!CovModVec.empty())
          {
            const size_t nparm = InvModel.size();
            const size_t ncovmod = CovModVec.size();

            if (nparm % ncovmod != 0)
              throw FatalException("Size of inversion model vector: "
                  + jiba::stringify(nparm)
                  + " is not a multiple of covariance model size: "
                  + jiba::stringify(ncovmod) + "!");
            rvec CovVec(nparm);
            const size_t nsections = nparm / ncovmod;
            for (size_t i = 0; i < nsections; ++i)
              {
                for (size_t j = 0; j < ncovmod; ++j)
                  {
                    CovVec(j + i * ncovmod) = std::abs(CovModVec(j));
                  }

              }
            Optimizer->SetModelCovDiag(CovVec);
          }
        return Optimizer;
      }
  }

