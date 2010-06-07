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
#include "../Tomo/ThreeDSeismicModel.h"
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
            "Use NLCG optimization")("covmod", po::value<std::string>(),
            "A file containing the model covariance");
        return desc;
      }

    void SetupInversion::ConfigureInversion(const po::variables_map &vm,
        boost::shared_ptr<jiba::GradientBasedOptimization> &Optimizer,
        boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction,
        const jiba::rvec &InvModel)
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
            Optimizer
                = boost::shared_ptr<jiba::GradientBasedOptimization>(
                    new jiba::LimitedMemoryQuasiNewton(ObjFunction,
                        correctionpairs));
          }
        if (vm.count("covmod"))
          {
            ThreeDSeismicModel CovModel;
            CovModel.ReadNetCDF(vm["covmod"].as<std::string> ());
            const size_t nparm = InvModel.size();
            const size_t ncovmod = CovModel.GetSlownesses().num_elements();

            if (ncovmod % nparm != 0)
              throw FatalException(
                  "Size of inversion model vector is not a multiple of covariance model size !");
            rvec CovVec(nparm);
            const size_t nsections = ncovmod / nparm;
            for (size_t i = 0; i < nsections; ++i)
              {
                for (size_t j = 0; j < ncovmod; ++j)
                  {
                    CovVec(j + i * ncovmod) = std::abs(
                        *(CovModel.GetSlownesses().origin() + j));
                  }

              }
            Optimizer->SetModelCovDiag(CovVec);
          }
      }
  }

