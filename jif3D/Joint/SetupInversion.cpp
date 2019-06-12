//============================================================================
// Name        : SetupInversion.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include <algorithm>
#include "SetupInversion.h"
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/NonLinearConjugateGradient.h"

namespace jif3D
  {

    SetupInversion::SetupInversion() :
        scalegrad(), corrpairs()
      {
      }

    SetupInversion::~SetupInversion()
      {
      }

    po::options_description SetupInversion::SetupOptions()
      {
        //setup the desctription object for the inversioj options
        po::options_description desc("Inversion options");
        desc.add_options()("corrpairs", po::value(&corrpairs)->default_value(5),
            "The number correction pairs for L-BFGS")("nlcg",
            "Use NLCG optimization, otherwise use L-BFGS")("scalegrad",
            po::value(&scalegrad)->default_value(true),
            "Scale the gradient for the first iteration of the L-BFGS algorithm");
        return desc;
      }

    boost::shared_ptr<jif3D::GradientBasedOptimization> SetupInversion::ConfigureInversion(
        const po::variables_map &vm,
        boost::shared_ptr<jif3D::ObjectiveFunction> ObjFunction,
        const jif3D::rvec &InvModel, boost::shared_ptr<jif3D::GeneralCovariance> CovObj)
      {
        //we can either use nlcg or L-BFGS for the optimizer
        boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer;
        if (vm.count("nlcg"))
          {
            Optimizer = boost::shared_ptr<jif3D::GradientBasedOptimization>(
                new jif3D::NonLinearConjugateGradient(ObjFunction));
          }
        else
          {

            //auto CovObj = boost::make_shared<jif3D::StochasticCovariance>(20,20,10,CovWidth,1.0,1.0);
            //auto CovObj = boost::make_shared<jif3D::DiagonalCovariance>(CovVec);
            //for L-BFGS we check whether the number of correction pairs is positive
            if (corrpairs < 0)
              throw jif3D::FatalException(
                  "Negative number of correction pairs specified !", __FILE__, __LINE__);
            Optimizer = boost::shared_ptr<jif3D::GradientBasedOptimization>(
                new jif3D::LimitedMemoryQuasiNewton(ObjFunction, CovObj, corrpairs,
                    scalegrad));
          }
        return Optimizer;
      }
  }

