//============================================================================
// Name        : NonLinearOptimization.h
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef NONLINEAROPTIMIZATION_H_
#define NONLINEAROPTIMIZATION_H_

#include <boost/shared_ptr.hpp>
#include "ObjectiveFunction.h"
#include "../Global/VecMat.h"

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    class NonLinearOptimization
      {
    private:
      void CalcMisfitandGradient(const jiba::rvec &CurrentModel);
      virtual void StepImplementation(jiba::rvec &CurrentModel) = 0;
      jiba::rvec ModelCovDiag;
      boost::shared_ptr<jiba::ObjectiveFunction> Objective;
      jiba::rvec LastModel;
    protected:
      jiba::rvec RawGrad;
      jiba::rvec CovGrad;
      jiba::rvec SearchDir;
      double Misfit;
      const boost::shared_ptr<jiba::ObjectiveFunction> GetObjective()
        {
          return Objective;
        }
    public:
      double GetMisfit()
        {
          return Misfit;
        }
      double GetGradNorm()
        {
          return ublas::norm_2(SearchDir);
        }
      const jiba::rvec &GetModelCovDiag() const
        {
          return ModelCovDiag;
        }

      void SetModelCovDiag(const jiba::rvec &Cov)
        {
          ModelCovDiag = Cov;
        }

      void MakeStep(jiba::rvec &CurrentModel);
      NonLinearOptimization(
          boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction);
      virtual ~NonLinearOptimization();
      };
  /* @} */
  }

#endif /* NONLINEAROPTIMIZATION_H_ */
