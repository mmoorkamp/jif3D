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
      jiba::rvec ModelCovDiag;
      boost::shared_ptr<jiba::ObjectiveFunction> Objective;
      virtual void StepImplementation(jiba::rvec &CurrentModel) = 0;
    protected:
      const boost::shared_ptr<jiba::ObjectiveFunction> GetObjective()
        {
          return Objective;
        }
    public:
      const jiba::rvec &GetModelCovDiag() const
        {
          return ModelCovDiag;
        }

      void SetModelCovDiag(const jiba::rvec &Cov)
        {
          ModelCovDiag = Cov;
        }

      void MakeStep(jiba::rvec &CurrentModel)
        {
          StepImplementation(CurrentModel);
        }
      NonLinearOptimization(boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction);
      virtual ~NonLinearOptimization();
      };
  /* @} */
  }

#endif /* NONLINEAROPTIMIZATION_H_ */
