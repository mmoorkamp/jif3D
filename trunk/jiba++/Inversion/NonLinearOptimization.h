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
    	 virtual void StepImplementation(jiba::rvec &CurrentModel) = 0;
      jiba::rvec ModelCovDiag;
      boost::shared_ptr<jiba::ObjectiveFunction> Objective;
    protected:
      jiba::rvec SearchDir;
      double Misfit;
      const boost::shared_ptr<jiba::ObjectiveFunction> GetObjective()
        {
          return Objective;
        }
    public:
      double GetMisfit() { return Misfit;}
      double GetGradNorm() {return ublas::norm_2(SearchDir);}
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
          if (ModelCovDiag.size() != CurrentModel.size())
            {
              ModelCovDiag.resize(CurrentModel.size());
              std::fill(ModelCovDiag.begin(),ModelCovDiag.end(),1.0);
            }
          if (SearchDir.size()!= CurrentModel.size() )
          {
        	  SearchDir.resize(CurrentModel.size());
          }
          StepImplementation(CurrentModel);
        }
      NonLinearOptimization(boost::shared_ptr<jiba::ObjectiveFunction> ObjFunction);
      virtual ~NonLinearOptimization();
      };
  /* @} */
  }

#endif /* NONLINEAROPTIMIZATION_H_ */
