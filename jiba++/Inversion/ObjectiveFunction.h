//============================================================================
// Name        : ObjectiveFunction.h
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef OBJECTIVEFUNCTION_H_
#define OBJECTIVEFUNCTION_H_

#include "../Global/VecMat.h"
#include "../Global/VectorTransform.h"
#include <boost/shared_ptr.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <cassert>

namespace jiba
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    class ObjectiveFunction
      {
    private:
      jiba::rvec DataDifference;
      jiba::rvec CovarDiag;
      boost::shared_ptr<VectorTransform> DataTransform;
      virtual void
          ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff) = 0;
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff) = 0;
      //we might have to do something in the derived class when we assign the data transform
      virtual void SetDataTransformAction()
        {

        }
    protected:
      boost::shared_ptr<VectorTransform> GetDataTransform()
        {
          return DataTransform;
        }
    public:
      void SetDataTransform(boost::shared_ptr<VectorTransform> Transform)
        {
          DataTransform = Transform;
          SetDataTransformAction();
        }
      void SetDataCovar(const jiba::rvec &Cov)
        {
          CovarDiag = Cov;
        }
      double CalcMisfit(const jiba::rvec &Model)
        {
          ImplDataDifference(Model, DataDifference);
          if (CovarDiag.size() != DataDifference.size())
            {
              CovarDiag.resize(DataDifference.size());
              std::fill(CovarDiag.begin(), CovarDiag.end(), 1.0);
            }
          DataDifference = ublas::element_div(DataDifference,CovarDiag);
          return ublas::inner_prod(DataDifference, DataDifference);
        }
      jiba::rvec CalcGradient(const jiba::rvec &Model)
        {
          assert(CovarDiag.size() == DataDifference.size());
          return ImplGradient(Model, DataDifference);
        }
      ObjectiveFunction();
      virtual ~ObjectiveFunction();
      };
  /* @} */
  }

#endif /* OBJECTIVEFUNCTION_H_ */
