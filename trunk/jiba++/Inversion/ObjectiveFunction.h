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
    //! The base class for all objective functions used in the inversion
    /*! This class defines the basic interface for all kinds of objective
     * functions including regularization etc. There are two main functions
     * that have to be implemented in the derived classes:
     * ImplDataDifference and ImplGradient, they return the difference between
     * observed and calculated data and the associated gradient, respectively.
     *
     * Other classes access these values through CalcMisfit and CalcGradient. Note
     * that CalcGradient is only guaranteed to return the correct result if CalcMisfit
     * has been called with the same model beforehand.
     *
     */
    class ObjectiveFunction
      {
    private:
      //How many forward plus gradient evaluations did we perform
      size_t nEval;
      jiba::rvec DataDifference;
      jiba::rvec CovarDiag;
      boost::shared_ptr<VectorTransform> DataTransform;
      //! The abstract interface for functions that implement the calculation  of the data difference
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff) = 0;
      //! The abstract interface for the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff) = 0;
      //we might have to do something in the derived class when we assign the data transform
      virtual void SetDataTransformAction()
        {

        }
    protected:
      //! Access to the data transform class for derived classes
      boost::shared_ptr<VectorTransform> GetDataTransform()
        {
          return DataTransform;
        }
    public:
      //! Access to the difference between observed and calculated data
      const jiba::rvec &GetDataDifference() const
        {
          return DataDifference;
        }
      //! Get the number of observed data
      size_t GetNData() const
        {
          return DataDifference.size();
        }
      //! Get the number of forward and gradient calculations
      size_t GetNEval() const
        {
          return nEval;
        }
      //! Assign a class that transforms the data from the forward calculation
      void SetDataTransform(boost::shared_ptr<VectorTransform> Transform)
        {
          DataTransform = Transform;
          SetDataTransformAction();
        }
      //! We assume that the data covariance is in diagonal form and store it as vector, if we do not assign a covariance it is assumed to be 1
      void SetDataCovar(const jiba::rvec &Cov)
        {
          CovarDiag = Cov;
        }
      //! Calculate the Chi-squared data misfit for the given model
      double CalcMisfit(const jiba::rvec &Model)
        {
          //calculate the data difference in the derived class
          ImplDataDifference(Model, DataDifference);
          //if we didn't assign a data covariance we assume it is 1
          const size_t ndata = DataDifference.size();
          if (CovarDiag.size() != ndata)
            {
              CovarDiag.resize(DataDifference.size());
              std::fill(CovarDiag.begin(), CovarDiag.end(), 1.0);
            }
          //go through the data and weigh each misfit by its covariance
          //add up to get the Chi-square misfit
          DataDifference = ublas::element_div(DataDifference, CovarDiag);
          ++nEval;
          return ublas::inner_prod(DataDifference, DataDifference);
        }
      //! Calculate the gradient associated with the last misfit calculation
      jiba::rvec CalcGradient(const jiba::rvec &Model)
        {
          assert(CovarDiag.size() == DataDifference.size());
          //for the gradient we need the difference between predicted and observed data
          //weighted by the squared covariance
          jiba::rvec GradDiff(ublas::element_div(DataDifference, CovarDiag));
          ++nEval;
          return ImplGradient(Model, GradDiff);
        }
      friend class JointObjective;
      ObjectiveFunction();
      virtual ~ObjectiveFunction();
      };
  /* @} */
  }

#endif /* OBJECTIVEFUNCTION_H_ */
