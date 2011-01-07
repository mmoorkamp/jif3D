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
      //! How many forward plus gradient evaluations did we perform
      size_t nEval;
      //! The difference between observed and synthetic data for the last forward calculation
      jiba::rvec DataDifference;
      //! The diagonal elements of the covariance matrix
      jiba::rvec CovarDiag;
      //! A possible storage for a data transformation object, whether and how this is used depends on the derived class
      boost::shared_ptr<VectorTransform> DataTransform;
      //! The abstract interface for functions that implement the calculation  of the data difference
      /*! This function has to be implemented in derived classes to calculate the data misfit from a given model vector.
       * Note that the model vector does not contain any information about model geometry,  etc. Therefore
       * derived classes have to obtain this information independently through appropriate initializations, see the
       * implementations of GravityObjective, TomoObjective or X3DObjective for examples. On exit the parameter
       * Diff has to contain the unweighted raw difference between synthetic and observed data \f$ d_{synth} - d_{obs} \f$
       * This class takes care of weighting and returning the squared sum at each call to CalcMisfit.
       * @param Model The vector of model parameters, e.g. velocity at each cell, the interpretation depends solely on the derived class
       * @param Diff The vector of unweighted differences between calculated and observed data, calculated in the derived class
       */
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff) = 0;
      //! The abstract interface for the gradient calculation
      /*! Here derived classes have to implement the calculation of the derivative of the objective function
       * with respect to the model parameters. It is important that the model vector Model and the data difference Diff correspond
       * to the same model, i.e. have been calculated from a call to ImplDataDifference. Therefore ideally we only calculate the
       * gradient after calculating the misfit for the same model.
       * @param Model The model vector for which we want to calculate the derivative of the objective function
       * @param Diff The difference between observed and calculated data that corresponds to the model
       * @return The gradient of the objective function with respect to the model parameters
       */
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff) = 0;
      //! We might have to do something in the derived class when we assign the data transform
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
      //! Some objective functions can reach a meaningful value that signals convergence, e.g. an RMS of 1 for data misfit, while regularization cannot
      /*! When we minimize data misfit we do not want to go significantly below a misfit that corresponds to the data error, i.e an RMS of 1.
       * When we reach this value the objective function has converged. For a regularization objective function we do not have such a value, the current
       * value for the regularization should not determine whether we stop the inversion. This function returns the threshold for which we assume convergence for
       * this objective and can stop the inversion. By default it is the number of data. For regularization classes it has to be overwritten.
       * @return The threshold for which we can consider this objective function to be converged
       */
      virtual double ConvergenceLimit() const
        {
          return boost::numeric_cast<double>(GetNData());
        }
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
      //! We assume that the data covariance is in diagonal form and store its square root as vector, if we do not assign a covariance it is assumed to be 1
      void SetDataCovar(const jiba::rvec &Cov)
        {
          CovarDiag = Cov;
        }
      //! Return a read-only version of the diagonal of the data covariance
      const jiba::rvec &GetDataCovar() const
        {
          return CovarDiag;
        }
      //! Calculate the Chi-squared data misfit for the given model
      /*! This function calculates the value of the objective function \f$ \Phi(m) = (d -f(m))C_D(d -f(m))\f$ for a given model vector m.
       * The details of the calculation are implemented in ImplDataDifference in the derived class. This function takes
       * care of weighting the misfit by the covariance and storing the data difference for further use in the gradient calculation.
       * @param Model The model vector. The interpretation of the values depends on the derived class.
       * @return The \f$ \chi^2 \f$ misfit of the given model
       */
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
          //we return the chi-squared misfit, not the RMS
          return ublas::inner_prod(DataDifference, DataDifference);
        }
      //! Calculate the gradient associated with the last misfit calculation
      /*! This function returns the gradient of the objective function with respect to the model
       * parameters for use in an inversion. The model vector passed to this function has to be
       * the same as for the last call to CalcMisfit. This is because we store the data difference
       * from the last call to CalcMisfit internally and use it in the gradient calculation. The
       * reason for this type of design is that we need the model vector in the inversion and therefore
       * only store it once in the inversion program and pass it by reference to CalcMisfit and CalcGradient to
       * save memory. The data difference is usually not needed and we store it internally to reduce
       * the number of parameters and the amount of variables in the main program.
       * @param Model The model vector, has to be the same as for the last call to CalcMisfit
       * @return The gradient of the objective function with respect to the model parameters
       */
      jiba::rvec CalcGradient(const jiba::rvec &Model)
        {
          assert(CovarDiag.size() == DataDifference.size());
          //for the gradient we need the difference between predicted and observed data
          //weighted by the squared covariance
          jiba::rvec GradDiff(ublas::element_div(DataDifference, CovarDiag));
          ++nEval;
          return ImplGradient(Model, GradDiff);
        }
      ObjectiveFunction();
      virtual ~ObjectiveFunction();
      };
  /* @} */
  }

#endif /* OBJECTIVEFUNCTION_H_ */
