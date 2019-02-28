//============================================================================
// Name        : ObjectiveFunction.h
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef OBJECTIVEFUNCTION_H_
#define OBJECTIVEFUNCTION_H_

#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/NumUtil.h"
#include "../Global/FatalException.h"
#include "../Global/Jif3DGlobal.h"

#include <cassert>
#include <cmath>
#include <boost/shared_ptr.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace jif3D
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
    class J3DEXPORT ObjectiveFunction
      {
    private:
      //! How many forward plus gradient evaluations did we perform
      size_t nEval;
      //! The difference between observed and synthetic data for the last forward calculation
      jif3D::rvec DataDifference;
      //! We store the misfit for each datum, as it we need it for the final evaluation of the inversion results
      jif3D::rvec IndividualMisfits;
      //! The inverse of the covariance matrix
      jif3D::comp_mat InvCovMat;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & nEval;
          ar & DataDifference;
          ar & InvCovMat;
        }
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
      ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff) = 0;
      //! The abstract interface for the gradient calculation
      /*! Here derived classes have to implement the calculation of the derivative of the objective function
       * with respect to the model parameters. It is important that the model vector Model and the data difference Diff correspond
       * to the same model, i.e. have been calculated from a call to ImplDataDifference. Therefore ideally we only calculate the
       * gradient after calculating the misfit for the same model.
       * @param Model The model vector for which we want to calculate the derivative of the objective function
       * @param Diff The difference between observed and calculated data that corresponds to the model
       * @return The gradient of the objective function with respect to the model parameters
       */
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model,
          const jif3D::rvec &Diff) = 0;
    public:
      //! We need a virtual constructor to create a new object from a pointer to a base class;
      /*! There are situations where we only have a pointer to the base class, but we need
       * a copy of the derived class without knowing what derived type it has. This virtual
       * constructor definition allows us to do this. Each derived class has to define this
       * function a return a pointer to a copy of itself.
       * @return A pointer to copy of the derived object
       */
      virtual ObjectiveFunction *clone() const = 0;
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
      const jif3D::rvec &GetDataDifference() const
        {
          return DataDifference;
        }
      //! Get the error weighted misfit for each datum
      const jif3D::rvec GetIndividualMisfit() const
        {
          return IndividualMisfits;
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
      //! We assume that the data covariance is in diagonal form and store its square root as vector, if we do not assign a covariance it is assumed to be 1
      void SetDataError(const jif3D::rvec &Cov)
        {
          //the last parameter false is important here so that ublas
          //does not try to preserve which is broken in boost 1.49
          InvCovMat.resize(Cov.size(), Cov.size(), false);
          for (size_t i = 0; i < Cov.size(); ++i)
            {
              InvCovMat(i, i) = 1.0 / jif3D::pow2(Cov(i));
            }
        }
      //! Return a read-only version of the diagonal of the data covariance
      const jif3D::rvec GetDataError() const
        {
          assert(InvCovMat.size1() == InvCovMat.size2());
          jif3D::rvec CovarDiag(InvCovMat.size1());
          for (size_t i = 0; i < InvCovMat.size1(); ++i)
            {
              CovarDiag(i) = 1.0 / sqrt(InvCovMat(i, i));
            }
          return CovarDiag;
        }
      //! Set the inverse of the data covariance matrix
      void SetInvCovMat(const jif3D::comp_mat &Mat)
        {
          if (Mat.size1() != Mat.size2())
            {
              throw jif3D::FatalException("Covariance matrix needs to be square ! ", __FILE__, __LINE__);
            }
          InvCovMat = Mat;
        }
      //! Access the inverse of the data covariance matrix
      const jif3D::comp_mat &GetInvCovMat() const
        {
          return InvCovMat;
        }
      //! Calculate the Chi-squared data misfit for the given model
      /*! This function calculates the value of the objective function \f$ \Phi(m) = (d -f(m))C_D(d -f(m))\f$ for a given model vector m.
       * The details of the calculation are implemented in ImplDataDifference in the derived class. This function takes
       * care of weighting the misfit by the covariance and storing the data difference for further use in the gradient calculation.
       * @param Model The model vector. The interpretation of the values depends on the derived class.
       * @return The \f$ \chi^2 \f$ misfit of the given model
       */
      double CalcMisfit(const jif3D::rvec &Model)
        {
          //calculate the data difference in the derived class
          ImplDataDifference(Model, DataDifference);
          //if we didn't assign a data covariance we assume it is 1
          const size_t ndata = DataDifference.size();
          if (InvCovMat.size1() != ndata || InvCovMat.size2() != ndata)
            {
              //the last parameter false is important here so that ublas
              //does not try to preserve which is broken in boost 1.49
              InvCovMat.resize(DataDifference.size(), DataDifference.size(), false);
              for (size_t i = 0; i < ndata; ++i)
                {
                  InvCovMat(i, i) = 1.0;
                }
            }
          //go through the data and weigh each misfit by its covariance
          //add up to get the Chi-square misfit
          jif3D::rvec TmpDiff(ndata);
          //it is essential here to store  the result of the vector matrix
          //multiplication in a temporary vector instead of reusing
          //DataDifference and saving some memory, the penalty
          //in terms of runtime is enormous
          ublas::axpy_prod(InvCovMat, DataDifference, TmpDiff, true);
          ++nEval;
          //we store the individual misfits for use at the end of the inversion
          //otherwise it is difficult to restore it from the covariance matrix and the
          //weighted datadifference
          IndividualMisfits = ublas::element_prod(DataDifference, TmpDiff);
          //we return the chi-squared misfit, not the RMS
          double Chi = ublas::sum(IndividualMisfits);
          for (size_t i = 0; i < ndata; ++i)
            {
              IndividualMisfits(i) = sqrt(IndividualMisfits(i))
                  * jif3D::sign(DataDifference(i));
            }
          //for the Gradient calculation we need the data difference
          //weighted by the covariance matrix, so we assign
          //the weighted values in TmpDiff to the object property DataDifference
          // for reuse in CalcGradient below.
          DataDifference = TmpDiff;
          return Chi;
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
      jif3D::rvec CalcGradient(const jif3D::rvec &Model)
        {
          assert(InvCovMat.size1() == DataDifference.size());
          assert(InvCovMat.size2() == DataDifference.size());
          //for the gradient we need the difference between predicted and observed data
          //weighted by the squared covariance
          //jif3D::rvec GradDiff(ublas::element_div(DataDifference, CovarDiag));
          ++nEval;
          return ImplGradient(Model, DataDifference);
        }
      ObjectiveFunction();
      virtual ~ObjectiveFunction();
      };
  /* @} */

  }
//BOOST_CLASS_EXPORT_KEY(jif3D::ObjectiveFunction)
#endif /* OBJECTIVEFUNCTION_H_ */
