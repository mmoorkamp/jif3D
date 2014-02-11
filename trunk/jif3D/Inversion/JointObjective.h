//============================================================================
// Name        : JointObjective.h
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef JOINTOBJECTIVE_H_
#define JOINTOBJECTIVE_H_

#include <boost/shared_ptr.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include "../Global/FatalException.h"
#include "ObjectiveFunction.h"
#include "ModelDistributor.h"

namespace jif3D
  {
    /** \addtogroup inversion General routines for inversion */
    /* @{ */
    //! A class that combines several objective function objects into a single objective function
    /*! For our joint inversion we construct the joint objective as a weighted sum of several
     * objective function to give us maximum flexibility. This class manages the calculation
     * of misfit and gradient information and provides the same interface as the individual
     * objective functions. This way the optimization routines do not have to know
     * whether they are dealing with a single objective or a joint objective problem.
     *
     * In addition, we use a ModelDistributor object to translate the model parameters
     * we use in the optimization to physical parameters that each individual objective
     * function expects. This allows us to invert for virtually arbitrary model parameterizations
     * as long as we can supply a transformation object that produces the appropriate physical
     * parameters for each method.
     */
    class JointObjective: public ObjectiveFunction
      {
      public:
    	//! For each objective function we keep track of its general purpose
    	enum ObjectiveType { unspecified, datafit, coupling, regularization };
    private:
      //! This vector holds pointers to the individual objective functions
      std::vector<boost::shared_ptr<ObjectiveFunction> > Objectives;
      //! The weight for each objective
      std::vector<double> Weights;
      //! The misfit of each objective for output and analysis purposes
      std::vector<double> IndividualFits;
      //! Stores the l2-norm of the gradient for each method
      std::vector<double> IndividualGradNorms;
      //! Store a name for each objective for display and output purposes
      std::vector<std::string> Names;
      //! The object that translates between the optimization parameters and the individual parameters for each objective function
      ModelDistributor Distributor;
      //! Do we want to print misfit and gradient information to the screen
      bool PrintMisfit;
      //! We can specify for each objective function what general purpose it has in the optimization
      std::vector<ObjectiveType> FunctionType;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ObjectiveFunction>(*this);
          ar & Objectives;
          ar & Weights;
          ar & IndividualFits;
          ar & IndividualGradNorms;
          ar & Names;
          ar & Distributor;
          ar & PrintMisfit;
        }
      //the implementation of the misfit calculation
      virtual void
      ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff);
      //the implementation of the gradient calculation
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff);
    public:
      //! The clone function provides a virtual constructor
      virtual JointObjective *clone() const
        {
          return new JointObjective(*this);
        }
      //! Read only access to the individual objective function objects, for output and maintenance purposes
      const ObjectiveFunction &GetObjective(const size_t i) const
        {
          return *Objectives.at(i);
        }
      //! Get the vector that contains the misfit for each individual objective at the current iteration
      const std::vector<double> &GetIndividualFits() const
        {
          return IndividualFits;
        }
      //! Get the RMS values for each data objective
      std::vector<double> GetRMS() const
		{
    	  std::vector<double> RMS;
    	  for (size_t i = 0; i < Weights.size(); ++i)
    	  {
    		  if (FunctionType.at(i) == datafit)
    		  {
    			  double r  = sqrt(IndividualFits.at(i) / Objectives.at(i)->GetNData());
    			  RMS.push_back(r);
    		  }
    	  }
    	  return RMS;
		}
      //! Get the vector of l2-norms of the gradient for each individual objective at the current iteration
      const std::vector<double> &GetIndividualGradNorms() const
        {
          return IndividualGradNorms;
        }
      //! Add an individual objective function
      /*! The joint objective delegates the misfit calculation to a collection of individual objective
       * function objects. Here we can add new objective functions.
       * @param Obj A shared pointer to the objective function object
       * @param Transform A shared pointer to a transformation object that calculates the correct physical parameters for this objective
       * @param lambda The weight of the objective function in the inversion
       * @param DisplayName A name for the objective function for display purposes
       */
      void AddObjective(boost::shared_ptr<ObjectiveFunction> Obj,
          boost::shared_ptr<GeneralModelTransform> Transform, const double lambda = 1.0,
          std::string DisplayName = "Objective", ObjectiveType Type = unspecified)
        {
          Objectives.push_back(Obj);
          Weights.push_back(lambda);
          Distributor.AddTransformer(Transform);
          Names.push_back(DisplayName);
          FunctionType.push_back(Type);
        }
      //! This is an experimental feature to adjust the weights for a certain type of objective function during the optimization
      /*! It generally makes sense to slowly reduce the weight for the regularization term during the inversion.
       * However, the theory for NLCG and L-BFGS does not allow this, see also the comments for SetWeights. This is
       * an experimental feature that allows to multiply the weight of objective functions
       * of a certain typeto be multiplied by a factor. USE AT YOUR OWN RISK.
       *
       * @param Which The type of objective function for which we want to adjust the weights
       * @param factor The factor by which we multiply the current value
       */
      void MultiplyWeights(ObjectiveType Which, double factor = 1.0)
      {
    	  for (size_t i = 0; i < Weights.size(); ++i)
    	  {
    		  if (FunctionType.at(i) == Which)
		      {
    			Weights.at(i) *= factor;
		       }
    	  }
      }
      //! Read only access to the weights of the objective functions
      const std::vector<double> &GetWeights() const
		{
    	  return Weights;
		}
      //! Change the weighting of the different methods
      /*! This function can be used to change the weight for
       * each method in the joint inversion after all objectives
       * have been added with AddObjective. The size of the argument
       * vector has to match the number of objective functions or an
       * exception will be thrown. Note that changing the weights during an inversion that stores previous gradient information,
       * i.e. NLCG and L-BFGS, will have unpredictable effects, as the approximation to the inverse of the Hessian will become
       * invalid.
       * @param W The vector of new weights for each method.
       */
      void SetWeights(const std::vector<double> &W)
        {
          if (W.size() != Weights.size())
            throw jif3D::FatalException(
                "Number of weights has to match the number of objective functions !");
          std::copy(W.begin(), W.end(), Weights.begin());
        }
      /*! When constructing a JointObjective object we can specify whether we want to print
       * information about the misfit of the individual objective functions to the screen or not.
       * This is part of the objective function so that we can see what happens inside the line search.
       * @param Verbose If true print misfit and gradient information
       */
      JointObjective(bool Verbose = false);
      virtual ~JointObjective();
      };

    //! Check whether we have reached the target misfit for one of the objective functions in the JointObjective object
    bool CheckConvergence(const jif3D::JointObjective &Objective);
  /* @} */
  }

#endif /* JOINTOBJECTIVE_H_ */
