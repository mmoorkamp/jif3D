//============================================================================
// Name        : JointObjective.h
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef JOINTOBJECTIVE_H_
#define JOINTOBJECTIVE_H_

#include "ObjectiveFunction.h"
#include "ModelDistributor.h"
#include "../Gravity/ThreeDGravityModel.h"
#include <boost/shared_ptr.hpp>
namespace jiba
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
    private:
      //! This vector holds pointers to the individual objective functions
      std::vector<boost::shared_ptr<ObjectiveFunction> > Objectives;
      //! The weight for each objective
      std::vector<double> Weights;
      //! The misfit of each objective for output and analysis purposes
      std::vector<double> IndividualFits;
      std::vector<double> IndividualGradNorms;
      //! The object that translates between the optimization parameters and the individual parameters for each objective function
      ModelDistributor Distributor;
      //the implementation of the misfit calculation
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      //the implementation of the gradient calculation
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff);
    public:
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
       */
      void AddObjective(boost::shared_ptr<ObjectiveFunction> Obj,
          boost::shared_ptr<GeneralModelTransform> Transform,
          const double lambda = 1.0)
        {
          Objectives.push_back(Obj);
          Weights.push_back(lambda);
          Distributor.AddTransformer(Transform);
        }
      JointObjective();
      virtual ~JointObjective();
      };
  /* @} */
  }

#endif /* JOINTOBJECTIVE_H_ */
