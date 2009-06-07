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

    class JointObjective: public ObjectiveFunction
      {
    private:
      std::vector<boost::shared_ptr<ObjectiveFunction> > Objectives;
      std::vector<double> Weights;
      std::vector<double> IndividualFits;
      ModelDistributor Distributor;
      virtual void
          ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model,
          const jiba::rvec &Diff);
    public:
      const ObjectiveFunction &GetObjective(const size_t i) const {return *Objectives.at(i);}
      const std::vector<double> &GetIndividualFits() const
        {
          return IndividualFits;
        }
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

  }

#endif /* JOINTOBJECTIVE_H_ */
