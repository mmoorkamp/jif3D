//============================================================================
// Name        : JointObjective.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "JointObjective.h"
#include "../Global/convert.h"
#include <iostream>
namespace jiba
  {

    JointObjective::JointObjective() :
      Objectives(), Weights(), IndividualFits(), Distributor()
      {

      }

    JointObjective::~JointObjective()
      {

      }

    void JointObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        size_t totaldata = 0;
        const size_t nobjective = Objectives.size();
        IndividualFits.resize(nobjective);
        //go through all the objective function objects and calculate the misfit
        //also count how much data points we have in total
        for (size_t i = 0; i < nobjective; ++i)
          {
            IndividualFits.at(i) = Objectives.at(i)->CalcMisfit(Distributor(
                Model, i));
            totaldata += Objectives.at(i)->GetDataDifference().size();
          }
        Diff.resize(totaldata);
        size_t currstart = 0;
        //now form a single misfit vector that combines all the individual misfits
        //the total misfit will be the squared sum, so we have to weight
        //each element by the square root of the weight
        for (size_t i = 0; i < nobjective; ++i)
          {
            const size_t ndata = Objectives.at(i)->GetDataDifference().size();
            ublas::vector_range<jiba::rvec>(Diff, ublas::range(currstart,
                currstart + ndata)) = sqrt(Weights.at(i))
                * Objectives.at(i)->GetDataDifference();
            currstart += ndata;
          }
      }

    jiba::rvec JointObjective::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        //initialize the gradient vector
        jiba::rvec Gradient(Model.size());
        Gradient.clear();
        const size_t nobjective = Objectives.size();
        //add up the contribution to the gradient from each method
        //considering the weighting and the transformation
        //that has been applied to the model parameters
        for (size_t i = 0; i < nobjective; ++i)
          {
            Gradient += Weights.at(i)
                * Distributor.TransformGradient(Model, Objectives.at(i)->CalcGradient(Model), i);
          }
        return Gradient;
      }
  }
