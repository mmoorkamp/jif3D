//============================================================================
// Name        : JointObjective.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "JointObjective.h"
#include <iostream>
namespace jiba
  {

    JointObjective::JointObjective()
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
        for (size_t i = 0; i < nobjective; ++i)
          {
        	IndividualFits.at(i) = Objectives.at(i)->CalcMisfit(Distributor(Model, i));
            std::cout << "Individual Misfits: " << i << " "
                << IndividualFits.at(i)
                << std::endl;
            totaldata += Objectives.at(i)->GetDataDifference().size();
          }
        Diff.resize(totaldata);
        size_t currstart = 0;
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
        jiba::rvec Gradient(Model.size());
        Gradient.clear();
        const size_t nobjective = Objectives.size();
        for (size_t i = 0; i < nobjective; ++i)
          {
            std::cout << "Gradient for objective: " << i << std::endl;
            Gradient += Weights.at(i) * Distributor.TransformGradient(Model,
                Objectives.at(i)->CalcGradient(), i);
          }
        return Gradient;
      }
  }
