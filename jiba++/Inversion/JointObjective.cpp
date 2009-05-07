//============================================================================
// Name        : JointObjective.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "JointObjective.h"

namespace jiba
  {

    JointObjective::JointObjective()
      {
        // TODO Auto-generated constructor stub

      }

    JointObjective::~JointObjective()
      {
        // TODO Auto-generated destructor stub
      }

    void JointObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        size_t totaldata = 0;
        const size_t nobjective = Objectives.size();
        for (size_t i = 0; i < nobjective; ++i)
          {
            Objectives.at(i)->CalcMisfit(Distributor(Model, i));
            totaldata += Objectives.at(i)->GetDataDifference().size();
          }
        Diff.resize(totaldata);
        size_t currstart = 0;
        for (size_t i = 0; i < nobjective; ++i)
          {
            const size_t ndata = Objectives.at(i)->GetDataDifference().size();
            ublas::vector_range<jiba::rvec>(Diff, ublas::range(currstart,
                currstart + ndata))
                = sqrt(Weights.at(i)) * Objectives.at(i)->GetDataDifference();
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
            Gradient += Weights.at(i) * Objectives.at(i)->CalcGradient();
          }
        return Gradient;
      }
  }
