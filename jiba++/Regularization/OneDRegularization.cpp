//============================================================================
// Name        : OneDRegularization.cpp
// Author      : 14 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "OneDRegularization.h"

namespace jiba
  {

    OneDRegularization::OneDRegularization(const size_t nlayers) :
        OperatorMatrix(nlayers, nlayers), RefMod(nlayers, 0.0)
      {
        const double eps = 1e-4;
        OperatorMatrix(0, 0) = 1.0 + eps;
        OperatorMatrix(0, 1) = -1.0;
        for (size_t i = 1; i < nlayers; ++i)
          {
            OperatorMatrix(i, i) = 1.0 + eps;
            OperatorMatrix(i, i - 1) = -1.0;
          }
      }

    OneDRegularization::~OneDRegularization()
      {
        // TODO Auto-generated destructor stub
      }

    void OneDRegularization::ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff)
      {
        Diff.resize(Model.size());
        if (RefMod.size() != Model.size())
          {
            RefMod.resize(Model.size(), false);
            RefMod.clear();
          }
        jiba::rvec x(Model - RefMod);
        Diff = ublas::prec_prod(OperatorMatrix, x);
      }

    jiba::rvec OneDRegularization::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        const size_t nlayers = Model.size();
        jiba::rvec Result(ublas::prec_prod(ublas::trans(OperatorMatrix), Diff));
        return 2.0 * Result;

      }
  } /* namespace jiba */
