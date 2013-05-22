//============================================================================
// Name        : OneDRegularization.cpp
// Author      : 14 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "OneDRegularization.h"

namespace jif3D
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

    void OneDRegularization::ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff)
      {
        Diff.resize(Model.size());
        if (RefMod.size() != Model.size())
          {
            RefMod.resize(Model.size(), false);
            RefMod.clear();
          }
        jif3D::rvec x(Model - RefMod);
        Diff = ublas::prec_prod(OperatorMatrix, x);
      }

    jif3D::rvec OneDRegularization::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {
        const size_t nlayers = Model.size();
        jif3D::rvec Result(ublas::prec_prod(ublas::trans(OperatorMatrix), Diff));
        return 2.0 * Result;

      }
  } /* namespace jif3D */
