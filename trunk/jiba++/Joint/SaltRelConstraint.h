//============================================================================
// Name        : SaltRelConstraint.h
// Author      : May 3, 2011
// Version     : 
// Copyright   : 2011, mmoorkamp
//============================================================================

#ifndef SALTRELCONSTRAINT_H_
#define SALTRELCONSTRAINT_H_

#include "ObjectiveFunction.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    static const double saltvel = 4500.0;
    static const double saltres = 500.0;
    static const double saltdens = 2.2;

    class SaltRelConstraint: public jiba::ObjectiveFunction
      {
    private:
      jiba::ThreeDSeismicModel Geometry;
      boost::shared_ptr<jiba::GeneralModelTransform> DensityTransform;
      boost::shared_ptr<jiba::GeneralModelTransform> ConductivityTransform;
    public:
      void SetExcludeCells(const jiba::ThreeDSeismicModel &Geo)
        {
          Geometry = Geo;
        }
      virtual void
      ImplDataDifference(const jiba::rvec &Model, jiba::rvec &Diff);
      virtual jiba::rvec ImplGradient(const jiba::rvec &Model, const jiba::rvec &Diff);
      SaltRelConstraint(boost::shared_ptr<jiba::GeneralModelTransform> DensTrans,
          boost::shared_ptr<jiba::GeneralModelTransform> CondTrans);
      virtual ~SaltRelConstraint();
      };

  }

#endif /* SALTRELCONSTRAINT_H_ */
