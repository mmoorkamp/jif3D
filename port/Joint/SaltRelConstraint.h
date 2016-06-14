//============================================================================
// Name        : SaltRelConstraint.h
// Author      : May 3, 2011
// Version     :
// Copyright   : 2011, mmoorkamp
//============================================================================

#ifndef SALTRELCONSTRAINT_H_
#define SALTRELCONSTRAINT_H_

#include "../Global/Jif3DGlobal.h"
#include "../Inversion/ObjectiveFunction.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include <boost/shared_ptr.hpp>

namespace jif3D
  {
    static const double saltvel = 4500.0;
    static const double saltres = 500.0;
    static const double saltdens = 2.2;

    class J3DEXPORT SaltRelConstraint : public jif3D::ObjectiveFunction
      {
    private:
      jif3D::ThreeDSeismicModel Geometry;
      boost::shared_ptr<jif3D::GeneralModelTransform> DensityTransform;
      boost::shared_ptr<jif3D::GeneralModelTransform> ConductivityTransform;
    public:
      //! The clone function provides a virtual constructor
      virtual SaltRelConstraint *clone() const
        {
          return new SaltRelConstraint(*this);
        }
      void SetExcludeCells(const jif3D::ThreeDSeismicModel &Geo)
        {
          Geometry = Geo;
        }
      virtual void
      ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff);
      virtual jif3D::rvec ImplGradient(const jif3D::rvec &Model, const jif3D::rvec &Diff);
      SaltRelConstraint(boost::shared_ptr<jif3D::GeneralModelTransform> DensTrans,
          boost::shared_ptr<jif3D::GeneralModelTransform> CondTrans);
      virtual ~SaltRelConstraint();
      };

  }

#endif /* SALTRELCONSTRAINT_H_ */
