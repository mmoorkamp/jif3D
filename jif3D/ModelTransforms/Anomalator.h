/*
 * Anomalator.h
 *
 *  Created on: Mar 23, 2021
 *      Author: max
 */

#ifndef MODELTRANSFORMS_ANOMALATOR_H_
#define MODELTRANSFORMS_ANOMALATOR_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include "GeneralModelTransform.h"

namespace jif3D
  {

    class Anomalator: public GeneralModelTransform
      {
    private:
      jif3D::rvec RefModel;
    public:
      Anomalator(const jif3D::rvec &Ref):
        RefModel(Ref)
        {

        }
      virtual ~Anomalator()
        {

        }
      virtual Anomalator* clone() const override
        {
          return new Anomalator(*this);
        }
      virtual jif3D::rvec GeneralizedToPhysical(const jif3D::rvec &FullModel) const override
        {
          return FullModel - RefModel;
        }
      virtual jif3D::rvec PhysicalToGeneralized(const jif3D::rvec &FullModel) const override
        {
          return FullModel + RefModel;
        }
      virtual jif3D::rvec Derivative(const jif3D::rvec &FullModel,
          const jif3D::rvec &Derivative) const override
        {
          return Derivative;
        }
      };

  } /* namespace jif3D */

#endif /* MODELTRANSFORMS_ANOMALATOR_H_ */
