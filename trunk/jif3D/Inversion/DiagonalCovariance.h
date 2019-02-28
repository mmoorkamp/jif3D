/*
 * DiagonalCovariance.h
 *
 *  Created on: 15 Jun 2018
 *      Author: mm489
 */

#ifndef INVERSION_DIAGONALCOVARIANCE_H_
#define INVERSION_DIAGONALCOVARIANCE_H_

#include "GeneralCovariance.h"
#include "../Global/VecMat.h"

namespace jif3D
  {

    class DiagonalCovariance: public GeneralCovariance
      {
    private:
      jif3D::rvec CovDiag;
    public:
      virtual jif3D::rvec ApplyCovar(const jif3D::rvec &vector) const override;
      virtual jif3D::rvec ApplyInvCovar(const jif3D::rvec &vector) const override;
      DiagonalCovariance(const jif3D::rvec &CD = jif3D::rvec());
      virtual ~DiagonalCovariance();
      };

  } /* namespace jif3D */

#endif /* INVERSION_DIAGONALCOVARIANCE_H_ */
