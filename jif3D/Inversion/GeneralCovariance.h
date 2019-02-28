/*
 * GeneralCovariance.h
 *
 *  Created on: 15 Jun 2018
 *      Author: mm489
 */

#ifndef INVERSION_GENERALCOVARIANCE_H_
#define INVERSION_GENERALCOVARIANCE_H_

#include "../Global/VecMat.h"

namespace jif3D
  {

    class GeneralCovariance
      {
    public:
      GeneralCovariance()
        {

        }
      virtual ~GeneralCovariance()
        {

        }
      virtual jif3D::rvec ApplyCovar(const jif3D::rvec &vector) const = 0;
      virtual jif3D::rvec ApplyInvCovar(const jif3D::rvec &vector) const = 0;
      };

  } /* namespace jif3D */

#endif /* INVERSION_GENERALCOVARIANCE_H_ */
