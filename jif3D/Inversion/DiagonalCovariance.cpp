/*
 * DiagonalCovariance.cpp
 *
 *  Created on: 15 Jun 2018
 *      Author: mm489
 */

#include "DiagonalCovariance.h"
#include <cassert>
namespace jif3D
  {

    jif3D::rvec DiagonalCovariance::ApplyCovar(const jif3D::rvec &vector) const
      {

        if (CovDiag.empty())
          return vector;
        assert(vector.size() == CovDiag.size());
        return ublas::element_prod(vector, CovDiag);
      }

    jif3D::rvec DiagonalCovariance::ApplyInvCovar(const jif3D::rvec &vector) const
      {
        if (CovDiag.empty())
          return vector;
        assert(vector.size() == CovDiag.size());
        return ublas::element_div(vector, CovDiag);
      }

    DiagonalCovariance::DiagonalCovariance(const jif3D::rvec &CD) :
        CovDiag(CD)
      {
        // TODO Auto-generated constructor stub

      }

    DiagonalCovariance::~DiagonalCovariance()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
