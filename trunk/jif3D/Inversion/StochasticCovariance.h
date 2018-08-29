/*
 * StochasticCovariance.h
 *
 *  Created on: 15 Jun 2018
 *      Author: mm489
 */

#ifndef INVERSION_STOCHASTICCOVARIANCE_H_
#define INVERSION_STOCHASTICCOVARIANCE_H_

#include "../Global/VecMat.h"
#include "GeneralCovariance.h"
#include "../ModelBase/ThreeDModelBase.h"
#include <Eigen/Dense>

namespace jif3D
  {

    class StochasticCovariance : public GeneralCovariance
      {
    public:
      typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
    private:
      size_t nx, ny, nz;
      double a;
      double nu;
      double sigma;
      MatrixXd Cm;
      bool HaveInv;
      Eigen::LDLT<MatrixXd > solver;
      jif3D::cvec Cmv, Cmf;
      void OffsetToIndex(int offset, int &xi, int &yi, int &zi) const
        {
          zi = offset % nz;
          xi = (offset - zi) / nz;
          yi = xi % ny;
          xi = (xi - yi) / ny;
        }
    public:
      const MatrixXd &GetCovar()
        {
          return Cm;
        }
      virtual jif3D::rvec ApplyCovar(const jif3D::rvec &vector) override;
      virtual jif3D::rvec ApplyInvCovar(const jif3D::rvec &vector) override;
      StochasticCovariance(size_t x, size_t y, size_t z, double ma, double mnu, double msigma);
      virtual ~StochasticCovariance();
      };

  } /* namespace jif3D */

#endif /* INVERSION_STOCHASTICCOVARIANCE_H_ */
