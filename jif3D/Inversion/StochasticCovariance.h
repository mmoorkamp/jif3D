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
#include "../Eigen/Sparse"
#include "../Eigen/SparseQR"
#include "../Eigen/OrderingMethods"
//#include <Eigen/Dense>
//#include <boost/math/interpolators/cubic_b_spline.hpp>

namespace jif3D
  {

    class StochasticCovariance: public GeneralCovariance
      {
    public:
      typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

      //typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
    private:
      Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper,
          Eigen::IncompleteCholesky<double>> cg;
      //Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> cg;
      SpMat A;
      int nx, ny, nz;
      double a;
      double nu;
      double sigma;
      mutable jif3D::rvec previous_vec;
      mutable jif3D::rvec previous_result;
      int distindex;
      //MatrixXd Cm;
      bool HaveInv;
      std::vector<double> values;
      jif3D::rvec CovDiag;
      //Make sure covariance sums to 1
      bool Normalize;
      //Eigen::LDLT<MatrixXd > solver;
      jif3D::cvec Cmv, Cmf;
      //boost::math::cubic_b_spline<double> spline;
      void OffsetToIndex(int offset, int &xi, int &yi, int &zi) const
        {
          zi = offset % nz;
          xi = (offset - zi) / nz;
          yi = xi % ny;
          xi = (xi - yi) / ny;
        }
      int IndexToOffset(int xi, int yi, int zi) const
        {
          return nz * (xi * ny + yi) + zi;
        }
    public:
      //const MatrixXd &GetCovar()
      //  {
      //    return Cm;
      //  }
      virtual jif3D::rvec ApplyCovar(const jif3D::rvec &vector) const override;
      virtual jif3D::rvec ApplyInvCovar(const jif3D::rvec &vector) const override;
      StochasticCovariance(const jif3D::rvec &CD, size_t x, size_t y, size_t z, double ma,
          double mnu, double msigma, bool Norm = false);
      virtual ~StochasticCovariance();
      };

  } /* namespace jif3D */

#endif /* INVERSION_STOCHASTICCOVARIANCE_H_ */
