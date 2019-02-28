//============================================================================
// Name        : test_optstep.cpp
// Author      : Apr 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE StochCov test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include "../Global/VecMat.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include <fstream>
#include "StochasticCovariance.h"
#include "DiagonalCovariance.h"
#include "../MT/X3DModel.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
//#include <unsupported/Eigen/IterativeSolvers>

class MatrixReplacement;
using Eigen::SparseMatrix;
namespace Eigen
  {
    namespace internal
      {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template<>
        struct traits<MatrixReplacement> : public Eigen::internal::traits<
            Eigen::SparseMatrix<double> >
          {
          };
      }
  }
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement: public Eigen::EigenBase<MatrixReplacement>
  {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum
    {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
    };
  Index rows() const
    {
      return nelem;
    }
  Index cols() const
    {
      return nelem;
    }
  template<typename Rhs>
  Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(
      const Eigen::MatrixBase<Rhs>& x) const
    {
      return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this,
          x.derived());
    }
  // Custom API:
  MatrixReplacement(size_t x, size_t y, size_t z, double ma, double mnu, double msigma) :
      Cov(x, y, z, ma, mnu, msigma)
    {
      nelem = x * y * z;
    }
  jif3D::StochasticCovariance Cov;

private:
  size_t nelem;
  };

// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen
  {
    namespace internal
      {
        template<typename Rhs>
        struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape,
            GemvProduct> // GEMV stands for matrix-vector
        : generic_product_impl_base<MatrixReplacement, Rhs,
            generic_product_impl<MatrixReplacement, Rhs> >
          {
          typedef typename Product<MatrixReplacement, Rhs>::Scalar Scalar;
          template<typename Dest>
          static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs,
              const Rhs& rhs, const Scalar& alpha)
            {
              // This method should implement "dst += alpha * lhs * rhs" inplace,
              // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
              assert(alpha == Scalar(1) && "scaling is not implemented"); EIGEN_ONLY_USED_FOR_DEBUG(alpha);
              jif3D::rvec in(rhs.size());
              std::copy(rhs.data(), rhs.data() + rhs.size(), in.begin());
              jif3D::rvec out = lhs.Cov.ApplyCovar(in);
              std::copy(out.begin(), out.end(), dst.data());

            }
          };
      }
  }

BOOST_AUTO_TEST_SUITE( StochCov_Test_Suite )

    BOOST_AUTO_TEST_CASE (StochCov_inv_test)
      {
        jif3D::X3DModel Model;
        const size_t nx = 20;
        const size_t ny = 20;
        const size_t nz = 10;
        Model.SetMeshSize(nx, ny, nz);
        Model.SetHorizontalCellSize(1.0, 1.0, nx, ny);
        jif3D::ThreeDModelBase::t3DModelDim ZCS(nz, 1.0);
        Model.SetZCellSizes(ZCS);

        double a = 1.0;
        double nu = 1.0;
        double sigma = 1.0;

        MatrixReplacement A(nx, ny, nz, a, nu, sigma);
        jif3D::rvec m(nx * ny * nz, 0.0);
        //std::generate(m.begin(), m.end(), drand48);
        //std::iota( m.begin(), m.end(),1);
        m(nx * ny * nz / 2 + nz / 2) = 1.0;
        Eigen::VectorXd b(nx * ny * nz), x;
        std::copy(m.begin(), m.end(), b.data());
        //b = Eigen::Map<Eigen::VectorXd>(m.data(),m.size());

        //Eigen::BiCGSTAB<MatrixReplacement, Eigen::IdentityPreconditioner> bicg;
        //bicg.compute(A);
        //x = bicg.solve(b);

        Eigen::VectorXd y = A * b;

        std::ofstream covfile("test_cov.out");
        std::copy(y.data(), y.data() + y.size(),
            std::ostream_iterator<double>(covfile, "\n"));

        Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper,
            Eigen::IdentityPreconditioner> cg;
        cg.compute(A);
        x = cg.solve(b);
        std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: "
            << cg.error() << std::endl;

        std::ofstream covinvfile("test_cov_inv.out");
        std::copy(x.data(), x.data() + x.size(),
            std::ostream_iterator<double>(covinvfile, "\n"));
      }

    BOOST_AUTO_TEST_CASE (StochCov_test)
      {
        jif3D::X3DModel Model;
        const size_t nx = 20;
        const size_t ny = 20;
        const size_t nz = 10;
        Model.SetMeshSize(nx, ny, nz);
        Model.SetHorizontalCellSize(1.0, 1.0, nx, ny);
        jif3D::ThreeDModelBase::t3DModelDim ZCS(nz, 1.0);
        Model.SetZCellSizes(ZCS);

        double a = 1.0;
        double nu = 1.0;
        double sigma = 1.0;

        jif3D::StochasticCovariance Cov(nx, ny, nz, a, nu, sigma);
        jif3D::rvec m(nx * ny * nz, 0.0);
        //std::generate(m.begin(), m.end(), drand48);
        //std::iota( m.begin(), m.end(),1);
        m(nx * ny * nz / 2 + nz / 2) = 1.0;
        jif3D::rvec mCm = Cov.ApplyCovar(m);
        //std::ofstream covfile("test_cov.out");
        //std::copy(mCm.begin(), mCm.end(),
        //		std::ostream_iterator<double>(covfile, "\n"));

        jif3D::rvec mCmi = Cov.ApplyInvCovar(mCm);

        for (size_t i = 0; i < nx * ny * nz; ++i)
          if (std::abs(m(i)) > 1e-6)
            {
              BOOST_CHECK_CLOSE(m(i), mCmi(i), 0.1);
            }

        const size_t ntries = 0;
        for (size_t i = 0; i < ntries; ++i)
          {
            std::generate(m.begin(), m.end(), drand48);
            jif3D::rvec mCm = Cov.ApplyCovar(m);
            jif3D::rvec mCmi = Cov.ApplyInvCovar(mCm);
          }
      }

    BOOST_AUTO_TEST_CASE (DiagCov_test)
      {
        const size_t nmod = 751;
        jif3D::rvec Covvec(nmod), m(nmod);
        for (double &c : Covvec)
          {
            c = std::abs(drand48()) + 0.01;
          }
        jif3D::DiagonalCovariance Cov(Covvec);

        std::generate(m.begin(), m.end(), drand48);
        jif3D::rvec cm = Cov.ApplyCovar(m);
        jif3D::rvec cmi = Cov.ApplyInvCovar(cm);
        for (size_t i = 0; i < nmod; ++i)
          {
            BOOST_CHECK_CLOSE(m(i), cmi(i), 0.1);
          }
      }

    BOOST_AUTO_TEST_SUITE_END()
