/*
 * StochasticCovariance.cpp
 *
 *  Created on: 15 Jun 2018
 *      Author: mm489
 */

#include "StochasticCovariance.h"
#include "../Global/NumUtil.h"
#include "../Global/kissfft.hh"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
#include <fstream>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
//#include <unsupported/Eigen/IterativeSolvers>

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

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
              assert(alpha == Scalar(1) && "scaling is not implemented");EIGEN_ONLY_USED_FOR_DEBUG(alpha);
              jif3D::rvec in(rhs.size());
              std::copy(rhs.data(), rhs.data() + rhs.size(), in.begin());
              jif3D::rvec out = lhs.Cov.ApplyCovar(in);
              std::copy(out.begin(), out.end(), dst.data());

            }
          };
      }
  }

namespace jif3D
  {

    StochasticCovariance::StochasticCovariance(size_t x, size_t y, size_t z, double ma,
        double mnu, double msigma) :
        A(x * y * z, x * y * z), nx(x), ny(y), nz(z), a(ma), nu(mnu), sigma(msigma), HaveInv(
            false)
      {
        distindex = 0;
        double factor = sigma * sigma / (std::pow(2, nu - 1) * boost::math::tgamma(nu));

        const size_t nmod = nx * ny * nz;
        //Cm.resize(nmod, nmod);
        //Cmf.resize(nmod * 2);
        //Cmv.resize(nmod * 2);
        double min = 0.0;
        double max = std::sqrt(jif3D::pow2(nx) + jif3D::pow2(ny) + jif3D::pow2(nz));
        double steps = 1000;
        double delta = (max - min) / steps;
        double aa = std::abs(a);
        for (double currval = 0; currval <= max; currval += delta)
          {
            double res =
                (currval == 0.0) ?
                    1.0 :
                    factor * std::pow(currval / aa, nu)
                        * boost::math::cyl_bessel_k(nu, currval / aa);
            values.push_back(res);
          }
        const double threshold = 0.001;
        double currval = 1.0;
        while (currval > threshold && distindex < values.size())
          {
            currval = values.at(distindex);
            distindex++;
          }
        distindex = std::ceil(distindex * delta);
        std::cout << "Index for correlation distance is " << distindex << std::endl;
        boost::posix_time::ptime starttime =
            boost::posix_time::microsec_clock::local_time();
        std::vector<T> coefficients; // list of non-zeros coefficients

//#pragma omp parallel for default(shared)
        for (size_t i = 0; i < nmod; ++i)
          {
            int xi, yi, zi;
            OffsetToIndex(i, xi, yi, zi);
            //double currx = Model.GetXCoordinates()[xi] + Model.GetXCellSizes()[xi] / 2.0;
            //double curry = Model.GetYCoordinates()[yi] + Model.GetYCellSizes()[yi] / 2.0;
            //double currz = Model.GetZCoordinates()[zi] + Model.GetZCellSizes()[zi] / 2.0;
            int startx = std::max(0, xi - distindex);
            int endx = std::min(nx, xi + distindex);
            int starty = std::max(0, yi - distindex);
            int endy = std::min(ny, yi + distindex);
            int startz = std::max(0, zi - distindex);
            int endz = std::min(nz, zi + distindex);
            for (size_t j = startx; j < endx; ++j)
              {
                for (size_t k = starty; k < endy; ++k)
                  {
                    for (size_t l = startz; l < endz; ++l)
                      {
                        int offset = IndexToOffset(j, k, l);
                        double r = std::sqrt(
                            jif3D::pow2(xi - j) + jif3D::pow2(yi - k)
                                + jif3D::pow2(zi - l));
                        int index = std::round((r - min) / delta);
                        double inter = values[index];
                        //+ (values[index + 1] - values[index]) / delta
                        //		* (r - index * delta);
                        if (a < 0.0)
                          {
                            coefficients.push_back(T(i, offset, inter));
                          }
                        //previous_result(i) += inter * vector(offset);
                      }
                  }
              }

          }
        if (a < 0.0)
          {
            std::cout << "Preparing inverse covariance " << std::endl;
            A.setFromTriplets(coefficients.begin(), coefficients.end());
            A.makeCompressed();
            cg.compute(A);

            boost::posix_time::ptime endtime =
                boost::posix_time::microsec_clock::local_time();
            double time = (endtime - starttime).total_seconds();

            std::cout << " took " << time << " seconds " << std::endl;
          }
        //previous_vec = vector;

//        boost::math::cubic_b_spline<double> localspline(values.begin(), values.end(), min, delta);
        //      spline = localspline;
        /*#pragma omp parallel for default(shared)
         for (size_t i = 0; i < nmod; ++i)
         {
         int xi, yi, zi;
         OffsetToIndex(i, xi, yi, zi);
         //double currx = Model.GetXCoordinates()[xi] + Model.GetXCellSizes()[xi] / 2.0;
         //double curry = Model.GetYCoordinates()[yi] + Model.GetYCellSizes()[yi] / 2.0;
         //double currz = Model.GetZCoordinates()[zi] + Model.GetZCellSizes()[zi] / 2.0;
         for (size_t j = 0; j < nmod; ++j)
         {
         int xj, yj, zj;
         OffsetToIndex(j, xj, yj, zj);
         //double x = Model.GetXCoordinates()[xj] + Model.GetXCellSizes()[xj] / 2.0;
         //double y = Model.GetYCoordinates()[yj] + Model.GetYCellSizes()[yj] / 2.0;
         //double z = Model.GetZCoordinates()[zj] + Model.GetZCellSizes()[zj] / 2.0;
         //double r = std::sqrt(
         //    jif3D::pow2(currx - x) + jif3D::pow2(curry - y)
         //        + jif3D::pow2(currz - z));
         double r = std::sqrt(
         jif3D::pow2(xi - xj) + jif3D::pow2(yi - yj) + jif3D::pow2(zi - zj));
         const double ra = std::abs(r / a);
         Cm(i, j) =
         (ra == 0.0) ?
         1.0 :
         factor * std::pow(ra, nu) * boost::math::cyl_bessel_k(nu, ra);

         }

         }*/
        /*Cmv(0) = Cm(0, 0);
         for (size_t i = 1; i < nmod; ++i)
         {
         Cmv(i) = Cm(i, 0);
         Cmv(2*nmod- i) = Cm(i, 0);
         }
         Cmv(nmod) = Cmv(nmod+1);
         kissfft<double> fft(nmod * 2, false);
         fft.transform(&Cmv.data()[0], &Cmf.data()[0]);*/
        // TODO Auto-generated constructor stub
        /*std::ofstream ovfile("cmv.out");
         for (auto val : Cmv)
         {
         ovfile << std::abs(val) << "\n";
         }

         std::ofstream offile("cmf.out");
         for (auto val : Cmf)
         {
         offile << std::abs(val) << "\n";
         }*/
      }

    jif3D::rvec StochasticCovariance::ApplyCovar(const jif3D::rvec &vector) const
      {

        if (previous_vec.empty() || previous_vec.size() != vector.size())
          {
            previous_vec.resize(vector.size(), 0.0);
            previous_result.resize(vector.size(), 0.0);
          }
        else
          {
            if (std::equal(vector.begin(), vector.end(), previous_vec.begin(),
                [](double a, double b)
                  {
                    return boost::math::relative_difference(a,b) < 0.001;
                  }))
              {
                std::cout << "Using previous covariance result " << std::endl;
                return previous_result;
              }
          }
        std::fill(previous_result.begin(), previous_result.end(), 0.0);
        //std::cout << " Applying covariance ";
        double factor = sigma * sigma / (std::pow(2, nu - 1) * boost::math::tgamma(nu));
        boost::posix_time::ptime starttime =
            boost::posix_time::microsec_clock::local_time();
        const size_t nmod = nx * ny * nz;

        double min = 0.0;
        double max = std::sqrt(jif3D::pow2(nx) + jif3D::pow2(ny) + jif3D::pow2(nz));
        double steps = 1000;
        double delta = (max - min) / steps;
        if (a > 0.0)
          {
#pragma omp parallel for default(shared)
            for (size_t i = 0; i < nmod; ++i)
              {
                int xi, yi, zi;
                OffsetToIndex(i, xi, yi, zi);
                //double currx = Model.GetXCoordinates()[xi] + Model.GetXCellSizes()[xi] / 2.0;
                //double curry = Model.GetYCoordinates()[yi] + Model.GetYCellSizes()[yi] / 2.0;
                //double currz = Model.GetZCoordinates()[zi] + Model.GetZCellSizes()[zi] / 2.0;
                int startx = std::max(0, xi - distindex);
                int endx = std::min(nx, xi + distindex);
                int starty = std::max(0, yi - distindex);
                int endy = std::min(ny, yi + distindex);
                int startz = std::max(0, zi - distindex);
                int endz = std::min(nz, zi + distindex);
                for (size_t j = startx; j < endx; ++j)
                  {
                    for (size_t k = starty; k < endy; ++k)
                      {
                        for (size_t l = startz; l < endz; ++l)
                          {
                            int offset = IndexToOffset(j, k, l);
                            double r = std::sqrt(
                                jif3D::pow2(xi - j) + jif3D::pow2(yi - k)
                                    + jif3D::pow2(zi - l));
                            int index = std::round((r - min) / delta);
                            double inter = values[index];
                            //+ (values[index + 1] - values[index]) / delta
                            //		* (r - index * delta);
                            previous_result(i) += inter * vector(offset);
                          }
                      }
                  }

              }
            previous_vec = vector;
          }
        else
          {
            typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXi;
            typedef Eigen::Map<const VectorXi> MapTypeConst;
            typedef Eigen::Map<VectorXi> MapType;

            MapTypeConst v(&(vector.data()[0]), vector.size(), 1);
            MapType r(&(previous_result.data()[0]), previous_result.size(), 1);
            r = A * v;

            boost::posix_time::ptime endtime =
                boost::posix_time::microsec_clock::local_time();
            //double time = (endtime - starttime).total_seconds();
          }
        //std::cout << " took " << time << " seconds " << std::endl;
        return previous_result;
        /*const size_t nmod = vector.size() * 2;
         jif3D::cvec inv(nmod, 0.0), invf(nmod, 0.0);
         jif3D::rvec result(vector.size(), 0.0);
         std::copy(vector.begin(), vector.end(), inv.begin());
         kissfft<double> fft(nmod, false);
         fft.transform(&inv.data()[0], &invf.data()[0]);
         invf = ublas::element_prod(invf, Cmf);
         kissfft<double> ifft(nmod, true);
         ifft.transform(&invf.data()[0], &inv.data()[0]);
         for (size_t i = 0; i < vector.size(); ++i)
         {
         result(i) = 1.0 / double(nmod) * inv(i).real();
         }
         return result;*/
      }

    jif3D::rvec StochasticCovariance::ApplyInvCovar(const jif3D::rvec &vector) const
      {

        //MatrixReplacement A(nx, ny, nz, a, nu, sigma);
        //jif3D::rvec m(nx * ny * nz, 0.0);
        //std::generate(m.begin(), m.end(), drand48);
        //std::iota( m.begin(), m.end(),1);
        //m(nx * ny * nz / 2 + nz / 2) = 1.0;
        Eigen::VectorXd b(nx * ny * nz), x;
        std::copy(vector.begin(), vector.end(), b.data());
        //b = Eigen::Map<Eigen::VectorXd>(m.data(),m.size());

        //Eigen::BiCGSTAB<MatrixReplacement, Eigen::IdentityPreconditioner> bicg;
        //bicg.compute(A);
        //x = bicg.solve(b);
        if (a < 0.0)
          {
            x = cg.solve(b);
            std::cout << "CG:       #iterations: " << cg.iterations()
                << ", estimated error: " << cg.error() << std::endl;

            jif3D::rvec Result(vector.size(), 0.0);
            std::copy(x.data(), x.data() + vector.size(), Result.begin());
            return Result;
          }
        // if (!HaveInv)
        //   {
        //     solver.compute(Cm);
        //  }
        //HaveInv = true;
        //const size_t nmod = vector.size();
        //jif3D::rvec result(nmod);
        //typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXi;
        //typedef Eigen::Map<const VectorXi> MapTypeConst;
        //typedef Eigen::Map<VectorXi> MapType;

        //MapTypeConst v(&(vector.data()[0]), nmod, 1);
        //MapType r(&(result.data()[0]), nmod, 1);
        //r = solver.solve(v);
        //return result;

        /*const size_t nmod = vector.size() * 2;
         jif3D::cvec inv(nmod, 0.0), invf(nmod, 0.0);
         jif3D::rvec result(vector.size(), 0.0);
         std::copy(vector.begin(), vector.end(), inv.begin());
         kissfft<double> fft(nmod, false);
         fft.transform(&inv.data()[0], &invf.data()[0]);
         double maxCmf = std::abs(*std::max_element(Cmf.begin(), Cmf.end(),
         [](std::complex<double> a, std::complex<double> b)
         { return std::abs(a) > std::abs(b);}));
         double threshold = 1e-6 * maxCmf;
         std::ofstream invfile("inv.out");
         for (size_t i = 0; i < nmod; ++i)
         {
         invfile << i << " " << std::abs(invf(i));
         if (std::abs(Cmf(i)) < threshold)
         {
         invf(i) = invf(i) / threshold;
         }
         else
         {
         invf(i) = invf(i) / Cmf(i);
         }
         invfile << " " << std::abs(invf(i)) << " " << std::abs(Cmf(i)) << "\n";
         }*/
        //invf = ublas::element_div(invf, Cmf);
        /*kissfft<double> ifft(nmod, true);
         std::fill_n(inv.begin(),nmod,0.0);
         ifft.transform(&invf.data()[0], &inv.data()[0]);
         for (size_t i = 0; i < vector.size(); ++i)
         {
         result(i) = 1.0 / double(nmod) * inv(i).real();
         }
         return result;*/
        return vector;
      }

    StochasticCovariance::~StochasticCovariance()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
