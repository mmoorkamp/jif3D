//============================================================================
// Name        : gravgrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/lapack/geev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include "LinearInversion.h"
#include "../Inversion/MatrixTools.h"
#include "../Inversion/GeneralizedInverse.h"

namespace jif3D
  {
    namespace lapack = boost::numeric::bindings::lapack;
    namespace atlas = boost::numeric::bindings::atlas;
    namespace ublas = boost::numeric::ublas;

    /*! Given the sensitivities, data, weights and data error, perform a linear dataspace inversion.
     * See Siripunvaraporn and Egbert, Geophysics 65, 2000 equations 9-11
     * @param Sensitivities The \f$ n \times m\f$ sensitivity matrix, after the call contains the filtered sensitivities
     * @param Data The n component vector containing the data or misfit, invalid after the call
     * @param WeightVector The m-component vector of model weights, the diagonal elements of the model covariance
     * @param DataError The n component vector of error estimates for the data, each element must be > 0
     * @param lambda The lagrangian multiplier to adjust the level of regularization
     * @param InvModel On input: The a priori model, On output The m component vector with the inversion result
     */
    void DataSpaceInversion::operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
        const rvec &DataError, const double lambda, rvec &InvModel)
      {
        //check that all lengths are consistent
        assert(Data.size() == Sensitivities.size1());
        assert(Data.size() == DataError.size());
        assert(WeightVector.size() == Sensitivities.size2());
        const size_t nmeas = Data.size();
        const size_t nparm = Sensitivities.size2();

        //calculate d^ (equation 6)
        atlas::gemv(1.0, Sensitivities, InvModel, 1.0, Data);
        //apply the model covariance to the sensitivity matrix
        for (size_t i = 0; i < nparm; ++i)
          {
            boost::numeric::ublas::matrix_column<jif3D::rmat> CurrentColumn(
                Sensitivities, i);
            CurrentColumn *= sqrt(WeightVector(i));
          }
        //calculate the data space matrix gamma
        jif3D::rmat Gamma(nmeas, nmeas);
        atlas::gemm(CblasNoTrans, CblasTrans, 1.0, Sensitivities,
            Sensitivities, 0.0, Gamma);
        //apply damping
        for (size_t i = 0; i < nmeas; ++i)
          {
            Gamma(i, i) += lambda * DataError(i);
          }
        //invert the data space matrix by solving
        //a linear system, the result is in Data
        lapack::gesv(Gamma, Data);
        //project the inverted matrix back into model space
        //equation 9 in Siripurnvaraporn and Egbert
        for (size_t i = 0; i < nparm; ++i)
          {
            boost::numeric::ublas::matrix_column<jif3D::rmat> CurrentColumn(
                Sensitivities, i);
            CurrentColumn *= sqrt(WeightVector(i));
          }
        //calculate the inversion model
        atlas::gemv(CblasTrans, 1.0, Sensitivities, Data, 0.0, InvModel);
      }

    /*! Given the sensitivities, data, weights and data error, perform a linear classic model space inversion.
     * See Siripurnvaraporn and Egbert, Geophysics 65, 2000 equation 7
     * @param Sensitivities The \f$ n \times m\f$ sensitivity matrix, after the call contains the filtered sensitivities
     * @param Data The n component vector containing the data, after the call conatins the error weighted data
     * @param WeightVector The m-component vector of model weights, the diagonal elements of the model covariance
     * @param DataError The n component vector of error estimates for the data, each element must be > 0
     * @param lambda The lagrangian multiplier to adjust the level of regularivation
     * @param InvModel The m component vector with the inversion result
     */
    void ModelSpaceInversion::operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
        const rvec &DataError,  const double lambda,
        rvec &InvModel)
      {
        const size_t nmeas = Data.size();
        const size_t nparm = Sensitivities.size2();
        atlas::gemv(1.0, Sensitivities, InvModel, 1.0, Data);
        for (size_t i = 0; i < nmeas; ++i)
          {
            boost::numeric::ublas::matrix_row<jif3D::rmat>(Sensitivities, i)
                /= sqrt(DataError(i));
            Data(i) /= sqrt(DataError(i));
          }

        jif3D::rmat Gamma(nparm, nparm);
        atlas::gemm(CblasTrans, CblasNoTrans, 1.0, Sensitivities,
            Sensitivities, 0.0, Gamma);
        for (size_t i = 0; i < nparm; ++i)
          {
            Gamma(i, i) += lambda * 1. / WeightVector(i);
          }
        atlas::gemv(CblasTrans, 1.0, Sensitivities, Data, 0.0, InvModel);
        lapack::gesv(Gamma, InvModel);
      }

    /*! Given the sensitivities, data, weights and data error, perform a Quasi-Newton Model update.
     * This procedure is very similar to ModelSpaceInversion, it differs however in the way the regularization
     * is applied. For details see Tarrantolla eq. 6.319, page 216
     * @param Sensitivities The \f$ n \times m\f$ sensitivity matrix, after the call contains the filtered sensitivities
     * @param Data The n component vector containing the data, after the call conatins the error weighted data
     * @param WeightVector The m-component vector of model weights, the diagonal elements of the model covariance
     * @param DataError The n component vector of error estimates for the data, each element must be > 0
     * @param lambda The lagrangian multiplier to adjust the level of regularivation
     * @param InvModel The m component vector with the inversion result
     */
    void QuasiNewtonInversion::operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
        const rvec &DataError, const double lambda, rvec &InvModel)
      {
        const size_t nmeas = Data.size();
        const size_t nparm = Sensitivities.size2();
        for (size_t i = 0; i < nmeas; ++i)
          {
            boost::numeric::ublas::matrix_row<jif3D::rmat>(Sensitivities, i)
                /= sqrt(DataError(i));
            Data(i) /= sqrt(DataError(i));
          }
        jif3D::rmat Gamma(nparm, nparm);
        atlas::gemm(CblasTrans, CblasNoTrans, 1.0, Sensitivities,
            Sensitivities, 0.0, Gamma);
        for (size_t i = 0; i < nparm; ++i)
          {
            Gamma(i, i) += lambda * 1. / WeightVector(i);
          }
        jif3D::rvec y(nparm);
        atlas::gemv(CblasTrans, 1.0, Sensitivities, Data, 0.0, y);
        for (size_t i = 0; i < nparm; ++i)
          {
            y(i) += lambda * 1. / WeightVector(i) * InvModel(i);
          }
        lapack::gesv(Gamma, y);
        InvModel = y;
      }
}
