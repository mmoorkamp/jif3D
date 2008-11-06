//============================================================================
// Name        : gravgrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include "LinearInversion.h"
#include "../Inversion/MatrixTools.h"

namespace jiba
  {
    namespace atlas = boost::numeric::bindings::atlas;
    namespace ublas = boost::numeric::ublas;

    /*! Given the sensitivities, data, weights and data error, perform a linear dataspace inversion.
     * See Siripurnvaraporn and Egbert, Geophysics 65, 2000
     * @param Sensitivities The \f$ n \times m\f$ sensitivity matrix, after the call contains the filtered sensitivities
     * @param Data The n component vector containing the data, after the call conatins the error weighted data
     * @param WeightVector The m-component vector of model weights, the diagonal elements of the model covariance
     * @param DataError The n component vector of error estimates for the data, each element must be > 0
     * @param evalthresh The relative threshold for the eigenvalues to be used for inverting the matrix
     * @param lambda The lagrangian multiplier to adjust the level of regularivation
     * @param InvModel The m component vector with the inversion result
     */
    void DataSpaceInversion::operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
        const rvec &DataError, const double evalthresh, const double lambda,
        rvec &InvModel)
      {
        //check that all lengths are consistent
        assert(Data.size() == Sensitivities.size1());
        assert(Data.size() == DataError.size());
        assert(WeightVector.size() == Sensitivities.size2());
        const size_t nmeas = Data.size();
        const size_t nparm = Sensitivities.size2();
        //weight the data and the sensitivity matrix by the error for each datum
        for (size_t i = 0; i < nmeas; ++i)
          {
            boost::numeric::ublas::matrix_row<jiba::rmat>(Sensitivities, i)
                /= DataError(i);
            Data(i) /= DataError(i);
          }
        //apply the model covariance to the sensitivity matrix
        for (size_t i = 0; i < nparm; ++i)
          {
            boost::numeric::ublas::matrix_column<jiba::rmat> CurrentColumn(
                Sensitivities, i);
            CurrentColumn *= sqrt(WeightVector(i));
          }
        //calculate the data space matrix gamma
        jiba::rmat Gamma(nmeas, nmeas);
        atlas::gemm(CblasNoTrans, CblasTrans, 1.0, Sensitivities,
            Sensitivities, 0.0, Gamma);
        //apply regularization
        for (size_t i = 0; i < nmeas; ++i)
          {
            Gamma(i, i) += lambda;
          }
        //invert the data space matrix
        jiba::rmat DataInverse(nmeas, nmeas);
        //jiba::InvertMatrix(Gamma, DataInverse);
        jiba::GeneralizedInverse()(Gamma, DataInverse, evalthresh, 1.0);

       //project the inverted matrix back into model space
        for (size_t i = 0; i < nparm; ++i)
          {
            boost::numeric::ublas::matrix_column<jiba::rmat> CurrentColumn(
                Sensitivities, i);
            CurrentColumn *= sqrt(WeightVector(i));
          }
        jiba::rmat ModelInverse(nparm, nmeas);
        atlas::gemm(CblasTrans, CblasNoTrans, 1.0, Sensitivities, DataInverse,
            0.0, ModelInverse);
        //calculate the inversion model
        InvModel.resize(nparm);
        atlas::gemv(ModelInverse, Data, InvModel);
      }

    /*! Given the sensitivities, data, weights and data error, perform a linear classic model space inversion.
     * @param Sensitivities The \f$ n \times m\F$ sensitivity matrix, after the call contains the filtered sensitivities
     * @param Data The n component vector containing the data, after the call conatins the error weighted data
     * @param WeightVector The m-component vector of model weights, the diagonal elements of the model covariance
     * @param DataError The n component vector of error estimates for the data, each element must be > 0
     * @param evalthresh The relative threshold for the eigenvalues to be used for inverting the matrix
     * @param lambda The lagrangian multiplier to adjust the level of regularivation
     * @param InvModel The m component vector with the inversion result
     */
    void ModelSpaceInversion::operator()(rmat &Sensitivities, rvec &Data, const rvec &WeightVector,
        const rvec &DataError, const double evalthresh, const double lambda,
        rvec &InvModel)
      {
        const size_t nmeas = Data.size();
        const size_t nparm = Sensitivities.size2();
        for (size_t i = 0; i < nmeas; ++i)
          {
            boost::numeric::ublas::matrix_row<jiba::rmat>(Sensitivities, i)
                /= sqrt(DataError(i));
            Data(i) /= sqrt(DataError(i));
          }

        jiba::rmat Gamma(nparm, nparm);
        atlas::gemm(CblasTrans, CblasNoTrans, 1.0, Sensitivities,
            Sensitivities, 0.0, Gamma);
        for (size_t i = 0; i < nparm; ++i)
          {
            Gamma(i, i) += lambda * 1. / WeightVector(i);
          }
        jiba::rmat Inverse(nparm, nparm);
        jiba::GeneralizedInverse()(Gamma, Inverse, evalthresh, 1.0);

        InvModel.resize(nparm);

        jiba::rmat ModelInverse(nparm, nmeas);
        atlas::gemm(CblasNoTrans, CblasTrans, 1.0, Inverse, Sensitivities, 0.0,
            ModelInverse);

        atlas::gemv(ModelInverse, Data, InvModel);
      }
  }
