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

    void DataSpaceInversion::operator()(rmat &Sensitivities,
        rvec &Data, const rvec &WeightVector, const rvec &DataError,
        const double evalthresh, const double lambda, rvec &InvModel)
      {
        const size_t nmeas = Data.size();
        const size_t nparm = Sensitivities.size2();
        for (size_t i = 0; i < nmeas; ++i)
          {
            boost::numeric::ublas::matrix_row<jiba::rmat>(
                Sensitivities, i) /= DataError(i);
            Data(i) /= DataError(i);
          }

        FilteredSens.resize(Sensitivities.size2(),Sensitivities.size1());
        FilteredSens.assign(trans(Sensitivities));
        for (size_t i = 0; i < nparm; ++i)
          {
            boost::numeric::ublas::matrix_row<jiba::rmat> CurrentRow(
                FilteredSens, i);

            CurrentRow *= WeightVector(i);
          }

        jiba::rmat Gamma(nmeas, nmeas);
        atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, Sensitivities,
            FilteredSens, 0.0, Gamma);
        for (size_t i = 0; i < nmeas; ++i)
          {
            Gamma(i,i) += lambda;
          }

        jiba::rmat DataInverse(nmeas, nmeas);
        jiba::GeneralizedInverse()(Gamma, DataInverse, evalthresh, 1.0);

        InvModel.resize(nparm);

        jiba::rmat ModelInverse(nparm, nmeas);
        atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, FilteredSens, DataInverse,
            0.0, ModelInverse);

        atlas::gemv(ModelInverse, Data, InvModel);
      }

    void ModelSpaceInversion::operator()(rmat &Sensitivities,
        rvec &Data, const rvec &WeightVector, const rvec &DataError,
        const double evalthresh, const double lambda, rvec &InvModel)
      {
        const size_t nmeas = Data.size();
        const size_t nparm = Sensitivities.size2();
        for (size_t i = 0; i < nmeas; ++i)
          {
            boost::numeric::ublas::matrix_row<jiba::rmat>(
                Sensitivities, i) /= sqrt(DataError(i));
            Data(i) /= sqrt(DataError(i));
          }

        //FilteredSens.resize(Sensitivities.size2(),Sensitivities.size1());
        //FilteredSens.assign(trans(Sensitivities));
        //for (size_t i = 0; i < nparm; ++i)
         // {
          //  boost::numeric::ublas::matrix_row<jiba::rmat> CurrentRow(
           //     FilteredSens, i);
//
  //          CurrentRow /= WeightVector(i);
    //      }

        jiba::rmat Gamma(nparm, nparm);
        atlas::gemm(CblasTrans, CblasNoTrans, 1.0,
            Sensitivities, Sensitivities, 0.0, Gamma);
        for (size_t i = 0; i < nparm; ++i)
          Gamma(i,i) += lambda * 1./WeightVector(i);

        jiba::rmat Inverse(nparm, nparm);
        jiba::GeneralizedInverse()(Gamma, Inverse, evalthresh, 1.0);

        InvModel.resize(nparm);

        jiba::rmat ModelInverse(nparm, nmeas);
        atlas::gemm(CblasNoTrans, CblasTrans, 1.0, Inverse, Sensitivities,
            0.0, ModelInverse);

        atlas::gemv(ModelInverse, Data, InvModel);
      }
  }
