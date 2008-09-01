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
        const double evalthresh, rvec &InvModel)
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

            CurrentRow /= WeightVector(i);
          }

        jiba::rmat Gamma(nmeas, nmeas);
        atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, Sensitivities,
            FilteredSens, 0.0, Gamma);
        //very crude data covariance treatment
        //for (size_t i = 0; i < nmeas; ++i)
        //  {
        //    Gamma(i,i) += 0.02 * 1e-3 * Data(i);
        //  }

        jiba::rmat DataInverse(nmeas, nmeas);
        jiba::GeneralizedInverse(Gamma, DataInverse, evalthresh, 1.0);

        InvModel.resize(nparm);

        jiba::rmat ModelInverse(nparm, nmeas);
        atlas::gemm(CblasNoTrans, CblasNoTrans, 1.0, FilteredSens, DataInverse,
            0.0, ModelInverse);

        atlas::gemv(ModelInverse, Data, InvModel);
      }
  }
