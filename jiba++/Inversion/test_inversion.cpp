//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Inversion test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "MatrixTools.h"
#include "LinearInversion.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Inversion_Test_Suite )

//test that the SVD fulfills the absolute minimum requirement
//calculate the SVD of an Identity matrix
BOOST_AUTO_TEST_CASE  (basic_svd_test)
    {
      srand(time(NULL));
      const size_t msize = 4;
      jiba::rmat TestMatrix(msize,msize);
      std::fill_n(TestMatrix.data().begin(),msize*msize,0.0);
      for (size_t i = 0; i < msize; ++i)
        {
          TestMatrix(i,i) = 1.0;
        }
      jiba::rvec s;
      jiba::rmat u,vt;
      jiba::SVD(TestMatrix,s,u,vt);
      for (size_t i = 0; i < msize; ++i)
        {
          BOOST_CHECK_CLOSE(s(i),1.0,std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(u(i,i),1.0,std::numeric_limits<float>::epsilon());
          BOOST_CHECK_CLOSE(vt(i,i),1.0,std::numeric_limits<float>::epsilon());
        }

    }
  //test that the determinant of an identity matrix is 1
  BOOST_AUTO_TEST_CASE (basic_determinant_test)
    {
      const size_t msize = 4;
      jiba::rmat TestMatrix(msize,msize);
      std::fill_n(TestMatrix.data().begin(),msize*msize,0.0);
      for (size_t i = 0; i < msize; ++i)
        {
          TestMatrix(i,i) = 1.0;
        }
      double det = jiba::Determinant(TestMatrix);
      BOOST_CHECK_CLOSE(det,1.0,std::numeric_limits<float>::epsilon());

    }

  //the determinant should be the product of the eigenvalues
  //because we use an SVD, it will always be positive,
  //so we take the absolute value
  BOOST_AUTO_TEST_CASE (determinant_test)
    {
      const size_t msize = 4;
      jiba::rmat TestMatrix(msize,msize);
      for (size_t i = 0; i < msize; ++i)
      for (size_t j = 0; j < msize; ++j)
      TestMatrix(i, j) = double(rand() % 10);

      double det = jiba::Determinant(TestMatrix);
      jiba::rvec s;
      jiba::rmat u,vt;
      jiba::SVD(TestMatrix,s,u,vt);
      double sdet = std::accumulate(s.begin(),s.end(),1.0,std::multiplies<double>());
      BOOST_CHECK_CLOSE(fabs(det), sdet,
          1e-4);
    }

  //check that matrix inversion works
  //by generating a random matrix M
  //and checking M^-1M=I
  BOOST_AUTO_TEST_CASE(invert_matrix_test)
    {
      const size_t msize = 4;
      jiba::rmat TestMatrix(msize, msize);
      jiba::rmat InverseMatrix(msize, msize);
      for (size_t i = 0; i < msize; ++i)
      for (size_t j = 0; j < msize; ++j)
      TestMatrix(i, j) = double(rand());

      jiba::rmat MatrixCopy(TestMatrix);
      jiba::InvertMatrix(MatrixCopy, InverseMatrix);
      jiba::rmat MultMatrix(boost::numeric::ublas::prec_prod(TestMatrix,
              InverseMatrix));
      for (size_t i = 0; i < msize; ++i)
        {
          BOOST_CHECK_CLOSE(MultMatrix(i, i), 1.0,
              std::numeric_limits<float>::epsilon());
        }
    }

  //the test for the generalized inverse is more tricky
  //so far we apply the same test as for the regular inverse
  BOOST_AUTO_TEST_CASE(generalized_inverse_test)
    {
      const size_t msize = 4;
      jiba::rmat TestMatrix(msize, msize);
      jiba::rmat InverseMatrix(msize, msize);
      srand(time(NULL));
      for (size_t i = 0; i < msize; ++i)
      for (size_t j = 0; j < msize; ++j)
      TestMatrix(i, j) = double(rand());

      jiba::GeneralizedInverse()(TestMatrix, InverseMatrix);
      jiba::rmat MultMatrix(boost::numeric::ublas::prec_prod(TestMatrix,
              InverseMatrix));
      for (size_t i = 0; i < msize; ++i)
        {
          BOOST_CHECK_CLOSE(MultMatrix(i, i), 1.0,
              std::numeric_limits<float>::epsilon());
        }
    }

  //compare the results of data space inversion and model space inversion
  //they should give identical results
  BOOST_AUTO_TEST_CASE(compare_inversions_test)
    {
      const size_t ndata =5;
      const size_t nparam = 5;
      jiba::rmat Sensitivities(ndata,nparam);
      jiba::rvec DataVec(ndata), DataError(ndata);
      jiba::rvec ModelWeight(nparam), DataSpaceInvModel(nparam), ModelSpaceInvModel(nparam), QuasiNewtonInvModel(nparam);
      std::fill_n(DataSpaceInvModel.begin(),nparam,1.0);
      std::fill_n(ModelSpaceInvModel.begin(),nparam,1.0);
      std::fill_n(QuasiNewtonInvModel.begin(),nparam,1.0);
      std::generate_n(DataVec.begin(),ndata,rand);
      std::generate_n(DataError.begin(),ndata,rand);
      std::generate_n(ModelWeight.begin(),nparam,rand);

      for (size_t i = 0; i < ndata* nparam; ++i)
        {
          Sensitivities.data()[i] = rand();
        }

      jiba::rmat OrigSens(Sensitivities);
      jiba::rvec OrigData(DataVec);
      jiba::rvec StartingModel(QuasiNewtonInvModel);
      jiba::DataSpaceInversion()(Sensitivities, DataVec, ModelWeight, DataError,1.0,
          DataSpaceInvModel);
      DataVec = OrigData;
      Sensitivities = OrigSens;
      jiba::ModelSpaceInversion()(Sensitivities, DataVec, ModelWeight, DataError,1.0,
          ModelSpaceInvModel);
      DataVec = OrigData;
      Sensitivities = OrigSens;
      jiba::QuasiNewtonInversion()(Sensitivities, DataVec, ModelWeight, DataError, 1.0,
          QuasiNewtonInvModel);
      for (size_t i = 0; i < nparam; ++i)
        {
          BOOST_CHECK_CLOSE(DataSpaceInvModel(i),ModelSpaceInvModel(i),0.01);
          BOOST_CHECK_CLOSE(DataSpaceInvModel(i),StartingModel(i) + QuasiNewtonInvModel(i),0.01);
        }
    }

  BOOST_AUTO_TEST_CASE(check_inversions_test)
    {
      const size_t ndata =5;
      const size_t nparam = 5;
      jiba::rmat Sensitivities(ndata,nparam);
      jiba::rvec DataVec(ndata), DataError(ndata);
      jiba::rvec TrueModel(nparam);
      jiba::rvec ModelWeight(nparam), DataSpaceInvModel(nparam), ModelSpaceInvModel(nparam), QuasiNewtonInvModel(nparam);
      std::fill_n(DataSpaceInvModel.begin(),nparam,1.0);
      std::fill_n(ModelSpaceInvModel.begin(),nparam,1.0);
      std::fill_n(QuasiNewtonInvModel.begin(),nparam,1.0);
      std::fill_n(DataError.begin(),ndata,1.0);
      std::fill_n(ModelWeight.begin(),nparam,1.0);
      std::generate_n(TrueModel.begin(),nparam,rand);
      for (size_t i = 0; i < ndata* nparam; ++i)
        {
          Sensitivities.data()[i] = rand();
        }
      DataVec = boost::numeric::ublas::prec_prod(Sensitivities,TrueModel);
      jiba::rvec OrigData(DataVec);
      jiba::rmat OrigSens(Sensitivities);
      jiba::rvec StartingModel(QuasiNewtonInvModel);
      jiba::QuasiNewtonInversion()(Sensitivities, DataVec, ModelWeight, DataError, 0.0,
          QuasiNewtonInvModel);

      DataVec = OrigData;
      Sensitivities = OrigSens;
      jiba::ModelSpaceInversion()(Sensitivities, DataVec, ModelWeight, DataError,0.0,
          ModelSpaceInvModel);

      DataVec = OrigData;
      Sensitivities = OrigSens;
      jiba::DataSpaceInversion()(Sensitivities, DataVec, ModelWeight, DataError,0.0,
          DataSpaceInvModel);
      for (size_t i = 0; i < nparam; ++i)
        {
          BOOST_CHECK_CLOSE(TrueModel(i),DataSpaceInvModel(i),0.01);
          BOOST_CHECK_CLOSE(TrueModel(i),ModelSpaceInvModel(i),0.01);
          BOOST_CHECK_CLOSE(TrueModel(i),StartingModel(i) + QuasiNewtonInvModel(i),0.01);
        }

    }
  BOOST_AUTO_TEST_SUITE_END()
