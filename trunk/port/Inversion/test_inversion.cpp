//============================================================================
// Name        : test_inversion.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#define BOOST_TEST_MODULE Inversion test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <stdlib.h>
#include <algorithm>
#include <numeric>
#include "MatrixTools.h"
#include "LinearInversion.h"
#include "GeneralizedInverse.h"


BOOST_AUTO_TEST_SUITE( Inversion_Test_Suite )

//test that the SVD fulfills the absolute minimum requirement
//calculate the SVD of an Identity matrix
BOOST_AUTO_TEST_CASE  (basic_svd_test)
    {
      srand(time(NULL));
      const size_t msize = 4;
      jif3D::rmat TestMatrix(msize,msize);
      std::fill_n(TestMatrix.data().begin(),msize*msize,0.0);
      for (size_t i = 0; i < msize; ++i)
        {
          TestMatrix(i,i) = 1.0;
        }
      jif3D::rvec s;
      jif3D::rmat u,vt;
      jif3D::SVD(TestMatrix,s,u,vt);
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
      jif3D::rmat TestMatrix(msize,msize);
      std::fill_n(TestMatrix.data().begin(),msize*msize,0.0);
      for (size_t i = 0; i < msize; ++i)
        {
          TestMatrix(i,i) = 1.0;
        }
      double det = jif3D::Determinant(TestMatrix);
      BOOST_CHECK_CLOSE(det,1.0,std::numeric_limits<float>::epsilon());

    }

  //the determinant should be the product of the eigenvalues
  //because we use an SVD, it will always be positive,
  //so we take the absolute value
  BOOST_AUTO_TEST_CASE (determinant_test)
    {
      const size_t msize = 4;
      jif3D::rmat TestMatrix(msize,msize);
      for (size_t i = 0; i < msize; ++i)
      for (size_t j = 0; j < msize; ++j)
      TestMatrix(i, j) = double(rand() % 10);

      double det = jif3D::Determinant(TestMatrix);
      jif3D::rvec s;
      jif3D::rmat u,vt;
      jif3D::SVD(TestMatrix,s,u,vt);
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
      jif3D::rmat TestMatrix(msize, msize);
      jif3D::rmat InverseMatrix(msize, msize);
      for (size_t i = 0; i < msize; ++i)
      for (size_t j = 0; j < msize; ++j)
      TestMatrix(i, j) = double(rand());

      jif3D::rmat MatrixCopy(TestMatrix);
      jif3D::InvertMatrix(MatrixCopy, InverseMatrix);
      jif3D::rmat MultMatrix(boost::numeric::ublas::prec_prod(TestMatrix,
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
      jif3D::rmat TestMatrix(msize, msize);
      jif3D::rmat InverseMatrix(msize, msize);
      srand(time(NULL));
      for (size_t i = 0; i < msize; ++i)
      for (size_t j = 0; j < msize; ++j)
      TestMatrix(i, j) = double(rand());

      jif3D::GeneralizedInverse()(TestMatrix, InverseMatrix);
      jif3D::rmat MultMatrix(boost::numeric::ublas::prec_prod(TestMatrix,
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
      jif3D::rmat Sensitivities(ndata,nparam);
      jif3D::rvec DataVec(ndata), DataError(ndata);
      jif3D::rvec ModelWeight(nparam), DataSpaceInvModel(nparam), ModelSpaceInvModel(nparam), QuasiNewtonInvModel(nparam);
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

      jif3D::rmat OrigSens(Sensitivities);
      jif3D::rvec OrigData(DataVec);
      jif3D::rvec StartingModel(QuasiNewtonInvModel);
      jif3D::DataSpaceInversion()(Sensitivities, DataVec, ModelWeight, DataError,1.0,
          DataSpaceInvModel);
      DataVec = OrigData;
      Sensitivities = OrigSens;
      jif3D::ModelSpaceInversion()(Sensitivities, DataVec, ModelWeight, DataError,1.0,
          ModelSpaceInvModel);
      DataVec = OrigData;
      Sensitivities = OrigSens;
      jif3D::QuasiNewtonInversion()(Sensitivities, DataVec, ModelWeight, DataError, 1.0,
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
      jif3D::rmat Sensitivities(ndata,nparam);
      jif3D::rvec DataVec(ndata), DataError(ndata);
      jif3D::rvec TrueModel(nparam);
      jif3D::rvec ModelWeight(nparam), DataSpaceInvModel(nparam), ModelSpaceInvModel(nparam), QuasiNewtonInvModel(nparam);
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
      jif3D::rvec OrigData(DataVec);
      jif3D::rmat OrigSens(Sensitivities);
      jif3D::rvec StartingModel(QuasiNewtonInvModel);
      jif3D::QuasiNewtonInversion()(Sensitivities, DataVec, ModelWeight, DataError, 0.0,
          QuasiNewtonInvModel);

      DataVec = OrigData;
      Sensitivities = OrigSens;
      jif3D::ModelSpaceInversion()(Sensitivities, DataVec, ModelWeight, DataError,0.0,
          ModelSpaceInvModel);

      DataVec = OrigData;
      Sensitivities = OrigSens;
      jif3D::DataSpaceInversion()(Sensitivities, DataVec, ModelWeight, DataError,0.0,
          DataSpaceInvModel);
      for (size_t i = 0; i < nparam; ++i)
        {
          BOOST_CHECK_CLOSE(TrueModel(i),DataSpaceInvModel(i),0.01);
          BOOST_CHECK_CLOSE(TrueModel(i),ModelSpaceInvModel(i),0.01);
          BOOST_CHECK_CLOSE(TrueModel(i),StartingModel(i) + QuasiNewtonInvModel(i),0.01);
        }

    }
  BOOST_AUTO_TEST_SUITE_END()
