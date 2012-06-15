//============================================================================
// Name        : VecMat.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef VECMAT_H_
#define VECMAT_H_
/*! \file VecMat.h
 * This file provides the basic includes for all matrix and vector operations. We store the matrices in column major format to
 * enable interfacing with fortran codes for blas, lapack etc.
 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <complex>

namespace ublas = boost::numeric::ublas;

namespace jiba
  {
    //! A matrix with complex entries
    typedef boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major>
        cmat;
    //! A matrix with real entries
    typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major>
        rmat;
    //! A sparse matrix with real valued entries
    typedef boost::numeric::ublas::mapped_matrix<double,boost::numeric::ublas::column_major>
        map_mat;
    //! The operator matrix for the spatial derivatives is sparse so we provide a typedef for convenience
    typedef boost::numeric::ublas::compressed_matrix<double,
        boost::numeric::ublas::column_major> comp_mat;
    //! A complex vector
    typedef boost::numeric::ublas::vector<std::complex<double> > cvec;
    //! A real vector
    typedef boost::numeric::ublas::vector<double> rvec;
  }
#endif /*VECMAT_H_*/
