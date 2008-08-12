//============================================================================
// Name        : VecMat.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef VECMAT_H_
#define VECMAT_H_
/*! \file This file provides the basic includes for all matrix and vector operations. We store the matrices in column major format to
 * enable interfacing with fortran codes for blas, lapack etc.
 */

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/lapack/geev.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/vector_traits.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <complex>

namespace jiba {
//! A matrix with complex entries
typedef boost::numeric::ublas::matrix<std::complex<double>,boost::numeric::ublas::column_major > cmat;
//! A matrix with real entries
typedef boost::numeric::ublas::matrix<double,boost::numeric::ublas::column_major > rmat;
//! A complex vector
typedef boost::numeric::ublas::vector<std::complex<double> > cvec;
//! A real vector
typedef boost::numeric::ublas::vector<double> rvec;
}
#endif /*VECMAT_H_*/
