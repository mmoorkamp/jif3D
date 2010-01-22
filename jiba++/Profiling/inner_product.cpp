//============================================================================
// Name        : time_vector.cpp
// Author      : Jun 24, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "cublas.h"
#include "../Global/VecMat.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <numeric>
int main()
  {
    const int nelements = 1e6;

    jiba::rvec a(nelements), b(nelements);

    std::generate(a.begin(), a.end(), drand48);
    std::generate(b.begin(), b.end(), drand48);

    boost::posix_time::ptime firststarttime =
        boost::posix_time::microsec_clock::local_time();
    double stdresult = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
    boost::posix_time::ptime firstendtime =
        boost::posix_time::microsec_clock::local_time();

    boost::posix_time::ptime secondstarttime =
        boost::posix_time::microsec_clock::local_time();
    double uresult = ublas::inner_prod(a, b);
    boost::posix_time::ptime secondendtime =
        boost::posix_time::microsec_clock::local_time();

    double *d_a, *d_b;

    cublasInit();
    cublasAlloc(nelements, sizeof(a[0]), (void**) &d_a);
    cublasAlloc(nelements, sizeof(a[0]), (void**) &d_b);

    cublasSetVector(nelements, sizeof(a[0]), (void**) &a[0], 1, (void**) &d_a,
        1);
    cublasSetVector(nelements, sizeof(a[0]), (void**) &b[0], 1, (void**) &d_b,
        1);
    boost::posix_time::ptime cudastarttime =
        boost::posix_time::microsec_clock::local_time();

    double curesult = cublasDdot(nelements, d_a, 1, d_b, 1);
    boost::posix_time::ptime cudaendtime =
        boost::posix_time::microsec_clock::local_time();
    cublasFree(d_a);
    cublasFree(d_b);

    cublasShutdown();

    std::cout << stdresult << " " << uresult << " " << curesult << std::endl;
    std::cout << (firstendtime - firststarttime).total_microseconds() << " "
        << (secondendtime - secondstarttime).total_microseconds() << " "
        << (cudaendtime - cudastarttime).total_microseconds() << std::endl;
  }
