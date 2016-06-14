//============================================================================
// Name        : time_vector.cpp
// Author      : Jun 24, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/VecMat.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <cstdlib>
int main()
  {
    const size_t nelements = 1e6;

    jif3D::rvec ModelCov(nelements), CovGrad(nelements), OldGradient(nelements);
    double omega = 0.0, alpha = 0.0;
    std::generate(ModelCov.begin(), ModelCov.end(), drand48);
    std::generate(CovGrad.begin(), CovGrad.end(), drand48);
    std::generate(OldGradient.begin(), OldGradient.end(), drand48);

    boost::posix_time::ptime firststarttime =
        boost::posix_time::microsec_clock::local_time();
    for (size_t i = 0; i < nelements; ++i)
      {
        const double factor = CovGrad(i) / ModelCov(i);
        omega += CovGrad(i) * factor;
        alpha += OldGradient(i) * factor;
      }
    boost::posix_time::ptime firstendtime =
        boost::posix_time::microsec_clock::local_time();
    std::cout << omega << " " << alpha << std::endl;
    boost::posix_time::ptime secondstarttime =
        boost::posix_time::microsec_clock::local_time();
    omega = ublas::inner_prod(CovGrad, ublas::element_div(CovGrad, ModelCov));
    alpha = ublas::inner_prod(OldGradient,
        ublas::element_div(CovGrad, ModelCov));
    boost::posix_time::ptime secondendtime =
        boost::posix_time::microsec_clock::local_time();

    std::cout << (firstendtime - firststarttime).total_microseconds() << " "
        << (secondendtime - secondstarttime).total_microseconds() << std::endl;
    std::cout << omega << " " << alpha << std::endl;
  }
