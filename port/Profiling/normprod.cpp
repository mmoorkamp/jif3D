//============================================================================
// Name        : normprod.cpp
// Author      : May 30, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Global/VecMat.h"
#include "../Global/NormProd.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <numeric>

int main()
  {
    const size_t nelements = 1e6;
    jif3D::rvec a(nelements), b(nelements), NormDiag(nelements);

    std::generate(a.begin(), a.end(), drand48);
    std::generate(b.begin(), b.end(), drand48);
    std::generate(NormDiag.begin(), NormDiag.end(), drand48);

    double firsttime = 0.0;
    double secondtime = 0.0;
    const size_t ntries = 50;
    for (size_t i = 0; i < ntries; ++i)
      {
        boost::posix_time::ptime firststarttime =
            boost::posix_time::microsec_clock::local_time();
        double myresult = jif3D::NormProd(a, b, NormDiag);
        boost::posix_time::ptime firstendtime =
            boost::posix_time::microsec_clock::local_time();

        boost::posix_time::ptime secondstarttime =
            boost::posix_time::microsec_clock::local_time();
        double uresult = ublas::inner_prod(a, ublas::element_div(b, NormDiag));
        boost::posix_time::ptime secondendtime =
            boost::posix_time::microsec_clock::local_time();
        firsttime += (firstendtime - firststarttime).total_microseconds();
        secondtime += (secondendtime - secondstarttime).total_microseconds();
      }
    std::cout << firsttime << " " << secondtime << std::endl;
  }
