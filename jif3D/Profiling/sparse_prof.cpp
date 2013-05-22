//============================================================================
// Name        : reg_prof.cpp
// Author      : 16 Aug 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include "../Global/VecMat.h"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <cstdlib>

int main()
  {
    for (size_t nelem = 10; nelem < 1000000; nelem *= 2)
      {
        jif3D::comp_mat Matrix(nelem, nelem, 2);
        jif3D::rvec Vector(nelem), Result(nelem);

        boost::posix_time::ptime firststarttime =
            boost::posix_time::microsec_clock::local_time();
        for (size_t i = 0; i < nelem - 1; ++i)
          {
            Matrix.push_back(i, i, 1.0);
            Matrix.push_back(i, nelem - 1 - i, 2.0);
          }

        boost::posix_time::ptime firstendtime =
            boost::posix_time::microsec_clock::local_time();

        boost::posix_time::ptime secondstarttime =
            boost::posix_time::microsec_clock::local_time();
        ublas::axpy_prod(trans(Matrix), Vector, Result);
        boost::posix_time::ptime secondendtime =
            boost::posix_time::microsec_clock::local_time();

        std::cout << nelem << " " << (firstendtime - firststarttime).total_microseconds()
            << " " << (secondendtime - secondstarttime).total_microseconds() << std::endl;
      }
  }
