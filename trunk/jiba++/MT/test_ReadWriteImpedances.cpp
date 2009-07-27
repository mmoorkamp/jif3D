//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE ReadWriteImpedances test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "ReadWriteImpedances.h"
#include "X3DModel.h"
BOOST_AUTO_TEST_SUITE( ReadWriteImpedances_Suite )

BOOST_AUTO_TEST_CASE  (read_write_netcdf_test)
    {
      const size_t nfreq = 5;
      const size_t nstat = 7;
      const std::string filename("imp.nc");
      std::vector<double> Frequencies(nfreq);
      std::vector<double> XCoord(nstat),YCoord(nstat),ZCoord(nstat);
      const size_t ndata = nfreq*nstat*8;
      jiba::rvec Impedances(ndata);

      std::generate_n(Frequencies.begin(),nfreq,drand48);
      sort(Frequencies.begin(),Frequencies.end());
      std::generate_n(XCoord.begin(),nstat,drand48);
      std::generate_n(YCoord.begin(),nstat,drand48);
      std::generate_n(ZCoord.begin(),nstat,drand48);
      std::generate_n(Impedances.begin(),ndata,drand48);
      jiba::WriteImpedancesToNetCDF(filename,Frequencies,XCoord,YCoord,ZCoord,Impedances);

      std::vector<double> ReadFrequencies;
      std::vector<double> ReadXCoord,ReadYCoord,ReadZCoord;
      jiba::rvec ReadImpedances;
      jiba::ReadImpedancesFromNetCDF(filename,ReadFrequencies,ReadXCoord,ReadYCoord,ReadZCoord,ReadImpedances);
      for (size_t i = 0; i < nfreq; ++i)
        {
          BOOST_CHECK_CLOSE(Frequencies[i],ReadFrequencies[i],0.001);
        }
      for (size_t i = 0; i < nstat; ++i)
        {
          BOOST_CHECK_CLOSE(XCoord[i],ReadXCoord[i],0.001);
          BOOST_CHECK_CLOSE(YCoord[i],ReadYCoord[i],0.001);
          BOOST_CHECK_CLOSE(ZCoord[i],ReadZCoord[i],0.001);
        }
      for (size_t i = 0; i < ndata; ++i)
        {
          BOOST_CHECK_CLOSE(Impedances(i),ReadImpedances(i),0.001);
        }
    }

  BOOST_AUTO_TEST_SUITE_END()
