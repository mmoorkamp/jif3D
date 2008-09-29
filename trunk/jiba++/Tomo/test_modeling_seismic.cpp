#define BOOST_TEST_MODULE SeismicModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include "Podvin.h"
#include "modeling_seismic.h"
#include <boost/test/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Seismic_Test_Suite )

BOOST_AUTO_TEST_CASE    (memory_test)
      {
        // check 0 length allocation
        BOOST_MESSAGE("We should see a warning now: ");
        BOOST_CHECK(jiba::memory(NULL,0,1,"test") == NULL);
        // check simple allocation
        size_t length = 100;
        char * mem_test = jiba::memory(NULL,length,1,"test");
        BOOST_CHECK( mem_test != NULL);
        for (size_t i =0; i < length; ++i) mem_test[i] = i;
        // check reallocation
        mem_test = jiba::memory(mem_test,2*length,1,"test");
        for (size_t i =0; i < length; ++i)
          {
            BOOST_CHECK( mem_test[i] == i);
          }
        std::vector<char> dummy;
        size_t failsize = dummy.max_size() * 2; // this should be enough to fail
        BOOST_CHECK_THROW(jiba::memory(NULL,failsize,1,"throw"),std::runtime_error);
        BOOST_CHECK_THROW(jiba::memory(mem_test,failsize,1,"throw"),std::runtime_error);
      }

    BOOST_AUTO_TEST_CASE(podvin3d_test)
      {
        const int nx =2;
        const int ny =2;
        const int nz =2;
        const int nparam = nx * ny *nz;
        float HS[nparam];
        float T[nparam];
        float XS = 0.5;
        float YS = 0.5;
        float ZS = 0.5;
        float HS_EPS_INIT = 0.001;
        int MSG = 2;
        std::fill_n(HS,nparam,1.0);
        std::fill_n(T,nparam,0.0);
        int status = time_3d(HS, T, nx,ny,nz,XS, YS, ZS, HS_EPS_INIT, MSG);

        std::cout << T[0] << " " <<  sqrt(XS*XS+YS*YS+ZS*ZS)  <<std::endl;
      }

    BOOST_AUTO_TEST_SUITE_END()
