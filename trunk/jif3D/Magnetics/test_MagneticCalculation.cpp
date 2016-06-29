//============================================================================
// Name        : test_ScalarCalculation.cpp
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../GravMag/MinMemGravMagCalculator.h"
#include "OMPMagneticImp.h"

#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/math/constants/constants.hpp>




using namespace boost::assign;

BOOST_AUTO_TEST_SUITE( Magnetic_Test_Suite )



    //check that for a box the results are independent of the discretization
    BOOST_AUTO_TEST_CASE(model_gravity_boxcomp_test)
      {
//        const double measx = 9.0; // unused
        const double measy = 10.0;
        const double measz = -1.0;

        jif3D::ThreeDMagneticModel MagneticTest;
        //create a model of 10x10x10 cells with 2m length in each dimension
        const size_t ncells = 10;
        const double cellsize = 2.0;
        MagneticTest.SetMeshSize(ncells,ncells,ncells);

        for (size_t i = 0; i < ncells; ++i)
          {
        	MagneticTest.SetXCellSizes()[i] = cellsize;
        	MagneticTest.SetYCellSizes()[i] = cellsize;
        	MagneticTest.SetZCellSizes()[i] = cellsize;
          }


        for (size_t i = 0; i < ncells; ++i)
          for (size_t j = 0; j < ncells; ++j)
            for (size_t k = 0; k < ncells; ++k)
              {
                MagneticTest.SetSusceptibilities()[i][j][k] = 1.0;

              }
        const size_t nmeas = 200;
        for (size_t i = 0; i < nmeas; ++i)
        {
        	MagneticTest.AddMeasurementPoint(-10.0 + i/5.0, measy, measz);
        }

        typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeDMagneticModel> CalculatorType;
        boost::shared_ptr<jif3D::OMPMagneticImp> Imp(new jif3D::OMPMagneticImp(acos(-1.0)/2.0));
        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Imp));
        jif3D::rvec magmeas(Calculator->Calculate(MagneticTest));

        BOOST_CHECK_EQUAL(nmeas * 3, magmeas.size());

        std::ofstream xoutfile("magx.out");
        std::ofstream youtfile("magy.out");
        std::ofstream zoutfile("magz.out");

        for (size_t i = 0; i < magmeas.size(); i+=3)
        {
        	xoutfile << MagneticTest.GetMeasPosX()[i/3] << " " <<  magmeas(i) << std::endl;
        	youtfile << MagneticTest.GetMeasPosX()[i/3] << " " <<magmeas(i+1)<< std::endl;
        	zoutfile << MagneticTest.GetMeasPosX()[i/3] << " " <<magmeas(i+2)<< std::endl;
        }

      }




    BOOST_AUTO_TEST_SUITE_END()
