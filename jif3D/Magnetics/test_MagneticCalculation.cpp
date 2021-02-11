//============================================================================
// Name        : test_ScalarCalculation.cpp
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../GravMag/MinMemGravMagCalculator.h"
#include "OMPMagneticImp.h"
#include "MagneticData.h"
#include "MagneticTransforms.h"
#include <boost/math/constants/constants.hpp>
#include <cmath>

#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>

BOOST_AUTO_TEST_SUITE( Magnetic_Test_Suite )

    double MagBlock(double xpos, double width, double depth, double height,
        double magnetization)
      {
        double m = width / 2.0;
        double factor = magnetization / (2 * boost::math::constants::pi<double>());
        double Bz = factor
            * (std::atan2(xpos + m, depth) - std::atan2(xpos - m, depth)
                - std::atan2(xpos + m, depth + height)
                + std::atan2(xpos - m, depth + height));
        return Bz;
      }

//compare calculation with text-book results
    BOOST_AUTO_TEST_CASE(model_magnetics_boxcomp_test)
      {
//        const double measx = 9.0; // unused
        const double measy = 10.0;
        const double measz = -1.0;

        jif3D::ThreeDMagneticModel MagneticTest;
        jif3D::MagneticData Data;
        //create a model of 10x10x10 cells with 2m length in each dimension
        const size_t ncells = 10;
        const double cellsize = 2.0;
        MagneticTest.SetMeshSize(ncells, ncells, ncells);

        jif3D::ThreeDModelBase::t3DModelDim XCS(ncells, cellsize);
        jif3D::ThreeDModelBase::t3DModelDim YCS(ncells, cellsize * 100.0);

        jif3D::ThreeDModelBase::t3DModelDim ZCS(ncells, cellsize / 4.0);

        double height = std::accumulate(ZCS.begin(), ZCS.end(), 0.0);
        double widthx = std::accumulate(XCS.begin(), XCS.end(), 0.0);
        double widthy = std::accumulate(YCS.begin(), YCS.end(), 0.0);
        std::cout << " Width x: " << widthx << " Width y: " << widthy << " Height: "
            << height << std::endl;

        MagneticTest.SetXCellSizes(XCS);
        MagneticTest.SetYCellSizes(YCS);
        MagneticTest.SetZCellSizes(ZCS);
        double susceptibility = 0.2;
        MagneticTest.SetOrigin(-widthx / 2.0, -widthy / 2.0, 0);
        for (size_t i = 0; i < ncells; ++i)
          for (size_t j = 0; j < ncells; ++j)
            for (size_t k = 0; k < ncells; ++k)
              {
                MagneticTest.SetSusceptibilities()[i][j][k] = susceptibility;

              }
        const size_t nmeas = 200;
        for (size_t i = 0; i < nmeas; ++i)
          {
            Data.AddMeasurementPoint(-20.0 + i / 5.0, measy, measz);
          }

        double inclination = 90 / 180.0 * boost::math::constants::pi<double>();
        double declination = 0.0 / 180.0 * boost::math::constants::pi<double>();
        double fieldstrength = 50000;
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::MagneticData> CalculatorType;
        boost::shared_ptr<jif3D::OMPMagneticImp> Imp(
            new jif3D::OMPMagneticImp(inclination, declination, fieldstrength));
        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Imp));
        jif3D::rvec magmeas(Calculator->Calculate(MagneticTest, Data));

        BOOST_CHECK_EQUAL(nmeas * 3, magmeas.size());
        Calculator->SetDataTransform(
            boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                new jif3D::TotalFieldAnomaly(inclination, declination)));
        jif3D::rvec magtotal(Calculator->Calculate(MagneticTest, Data));
        std::ofstream toutfile("magt.out");
        for (size_t i = 0; i < magtotal.size(); ++i)
          {

            double Bzref = MagBlock(Data.GetMeasPosX()[i], widthx, -measz, height,
                fieldstrength * susceptibility);
            toutfile << Data.GetMeasPosX()[i] << " " << magtotal(i) << " " << Bzref
                << std::endl;
            BOOST_CHECK_CLOSE(magtotal(i), Bzref, 0.1);
          }
      }

    BOOST_AUTO_TEST_CASE(model_magnetics_simpeg_test)
      {
//        const double measx = 9.0; // unused
        const double measy = 10.0;
        const double measz = -1.0;

        jif3D::ThreeDMagneticModel MagneticTest;
        jif3D::MagneticData Data;
        //create a model of 10x10x10 cells with 2m length in each dimension
        const size_t ncells = 10;
        const double cellsize = 2.0;
        MagneticTest.SetMeshSize(ncells, ncells, ncells);

        jif3D::ThreeDModelBase::t3DModelDim XCS(ncells, cellsize);
        jif3D::ThreeDModelBase::t3DModelDim YCS(ncells, cellsize * 2.0);

        jif3D::ThreeDModelBase::t3DModelDim ZCS(ncells, cellsize / 4.0);

        double height = std::accumulate(ZCS.begin(), ZCS.end(), 0.0);
        double widthx = std::accumulate(XCS.begin(), XCS.end(), 0.0);
        double widthy = std::accumulate(YCS.begin(), YCS.end(), 0.0);
        std::cout << " Width x: " << widthx << " Width y: " << widthy << " Height: "
            << height << std::endl;

        MagneticTest.SetXCellSizes(XCS);
        MagneticTest.SetYCellSizes(YCS);
        MagneticTest.SetZCellSizes(ZCS);
        double susceptibility = 0.2;
        MagneticTest.SetOrigin(-widthx / 2.0, -widthy / 2.0, 0);
        for (size_t i = 0; i < ncells; ++i)
          for (size_t j = 0; j < ncells; ++j)
            for (size_t k = 0; k < ncells; ++k)
              {
                MagneticTest.SetSusceptibilities()[i][j][k] = 0.2;

              }
        const size_t nmeas = 200;
        for (size_t i = 0; i < nmeas; ++i)
          {
            Data.AddMeasurementPoint(-20.0 + i / 5.0, measy, measz);
          }

        double inclination = 15 / 180.0 * boost::math::constants::pi<double>();
        double declination = 10.0 / 180.0 * boost::math::constants::pi<double>();
        double fieldstrength = 50000;
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::MagneticData> CalculatorType;
        boost::shared_ptr<jif3D::OMPMagneticImp> Imp(
            new jif3D::OMPMagneticImp(inclination, declination, fieldstrength));
        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Imp));
        jif3D::rvec magmeas(Calculator->Calculate(MagneticTest, Data));

        BOOST_CHECK_EQUAL(nmeas * 3, magmeas.size());

        std::ifstream xinfile("magx.simpeg");
        std::ifstream yinfile("magy.simpeg");
        std::ifstream zinfile("magz.simpeg");
        std::ifstream tinfile("magt.simpeg");

        std::vector<double> magx, magy, magz,magt;
        std::copy(std::istream_iterator<double>(xinfile), std::istream_iterator<double>(),
            std::back_inserter(magx));
        std::copy(std::istream_iterator<double>(yinfile), std::istream_iterator<double>(),
            std::back_inserter(magy));
        std::copy(std::istream_iterator<double>(zinfile), std::istream_iterator<double>(),
            std::back_inserter(magz));
        std::copy(std::istream_iterator<double>(tinfile), std::istream_iterator<double>(),
            std::back_inserter(magt));
        for (size_t i = 0; i < magx.size(); ++i)
          {

            BOOST_CHECK_CLOSE(magmeas(i * 3), magx.at(i), 0.2);
            BOOST_CHECK_CLOSE(magmeas(i * 3 + 1), magy.at(i), 0.2);
            BOOST_CHECK_CLOSE(magmeas(i * 3 + 2), magz.at(i), 0.2);
          }

        Calculator->SetDataTransform(
            boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                new jif3D::TotalFieldAnomaly(inclination, declination)));
        jif3D::rvec magtotal(Calculator->Calculate(MagneticTest, Data));
        for (int i = 0; i < magtotal.size(); ++i)
          {
            BOOST_CHECK_CLOSE(magtotal(i), magt.at(i), 0.2);
          }

      }

    BOOST_AUTO_TEST_SUITE_END()
