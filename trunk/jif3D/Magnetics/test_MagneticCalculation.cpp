//============================================================================
// Name        : test_ScalarCalculation.cpp
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "MagneticTransforms.h"
#include <boost/math/constants/constants.hpp>
#include <cmath>

#include "OMPMagneticSusceptibilityImp.h"
#include "CudaMagneticSusceptibilityImp.h"
#include "TotalFieldMagneticData.h"
#include "OMPMagnetizationImp.h"
#include "ThreeDMagnetizationModel.h"

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

        jif3D::ThreeDSusceptibilityModel MagneticTest;
        jif3D::TotalFieldMagneticData Data;
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
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TotalFieldMagneticData> CalculatorType;
        boost::shared_ptr<jif3D::OMPMagneticSusceptibilityImp> Imp(
            new jif3D::OMPMagneticSusceptibilityImp(inclination, declination,
                fieldstrength));
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

        jif3D::ThreeDSusceptibilityModel MagneticTest;
        jif3D::TotalFieldMagneticData Data;
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
                MagneticTest.SetSusceptibilities()[i][j][k] = susceptibility;

              }
        const size_t nmeas = 200;
        for (size_t i = 0; i < nmeas; ++i)
          {
            Data.AddMeasurementPoint(-20.0 + i / 5.0, measy, measz);
          }

        double inclination = 15 / 180.0 * boost::math::constants::pi<double>();
        double declination = 10.0 / 180.0 * boost::math::constants::pi<double>();
        double fieldstrength = 50000;
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TotalFieldMagneticData> CalculatorType;
        boost::shared_ptr<jif3D::OMPMagneticSusceptibilityImp> Imp(
            new jif3D::OMPMagneticSusceptibilityImp(inclination, declination,
                fieldstrength));
        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Imp));
        jif3D::rvec magmeas(Calculator->Calculate(MagneticTest, Data));

        BOOST_CHECK_EQUAL(nmeas * 3, magmeas.size());

        std::ifstream xinfile("magx.simpeg");
        std::ifstream yinfile("magy.simpeg");
        std::ifstream zinfile("magz.simpeg");
        std::ifstream tinfile("magt.simpeg");

        std::vector<double> magx, magy, magz, magt;
        std::copy(std::istream_iterator<double>(xinfile), std::istream_iterator<double>(),
            std::back_inserter(magx));
        std::copy(std::istream_iterator<double>(yinfile), std::istream_iterator<double>(),
            std::back_inserter(magy));
        std::copy(std::istream_iterator<double>(zinfile), std::istream_iterator<double>(),
            std::back_inserter(magz));
        std::copy(std::istream_iterator<double>(tinfile), std::istream_iterator<double>(),
            std::back_inserter(magt));

        std::ofstream xoutfile("magx.jif3d");
        std::ofstream youtfile("magy.jif3d");
        std::ofstream zoutfile("magz.jif3d");

        for (size_t i = 0; i < magx.size(); ++i)
          {
            xoutfile << magmeas(i * 3) << "\n";
            youtfile << magmeas(i * 3 + 1) << "\n";
            zoutfile << magmeas(i * 3 + 2) << "\n";
            BOOST_CHECK_CLOSE(magmeas(i * 3), magx.at(i), 0.2);
            BOOST_CHECK_CLOSE(magmeas(i * 3 + 1), magy.at(i), 0.2);
            BOOST_CHECK_CLOSE(magmeas(i * 3 + 2), magz.at(i), 0.2);
          }
#ifdef HAVEGPU
        std::cout << "Running CUDA" << std::endl;
        typedef typename jif3D::MinMemGravMagCalculator<jif3D::TotalFieldMagneticData> CalculatorType;
        boost::shared_ptr<jif3D::CudaMagneticSusceptibilityImp> CudaImp(
            new jif3D::CudaMagneticSusceptibilityImp(inclination, declination,
                fieldstrength));
        boost::shared_ptr<CalculatorType> CudaCalculator(new CalculatorType(CudaImp));
        jif3D::rvec cudameas(CudaCalculator->Calculate(MagneticTest, Data));
        for (size_t i = 0; i < magx.size(); ++i)
                 {
            std::cout << magmeas(i * 3) << " " << cudameas(i * 3) << std::endl;
                   BOOST_CHECK_CLOSE(magmeas(i * 3), cudameas(i * 3), 0.2);
                   BOOST_CHECK_CLOSE(magmeas(i * 3 + 1), magmeas(i * 3 +1), 0.2);
                   BOOST_CHECK_CLOSE(magmeas(i * 3 + 2), magmeas(i * 3 +2), 0.2);
                 }
#endif
        Calculator->SetDataTransform(
            boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                new jif3D::TotalFieldAnomaly(inclination, declination)));
        jif3D::rvec magtotal(Calculator->Calculate(MagneticTest, Data));
        for (size_t i = 0; i < magtotal.size(); ++i)
          {
            BOOST_CHECK_CLOSE(magtotal(i), magt.at(i), 0.2);
          }

      }

    BOOST_AUTO_TEST_CASE(model_magnetization_simpeg_test)
      {
//        const double measx = 9.0; // unused
        const double measy = 10.0;
        const double measz = -1.0;

        jif3D::ThreeDMagnetizationModel MagneticTest;
        jif3D::ThreeComponentMagneticData Data;
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

        double inclination = 15 / 180.0 * boost::math::constants::pi<double>();
        double declination = 10.0 / 180.0 * boost::math::constants::pi<double>();
        double fieldstrength = 50000;
        double susceptibility = 0.2;
        const double BxComp = cos(inclination) * cos(declination) * fieldstrength
            * susceptibility;
        const double ByComp = cos(inclination) * sin(declination) * fieldstrength
            * susceptibility;
        const double BzComp = sin(inclination) * fieldstrength * susceptibility;

        MagneticTest.SetXCellSizes(XCS);
        MagneticTest.SetYCellSizes(YCS);
        MagneticTest.SetZCellSizes(ZCS);
        MagneticTest.SetOrigin(-widthx / 2.0, -widthy / 2.0, 0);
        for (size_t i = 0; i < ncells; ++i)
          for (size_t j = 0; j < ncells; ++j)
            for (size_t k = 0; k < ncells; ++k)
              {
                MagneticTest.SetMagnetization_X()[i][j][k] = BxComp;
                MagneticTest.SetMagnetization_Y()[i][j][k] = ByComp;
                MagneticTest.SetMagnetization_Z()[i][j][k] = BzComp;

              }
        const size_t nmeas = 200;
        for (size_t i = 0; i < nmeas; ++i)
          {
            Data.AddMeasurementPoint(-20.0 + i / 5.0, measy, measz);
          }

        typedef typename jif3D::DiskGravMagCalculator<jif3D::ThreeComponentMagneticData> CalculatorType;
        boost::shared_ptr<jif3D::OMPMagnetizationImp> Imp(
            new jif3D::OMPMagnetizationImp());
        boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Imp));
        jif3D::rvec magmeas(Calculator->Calculate(MagneticTest, Data));
        //! This should come from cache
        jif3D::rvec magmeas2(Calculator->Calculate(MagneticTest, Data));
        for (size_t i = 0; i < magmeas.size(); ++i)
          {
            BOOST_CHECK_CLOSE(magmeas(i), magmeas2(i), 0.01);
          }
        BOOST_CHECK_EQUAL(nmeas * 3, magmeas.size());

        std::ifstream xinfile("magx.simpeg");
        std::ifstream yinfile("magy.simpeg");
        std::ifstream zinfile("magz.simpeg");
        std::ifstream tinfile("magt.simpeg");

        std::vector<double> magx, magy, magz, magt;
        std::copy(std::istream_iterator<double>(xinfile), std::istream_iterator<double>(),
            std::back_inserter(magx));
        std::copy(std::istream_iterator<double>(yinfile), std::istream_iterator<double>(),
            std::back_inserter(magy));
        std::copy(std::istream_iterator<double>(zinfile), std::istream_iterator<double>(),
            std::back_inserter(magz));
        std::copy(std::istream_iterator<double>(tinfile), std::istream_iterator<double>(),
            std::back_inserter(magt));

        std::ofstream xoutfile("magx2.jif3d");
        std::ofstream youtfile("magy2.jif3d");
        std::ofstream zoutfile("magz2.jif3d");

        for (size_t i = 0; i < magx.size(); ++i)
          {
            xoutfile << magmeas(i * 3) << "\n";
            youtfile << magmeas(i * 3 + 1) << "\n";
            zoutfile << magmeas(i * 3 + 2) << "\n";
            BOOST_CHECK_CLOSE(magmeas(i * 3), magx.at(i), 0.2);
            BOOST_CHECK_CLOSE(magmeas(i * 3 + 1), magy.at(i), 0.2);
            BOOST_CHECK_CLOSE(magmeas(i * 3 + 2), magz.at(i), 0.2);
          }

        Calculator->SetDataTransform(
            boost::shared_ptr<jif3D::TotalFieldAnomaly>(
                new jif3D::TotalFieldAnomaly(inclination, declination)));
        jif3D::rvec magtotal(Calculator->Calculate(MagneticTest, Data));
        for (size_t i = 0; i < magtotal.size(); ++i)
          {
            BOOST_CHECK_CLOSE(magtotal(i), magt.at(i), 0.2);
          }

      }
    BOOST_AUTO_TEST_SUITE_END()
