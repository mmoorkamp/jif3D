//============================================================================
// Name        : test_Rinterface.cpp
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE ThreeDGravityModel test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/numeric/conversion/cast.hpp>


#include "test_common.h"
#include "BasicGravElements.h"
#include "MinMemGravityCalculator.h"
#include "FullSensitivityGravityCalculator.h"
#include "ScalarOMPGravityImp.h"
#include "TensorOMPGravityImp.h"
#include "ThreeDGravityFactory.h"

BOOST_AUTO_TEST_SUITE( TensorGravity_Test_Suite )



//write a C++ vectorial quantity into a file so that R can understand the values
template<typename VectorType>
void WriteVectorToScript(std::ofstream &file, VectorType thevector,
    std::string Name)
  {
    file << Name << "<-c(";
    copy(thevector.begin(), thevector.end() - 1, std::ostream_iterator<double>(
            file, ","));
    file << thevector[thevector.size() - 1];
    file << ")\n";
  }

//Make a random model without background layers and write the information
//about the model into a file so that R can understand it
void PrepareModelForR(jif3D::ThreeDGravityModel &GravityTest,
    std::ofstream &scriptfile)
  {
    //first we setup our local model with some test measurements
    const size_t nmeas = 10;
    MakeRandomModel(GravityTest, nmeas);
    //we generate a script file for R that produces the same model
    //the build environment has copied these libraries in the right path
    scriptfile << "dyn.load(\"libmodelbase.so\") \n";
    scriptfile << "dyn.load(\"libgravity.so\") \n";
    scriptfile << "source(\"Gravity/GravForward.R\") \n";
    scriptfile << "alloc<-AllocateModel(1,1) \n";
    //write the different model quantities into the script
    //the cell size coordinates
    WriteVectorToScript(scriptfile, GravityTest.GetXCellSizes(), "XSizes");
    WriteVectorToScript(scriptfile, GravityTest.GetYCellSizes(), "YSizes");
    WriteVectorToScript(scriptfile, GravityTest.GetZCellSizes(), "ZSizes");
    //and the measurement coordinates
    WriteVectorToScript(scriptfile, GravityTest.GetMeasPosX(), "XMeasPos");
    WriteVectorToScript(scriptfile, GravityTest.GetMeasPosY(), "YMeasPos");
    WriteVectorToScript(scriptfile, GravityTest.GetMeasPosZ(), "ZMeasPos");
    //copy the 3D density model into a vector
    jif3D::rvec DensityVector(GravityTest.GetDensities().num_elements());
    copy(GravityTest.GetDensities().origin(),
        GravityTest.GetDensities().origin()
        + GravityTest.GetDensities().num_elements(), DensityVector.begin());
    //write the density vector into the script as well
    WriteVectorToScript(scriptfile, DensityVector, "Densities");
    std::vector<double> DummyBG(1,0.0);
    WriteVectorToScript(scriptfile, DummyBG, "BGDensities");
    WriteVectorToScript(scriptfile, DummyBG, "BGThicknesses");
  }

//just a dummy to shut up the test system in case we do not compile the proper R tests
BOOST_AUTO_TEST_CASE(dummy_test)
  {
   BOOST_CHECK(true);
  }


//we only compile and run this if we want to test the R interface
#ifdef TESTR
//test the scalar forward interface for R
//this needs the current svn of boost test, as boost test in 1.35.0 has problems with the system call
BOOST_AUTO_TEST_CASE(R_scalar_interface_test)
  {
    //create a 3D Gravity object
    jif3D::ThreeDGravityModel GravityTest;
    //create a random model and write the information into the R-script file
    std::ofstream Rscript("scalar_test.R");
    PrepareModelForR(GravityTest, Rscript);
    //calculate our results
    boost::shared_ptr<jif3D::MinMemGravityCalculator> ScalarCalculator(jif3D::CreateGravityCalculator<jif3D::MinMemGravityCalculator>::MakeScalar());
    jif3D::rvec scalarmeas(
        ScalarCalculator->Calculate(GravityTest));

    //finish the R script
    //call the R interface function
    Rscript
    << " raw<-system.time(result<-GravScalarForward(XSizes,YSizes,ZSizes,Densities,BGDensities,BGThicknesses,XMeasPos,YMeasPos,ZMeasPos))\n";
    Rscript
    << " cached<-system.time(result2<-GravScalarForward(XSizes,YSizes,ZSizes,Densities,BGDensities,BGThicknesses,XMeasPos,YMeasPos,ZMeasPos))\n";
    Rscript << " sink(\"scalar_output\")\n";
    Rscript << " cat(result$GravAcceleration)\n";
    Rscript << " sink(\"r_timing\")\n";
    Rscript << " cat(raw)\n";
    Rscript << " cat(cached)\n";
    Rscript << " q()\n";
    Rscript << std::flush;
    //execute R with the script
    system("R  --slave -f scalar_test.R");
    //read in the output R has generated
    std::ifstream routput("scalar_output");
    std::vector<double> rvalues;
    //this is simply an ascii file with a bunch of numbers
    copy(std::istream_iterator<double>(routput), std::istream_iterator<double>(),
        back_inserter(rvalues));
    //first of all we should have values for each measurement
    BOOST_CHECK_EQUAL(rvalues.size(), scalarmeas.size());
    //and they should be equal, the tolerance is 0.01% as writing to the file truncates the numbers
    for (size_t i = 0; i < scalarmeas.size(); ++i)
      {
        BOOST_CHECK_CLOSE(scalarmeas(i), rvalues.at(i), 0.01);
      }

  }

//test the tensor forward interface for R
//this needs the current svn of boost test, as boost test in 1.35.0 has problems with the system call
BOOST_AUTO_TEST_CASE(R_tensor_interface_test)
  {
    jif3D::ThreeDGravityModel GravityModel;

    std::ofstream Rscript("tensor_test.R");
    PrepareModelForR(GravityModel, Rscript);
    boost::shared_ptr<jif3D::MinMemGravityCalculator> TensorCalculator(jif3D::CreateGravityCalculator<jif3D::MinMemGravityCalculator>::MakeTensor());
    jif3D::rvec tensormeas(
        TensorCalculator->Calculate(GravityModel));
    Rscript << " result<-GravTensorForward(XSizes,YSizes,ZSizes,Densities,XMeasPos,YMeasPos,ZMeasPos)\n";
    Rscript << " sink(\"tensor_output\")\n";
    Rscript << " cat(result$GravAcceleration)\n";
    Rscript << " q()\n";
    Rscript << std::flush;
    //execute R with the script
    system("R  --slave -f tensor_test.R");
    //read in the output R has generated
    std::ifstream routput("tensor_output");
    std::vector<double> rvalues;
    //this is simply an ascii file with a bunch of numbers
    copy(std::istream_iterator<double>(routput), std::istream_iterator<double>(),
        back_inserter(rvalues));
    //first of all we should have values for each measurement
    //tensormeas is a vector of matrices while rvalues is just a vector where 9 consecutive elements
    //correspond to 1 FTG matrix
    BOOST_CHECK_EQUAL(rvalues.size(), tensormeas.size());
    //and they should be equal, the tolerance is 0.01% as writing to the file truncates the numbers
    for (size_t i = 0; i < tensormeas.size(); ++i)
      {
        BOOST_CHECK_CLOSE(tensormeas(i), rvalues.at(i), 0.01);

      }

  }
#endif


BOOST_AUTO_TEST_SUITE_END()
