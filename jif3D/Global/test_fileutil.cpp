//============================================================================
// Name        : test_fileutil.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2010, MM
//============================================================================


//test the file utility functions
#define BOOST_TEST_MODULE FileUtil test
#define BOOST_TEST_MAIN ...
#include "Jif3DTesting.h"
#include <boost/test/floating_point_comparison.hpp>
#include <boost/algorithm/string.hpp>
#include <stdlib.h>
#include "FileUtil.h"

BOOST_AUTO_TEST_SUITE( FileUtil_Test_Suite )

//test whether we recognize the file extension under different circumstances
BOOST_AUTO_TEST_CASE  (test_fileextension)
    {
     //first hte simplest case with a single dot
      std::string OneDot("test.ext");
      BOOST_CHECK(boost::algorithm::equals( jif3D::GetFileExtension(OneDot),".ext"));
      //when there are several dots we consider the part after the last dot as the extension
      std::string TwoDots("test.false.ext");
      BOOST_CHECK(boost::algorithm::equals( jif3D::GetFileExtension(TwoDots),".ext"));
      //no dot measn no extension
      std::string NoDots("test");
      BOOST_CHECK(boost::algorithm::equals( jif3D::GetFileExtension(NoDots),""));
    }

//check whether we correctly find a token if it is in a file
  BOOST_AUTO_TEST_CASE (test_findtoken_positive)
    {
      std::ofstream outfile;
      //we create a file that contains random characters
      //then two tokens and then more random characters
      const std::string filename("token.test");
      outfile.open(filename.c_str());
      srand(time(0));
      //first create up to 500 lines of character noise
      size_t tokenline = rand() % 500;
      for (size_t i = 0; i < tokenline; ++i)
        {
          outfile << rand();
          char c = (rand() % 26) + 'a';
          outfile << c << "\n";
        }
      //then write two lines with tokens that we want to find
      std::string Token("findme");
      outfile << Token << "\n";
      std::string Afterline("after");
      outfile << Afterline << "\n";
      //and then another up to 500 lines of noise
      size_t garbagelines = rand() % 500;
      for (size_t i = 0; i < garbagelines; ++i)
        {
          outfile << rand();
          char c = (rand() % 26) + 'a';
          outfile << c << rand() << "\n";
        }
      outfile.close();
      //now open the file for reading and find the token
      std::ifstream infile;
      infile.open(filename.c_str());
      std::string found = jif3D::FindToken(infile,Token);
      //we also make sure that the next thing we read from the file
      //is the second line we wrote
      std::string AfterComp;
      infile >> AfterComp;
      infile.close();
      //check whether we read what we expected
      BOOST_CHECK(boost::algorithm::equals(Token,found));
      BOOST_CHECK(boost::algorithm::equals(Afterline,AfterComp));
      //and clean up
      boost::filesystem::remove(filename.c_str());
    }

//check that FindToken throws if the required token does not exist in the file
  BOOST_AUTO_TEST_CASE (test_findtoken_negative)
    {
      std::ofstream outfile;
      const std::string filename("token.test");
      outfile.open(filename.c_str());
      srand(time(0));
      //produce up to 500 lines of character noise
      //in theory this could accidently produce the token
      //however it is not likely
      size_t tokenline = rand() % 500;
      for (size_t i = 0; i < tokenline; ++i)
        {
          outfile << rand();
          char c = (rand() % 26) + 'a';
          outfile << c << "\n";
        }
      outfile.close();
      //this token does not exist in the file
      std::string Token("findme");
      std::ifstream infile;
      infile.open(filename.c_str());
      //make sure FindToken throws an exception
      BOOST_CHECK_THROW(jif3D::FindToken(infile,Token),jif3D::FatalException);


      boost::filesystem::remove(filename.c_str());
    }
  BOOST_AUTO_TEST_SUITE_END()
