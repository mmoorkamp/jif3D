//============================================================================
// Name        : test_fileutil.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2010, MM
//============================================================================


//test the file utility functions
#define BOOST_TEST_MODULE FileUtil test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/algorithm/string.hpp>
#include <stdlib.h>
#include "FileUtil.h"

BOOST_AUTO_TEST_SUITE( FileUtil_Test_Suite )

BOOST_AUTO_TEST_CASE  (test_fileextension)
    {
      std::string OneDot("test.ext");
      BOOST_CHECK(boost::algorithm::equals( jiba::GetFileExtension(OneDot),".ext"));
      std::string TwoDots("test.false.ext");
      BOOST_CHECK(boost::algorithm::equals( jiba::GetFileExtension(TwoDots),".ext"));
      std::string NoDots("test");
      BOOST_CHECK(boost::algorithm::equals( jiba::GetFileExtension(NoDots),""));
    }

  BOOST_AUTO_TEST_CASE (test_findtoken_positive)
    {
      std::ofstream outfile;
      const std::string filename("token.test");
      outfile.open(filename.c_str());
      srand(time(0));
      size_t tokenline = rand() % 500;
      for (size_t i = 0; i < tokenline; ++i)
        {
          outfile << rand();
          char c = (rand() % 26) + 'a';
          outfile << c << "\n";
        }
      std::string Token("findme");
      outfile << Token << "\n";
      std::string Afterline("after");
      outfile << Afterline << "\n";
      size_t garbagelines = rand() % 500;
      for (size_t i = 0; i < garbagelines; ++i)
        {
          outfile << rand();
          char c = (rand() % 26) + 'a';
          outfile << c << rand() << "\n";
        }
      outfile.close();
      std::ifstream infile;
      infile.open(filename.c_str());
      std::string found = jiba::FindToken(infile,Token);
      std::string AfterComp;
      infile >> AfterComp;
      infile.close();
      BOOST_CHECK(boost::algorithm::equals(Token,found));
      BOOST_CHECK(boost::algorithm::equals(Afterline,AfterComp));

      boost::filesystem::remove(filename.c_str());
    }


  BOOST_AUTO_TEST_CASE (test_findtoken_negative)
    {
      std::ofstream outfile;
      const std::string filename("token.test");
      outfile.open(filename.c_str());
      srand(time(0));
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
      BOOST_CHECK_THROW(jiba::FindToken(infile,Token),jiba::FatalException);


      boost::filesystem::remove(filename.c_str());
    }
  BOOST_AUTO_TEST_SUITE_END()
