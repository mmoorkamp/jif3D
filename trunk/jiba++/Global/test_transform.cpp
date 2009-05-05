//============================================================================
// Name        : test_transform.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE VectorTransform test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "VectorTransform.h"

BOOST_AUTO_TEST_SUITE( VectorTransform_Test_Suite )

//This class is purely to test the functionality of ApplyTransform
//We take a random number of input and output parameters
//and copy the elements from the input that fit into the output
//if we have more output parameters we fill with zeros
class  TestTransform : public jiba::VectorTransform
    {
    private:
      size_t nin;
      size_t nout;
    public:
      virtual size_t GetInputSize()
        {
          return nin;
        }
      virtual size_t GetOutputSize()
        {
          return nout;
        }
      virtual jiba::rvec Transform(const jiba::rvec &InputVector)
        {
          jiba::rvec Out(nout);
          std::fill(Out.begin(),Out.end(),0.0);
          size_t end = std::min(nin,nout);
          for (size_t i = 0; i < end; ++i)
            {
              Out(i) = InputVector(i);
            }
          return Out;
        }
      //we are not interested in the derivative here
      virtual jiba::rmat Derivative(const jiba::rvec &InputVector)
        {
          return jiba::rmat();

        }
      TestTransform(const size_t inputsize,const size_t outputsize):
      nin(inputsize), nout(outputsize)
        {}
      virtual ~TestTransform()
        {}
    };

  BOOST_AUTO_TEST_CASE(test_apply)
    {
      // we test a few random combinations of input
      // and output size
      const size_t ntest = 10;
      for (size_t i = 0; i < ntest; ++i)
        {
          const size_t nin = rand() %10 + 2;
          const size_t nout = rand() %10 + 1;
          TestTransform Transform(nin,nout);
          const size_t nelements = nin * ((rand() % 10) + 1);
          jiba::rvec InVector(nelements);
          jiba::rvec OutVector(jiba::ApplyTransform(InVector,Transform));
          if (nin <= nout)
            {
              size_t outindex = 0;
              for (size_t j = 0; j < nelements; j+= nin)
                {
                  for (size_t k = 0; k < nin; ++k)
                    {
                      BOOST_CHECK(InVector(j+k) == OutVector(outindex));
                      ++outindex;
                    }
                  for (size_t k = 0; k < nout - nin; ++k)
                    {
                      BOOST_CHECK(OutVector(outindex) == 0.0);
                      ++outindex;
                    }
                }
            }
          else
            {
              size_t inindex = 0;
              for (size_t j =0; j < OutVector.size(); j+=nout)
                {
                  for (size_t k = 0; k < nout; ++k)
                    {
                      BOOST_CHECK(OutVector(j+k) == InVector(inindex));
                      inindex++;
                    }
                  inindex += nin-nout;
                }
            }
          //end of one random test run
        }
    }

  BOOST_AUTO_TEST_SUITE_END()
