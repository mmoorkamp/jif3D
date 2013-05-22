//============================================================================
// Name        : test_transform.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

//test the vector transforms
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
class  TestTransform : public jif3D::VectorTransform
    {
    private:
      //the number of input parameters
      size_t nin;
      //number of output parameters
      size_t nout;
    public:
      //first access functions that need to be implemented
      virtual size_t GetInputSize()
        {
          return nin;
        }
      virtual size_t GetOutputSize()
        {
          return nout;
        }
      //the transform that tests the functionality
      virtual jif3D::rvec Transform(const jif3D::rvec &InputVector)
        {
          jif3D::rvec Out(nout);
          std::fill(Out.begin(),Out.end(),0.0);
          size_t end = std::min(nin,nout);
          for (size_t i = 0; i < end; ++i)
            {
              Out(i) = InputVector(i);
            }
          return Out;
        }
      //we are not interested in the derivative here
      virtual jif3D::rmat Derivative(const jif3D::rvec &InputVector)
        {
          return jif3D::rmat();

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
          //make a random number for the input and output sizes
          const size_t nin = rand() %10 + 2;
          const size_t nout = rand() %10 + 1;
          TestTransform Transform(nin,nout);
          //the number of elements in the vector can be different again
          const size_t nelements = nin * ((rand() % 10) + 1);
          //make a random input vector
          jif3D::rvec InVector(nelements);
          std::generate(InVector.begin(),InVector.end(),drand48);
          //transform the vector
          jif3D::rvec OutVector(jif3D::ApplyTransform(InVector,Transform));
          //check that the transform has been applied correctly
          //regardless of the size of the input/output
          if (nin <= nout)
            {
              //if we have more output values, the first elements have to match
              //and the others have to be zero
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
             //if we nout is less than nin we "calculate" the output
            //from a block of nin values, for each of these blocks
            //we take the first nout values, this is like an invariant
            //for MT or FTG
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
