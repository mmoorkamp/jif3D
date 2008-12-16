//============================================================================
// Name        : ThreeDGravityCalculator.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#ifndef THREEDGRAVITYCALCULATOR_H_
#define THREEDGRAVITYCALCULATOR_H_

#include "ThreeDGravityModel.h"
#include "ThreeDGravityImplementation.h"
#include "ScalarOMPGravityImp.h"
#include "ScalarCudaGravityImp.h"
#include "TensorOMPGravityImp.h"
#include <boost/shared_ptr.hpp>

namespace jiba
  {

    template<class CalculatorClass>
    class CreateGravityCalculator
      {
    public:
      typedef boost::shared_ptr<CalculatorClass> sp_CalculatorClass;
      static boost::shared_ptr<CalculatorClass> MakeScalar(bool wantcuda =
          false);
      static boost::shared_ptr<CalculatorClass> MakeTensor();
      };

    template<class CalculatorClass>
    boost::shared_ptr<CalculatorClass> CreateGravityCalculator<CalculatorClass>::MakeScalar(
        bool wantcuda)
      {
        boost::shared_ptr<ThreeDGravityImplementation> Imp;
        if (wantcuda)
          {
            Imp = boost::shared_ptr<ThreeDGravityImplementation>(
                new ScalarCudaGravityImp);
          }
        else
          {
            Imp = boost::shared_ptr<ThreeDGravityImplementation>(
                new ScalarOMPGravityImp);
          }
        return boost::shared_ptr<CalculatorClass>(new CalculatorClass(Imp));
      }

    template<class CalculatorClass>
    boost::shared_ptr<CalculatorClass> CreateGravityCalculator<CalculatorClass>::MakeTensor()
      {
        boost::shared_ptr<ThreeDGravityImplementation> Imp(
            new TensorOMPGravityImp);
        return boost::shared_ptr<CalculatorClass>(new CalculatorClass(Imp));
      }

    class ThreeDGravityCalculator
      {
    private:
      rmat CurrentSensitivities;
    protected:
      boost::shared_ptr<ThreeDGravityImplementation> Imp;

      void CheckModelConsistency(const ThreeDGravityModel &Model);
    public:

      virtual rvec Calculate(const ThreeDGravityModel &Model) = 0;
      ThreeDGravityCalculator(
          boost::shared_ptr<ThreeDGravityImplementation> TheImp);
      rmat &SetCurrentSensitivities()
        {
          return CurrentSensitivities;
        }
      virtual void HandleSensitivities(const size_t measindex) = 0;
      ThreeDGravityCalculator();
      virtual ~ThreeDGravityCalculator();
      };

  }

#endif /* THREEDGRAVITYCALCULATOR_H_ */
