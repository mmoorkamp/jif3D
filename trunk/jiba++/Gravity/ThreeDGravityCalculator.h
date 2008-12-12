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
namespace jiba
  {

    class ThreeDGravityCalculator
      {
      private:
        rmat CurrentSensitivities;
    protected:
      ThreeDGravityImplementation &Imp;


      void CheckModelConsistency(const ThreeDGravityModel &Model);
    public:

      virtual rvec Calculate(const ThreeDGravityModel &Model) = 0;
      ThreeDGravityCalculator(ThreeDGravityImplementation &TheImp);
      rmat &SetCurrentSensitivities(){return CurrentSensitivities;}
      virtual void HandleSensitivities(const size_t measindex) = 0;
      ThreeDGravityCalculator();
      virtual ~ThreeDGravityCalculator();
      };

  }

#endif /* THREEDGRAVITYCALCULATOR_H_ */
