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
    protected:
<<<<<<< .mine
      ThreeDGravityImplementation &Imp;
=======
      rmat CurrentSensitivities;
>>>>>>> .r123
      void CheckModelConsistency(const ThreeDGravityModel &Model);
    public:
<<<<<<< .mine
      virtual rvec Calculate(const ThreeDGravityModel &Model) = 0;
      ThreeDGravityCalculator(ThreeDGravityImplementation &TheImp);
=======
      rmat &SetCurrentSensitivities(){return CurrentSensitivities;}
      virtual void HandleSensitivities() = 0;
      virtual rvec Calculate(const ThreeDGravityModel &Model,
          ThreeDGravityImplementation &Imp) = 0;
      ThreeDGravityCalculator();
>>>>>>> .r123
      virtual ~ThreeDGravityCalculator();
      };

  }

#endif /* THREEDGRAVITYCALCULATOR_H_ */
