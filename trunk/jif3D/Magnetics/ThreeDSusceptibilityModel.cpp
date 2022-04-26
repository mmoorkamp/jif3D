//============================================================================
// Name        : ThreeDMagneticModel.cpp
// Author      : 11 Oct 2013
// Version     :
// Copyright   : 2013, mm489
//============================================================================

#include "ThreeDSusceptibilityModel.h"



namespace jif3D
  {

    ThreeDSusceptibilityModel::ThreeDSusceptibilityModel()
      {
        // TODO Auto-generated constructor stub

      }


    ThreeDSusceptibilityModel::~ThreeDSusceptibilityModel()
      {
        // TODO Auto-generated destructor stub
      }


    ThreeDSusceptibilityModel& ThreeDSusceptibilityModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
      }

  } /* namespace jif3D */
