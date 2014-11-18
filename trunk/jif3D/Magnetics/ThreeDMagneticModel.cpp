//============================================================================
// Name        : ThreeDMagneticModel.cpp
// Author      : 11 Oct 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include "ThreeDMagneticModel.h"
#include <boost/bind.hpp>
#include <cassert>
#include <netcdfcpp.h>
#include <fstream>
#include <iomanip>
#include <iterator>


namespace jif3D
  {

    ThreeDMagneticModel::ThreeDMagneticModel()
      {
        // TODO Auto-generated constructor stub

      }


    ThreeDMagneticModel::~ThreeDMagneticModel()
      {
        // TODO Auto-generated destructor stub
      }

    void ThreeDMagneticModel::WriteNetCDF(const std::string filename) const
      {

      }

    ThreeDMagneticModel& ThreeDMagneticModel::operator=(const ThreeDModelBase& source)
      {
        if (&source != this)
          {
            ThreeDModelBase::operator=(source);
          }
        return *this;
      }

  } /* namespace jif3D */
