//============================================================================
// Name        : ThreeDGravityCalculator.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "ThreeDGravityCalculator.h"

namespace jiba
  {

    ThreeDGravityCalculator::ThreeDGravityCalculator()
      {
        // TODO Auto-generated constructor stub

      }

    ThreeDGravityCalculator::~ThreeDGravityCalculator()
      {
        // TODO Auto-generated destructor stub
      }

    void ThreeDGravityCalculator::CheckModelConsistency(const ThreeDGravityModel &Model)
      {
        //get the amount of cells in each direction
        const size_t xsize = Model.GetDensities().shape()[0];
        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];
        //do some sanity checks
        assert(xsize == Model.GetXCoordinates().shape()[0]);
        assert(ysize == Model.GetYCoordinates().shape()[0]);
        assert(zsize == Model.GetZCoordinates().shape()[0]);
        assert(xsize == Model.GetXCellSizes().shape()[0]);
        assert(ysize == Model.GetYCellSizes().shape()[0]);
        assert(zsize == Model.GetZCellSizes().shape()[0]);

        // make sure we have coordinates for all sites
        const size_t nmeas = Model.GetMeasPosX().size();
        assert(nmeas == Model.GetMeasPosY().size());
        assert(nmeas == Model.GetMeasPosZ().size());
      }
  }
