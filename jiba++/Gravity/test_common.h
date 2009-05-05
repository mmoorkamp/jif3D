//============================================================================
// Name        : test_common.h
// Author      : Feb 11, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef TEST_COMMON_H_
#define TEST_COMMON_H_

#include <stdlib.h>
#include "ThreeDGravityModel.h"
#include <boost/numeric/conversion/cast.hpp>

//a helper function to create a model dimension of random size
jiba  ::ThreeDModelBase::t3DModelDim GenerateDimension()
    {

      const size_t DimLength = rand() % 11 +15;
      jiba::ThreeDModelBase::t3DModelDim TestDim(boost::extents[DimLength]);
      for (size_t i = 0; i < DimLength; ++i)
        {
          TestDim[i] = rand() % 25 + 10;
        }
      return TestDim;
    }

  //create a random density model
  void MakeRandomModel(jiba::ThreeDGravityModel &Model, const size_t nmeas, const bool withbackground = true)
    {
      srand(time(NULL));

      jiba::ThreeDModelBase::t3DModelDim XDim = GenerateDimension();
      jiba::ThreeDModelBase::t3DModelDim YDim = GenerateDimension();
      jiba::ThreeDModelBase::t3DModelDim ZDim = GenerateDimension();
      const size_t xsize = XDim.size();
      const size_t ysize = YDim.size();
      const size_t zsize = ZDim.size();
      Model.SetXCellSizes().resize(boost::extents[xsize]);
      Model.SetXCellSizes() = XDim;
      Model.SetYCellSizes().resize(boost::extents[ysize]);
      Model.SetYCellSizes() = YDim;
      Model.SetZCellSizes().resize(boost::extents[zsize]);
      Model.SetZCellSizes() = ZDim;
      int xlength = boost::numeric_cast<int>(floor(std::accumulate(XDim.begin(), XDim.end(), 0.0)));
      int ylength = boost::numeric_cast<int>(floor(std::accumulate(YDim.begin(), YDim.end(), 0.0)));

      Model.SetDensities().resize(boost::extents[xsize][ysize][zsize]);

      for (size_t i = 0; i < xsize; ++i)
      for (size_t j = 0; j < ysize; ++j)
      for (size_t k = 0; k < zsize; ++k)
        {
          if (i < xsize / 2)
          Model.SetDensities()[i][j][k] = double(rand() % 50) / 10.0 + 1.0;
          else
          Model.SetDensities()[i][j][k] = -double(rand() % 50) / 10.0 + 1.0;
        }
      //generate measurement  points
      //the z-axis is positive down, so we choose negative z-coordinates
      // => we are measuring above the surface
      for (size_t i = 0; i < nmeas; ++i)
      Model.AddMeasurementPoint(rand() % xlength + 1, rand() % ylength + 1,
          -(rand() % 100 + 1.0));
      if (withbackground)
        {
          //generate a random background
          const size_t nbglayers = rand() % 20;
          std::vector<double> bg_dens(nbglayers,0.0), bg_thick(nbglayers,0.0);
          for (size_t k = 0; k < nbglayers; ++k)
            {

              bg_dens.at(k) = double(rand() % 50) / 10.0 + 1.0;
              bg_thick.at(k) = rand() % 50 + 1.0;
            }
          Model.SetBackgroundDensities(bg_dens);
          Model.SetBackgroundThicknesses(bg_thick);
        }
    }


#endif /* TEST_COMMON_H_ */
