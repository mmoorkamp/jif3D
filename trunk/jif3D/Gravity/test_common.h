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
#include <boost/lambda/lambda.hpp>

//a helper function to create a model dimension of random size
jif3D::ThreeDModelBase::t3DModelDim GenerateDimension(const size_t maxcells)
  {
    //we want between 1 and maxcells+1 cells
    const size_t DimLength = rand() % maxcells + 1;
    //allocate memory
    jif3D::ThreeDModelBase::t3DModelDim TestDim(boost::extents[DimLength]);
    //and assign random sizes between 10 and 35
    for (size_t i = 0; i < DimLength; ++i)
      {
        TestDim[i] = rand() % 25 + 10;
      }
    return TestDim;
  }

//create a random density model
template <typename ModelType>
void MakeRandomModel(ModelType &Model, const size_t maxcells,
    const size_t nmeas = 10, const bool withbackground = true)
  {
    srand( time(NULL));
    //make three random axis, each with possibly different lengths and cell sizes
    jif3D::ThreeDModelBase::t3DModelDim XDim = GenerateDimension(maxcells);
    jif3D::ThreeDModelBase::t3DModelDim YDim = GenerateDimension(maxcells);
    jif3D::ThreeDModelBase::t3DModelDim ZDim = GenerateDimension(maxcells);
    const size_t xsize = XDim.size();
    const size_t ysize = YDim.size();
    const size_t zsize = ZDim.size();
    //copy the generated sizes to the model object
    Model.SetXCellSizes().resize(boost::extents[xsize]);
    Model.SetXCellSizes() = XDim;
    Model.SetYCellSizes().resize(boost::extents[ysize]);
    Model.SetYCellSizes() = YDim;
    Model.SetZCellSizes().resize(boost::extents[zsize]);
    Model.SetZCellSizes() = ZDim;
    //get the horizontal dimensions of the model
    int xlength = boost::numeric_cast<int>(floor(std::accumulate(XDim.begin(),
        XDim.end(), 0.0)));
    int ylength = boost::numeric_cast<int>(floor(std::accumulate(YDim.begin(),
        YDim.end(), 0.0)));
    //allocate the grid
    Model.SetData().resize(boost::extents[xsize][ysize][zsize]);
    //and fill the grid with random values
    std::generate_n(Model.SetData().origin(), xsize * ysize * zsize,
        drand48() * boost::lambda::constant(10) + 1.0);
    //generate measurement  points
    //the z-axis is positive down, so we choose negative z-coordinates
    // => we are measuring above the surface
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(rand() % xlength + 1, rand() % ylength + 1,
            -(rand() % 100 + 1.0));
      }
//    //if we also need to generate a background
//    if (withbackground)
//      {
//        //generate a random background
//        const size_t nbglayers = rand() % 20;
//        std::vector<double> bg_dens(nbglayers, 0.0), bg_thick(nbglayers, 0.0);
//        for (size_t k = 0; k < nbglayers; ++k)
//          {
//            bg_dens.at(k) = double(rand() % 50) / 10.0 + 1.0;
//            bg_thick.at(k) = rand() % 50 + 1.0;
//          }
//        Model.SetBackgroundDensities(bg_dens);
//        Model.SetBackgroundThicknesses(bg_thick);
//      }
  }

#endif /* TEST_COMMON_H_ */
