//============================================================================
// Name        : ThreeDGravityImplementation.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "ThreeDGravityImplementation.h"
#include "ThreeDGravityCalculator.h"

namespace jiba
  {

    ThreeDGravityImplementation::ThreeDGravityImplementation()
      {

      }

    ThreeDGravityImplementation::~ThreeDGravityImplementation()
      {
      }

    rvec ThreeDGravityImplementation::Calculate(
        const ThreeDGravityModel &Model, ThreeDGravityCalculator &Calculator)
      {

        // get the number of measurements
        // this class is only called from a calculator object
        // which performs consistency check
        const size_t nmeas = Model.GetMeasPosX().size();
        //Cache the coordinate and size information
        //to avoid problems with thread safety
        XCoord.resize(boost::extents[Model.GetXCoordinates().size()]);
        YCoord.resize(boost::extents[Model.GetYCoordinates().size()]);
        ZCoord.resize(boost::extents[Model.GetZCoordinates().size()]);
        XSizes.resize(boost::extents[Model.GetXCellSizes().size()]);
        YSizes.resize(boost::extents[Model.GetYCellSizes().size()]);
        ZSizes.resize(boost::extents[Model.GetZCellSizes().size()]);
        std::copy(Model.GetXCoordinates().begin(),Model.GetXCoordinates().end(),XCoord.begin());
        std::copy(Model.GetYCoordinates().begin(),Model.GetYCoordinates().end(),YCoord.begin());
        std::copy(Model.GetZCoordinates().begin(),Model.GetZCoordinates().end(),ZCoord.begin());
        std::copy(Model.GetXCellSizes().begin(),Model.GetXCellSizes().end(),XSizes.begin());
        std::copy(Model.GetYCellSizes().begin(),Model.GetYCellSizes().end(),YSizes.begin());
        std::copy(Model.GetZCellSizes().begin(),Model.GetZCellSizes().end(),ZSizes.begin());

        //allocate enough memory for all datapoints in the result vector
        rvec result(nmeas * GetDataPerMeasurement());

        // calculate the size of the modelling domain for the background adjustment
        // this way we do not have to recalculate within the loop
        const double modelxwidth = std::accumulate(
            Model.GetXCellSizes().begin(), Model.GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(
            Model.GetYCellSizes().begin(), Model.GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(
            Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(), 0.0);
        //the total number of model parameters, gridded domain + background
        const size_t nmod = Model.GetDensities().num_elements() + Model.GetBackgroundDensities().size();

        // for all measurement points add the responses of the discretized part and the 1D background
        for (size_t i = 0; i < nmeas; ++i)
          {
            // the vector to hold the result for the current measurement
            rvec currresult(GetDataPerMeasurement());

            currresult = CalcGridded(i, Model,
                Calculator.SetCurrentSensitivities());

            //adjust for the effects of finite extents of the grid
            currresult += CalcBackground(i, modelxwidth,
                modelywidth, modelzwidth, Model,
                Calculator.SetCurrentSensitivities());
            //give the calculator object the chance to handle the sensitivities
            //for the current measurement
            Calculator.HandleSensitivities(i);
            std::copy(currresult.begin(), currresult.end(), result.begin() + (i
                * GetDataPerMeasurement()));
          }
        return result;
      }
  }
