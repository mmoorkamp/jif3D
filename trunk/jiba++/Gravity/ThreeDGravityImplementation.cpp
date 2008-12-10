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
        // TODO Auto-generated constructor stub

      }

    ThreeDGravityImplementation::~ThreeDGravityImplementation()
      {
        // TODO Auto-generated destructor stub
      }

    rvec ThreeDGravityImplementation::Calculate(const ThreeDGravityModel &Model,ThreeDGravityCalculator &Calculator)
    {
    	 //get the amount of cells in each direction


    	        // make sure we have coordinates for all sites
    	        const size_t nmeas = Model.GetMeasPosX().size();

    	        rvec result(nmeas * GetDataPerMeasurement());

    	        // calculate the size of the modelling domain for the background adjustment
    	        const double modelxwidth = std::accumulate(Model.GetXCellSizes().begin(),
    	            Model.GetXCellSizes().end(), 0.0);
    	        const double modelywidth = std::accumulate(Model.GetYCellSizes().begin(),
    	            Model.GetYCellSizes().end(), 0.0);
    	        const double modelzwidth = std::accumulate(Model.GetZCellSizes().begin(),
    	            Model.GetZCellSizes().end(), 0.0);
    	        const size_t nmod = Model.GetDensities().shape()[0] * Model.GetDensities().shape()[1]
    	                       * Model.GetDensities().shape()[2] + Model.GetBackgroundDensities().size();

    	        // for all measurement points add the responses of the discretized part and the 1D background
    	        for (size_t i = 0; i < nmeas; ++i)
    	          {

    	            rvec currresult(GetDataPerMeasurement());

    	            currresult = CalcGridded(Model.GetMeasPosX()[i],
    	                Model.GetMeasPosY()[i], Model.GetMeasPosZ()[i], Model,Calculator.SetCurrentSensitivities());

    	            //adjust for the effects of finite extents of the grid
    	            currresult += CalcBackground(Model.GetMeasPosX()[i],
    	                Model.GetMeasPosY()[i], Model.GetMeasPosZ()[i], modelxwidth, modelywidth,
    	                modelzwidth, Model, Calculator.SetCurrentSensitivities());
    	            Calculator.HandleSensitivities();
    	            std::copy(currresult.begin(), currresult.end(), result.begin() + (i* GetDataPerMeasurement()));

    	          }
    	        return result;
    }
  }
