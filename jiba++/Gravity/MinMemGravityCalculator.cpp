//============================================================================
// Name        : MinMemGravityCalculator.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "MinMemGravityCalculator.h"

namespace jiba
  {

    MinMemGravityCalculator::MinMemGravityCalculator()
      {
        // TODO Auto-generated constructor stub

      }

    MinMemGravityCalculator::~MinMemGravityCalculator()
      {
        // TODO Auto-generated destructor stub
      }

    rvec MinMemGravityCalculator::Calculate(const ThreeDGravityModel &Model,
        ThreeDGravityImplementation &Imp)
      {
        CheckModelConsistency(Model);
        //get the amount of cells in each direction
        const size_t xsize = Model.GetDensities().shape()[0];
        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];

        // make sure we have coordinates for all sites
        const size_t nmeas = Model.GetMeasPosX().size();

        rvec result(nmeas * Imp.GetDataPerMeasurement());

        // calculate the size of the modelling domain for the background adjustment
        const double modelxwidth = std::accumulate(Model.GetXCellSizes().begin(),
            Model.GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(Model.GetYCellSizes().begin(),
            Model.GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(Model.GetZCellSizes().begin(),
            Model.GetZCellSizes().end(), 0.0);

        rmat Sensitivities(0, 0);
        ublas::matrix_range<rmat> mr(Sensitivities, ublas::range(0, 0),
            ublas::range(0, 0));
        // for all measurement points add the responses of the discretized part and the 1D background
        for (size_t i = 0; i < nmeas; ++i)
          {

            rvec currresult(Imp.GetDataPerMeasurement());

            currresult = Imp.CalcGridded(Model.GetMeasPosX()[i],
                Model.GetMeasPosY()[i], Model.GetMeasPosZ()[i], Model, mr);

            //adjust for the effects of finite extents of the grid
            currresult += Imp.CalcBackground(Model.GetMeasPosX()[i],
                Model.GetMeasPosY()[i], Model.GetMeasPosZ()[i], modelxwidth, modelywidth,
                modelzwidth, Model, mr);

            std::copy(currresult.begin(), currresult.end(), result.begin() + (i* Imp.GetDataPerMeasurement()));

          }
        return result;

      }

  }
