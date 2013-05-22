//============================================================================
// Name        : ScalarOMPGravityImp.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "ScalarOMPGravityImp.h"
#include "BasicGravElements.h"
#include "GravityBackground.h"
namespace jif3D
  {

    ScalarOMPGravityImp::ScalarOMPGravityImp()
      {

      }

    ScalarOMPGravityImp::~ScalarOMPGravityImp()
      {

      }

    /*!  Calculate the contribution of a layered background to a scalar gravity measurement.
     * @param measindex the index of the measurement
     * @param xwidth The total width of the discretized model area in x-direction in m
     * @param ywidth The total width of the discretized model area in y-direction in m
     * @param zwidth The total width of the discretized model area in z-direction in m
     * @param Model The gravity model
     * @param Sensitivities If the matrix hold enough elements the sensitivities for the background are stored for a single measurement
     * @return The gravitational acceleration in m/s^2 due to the background
     */
    rvec ScalarOMPGravityImp::CalcBackground(const size_t measindex,
        const double xwidth, const double ywidth, const double zwidth,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
        return CalcScalarBackground(measindex, xwidth, ywidth, zwidth, Model,
            Sensitivities);
      }

    /*! Calculate the gravitational effect of the 3D model at a single measurement site.
     * @param measindex The index of the measurement
     * @param Model The gravity model
     * @param Sensitivities If the matrix hold enough elements the sensitivities for the background are stored for a single measurement
     * @return The gravitational acceleration in m/s^2 due to the model at this site
     */
    rvec ScalarOMPGravityImp::CalcGridded(const size_t measindex,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
        //get the dimensions of the model
        const size_t xsize = Model.GetDensities().shape()[0];
        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];
        const double x_meas = Model.GetMeasPosX()[measindex];
        const double y_meas = Model.GetMeasPosY()[measindex];
        const double z_meas = Model.GetMeasPosZ()[measindex];
        const int nmod = xsize * ysize * zsize;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= size_t(nmod));
        double returnvalue = 0.0;
        double currvalue = 0.0;

        //sum up the contributions of all prisms in an openmp parallel loop
#pragma omp parallel default(shared) private(currvalue) reduction(+:returnvalue)
          {
            //instead of nested loops over each dimension, we have one big
            //loop over all elements, this allows for a nearly infinite number
            //of parallel processors
#pragma omp for
            for (int offset = 0; offset < nmod; ++offset)
              {

                //we store the current value for possible sensitivity calculations
                //currvalue contains the geometric term, i.e. the sensitivity
                int xindex, yindex, zindex;
                Model.OffsetToIndex(offset, xindex, yindex, zindex);
                currvalue = CalcGravBoxTerm(x_meas, y_meas, z_meas,
                    XCoord[xindex], YCoord[yindex], ZCoord[zindex],
                    XSizes[xindex], YSizes[yindex], ZSizes[zindex]);
                returnvalue += currvalue
                    * Model.GetDensities()[xindex][yindex][zindex];
                if (storesens)
                  {
                    Sensitivities(0, offset) = currvalue;
                  }
              }

          }//end of parallel section
        rvec returnvector(1);
        returnvector(0) = returnvalue;
        return returnvector;

      }
  }
