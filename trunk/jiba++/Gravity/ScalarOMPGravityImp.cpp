//============================================================================
// Name        : ScalarOMPGravityImp.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "ScalarOMPGravityImp.h"
#include "BasicGravElements.h"

namespace jiba
  {

    ScalarOMPGravityImp::ScalarOMPGravityImp()
      {
        // TODO Auto-generated constructor stub

      }

    ScalarOMPGravityImp::~ScalarOMPGravityImp()
      {
        // TODO Auto-generated destructor stub
      }

    /*!  Calculate the contribution of a layered background to a scalar gravity measurement.
     * @param xmeas The x-coordinate of the measurement point in m
     * @param ymeas The y-coordinate of the measurement point in m
     * @param zmeas The z-coordinate of the measurement point in m
     * @param xwidth The total width of the discretized model area in x-direction in m
     * @param ywidth The total width of the discretized model area in y-direction in m
     * @param zwidth The total width of the discretized model area in z-direction in m
     * @param meas_index The index of the measurement
     * @return The gravitational acceleration in m/s^2 due to the background
     */
    rvec ScalarOMPGravityImp::CalcBackground(const double xmeas,
        const double ymeas, const double zmeas, const double xwidth,
        const double ywidth, const double zwidth,
        const ThreeDGravityModel &Model,
        rmat &Sensitivities)
      {
        //make sure we have thicknesses and densities for all layers
        assert(Model.GetBackgroundDensities().size() == Model.GetBackgroundThicknesses().size());
        const size_t nbglayers = Model.GetBackgroundThicknesses().size();
        double result = 0.0;
        double currtop = 0.0;
        double currvalue = 0.0;
        double currbottom = 0.0;
        const size_t modelsize = Model.GetDensities().shape()[0]
            * Model.GetDensities().shape()[1] * Model.GetDensities().shape()[2];
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= modelsize + nbglayers);
        // for all layers of the background
        for (size_t j = 0; j < nbglayers; ++j)
          {
            const double currthick = Model.GetBackgroundThicknesses()[j];
            // first assume an infinite sheet for the current layer
            currvalue = CalcInfSheetTerm(currthick);
            currbottom = currtop + currthick;
            // and then subtract the value for the modelling domain, as this is already calculated in the discretized routine
            // if the background layer complete coincides with the discretized area
            if (currtop < zwidth && (currbottom <= zwidth))

              {
                currvalue -= CalcGravBoxTerm(xmeas, ymeas, zmeas, 0.0, 0.0,
                    currtop, xwidth, ywidth, currthick);
              }
            //if some of the background coincides and some is below
            if (currtop < zwidth && currbottom > zwidth)

              {
                currvalue -= CalcGravBoxTerm(xmeas, ymeas, zmeas, 0.0, 0.0,
                    currtop, xwidth, ywidth, (zwidth - currtop));
              }
            if (storesens)
              {
                Sensitivities(0, modelsize + j) = currvalue;
              }
            result += currvalue * Model.GetBackgroundDensities()[j];
            currtop += currthick;
          }
        rvec returnvector(1);
        returnvector(0) = result;
        return returnvector;
      }

    /*! Calculate the gravitational effect of the 3D model at a single measurement site. The way we calculate
     * the sensitivity matrix at the moment, the model cannot contain densities of 0 if we
     * store the sensitivity matrix
     * @param x_meas x-coordinate of the measurement site in m
     * @param y_meas y-coordinate of the measurement site in m
     * @param z_meas z-coordinate of the measurement site in m
     * @param meas_index index of the measurement site among all measurements, this is only for storing sensitivities
     * @return The gravitational acceleration in m/s^2 due to the model at this site
     */
    rvec ScalarOMPGravityImp::CalcGridded(const double x_meas,
        const double y_meas, const double z_meas,
        const ThreeDGravityModel &Model,
        rmat &Sensitivities)
      {
        //get the dimensions of the model
        const size_t xsize = Model.GetDensities().shape()[0];
        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];
        const size_t nmod = xsize * ysize * zsize;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= nmod);
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
