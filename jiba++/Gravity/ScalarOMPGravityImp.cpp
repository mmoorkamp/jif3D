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
     * @param measindex the index of the measurement
     * @param xwidth The total width of the discretized model area in x-direction in m
     * @param ywidth The total width of the discretized model area in y-direction in m
     * @param zwidth The total width of the discretized model area in z-direction in m
     * @param Model The gravity model
     * @param Sensitivities If the matrix hold enough elements the sensitivities for the background are stored for a single measurement
     * @return The gravitational acceleration in m/s^2 due to the background
     */
    rvec ScalarOMPGravityImp::CalcBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth,
        const ThreeDGravityModel &Model,
        rmat &Sensitivities)
      {
        //make sure we have thicknesses and densities for all layers
        assert(Model.GetBackgroundDensities().size() == Model.GetBackgroundThicknesses().size());
        const double x_meas = Model.GetMeasPosX()[measindex];
        const double y_meas = Model.GetMeasPosY()[measindex];
        const double z_meas = Model.GetMeasPosZ()[measindex];
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
                currvalue -= CalcGravBoxTerm(x_meas, y_meas, z_meas, 0.0, 0.0,
                    currtop, xwidth, ywidth, currthick);
              }
            //if some of the background coincides and some is below
            if (currtop < zwidth && currbottom > zwidth)

              {
                currvalue -= CalcGravBoxTerm(x_meas, y_meas, z_meas, 0.0, 0.0,
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

    /*! Calculate the gravitational effect of the 3D model at a single measurement site.
     * @param measindex The index of the measurement
     * @param Model The gravity model
     * @param Sensitivities If the matrix hold enough elements the sensitivities for the background are stored for a single measurement
     * @return The gravitational acceleration in m/s^2 due to the model at this site
     */
    rvec ScalarOMPGravityImp::CalcGridded(const size_t measindex,
        const ThreeDGravityModel &Model,
        rmat &Sensitivities)
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
