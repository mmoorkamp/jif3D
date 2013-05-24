//============================================================================
// Name        : GravityBackground.cpp
// Author      : Jun 18, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "BasicGravElements.h"
#include "ThreeDGravityModel.h"
#include "../Global/VecMat.h"
namespace jif3D
  {

    /*! Calculate the effect of the background layers for a single scalar measurement
     * @param measindex The index of the measurement in the Model
     * @param xwidth The total width of the gridded domain in x-direction in m
     * @param ywidth The total width of the gridded domain in y-direction in m
     * @param zwidth The total width of the gridded domain in z-direction in m
     * @param Model The Gravity model
     * @param Sensitivities If the matrix passed here holds \f$ 1 \times nbg+ngrid \f$ or more elements, store sensitivity information in the right fields
     * @return A vector with a single component that contains the gravitational effect of the background
     */
    rvec CalcScalarBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
        const size_t ndatapermeas = 1;
        const double x_meas = Model.GetMeasPosX()[measindex];
        const double y_meas = Model.GetMeasPosY()[measindex];
        const double z_meas = Model.GetMeasPosZ()[measindex];
        const size_t nbglayers = Model.GetBackgroundThicknesses().size();
        double result = 0.0;
        double currtop = 0.0;
        double currvalue = 0.0;
        double currbottom = 0.0;
        const size_t modelsize = Model.GetDensities().num_elements();
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= modelsize + nbglayers);
        // for all layers of the background
        for (size_t j = 0; j < nbglayers; ++j)
          {
            const double currthick = Model.GetBackgroundThicknesses()[j];
            currbottom = currtop + currthick;
            // first assume an infinite sheet for the current layer
            currvalue = CalcInfSheetTerm(z_meas, currtop, currbottom);
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
        rvec returnvector(ndatapermeas);
        returnvector(0) = result;
        return returnvector;
      }


    /*! Calculate the effect of the background layers for a single scalar measurement
     * the general structure of this function is very similar to CalcScalarBackground
     * however there are enough subtle differences in types and functions being called
     * that it does not make too much sense to try to combine the two, particularly because
     * the function is relatively short and simple
     * @param measindex The index of the measurement in the Model
     * @param xwidth The total width of the gridded domain in x-direction in m
     * @param ywidth The total width of the gridded domain in y-direction in m
     * @param zwidth The total width of the gridded domain in z-direction in m
     * @param Model The Gravity model
     * @param Sensitivities If the matrix passed here holds \f$ 1 \times nbg+ngrid \f$ or more elements, store sensitivity information in the right fields
     * @return A vector with a 9 components that contains the gravitational effect of the background on the different tensor elements
     */
    rvec CalcTensorBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
        const size_t ndatapermeas = 9;
        //make sure we have thicknesses and densities for all layers
        assert(Model.GetBackgroundDensities().size() == Model.GetBackgroundThicknesses().size());
        const size_t nbglayers = Model.GetBackgroundDensities().size();
        const double x_meas = Model.GetMeasPosX()[measindex];
        const double y_meas = Model.GetMeasPosY()[measindex];
        const double z_meas = Model.GetMeasPosZ()[measindex];
        rmat result(3, 3);
        rmat currvalue(3, 3);
        std::fill_n(result.data().begin(), ndatapermeas, 0.0);
        std::fill_n(currvalue.data().begin(), ndatapermeas, 0.0);
        double currtop = 0.0;
        double currbottom = 0.0;
        const size_t nmod = Model.GetDensities().num_elements();
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= nmod);
        // for all layers of the background
        for (size_t j = 0; j < nbglayers; ++j)
          {
            std::fill_n(currvalue.data().begin(), ndatapermeas, 0.0);
            const double currthick = Model.GetBackgroundThicknesses()[j];
            currbottom = currtop + currthick;
            currvalue(2, 2) = CalcUzzInfSheetTerm(z_meas, currtop, currbottom);
            if (currtop < zwidth && (currbottom <= zwidth)) // if the background layer complete coincides with the discretized area
              {
                currvalue -= CalcTensorBoxTerm(x_meas, y_meas, z_meas, 0.0,
                    0.0, currtop, xwidth, ywidth, currthick);
              }
            if (currtop < zwidth && currbottom > zwidth) //if some of the background coincides and some is below
              {
                currvalue -= CalcTensorBoxTerm(x_meas, y_meas, z_meas, 0.0,
                    0.0, currtop, xwidth, ywidth, (zwidth - currtop));
              }
            if (storesens)
              {
                for (size_t i = 0; i < ndatapermeas; ++i)
                  Sensitivities(i, nmod + j) = currvalue.data()[i];
              }
            result += currvalue * Model.GetBackgroundDensities()[j];
            currtop += currthick;
          }
        rvec resultvector(ndatapermeas);
        std::copy(result.data().begin(), result.data().end(),
            resultvector.begin());
        return resultvector;
      }
  }
