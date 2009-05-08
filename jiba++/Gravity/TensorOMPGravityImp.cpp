//============================================================================
// Name        : TensorOMPGravityImp.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "TensorOMPGravityImp.h"
#include "BasicGravElements.h"
namespace jiba
  {

    TensorOMPGravityImp::TensorOMPGravityImp()
      {


      }

    TensorOMPGravityImp::~TensorOMPGravityImp()
      {

      }

    /*!  Calculate the contribution of a layered background to a tensor gravity measurement.
     * @param measindex The index of the measurement
     * @param xwidth The total width of the discretized model area in x-direction in m
     * @param ywidth The total width of the discretized model area in y-direction in m
     * @param zwidth The total width of the discretized model area in z-direction in m
     * @param Model The gravity model
     * @param Sensitivities The \f$ 9 \times m\f$ matrix of sensitivities for the current measurement
     * @return The gravitational tensor due to the background
     */
    rvec TensorOMPGravityImp::CalcBackground(const size_t measindex,
        const double xwidth, const double ywidth, const double zwidth,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
        //make sure we have thicknesses and densities for all layers
        assert(Model.GetBackgroundDensities().size() == Model.GetBackgroundThicknesses().size());
        const size_t nbglayers = Model.GetBackgroundDensities().size();
        const double x_meas = Model.GetMeasPosX()[measindex];
        const double y_meas = Model.GetMeasPosY()[measindex];
        const double z_meas = Model.GetMeasPosZ()[measindex];
        GravimetryMatrix result(3, 3);
        GravimetryMatrix currvalue(3, 3);
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
    /*! Calculate the FTG response of the gridded domain.
     * @param measindex The index of the measurement
     * @param Model The gravity model
     * @param Sensitivities The \f$ 9 \times m\f$ matrix of sensitivities for the current measurement
     * @return A 9 component vector with the FTG matrix components
     */
    rvec TensorOMPGravityImp::CalcGridded(const size_t measindex,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
        const size_t xsize = Model.GetDensities().shape()[0];
        const size_t ysize = Model.GetDensities().shape()[1];
        const size_t zsize = Model.GetDensities().shape()[2];
        const double x_meas = Model.GetMeasPosX()[measindex];
        const double y_meas = Model.GetMeasPosY()[measindex];
        const double z_meas = Model.GetMeasPosZ()[measindex];
        const int nmod = xsize * ysize * zsize;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= size_t(nmod));
        GravimetryMatrix currvalue(3, 3);

        //we cannot add up a user defined quantity in parallel
        //so break up the tensor into its component with different variables
        //and assign the results after the parallel loop
        double U0 = 0.0, U1 = 0.0, U2 = 0.0, U3 = 0.0, U4 = 0.0, U5 = 0.0, U6 =
            0.0, U7 = 0.0, U8 = 0.0;
        //sum up the contributions of all prisms
#pragma omp parallel default(shared) private(currvalue) reduction(+:U0,U1,U2,U3,U4,U5,U6,U7,U8)
          {
            //instead of nested loops over each dimension, we have one big
            //loop over all elements, this allows for a nearly infinite number
            //of parallel processors
#pragma omp for
            for (int offset = 0; offset < nmod; ++offset)
              {
                int xindex, yindex, zindex;
                //we still need the indices for each dimension
                //so we have to convert our loop variable
                Model.OffsetToIndex(offset, xindex, yindex, zindex);
                //currvalue contains only the geometric term
                currvalue = CalcTensorBoxTerm(x_meas, y_meas, z_meas,
                    XCoord[xindex], YCoord[yindex], ZCoord[zindex],
                    XSizes[xindex], YSizes[yindex], ZSizes[zindex]);
                //to we have to multiply each element by the density
                const double Density =
                    Model.GetDensities()[xindex][yindex][zindex];
                U0 += currvalue(0, 0) * Density;
                U1 += currvalue(0, 1) * Density;
                U2 += currvalue(0, 2) * Density;
                U3 += currvalue(1, 0) * Density;
                U4 += currvalue(1, 1) * Density;
                U5 += currvalue(1, 2) * Density;
                U6 += currvalue(2, 0) * Density;
                U7 += currvalue(2, 1) * Density;
                U8 += currvalue(2, 2) * Density;
                if (storesens)
                  {
                    for (size_t i = 0; i < ndatapermeas; ++i)
                      Sensitivities(i, offset) = currvalue.data()[i];
                  }
              }
          }//end of parallel region

        rvec returnvalue(ndatapermeas);
        returnvalue(0) = U0;
        returnvalue(1) = U1;
        returnvalue(2) = U2;
        returnvalue(3) = U3;
        returnvalue(4) = U4;
        returnvalue(5) = U5;
        returnvalue(6) = U6;
        returnvalue(7) = U7;
        returnvalue(8) = U8;
        return returnvalue;

      }

  }
