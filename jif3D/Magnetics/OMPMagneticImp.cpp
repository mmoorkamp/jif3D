//============================================================================
// Name        : OMPMagneticImp.cpp
// Author      : 6 Nov 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include "OMPMagneticImp.h"
#include "../Gravity/BasicGravElements.h"
#include "../Gravity/GravityBackground.h"
#include <boost/math/constants/constants.hpp>

namespace jif3D
  {


  /*! Calculate the effect of the background layers for a single magnetic measurement.
   * tThe general structure of this function is very similar to CalcTensorBackground
   * in GravityBackground.cpp however there are enough subtle differences in types and functions being called
   * that it does not make too much sense to try to combine the two, particularly because
   * the function is relatively short and simple
   * @param measindex The index of the measurement in the Model
   * @param xwidth The total width of the gridded domain in x-direction in m
   * @param ywidth The total width of the gridded domain in y-direction in m
   * @param zwidth The total width of the gridded domain in z-direction in m
   * @param Model The Susceptibility model
   * @param Sensitivities If the matrix passed here holds \f$ 1 \times nbg+ngrid \f$ or more elements, store sensitivity information in the right fields
   * @return A vector with a 3 components that contains the magnetic field
   */
  rvec OMPMagneticImp::CalcMagneticBackground(const size_t measindex, const double xwidth,
      const double ywidth, const double zwidth, const ThreeDSusceptibilityModel &Model,
      const TotalFieldMagneticData &Data, rmat &Sensitivities)
    {
      //make sure we have thicknesses and susceptibilities for all layers
      assert(
          Model.GetBackgroundSusceptibilities().size()
              == Model.GetBackgroundThicknesses().size());
      const size_t nbglayers = Model.GetBackgroundSusceptibilities().size();
      const double x_meas = Data.GetMeasPosX()[measindex];
      const double y_meas = Data.GetMeasPosY()[measindex];
      const double z_meas = Data.GetMeasPosZ()[measindex];
      //calculate the directional angles
      const double BxComp = cos(Inclination) * cos(Declination);
      const double ByComp = cos(Inclination) * sin(Declination);
      const double BzComp = sin(Inclination);
      //the conversion factor to go from tensor gravity sensitivities to magnetics
      const double factor = 1.0 / ( 4 * boost::math::constants::pi<double>()   * jif3D::Grav_const);

      rvec result(3, 0.0);
      rmat currvalue(3, 3, 0.0);
      double currtop = Model.GetZOrigin();
      double currbottom = currtop;
      const size_t nmod = Model.GetSusceptibilities().num_elements();
      const bool storesens = (Sensitivities.size1() >= ndatapermeas)
          && (Sensitivities.size2() >= nmod);
      // for all layers of the background
      for (size_t j = 0; j < nbglayers; ++j)
        {
          //reset current sensitivity values
          std::fill(currvalue.data().begin(), currvalue.data().end(), 0.0);
          const double currthick = Model.GetBackgroundThicknesses()[j];
          currbottom = currtop + currthick;
          //if we are inside the layer we might have a contribution to the zz-component
          currvalue(2, 2) = CalcUzzInfSheetTerm(z_meas, currtop, currbottom);
          //most of the times and possible contribution comes from the interface between the grid and the background
          if (currtop < zwidth && (currbottom <= zwidth)) // if the background layer complete coincides with the discretized area
            {
              currvalue -= CalcTensorBoxTerm(x_meas, y_meas, z_meas,  Model.GetXOrigin(), Model.GetYOrigin(), currtop,
                  xwidth, ywidth, currthick);
            }
          if (currtop < zwidth && currbottom > zwidth) //if some of the background coincides and some is below
            {
              currvalue -= CalcTensorBoxTerm(x_meas, y_meas, z_meas,  Model.GetXOrigin(), Model.GetYOrigin(), currtop,
                  xwidth, ywidth, (zwidth - currtop));
            }
          //project the tensor sensitivities to magnetic field components
          const double BxSens = (currvalue(0, 0) * BxComp + currvalue(0, 1) * ByComp
              + currvalue(0, 2) * BzComp) * FieldStrength * factor;
          const double BySens = (currvalue(1, 0) * BxComp + currvalue(1, 1) * ByComp
              + currvalue(1, 2) * BzComp) * FieldStrength * factor;
          const double BzSens = (currvalue(2, 0) * BxComp + currvalue(2, 1) * ByComp
              + currvalue(2, 2) * BzComp) * FieldStrength * factor;
          // store in the sensitivity matrix if necessary
          if (storesens)
            {
              Sensitivities(0, nmod + j) = BxSens;
              Sensitivities(1, nmod + j) = BySens;
              Sensitivities(2, nmod + j) = BzSens;

            }
          //calculate resultant field components
          const double bgsus = Model.GetBackgroundSusceptibilities()[j];
          result(0) += BxSens * bgsus;
          result(1) += BySens * bgsus;
          result(2) += BzSens * bgsus;
          currtop += currthick;
        }

      return result;
    }

 rvec OMPMagneticImp::CalcBackground(const size_t measindex, const double xwidth,
      const double ywidth, const double zwidth, const ThreeDSusceptibilityModel &Model,
      const TotalFieldMagneticData &Data,
      rmat &Sensitivities)
    {

      return CalcMagneticBackground(measindex, xwidth, ywidth, zwidth, Model, Data,
          Sensitivities);
    }

    rvec OMPMagneticImp::CalcGridded(const size_t measindex,
        const ThreeDSusceptibilityModel &Model, const TotalFieldMagneticData &Data, rmat &Sensitivities)
      {
        const double BxComp = cos(Inclination) * cos(Declination);
        const double ByComp = cos(Inclination) * sin(Declination);
        const double BzComp = sin(Inclination);

        const size_t xsize = Model.GetSusceptibilities().shape()[0];
        const size_t ysize = Model.GetSusceptibilities().shape()[1];
        const size_t zsize = Model.GetSusceptibilities().shape()[2];
        const double x_meas = Data.GetMeasPosX()[measindex];
        const double y_meas = Data.GetMeasPosY()[measindex];
        const double z_meas = Data.GetMeasPosZ()[measindex];
        const int nmod = xsize * ysize * zsize;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= size_t(nmod));

        rmat currvalue(3, 3);

        //we cannot add up a user defined quantity in parallel
        //so break up the tensor into its component with different variables
        //and assign the results after the parallel loop
        double Bx = 0.0, By = 0.0, Bz = 0.0;

        //sum up the contributions of all prisms
#pragma omp parallel default(shared) private(currvalue) reduction(+:Bx,By,Bz)
          {
            //instead of nested loops over each dimension, we have one big
            //loop over all elements, this allows for a large number
            //of parallel processors
#pragma omp for
            for (int offset = 0; offset < nmod; ++offset)
              {
                int xindex, yindex, zindex;
                //we still need the indices for each dimension
                //so we have to convert our loop variable
                Model.OffsetToIndex(offset, xindex, yindex, zindex);
                // we reuse the calculation for the FTG matrix, as the equations are
                //identical
                //currvalue contains the geometric term times the gravitational constant
                currvalue = CalcTensorBoxTerm(x_meas, y_meas, z_meas,
                    Model.GetXCoordinates()[xindex], Model.GetYCoordinates()[yindex],
                    Model.GetZCoordinates()[zindex], Model.GetXCellSizes()[xindex],
                    Model.GetYCellSizes()[yindex], Model.GetZCellSizes()[zindex]);
                //we have to multiply each element by the susceptibility
                const double Susceptibility =
                    Model.GetSusceptibilities()[xindex][yindex][zindex];
                // we have to convert the units of the FTG calculation to magnetics
                //using poisson's relationship 1/(4 pi gamma)
                const double factor = 1.0 / ( 4 * boost::math::constants::pi<double>()   * jif3D::Grav_const);
                //the sensitivity element for the current grid cell and the x-component
                //of the magnetic field is the vector product of the FTG geometric terms
                //with the direction of the inducing magnetic field and the field strength
                const double BxSens = (currvalue(0, 0) * BxComp + currvalue(0, 1) * ByComp
                    + currvalue(0, 2) * BzComp) * FieldStrength * factor;
                //and the same for By and Bz
                const double BySens = (currvalue(1, 0) * BxComp + currvalue(1, 1) * ByComp
                    + currvalue(1, 2) * BzComp) * FieldStrength * factor;
                const double BzSens = (currvalue(2, 0) * BxComp + currvalue(2, 1) * ByComp
                    + currvalue(2, 2) * BzComp) * FieldStrength * factor;
                //the field strength due to a single cell is the Sensitivity for the cell
                //times the Susceptibility
                Bx += BxSens * Susceptibility;
                By += BySens * Susceptibility;
                Bz += BzSens * Susceptibility;
                //if we need the sensitivity values for each cell in the inversion
                //we store it in the appropriate component
                if (storesens)
                  {
                    Sensitivities(0, offset) = BxSens;
                    Sensitivities(1, offset) = BySens;
                    Sensitivities(2, offset) = BzSens;
                  }
              }
          } //end of parallel region

        rvec returnvalue(ndatapermeas);
        returnvalue(0) = Bx;
        returnvalue(1) = By;
        returnvalue(2) = Bz;

        return returnvalue;

      }

    OMPMagneticImp::OMPMagneticImp(double Inc, double Dec, double Fs) :
        Inclination(Inc), Declination(Dec), FieldStrength(Fs)
      {
        // TODO Auto-generated constructor stub

      }

    OMPMagneticImp::~OMPMagneticImp()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
