/*
 * OMPMagnetizationImp.cpp
 *
 *  Created on: May 19, 2022
 *      Author: max
 */

#include "OMPMagnetizationImp.h"
#include "../Gravity/BasicGravElements.h"
#include "../Gravity/GravityBackground.h"
#include <boost/math/constants/constants.hpp>

namespace jif3D
  {

    rvec OMPMagnetizationImp::CalcBackground(const size_t measindex, const double xwidth,
        const double ywidth, const double zwidth, const ThreeDMagnetizationModel &Model,
        const ThreeComponentMagneticData &Data, rmat &Sensitivities)
      {
        rvec result(3, 0.0);
        return result;
      }

    //! Calculate the response of the gridded part
    rvec OMPMagnetizationImp::CalcGridded(const size_t measindex,
        const ThreeDMagnetizationModel &Model, const ThreeComponentMagneticData &Data,
        rmat &Sensitivities)
      {

        const size_t xsize = Model.GetMagnetization_X().shape()[0];
        const size_t ysize = Model.GetMagnetization_X().shape()[1];
        const size_t zsize = Model.GetMagnetization_X().shape()[2];
        const double x_meas = Data.GetMeasPosX()[measindex];
        const double y_meas = Data.GetMeasPosY()[measindex];
        const double z_meas = Data.GetMeasPosZ()[measindex];
        const int nmod = xsize * ysize * zsize;
        const size_t nsens = 3 * nmod;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= size_t(nsens));

        //we cannot add up a user defined quantity in parallel
        //so break up the tensor into its component with different variables
        //and assign the results after the parallel loop
        double Bx = 0.0, By = 0.0, Bz = 0.0;

        //sum up the contributions of all prisms
#pragma omp parallel default(shared) reduction(+:Bx,By,Bz)
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

                // we have to convert the units of the FTG calculation to magnetics
                //using poisson's relationship 1/(4 pi gamma)
                const double factor = 1.0
                    / (4 * boost::math::constants::pi<double>() * jif3D::Grav_const);
                // we reuse the calculation for the FTG matrix, as the equations are
                //identical
                //currvalue contains the geometric term times the gravitational constant
                rmat currvalue = factor
                    * CalcTensorBoxTerm(x_meas, y_meas, z_meas,
                        Model.GetXCoordinates()[xindex], Model.GetYCoordinates()[yindex],
                        Model.GetZCoordinates()[zindex], Model.GetXCellSizes()[xindex],
                        Model.GetYCellSizes()[yindex], Model.GetZCellSizes()[zindex]);
                //we have to multiply each element by the susceptibility
                const double magx = Model.GetMagnetization_X()[xindex][yindex][zindex];
                const double magy = Model.GetMagnetization_Y()[xindex][yindex][zindex];
                const double magz = Model.GetMagnetization_Z()[xindex][yindex][zindex];

                //the sensitivity element for the current grid cell and the x-component
                //of the magnetic field is the vector product of the FTG geometric terms
                //with the direction of the inducing magnetic field and the field strength
                const double BxSens = (currvalue(0, 0) * magx + currvalue(0, 1) * magy
                    + currvalue(0, 2) * magz);
                //and the same for By and Bz
                const double BySens = (currvalue(1, 0) * magx + currvalue(1, 1) * magy
                    + currvalue(1, 2) * magz);
                const double BzSens = (currvalue(2, 0) * magx + currvalue(2, 1) * magy
                    + currvalue(2, 2) * magz);
                //the field strength due to a single cell is the Sensitivity for the cell
                //times the Susceptibility
                Bx += BxSens;
                By += BySens;
                Bz += BzSens;
                //if we need the sensitivity values for each cell in the inversion
                //we store it in the appropriate component
                if (storesens)
                  {
                    Sensitivities(0, offset) = currvalue(0, 0);
                    Sensitivities(0, offset + nmod) = currvalue(0, 1);
                    Sensitivities(0, offset + 2 * nmod) = currvalue(0, 2);
                    Sensitivities(1, offset) = currvalue(1, 0);
                    Sensitivities(1, offset + nmod) = currvalue(1, 1);
                    Sensitivities(1, offset + 2 * nmod) = currvalue(1, 2);
                    Sensitivities(2, offset) = currvalue(2, 0);
                    Sensitivities(2, offset + nmod) = currvalue(2, 1);
                    Sensitivities(2, offset + 2 * nmod) = currvalue(2, 2);
                  }
              }
          } //end of parallel region

        rvec returnvalue(ndatapermeas);
        returnvalue(0) = Bx;
        returnvalue(1) = By;
        returnvalue(2) = Bz;

        return returnvalue;
      }

    OMPMagnetizationImp::OMPMagnetizationImp()
      {
        // TODO Auto-generated constructor stub

      }

    OMPMagnetizationImp::~OMPMagnetizationImp()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
