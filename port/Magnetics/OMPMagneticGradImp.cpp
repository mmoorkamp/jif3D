//============================================================================
// Name        : OMPMagneticGradImp.cpp
// Author      : 17 July 2014
// Version     : 
// Copyright   : 2014, zhanjie
//============================================================================

#include "OMPMagneticGradImp.h"
#include "../Gravity/BasicGravElements.h"

namespace jif3D
  {
    rvec OMPMagneticGradImp::CalcGridded(const size_t measindex,
        const ThreeDMagneticModel &Model, rmat &Sensitivities)
      {
        //we want to return vertical gradient of z component of magnetic field.
    	//firstly, calculate the z magnetic component of the two sensor with vertical distance of 1m.
        //then get the difference of the two sensor by subtracting z component of Top Sensor using that from Bottom Sensor.
        //And the objective function gradient is calculated using the difference of the two Sensitivity matrix of the two sensor.
        const double BxComp = cos(Inclination) * cos(Declination);
        const double ByComp = cos(Inclination) * sin(Declination);
        const double BzComp = sin(Inclination);

        const size_t xsize = Model.GetSusceptibilities().shape()[0];
        const size_t ysize = Model.GetSusceptibilities().shape()[1];
        const size_t zsize = Model.GetSusceptibilities().shape()[2];
        const double x_meas = Model.GetMeasPosX()[measindex];
        const double y_meas = Model.GetMeasPosY()[measindex];
        const double z_meas = Model.GetMeasPosZ()[measindex];
        const int nmod = xsize * ysize * zsize;
        const bool storesens = (Sensitivities.size1() >= ndatapermeas)
            && (Sensitivities.size2() >= size_t(nmod));

        rmat currvalueBot(3, 3);//currvalue matrix used to calculate z component of Bottom Sensor
        rmat currvalueTop(3, 3);//currvalue matrix used to calculate z component of Top Sensor

        //we cannot add up a user defined quantity in parallel
        //so break up the tensor into its component with different variables
        //and assign the results after the parallel loop
        double BzBot = 0.0;
        double BzTop = 0.0;

        //sum up the contributions of all prisms
#pragma omp parallel default(shared) private(currvalueBot, currvalueTop) reduction(+:BzBot,BzTop)
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
                // we reuse the calculation for the FTG matrix, as the equations are
                //identical
                //currvalue contains the geometric term times the gravitational constant
                currvalueBot = CalcTensorBoxTerm(x_meas, y_meas, z_meas, XCoord[xindex],
                    YCoord[yindex], ZCoord[zindex], XSizes[xindex], YSizes[yindex],
                    ZSizes[zindex]);
                //we have to multiply each element by the susceptibility
                const double Susceptibility =
                    Model.GetSusceptibilities()[xindex][yindex][zindex];
                // we have to convert the units of the FTG calculation to magnetics
                //using poisson's relation ship mu_0/4 pi = 1e-7
                const double factor = 1e-7 /jif3D::Grav_const;
                //the sensitivity element for the current grid cell and the z-component
                //of the magnetic field is the vector product of the FTG geometric terms
                //with the direction of the inducing magnetic field and the field strength
                const double BzSensBot = (currvalueBot(2, 0) * BxComp + currvalueBot(2, 1) * ByComp
                    + currvalueBot(2, 2) * BzComp) * FieldStrength * factor;
                //the field strength due to a single cell is the Sensitivity for the cell
                //times the Susceptibility
                BzBot += BzSensBot * Susceptibility;

                //currvalue of Top Sensor contains the geometric term times the gravitational constant
                currvalueTop = CalcTensorBoxTerm(x_meas, y_meas, z_meas+1, XCoord[xindex],
                    YCoord[yindex], ZCoord[zindex], XSizes[xindex], YSizes[yindex],
                    ZSizes[zindex]);
                //the sensitivity element for the current grid cell and the z-component
                //of the magnetic field is the vector product of the FTG geometric terms
                //with the direction of the inducing magnetic field and the field strength
                const double BzSensTop = (currvalueTop(2, 0) * BxComp + currvalueTop(2, 1) * ByComp
                    + currvalueTop(2, 2) * BzComp) * FieldStrength * factor;
                //the field strength due to a single cell is the Sensitivity for the cell
                //times the Susceptibility
                BzTop += BzSensTop * Susceptibility;

                //if we need the sensitivity values for each cell in the inversion
                //we store it in the appropriate component
                if (storesens)
                  {
                    Sensitivities(0, offset) = BzSensBot - BzSensTop;
                  }
              }
          } //end of parallel region

        rvec returnvalue(ndatapermeas);
        returnvalue(0) = BzBot - BzTop;
        return returnvalue;
      }

    OMPMagneticGradImp::OMPMagneticGradImp(double Inc, double Dec, double Fs) :
        Inclination(Inc), Declination(Dec), FieldStrength(Fs)
      {
        // TODO Auto-generated constructor stub

      }

    OMPMagneticGradImp::~OMPMagneticGradImp()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
