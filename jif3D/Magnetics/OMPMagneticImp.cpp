//============================================================================
// Name        : OMPMagneticImp.cpp
// Author      : 6 Nov 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include "OMPMagneticImp.h"

namespace jif3D
  {

    rvec OMPMagneticImp::CalcGridded(const size_t measindex,
        const ThreeDMagneticModel &Model, rmat &Sensitivities)
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
        // we reuse the calculation for the FTG matrix, as the equations are
        //identical
        GravimetryMatrix currvalue(3, 3);

        //we cannot add up a user defined quantity in parallel
        //so break up the tensor into its component with different variables
        //and assign the results after the parallel loop
        double U0 = 0.0, U1 = 0.0, U2 = 0.0, U3 = 0.0, U4 = 0.0, U5 = 0.0, U6 = 0.0, U7 =
            0.0, U8 = 0.0;
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
                currvalue = CalcTensorBoxTerm(x_meas, y_meas, z_meas, XCoord[xindex],
                    YCoord[yindex], ZCoord[zindex], XSizes[xindex], YSizes[yindex],
                    ZSizes[zindex]);
                //to we have to multiply each element by the density
                const double Susceptibility = Model.GetSusceptibilities()[xindex][yindex][zindex];
                U0 += currvalue(0, 0) * Susceptibility;
                U1 += currvalue(0, 1) * Susceptibility;
                U2 += currvalue(0, 2) * Susceptibility;
                U3 += currvalue(1, 0) * Susceptibility;
                U4 += currvalue(1, 1) * Susceptibility;
                U5 += currvalue(1, 2) * Susceptibility;
                U6 += currvalue(2, 0) * Susceptibility;
                U7 += currvalue(2, 1) * Susceptibility;
                U8 += currvalue(2, 2) * Susceptibility;
                if (storesens)
                  {
                    for (size_t i = 0; i < ndatapermeas; ++i)
                      Sensitivities(i, offset) = currvalue.data()[i];
                  }
              }
          } //end of parallel region

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

OMPMagneticImp::OMPMagneticImp(double Inc) :
    Inclination(Inc), Declination(Dec)
  {
    // TODO Auto-generated constructor stub

  }

OMPMagneticImp::~OMPMagneticImp()
  {
    // TODO Auto-generated destructor stub
  }

} /* namespace jif3D */
