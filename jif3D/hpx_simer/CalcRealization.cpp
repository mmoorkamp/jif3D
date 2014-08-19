/*
 * CalcRealization.cpp
 *
 *  Created on: 19 Aug 2014
 *      Author: mm489
 */


#include "CalcRealization.h"
#include "../MT/X3DMTCalculator.h"
#include "../Global/VecMat.h"
jif3D::rvec CalcRealization(realinfo Info)
  {
    const size_t nx = Info.Model.GetXCoordinates().num_elements();
    const size_t ny = Info.Model.GetYCoordinates().num_elements();
    const size_t nz = Info.Model.GetZCoordinates().num_elements();

    for (size_t i = 0; i < nx; ++i)
      {
        for (size_t j = 0; j < ny; ++j)
          {
            Info.Model.SetConductivities()[i][j][0] = Info.bg_conductivity;
            for (size_t k = 1; k < nz; ++k)
              {
                if (drand48() < Info.phase1frac)
                  {
                    Info.Model.SetConductivities()[i][j][k] = Info.phase1cond;
                  }
                else
                  {
                    Info.Model.SetConductivities()[i][j][k] = Info.phase2cond;
                  }
              }
          }
      }

    jif3D::X3DMTCalculator Calculator(Info.tempdir, Info.x3dname);
    jif3D::rvec Impedances(Calculator.Calculate(Info.Model));
    return Impedances;
  }


#ifdef HAVEHPX

HPX_REGISTER_PLAIN_ACTION(Calc_action)
#endif
