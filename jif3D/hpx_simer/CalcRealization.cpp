/*
 * CalcRealization.cpp
 *
 *  Created on: 19 Aug 2014
 *      Author: mm489
 */


#include "CalcRealization.h"
#include "../MT/X3DMTCalculator.h"
#include "../Global/VecMat.h"



std::vector<double>  CalcRealization(realinfo Info)
  {
    const size_t nx = Info.Model.GetData().shape()[0];
    const size_t ny = Info.Model.GetData().shape()[1];
    const size_t nz = Info.Model.GetData().shape()[2];

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
    return std::vector<double>(Impedances.begin(),Impedances.end());
  }






#ifdef HAVEHPX

HPX_REGISTER_ACTION(Calc_action)
#endif
