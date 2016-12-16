//============================================================================
// Name        : ScalarThrustGravityImp.cpp
// Author      : 15 Dec 2016
// Version     : 
// Copyright   : 2016, mm489
//============================================================================

#include "ScalarThrustGravityImp.h"
#include "../GravMag/ThreeDGravMagCalculator.h"
#include "BasicGravElements.h"
#include "GravityBackground.h"


namespace jif3D
  {

    ScalarThrustGravityImp::ScalarThrustGravityImp()
      {
        // TODO Auto-generated constructor stub

      }

    ScalarThrustGravityImp::~ScalarThrustGravityImp()
      {
        // TODO Auto-generated destructor stub
      }

    //! Calculate the response of the background, currently this is done on the CPU
    rvec ScalarThrustGravityImp::CalcBackground(const size_t measindex,
        const double xwidth, const double ywidth, const double zwidth,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
        return CalcScalarBackground(measindex, xwidth, ywidth, zwidth, Model,
            Sensitivities);
      }
    //! Calculate the response of the gridded part, this is done on the GPU with CUDA
    rvec ScalarThrustGravityImp::CalcGridded(const size_t measindex,
        const ThreeDGravityModel &Model, rmat &Sensitivities)
      {
      }
    rvec ScalarThrustGravityImp::Calculate(const ThreeDGravityModel &Model,
        ThreeDGravMagCalculator<ThreeDGravityModel> &Calculator)
      {

        const unsigned int nx = Model.GetDensities().shape()[0];
        const unsigned int ny = Model.GetDensities().shape()[1];
        const unsigned int nz = Model.GetDensities().shape()[2];
        const unsigned int nmod = nx * ny * nz;
        std::copy(Model.GetXCoordinates().data().origin(),
            Model.GetXCoordinates().data().origin() + nx, d_xcoord.begin());
        std::copy(Model.GetYCoordinates().data().origin(),
            Model.GetYCoordinates().data().origin() + ny, d_ycoord.begin());
        std::copy(Model.GetZCoordinates().data().origin(),
            Model.GetZCoordinates().data().origin() + nz, d_zcoord.begin());
        std::copy(Model.GetXCellSizes().data().origin(),
            Model.GetXCellSizes().data().origin() + nx, d_xsize.begin());
        std::copy(Model.GetYCellSizes().data().origin(),
            Model.GetYCellSizes().data().origin() + ny, d_ysize.begin());
        std::copy(Model.GetZCellSizes().data().origin(),
            Model.GetZCellSizes().data().origin() + nz, d_zsize.begin());
        d_result.resize(nmod);

        // call the base class that coordinates the calculation of gridded and background parts
        rvec result(ThreeDGravMagImplementation::Calculate(Model, Calculator));
        // free memory

        return result;
      }
  }

} /* namespace jif3D */
