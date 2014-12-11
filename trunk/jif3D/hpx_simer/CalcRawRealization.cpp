//============================================================================
// Name        : CalcRawRealization.cpp
// Author      : 5 Dec 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#ifndef HPX_SIMER_CALCRAWREALIZATION_CPP_
#define HPX_SIMER_CALCRAWREALIZATION_CPP_


#include "CalcRawRealization.h"
#include <boost/python.hpp>
#include "../MT/X3DMTCalculator.h"



jif3D::rvec CalcRawRealization(int nx, int ny, int nz, double delta, double topthick,
    double freq, double bgcond, std::vector<double> Conductivities, std::string tempdir,
    std::string x3dname)
  {
    jif3D::X3DModel Model;
    Model.SetMeshSize(nx, ny, nz);
    Model.SetHorizontalCellSize(delta, delta, nx, ny);
    std::fill_n(Model.SetZCellSizes().begin(), nz, delta);
    if (topthick > 0.0)
      {
        Model.SetZCellSizes()[0] = topthick;
      }

    double posx, posy, posz = 0;
    posx = (delta * nx) / 2.0;
    posy = (delta * ny) / 2.0;
    Model.SetFrequencies().assign(1, freq);
    Model.AddMeasurementPoint(posx, posy, posz);

    std::vector<double> bg_thicknesses(Model.GetZCellSizes().size(), delta);
    std::vector<double> bg_conductivities(Model.GetZCellSizes().size(), bgcond);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);

    std::copy(Conductivities.begin(), Conductivities.end(),
        Model.SetConductivities().origin());
    jif3D::X3DMTCalculator Calculator(tempdir, x3dname);
    jif3D::rvec Impedances(Calculator.Calculate(Model));
    return Impedances;

  }


BOOST_PYTHON_MODULE(simerreal_ext)
{
    using namespace boost::python;
    def("CalcRawRealization", CalcRawRealization);
}


#endif /* HPX_SIMER_CALCRAWREALIZATION_CPP_ */
