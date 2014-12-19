//============================================================================
// Name        : CalcRawRealization.cpp
// Author      : 5 Dec 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#ifndef HPX_SIMER_CALCRAWREALIZATION_CPP_
#define HPX_SIMER_CALCRAWREALIZATION_CPP_


#include <boost/python.hpp>
#include "../Global/FatalException.h"
#include "../MT/X3DMTCalculator.h"
#include "CalcRawRealization.h"
#include "InvertBlock.h"




double CalcRawRealization(int nx, int ny, int nz, double delta, double topthick,
    double freq, double bgcond, boost::python::list& Conductivities, std::string tempdir,
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

    double cubethick = topthick + (nz-1) * delta;
    std::vector<double> bg_thicknesses(2, cubethick);
    std::vector<double> bg_conductivities(2, bgcond);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);
    const size_t nc = len(Conductivities);
    if (nx * ny * nz != nc)
      throw jif3D::FatalException(
          "Number of conducitvity values does not match mesh size !");
    for (size_t i = 0; i < nc; ++i)
      {
        *(Model.SetConductivities().origin() + i) = boost::python::extract<double>(Conductivities[i]);
      }

    jif3D::X3DMTCalculator Calculator(tempdir, x3dname);
    jif3D::rvec Impedances(Calculator.Calculate(Model));
    double res = InvertBlock(Model,Impedances, tempdir, x3dname);

    return res;

  }

BOOST_PYTHON_MODULE(simerreal_ext)
  {
    using namespace boost::python;
    def("CalcRawRealization", CalcRawRealization);
  }

#endif /* HPX_SIMER_CALCRAWREALIZATION_CPP_ */
