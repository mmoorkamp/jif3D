//============================================================================
// Name        : MTTests.h
// Author      : 10 Mar 2014
// Version     :
// Copyright   : 2014, mm489
//============================================================================

#ifndef MTTESTS_H_
#define MTTESTS_H_

#include "X3DModel.h"
#include "MTData.h"
#include <vector>
#include <algorithm>

inline void MakeMTModel(jif3D::X3DModel &Model, jif3D::MTData &Data)
  {
    const size_t xsize = 4;
    const size_t ysize = 3;
    const size_t zsize = 2;
    const size_t nbglayers = 2;
    const size_t nmod = xsize * ysize * zsize;

    Model.SetMeshSize(xsize, ysize, zsize);

    std::vector<double> bg_thicknesses(nbglayers), bg_conductivities(nbglayers);

    const double deltax = 100.0;
    const double deltay = 100.0;
    const double deltaz = 100.0;
    Model.SetHorizontalCellSize(deltax, deltay, xsize, ysize);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize, deltaz);
    ZCS[1] = deltaz * 2.0;
    Model.SetZCellSizes(ZCS);
    std::fill_n(Model.SetConductivities().origin(), nmod, 0.02);

    Model.SetConductivities()[0][0][0] = 0.025;
    std::fill_n(bg_conductivities.begin(), nbglayers, 0.02);
    for (size_t i = 0; i < nbglayers; ++i)
      {
        bg_conductivities[i] *= 1.05 + i * 0.1;
      }
    bg_thicknesses = ZCS;

    Model.SetBackgroundConductivities(bg_conductivities);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Data.SetFrequencies(
      { 1.0, 2.0, 5.0, 10.0 });
    for (size_t i = 1; i < xsize - 1; ++i)
      {
        for (size_t j = 1; j < ysize - 1; ++j)
          {
            double currx = Model.GetXCoordinates()[i] + deltax / 3.0;
            double curry = Model.GetYCoordinates()[j] + deltay / 4.0;
            double currz = (j - 1) * deltaz;
            Data.AddMeasurementPoint(currx, curry, currz);
          }
      }
    Data.CompleteObject();
  }

#endif
