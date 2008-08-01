//============================================================================
// Name        : GravityInterface.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================
#include "GravityInterface.h"
#include "ThreeDGravityModel.h"
#include <algorithm>
static jiba::ThreeDGravityModel GravForward;

void CalcScalarForward(const double *XSizes, const unsigned int nx,
    const double *YSizes, const unsigned int ny, const double *ZSizes,
    const unsigned int nz, const double *Densities, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int nmeas,
    double *GravAcceleration)
  {
    GravForward.ClearMeasurementPoints();
    GravForward.SetXCellSizes().resize(boost::extents[nx]);
    GravForward.SetYCellSizes().resize(boost::extents[ny]);
    GravForward.SetZCellSizes().resize(boost::extents[nz]);
    std::copy(XSizes, XSizes + nx, GravForward.SetXCellSizes().begin());
    std::copy(YSizes, YSizes + ny, GravForward.SetYCellSizes().begin());
    std::copy(ZSizes, ZSizes + nz, GravForward.SetZCellSizes().begin());
    GravForward.SetDensities().resize(boost::extents[nx][ny][nz]);
    std::copy(Densities, Densities + nx * ny * nz,
        GravForward.SetDensities().origin());
    for (size_t i = 0; i < nmeas; ++i)
      GravForward.AddMeasurementPoint(XMeasPos[i], YMeasPos[i], ZMeasPos[i]);
    jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas(
        GravForward.CalcGravity());
    std::copy(scalarmeas.begin(), scalarmeas.end(), GravAcceleration);
  }

void CalcTensorForward()
  {

  }

