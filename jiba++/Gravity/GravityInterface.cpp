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

//This is a helper function that summarizes the steps
//that are identical for scalar and tensorial calculations
void SetupModel(const double *XSizes, const unsigned int *nx,
    const double *YSizes, const unsigned int *ny, const double *ZSizes,
    const unsigned int *nz, const double *Densities, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int *nmeas)
  {
    //we are using a static object, so we have to clear possible old measurement points
    GravForward.ClearMeasurementPoints();
    //make enough space for the coordinates
    GravForward.SetXCellSizes().resize(boost::extents[*nx]);
    GravForward.SetYCellSizes().resize(boost::extents[*ny]);
    GravForward.SetZCellSizes().resize(boost::extents[*nz]);
    //and copy the values into the model object
    std::copy(XSizes, XSizes + *nx, GravForward.SetXCellSizes().begin());
    std::copy(YSizes, YSizes + *ny, GravForward.SetYCellSizes().begin());
    std::copy(ZSizes, ZSizes + *nz, GravForward.SetZCellSizes().begin());
    //allocates space for densities
    GravForward.SetDensities().resize(boost::extents[*nx][*ny][*nz]);
    //and copy as well
    std::copy(Densities, Densities + *nx * *ny * *nz,
        GravForward.SetDensities().origin());
    for (size_t i = 0; i < *nmeas; ++i)
      GravForward.AddMeasurementPoint(XMeasPos[i], YMeasPos[i], ZMeasPos[i]);
  }

void CalcScalarForward(const double *XSizes, const unsigned int *nx,
    const double *YSizes, const unsigned int *ny, const double *ZSizes,
    const unsigned int *nz, const double *Densities, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int *nmeas,
    double *GravAcceleration)
  {
    SetupModel(XSizes,nx,YSizes,ny,ZSizes,nz,Densities,XMeasPos,YMeasPos,ZMeasPos,nmeas);
    jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas(
        GravForward.CalcGravity());
    std::copy(scalarmeas.begin(), scalarmeas.end(), GravAcceleration);
  }

void CalcTensorForward(const double *XSizes, const unsigned int *nx,
    const double *YSizes, const unsigned int *ny, const double *ZSizes,
    const unsigned int *nz, const double *Densities, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int *nmeas,
    double *GravTensor)
  {
    SetupModel(XSizes,nx,YSizes,ny,ZSizes,nz,Densities,XMeasPos,YMeasPos,ZMeasPos,nmeas);
        jiba::ThreeDGravityModel::tTensorMeasVec tensormeas(
            GravForward.CalcTensorGravity());
    for (size_t i = 0; i < *nmeas; ++i)
      std::copy(tensormeas.at(i).data().begin(), tensormeas.at(i).data().end(), GravTensor + i * 9);
  }

