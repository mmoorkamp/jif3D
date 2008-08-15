//============================================================================
// Name        : GravityInterface.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================
#include "GravityInterface.h"
#include "ThreeDGravityModel.h"
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <iostream>

static boost::shared_ptr<jiba::ThreeDGravityModel> GravForward;

typedef boost::function<jiba::ThreeDGravityModel::t3DModelDim & ()> tcoordinatefunc;

void AllocateModel(const int *storescalar, const int *storetensor)
  {
    bool wantscalar = (storescalar > 0);
    bool wanttensor = (storetensor > 0);
    GravForward = boost::shared_ptr<jiba::ThreeDGravityModel>(
        new jiba::ThreeDGravityModel(wantscalar, wanttensor));
  }

//Check whether a mesh coordinate has changed and copy values if necessary
void CheckGridCoordinate(const double *Sizes, const unsigned int *n,
    const jiba::ThreeDGravityModel::t3DModelDim &Dim,
    tcoordinatefunc ChangeFunc)
  {
    //first we have to check whether the vectors have the same length
    bool samesize = (*n == Dim.size());
    //if they do, check whether the elements are equal
    if (samesize)
      {
        samesize = std::equal(Sizes, Sizes + *n, Dim.begin());
      }
    //if something has been changed reallocate and copy
    if (!samesize)
      {
        ChangeFunc().resize(boost::extents[*n]);
        std::copy(Sizes, Sizes + *n, ChangeFunc().begin());
      }
  }

bool CheckMeasPos(const double *Sizes, const unsigned int *n,
    const jiba::ThreeDGravityModel::tMeasPosVec &Meas)
  {
    bool same = (*n == Meas.size());
    if (same)
      {
        same = std::equal(Sizes, Sizes + *n, Meas.begin());
      }
    return !same;
  }

//This is a helper function that summarizes the steps
//that are identical for scalar and tensorial calculations
void SetupModel(const double *XSizes, const unsigned int *nx,
    const double *YSizes, const unsigned int *ny, const double *ZSizes,
    const unsigned int *nz, const double *Densities, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int *nmeas)
  {

    //we check whether the geometry has changed to take advantage of possible accelerations
    //through sensitivity matrix storage
    CheckGridCoordinate(XSizes, nx, GravForward->GetXCellSizes(), boost::bind(
        &jiba::ThreeDGravityModel::SetXCellSizes, GravForward));
    CheckGridCoordinate(YSizes, ny, GravForward->GetYCellSizes(), boost::bind(
        &jiba::ThreeDGravityModel::SetYCellSizes, GravForward));
    CheckGridCoordinate(ZSizes, nz, GravForward->GetZCellSizes(), boost::bind(
        &jiba::ThreeDGravityModel::SetZCellSizes, GravForward));
    //allocates space for densities
    GravForward->SetDensities().resize(boost::extents[*nx][*ny][*nz]);
    //and copy as well
    std::copy(Densities, Densities + *nx * *ny * *nz,
        GravForward->SetDensities().origin());
    //we have to check whether the measurement points changed
    bool changed = CheckMeasPos(XMeasPos, nmeas, GravForward->GetMeasPosX())
        || CheckMeasPos(YMeasPos, nmeas, GravForward->GetMeasPosY())
        || CheckMeasPos(ZMeasPos, nmeas, GravForward->GetMeasPosZ());
    if (changed)
      {
        //we are using a static object, so we have to clear possible old measurement points
        GravForward->ClearMeasurementPoints();
        for (size_t i = 0; i < *nmeas; ++i)
          {
            GravForward->AddMeasurementPoint(XMeasPos[i], YMeasPos[i],
                ZMeasPos[i]);
          }
      }
  }

void CalcScalarForward(const double *XSizes, const unsigned int *nx,
    const double *YSizes, const unsigned int *ny, const double *ZSizes,
    const unsigned int *nz, const double *Densities, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int *nmeas,
    double *GravAcceleration)
  {
    SetupModel(XSizes, nx, YSizes, ny, ZSizes, nz, Densities, XMeasPos,
        YMeasPos, ZMeasPos, nmeas);
    jiba::ThreeDGravityModel::tScalarMeasVec scalarmeas(
        GravForward->CalcGravity());
    std::copy(scalarmeas.begin(), scalarmeas.end(), GravAcceleration);
  }

void CalcTensorForward(const double *XSizes, const unsigned int *nx,
    const double *YSizes, const unsigned int *ny, const double *ZSizes,
    const unsigned int *nz, const double *Densities, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int *nmeas,
    double *GravTensor)
  {
    SetupModel(XSizes, nx, YSizes, ny, ZSizes, nz, Densities, XMeasPos,
        YMeasPos, ZMeasPos, nmeas);
    jiba::ThreeDGravityModel::tTensorMeasVec tensormeas(
        GravForward->CalcTensorGravity());
    for (size_t i = 0; i < *nmeas; ++i)
      std::copy(tensormeas.at(i).data().begin(), tensormeas.at(i).data().end(),
          GravTensor + i * 9);
  }

