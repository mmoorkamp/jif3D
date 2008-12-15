//============================================================================
// Name        : GravityInterface.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================
#include "GravityInterface.h"
#include "ThreeDGravityModel.h"
#include "../Gravity/MinMemGravityCalculator.h"
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <iostream>

static boost::shared_ptr<jiba::ThreeDGravityModel> GravModel;
static boost::shared_ptr<jiba::ThreeDGravityCalculator> ScalarGravCalculator;
static boost::shared_ptr<jiba::ThreeDGravityCalculator> TensorGravCalculator;
typedef boost::function<jiba::ThreeDGravityModel::t3DModelDim & ()>
    tcoordinatefunc;

void AllocateModel(const int *storescalar, const int *storetensor)
  {
    bool wantscalar = (storescalar > 0);
    bool wanttensor = (storetensor > 0);
    GravModel = boost::shared_ptr<jiba::ThreeDGravityModel>(
        new jiba::ThreeDGravityModel(wantscalar, wanttensor));
    ScalarGravCalculator = jiba::CreateGravityCalculator<
        jiba::MinMemGravityCalculator>::MakeScalar();
    TensorGravCalculator = jiba::CreateGravityCalculator<
        jiba::MinMemGravityCalculator>::MakeTensor();
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
    CheckGridCoordinate(XSizes, nx, GravModel->GetXCellSizes(), boost::bind(
        &jiba::ThreeDGravityModel::SetXCellSizes, GravModel));
    CheckGridCoordinate(YSizes, ny, GravModel->GetYCellSizes(), boost::bind(
        &jiba::ThreeDGravityModel::SetYCellSizes, GravModel));
    CheckGridCoordinate(ZSizes, nz, GravModel->GetZCellSizes(), boost::bind(
        &jiba::ThreeDGravityModel::SetZCellSizes, GravModel));
    //allocates space for densities
    GravModel->SetDensities().resize(boost::extents[*nx][*ny][*nz]);
    //and copy as well
    std::copy(Densities, Densities + *nx * *ny * *nz,
        GravModel->SetDensities().origin());
    //we have to check whether the measurement points changed
    bool changed = CheckMeasPos(XMeasPos, nmeas, GravModel->GetMeasPosX())
        || CheckMeasPos(YMeasPos, nmeas, GravModel->GetMeasPosY())
        || CheckMeasPos(ZMeasPos, nmeas, GravModel->GetMeasPosZ());
    if (changed)
      {
        //we are using a static object, so we have to clear possible old measurement points
        GravModel->ClearMeasurementPoints();
        for (size_t i = 0; i < *nmeas; ++i)
          {
            GravModel->AddMeasurementPoint(XMeasPos[i], YMeasPos[i],
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
    jiba::rvec scalarmeas(ScalarGravCalculator->Calculate(*GravModel.get()));
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
    jiba::rvec tensormeas(TensorGravCalculator->Calculate(*GravModel.get()));

    std::copy(tensormeas.begin(), tensormeas.end(), GravTensor);
  }

