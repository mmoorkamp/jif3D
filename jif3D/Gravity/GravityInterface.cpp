//============================================================================
// Name        : GravityInterface.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#include "GravityInterface.h"
#include "ThreeDGravityModel.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "ThreeDGravityFactory.h"
#include <algorithm>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <iostream>

//these static declarations make the code not thread-safe
//however the forward calculations are parallelized and can be used
//we just cannot have parallel calls through this interface

static boost::shared_ptr<jif3D::ThreeDGravityModel> GravModel;
static boost::shared_ptr<jif3D::ThreeDGravMagCalculator<jif3D::ThreeDGravityModel> > ScalarGravCalculator;
static boost::shared_ptr<jif3D::ThreeDGravMagCalculator<jif3D::ThreeDGravityModel> > TensorGravCalculator;
typedef boost::function0<jif3D::ThreeDGravityModel::t3DModelDim &>
    tcoordinatefunc;

void AllocateModel(const int *storescalar, const int *storetensor)
  {
    //c doesn't have bool so we interpret every positive value as true
    bool cachescalar = (storescalar > 0);
    bool cachetensor = (storetensor > 0);
    // we create a new model object that we can use for the calculations
    GravModel = boost::make_shared<jif3D::ThreeDGravityModel>();
    //depending on whether we want to save the sensitivities
    //or not we allocate different forward calculation objects
    //for FTG and scalar
    if (cachescalar)
      {
        ScalarGravCalculator = jif3D::CreateGravityCalculator<
            jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeScalar();
      }
    else
      {
        ScalarGravCalculator = jif3D::CreateGravityCalculator<
            jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeScalar();
      }
    if (cachetensor)
      {
        TensorGravCalculator = jif3D::CreateGravityCalculator<
            jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeTensor();
      }
    else
      {
        TensorGravCalculator = jif3D::CreateGravityCalculator<
            jif3D::MinMemGravMagCalculator<jif3D::ThreeDGravityModel> >::MakeTensor();
      }

  }

//Check whether a mesh coordinate has changed and copy values if necessary
void CheckGridCoordinate(const double *Sizes, const unsigned int *n,
    const jif3D::ThreeDGravityModel::t3DModelDim &Dim,
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
        ChangeFunc().resize(*n);
        std::copy(Sizes, Sizes + *n, ChangeFunc().begin());
      }
  }

//check whether a measurement position has changed
bool CheckMeasPos(const double *Sizes, const unsigned int *n,
    const jif3D::ThreeDGravityModel::tMeasPosVec &Meas)
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
        &jif3D::ThreeDGravityModel::SetXCellSizes, GravModel));
    CheckGridCoordinate(YSizes, ny, GravModel->GetYCellSizes(), boost::bind(
        &jif3D::ThreeDGravityModel::SetYCellSizes, GravModel));
    CheckGridCoordinate(ZSizes, nz, GravModel->GetZCellSizes(), boost::bind(
        &jif3D::ThreeDGravityModel::SetZCellSizes, GravModel));
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

//Another helper function to copy values for the background thicknesses and densities
void SetupBackground(const double *Thicknesses, const double *Densities,
    const unsigned int nlayers)
  {
    std::vector<double> Thick, Dens;
    std::copy(Thicknesses, Thicknesses + nlayers, std::back_inserter(Thick));
    std::copy(Densities, Densities + nlayers, std::back_inserter(Dens));
    GravModel->SetBackgroundDensities(Dens);
    GravModel->SetBackgroundThicknesses(Thick);
  }

//perform a forward calculation for scalar gravity data
void CalcScalarForward(const double *XSizes, const unsigned int *nx,
    const double *YSizes, const unsigned int *ny, const double *ZSizes,
    const unsigned int *nz, const double *Densities,
    const double *BG_Thicknesses, const double *BG_Densities,
    const unsigned int *nbg_layers, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int *nmeas,
    double *GravAcceleration)
  {
    SetupModel(XSizes, nx, YSizes, ny, ZSizes, nz, Densities, XMeasPos,
        YMeasPos, ZMeasPos, nmeas);
    SetupBackground(BG_Thicknesses, BG_Densities, *nbg_layers);
    jif3D::rvec scalarmeas(ScalarGravCalculator->Calculate(*GravModel.get()));
    std::copy(scalarmeas.begin(), scalarmeas.end(), GravAcceleration);
  }

//perform a forward calculation for tensor gravity data
void CalcTensorForward(const double *XSizes, const unsigned int *nx,
    const double *YSizes, const unsigned int *ny, const double *ZSizes,
    const unsigned int *nz, const double *Densities,
    const double *BG_Thicknesses, const double *BG_Densities,
    const unsigned int *nbg_layers, const double *XMeasPos,
    const double *YMeasPos, const double *ZMeasPos, const unsigned int *nmeas,
    double *GravTensor)
  {
    SetupModel(XSizes, nx, YSizes, ny, ZSizes, nz, Densities, XMeasPos,
        YMeasPos, ZMeasPos, nmeas);
    SetupBackground(BG_Thicknesses, BG_Densities, *nbg_layers);
    jif3D::rvec tensormeas(TensorGravCalculator->Calculate(*GravModel.get()));

    std::copy(tensormeas.begin(), tensormeas.end(), GravTensor);
  }

