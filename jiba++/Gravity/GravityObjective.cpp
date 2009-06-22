//============================================================================
// Name        : GravityObjective.cpp
// Author      : Apr 16, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "GravityObjective.h"
#include "MinMemGravityCalculator.h"
#include "DiskGravityCalculator.h"
#include "ThreeDGravityFactory.h"
#include <cassert>
namespace jiba
  {
    //copy a vector of densities into a ThreeDGravityModel object
    //this function checks whether the vector only contains elements
    //corresponding to the gridded part or also to the background
    //the background is assumed to be at the end
    void VectorToModel(const jiba::rvec &InVector,
        jiba::ThreeDGravityModel &Model)
      {
        const size_t ngrid = Model.GetDensities().num_elements();
        assert(InVector.size() == ngrid ||
            InVector.size() == ngrid+ Model.GetBackgroundDensities().size());
        std::copy(InVector.begin(), InVector.begin() + ngrid,
            Model.SetDensities().origin());
        //if we also have values for the background
        //copy them an set the values appropriately
        if (InVector.size() == ngrid + Model.GetBackgroundDensities().size())
          {
            std::vector<double> Background(InVector.begin() + ngrid,
                InVector.end());
            Model.SetBackgroundDensities(Background);
          }
      }

    GravityObjective::GravityObjective(bool ftg, bool cuda) :
      Calculator(), ObservedData(), DensityModel()
      {
        //check whether we want to minimize FTG data and allocate the appropriate object
        //for forward calculation
        if (ftg)
          {
            Calculator
                = boost::shared_ptr<jiba::ThreeDGravityCalculator>(
                    jiba::CreateGravityCalculator<jiba::DiskGravityCalculator>::MakeTensor(
                        cuda));
          }
        else
          {
            Calculator
                = boost::shared_ptr<jiba::ThreeDGravityCalculator>(
                    jiba::CreateGravityCalculator<jiba::DiskGravityCalculator>::MakeScalar(
                        cuda));
          }
      }

    GravityObjective::~GravityObjective()
      {

      }

    void GravityObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        assert(DensityModel.GetMeasPosX().size() * Calculator->GetDataPerMeasurement() == ObservedData.size() );
        //copy the model vector to the 3D object that contains the geometry information
        VectorToModel(Model, DensityModel);
        //calculate the synthetic data
        jiba::rvec SynthData(Calculator->Calculate(DensityModel));
        Diff.resize(ObservedData.size());
        //calculate the difference between observed and synthetic
        std::transform(SynthData.begin(), SynthData.end(),
            ObservedData.begin(), Diff.begin(), std::minus<double>());
      }

    jiba::rvec GravityObjective::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        assert(DensityModel.GetMeasPosX().size()* Calculator->GetDataPerMeasurement() == Diff.size() );
        //copy the model vector to the 3D object that contains the geometry information
        VectorToModel(Model, DensityModel);
        //calculate and return the gradient
        return Calculator->LQDerivative(DensityModel, Diff);
      }
  }
