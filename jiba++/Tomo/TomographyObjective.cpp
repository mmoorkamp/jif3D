//============================================================================
// Name        : TomographyObjective.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#include "TomographyObjective.h"

namespace jiba
  {

    TomographyObjective::TomographyObjective() :
      FineSlownessModel(), CoarseSlownessModel(), ObservedData(), Calculator()
      {
      }

    TomographyObjective::~TomographyObjective()
      {
      }

    void TomographyObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {

        assert(Model.size() == CoarseSlownessModel.GetSlownesses().num_elements());
        //Copy the model vector into the object with the geometry information
        std::copy(Model.begin(), Model.end(),
            CoarseSlownessModel.SetSlownesses().origin());

        Refiner.SetXCoordinates(FineSlownessModel.GetXCoordinates());
        Refiner.SetYCoordinates(FineSlownessModel.GetYCoordinates());
        Refiner.SetZCoordinates(FineSlownessModel.GetZCoordinates());
        Refiner.RefineModel(CoarseSlownessModel, FineSlownessModel);

        //Calculate the travel times for the 3D model
        jiba::rvec SynthData(Calculator.Calculate(FineSlownessModel));
        Diff.resize(ObservedData.size());
        assert(SynthData.size() == ObservedData.size());
        //calculate the difference between observed and synthetic
        std::transform(SynthData.begin(), SynthData.end(),
            ObservedData.begin(), Diff.begin(), std::minus<double>());
      }

    jiba::rvec TomographyObjective::ImplGradient(const jiba::rvec &Model,
        const jiba::rvec &Diff)
      {
        assert(Model.size() == CoarseSlownessModel.GetSlownesses().num_elements());
        //Copy the model vector into the object with the geometry information
        std::copy(Model.begin(), Model.end(),
            CoarseSlownessModel.SetSlownesses().origin());

        jiba::ThreeDSeismicModel FineModel;
        Refiner.RefineModel(CoarseSlownessModel, FineSlownessModel);
        //calculate the gradient
        jiba::rvec FineGradient(
            Calculator.LQDerivative(FineSlownessModel, Diff));
        return Refiner.CombineGradient(FineGradient, CoarseSlownessModel,
            FineSlownessModel);
      }
  }
