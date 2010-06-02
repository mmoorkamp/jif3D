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
      FineSlownessModel(), CoarseSlownessModel(), Refiner(), wantrefinement(
          false), ObservedData(), Calculator()
      {
      }

    TomographyObjective::~TomographyObjective()
      {
      }

    void TomographyObjective::ImplDataDifference(const jiba::rvec &Model,
        jiba::rvec &Diff)
      {
        //make sure that our coarse model has the right number of grid elements for the model vector
        assert(Model.size() == CoarseSlownessModel.GetSlownesses().num_elements());
        //Copy the model vector into the object with the geometry information
        std::copy(Model.begin(), Model.end(),
            CoarseSlownessModel.SetSlownesses().origin());

        jiba::rvec SynthData;
        //project the values from the coarse model onto the fine model
        if (wantrefinement)
          {
            Refiner.RefineModel(CoarseSlownessModel, FineSlownessModel);
            //Calculate the travel times for the 3D model
            SynthData = Calculator.Calculate(FineSlownessModel);
          }
        else
          {
            SynthData = Calculator.Calculate(CoarseSlownessModel);
          }
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

        if (wantrefinement)
          {
            Refiner.RefineModel(CoarseSlownessModel, FineSlownessModel);
            //calculate the gradient for the fine model
            jiba::rvec FineGradient(Calculator.LQDerivative(FineSlownessModel,
                Diff));
            //and return the projection of the fine gradient onto the coarse model
            return Refiner.CombineGradient(FineGradient, CoarseSlownessModel,
                FineSlownessModel);
          }
        //we only get here if we do not do any refinement
        //omitting the else saves us a compiler warning
        return Calculator.LQDerivative(CoarseSlownessModel, Diff);
      }
  }
