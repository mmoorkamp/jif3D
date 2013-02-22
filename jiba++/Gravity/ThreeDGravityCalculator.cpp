//============================================================================
// Name        : ThreeDGravityCalculator.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#include "../Global/FatalException.h"
#include "ThreeDGravityCalculator.h"

namespace jiba
  {

    ThreeDGravityCalculator::ThreeDGravityCalculator(
        boost::shared_ptr<ThreeDGravityImplementation> TheImp) :
        CurrentSensitivities(), Transform(), Imp(TheImp)
      {

      }

    ThreeDGravityCalculator::~ThreeDGravityCalculator()
      {
      }

    void ThreeDGravityCalculator::CheckModelConsistency(const ThreeDGravityModel &Model)
      {
        //do some sanity checks
        // we can assign cell sizes and model grid independently
        // so we have to check every time
        if (Model.GetDensities().shape()[0] != Model.GetXCellSizes().shape()[0])
          {
            throw jiba::FatalException(
                "Model x-dimension does not match size for specification of cell sizes.");
          }
        if (Model.GetDensities().shape()[1] != Model.GetYCellSizes().shape()[0])
          {
            throw jiba::FatalException(
                "Model y-dimension does not match size for specification of cell sizes.");
          }
        if (Model.GetDensities().shape()[2] != Model.GetZCellSizes().shape()[0])
          {
            throw jiba::FatalException(
                "Model x-dimension does not match size for specification of cell sizes.");
          }

        // make sure we have coordinates for all sites
        //these should always be equal, so we use an assertion
        //to catch strange cases that should not occur
        const size_t nmeas = Model.GetMeasPosX().size();
        assert(nmeas == Model.GetMeasPosY().size());
        assert(nmeas == Model.GetMeasPosZ().size());
      }

    /*! The least squares derivative is the building block for most types of objective functions, here we define
     * the abstract interface. The implementation depends on the type of data and whether we have sensitivity information or not
     * @param Model The 3D gravity model
     * @param Misfit The misfit at which we need the derivative, has to match the type of data in the derived class
     * @return The partial derivative of the objective function, size and storage order depends on the type of data
     */
    rvec ThreeDGravityCalculator::LQDerivative(const ThreeDGravityModel &Model,
        const rvec &Misfit)
      {
        CheckModelConsistency(Model);
        return Imp->LQDerivative(Model, Misfit);
      }

    /*! Given a 3D model this routine calculates the forward response. The type of data is determined
     * by the derived class, e.g. scalar gravity, FTG
     * @param Model The model for which we want the response
     * @return The calculated data, the length of the vector and the order of the data depends on the derived class
     */
    rvec ThreeDGravityCalculator::Calculate(const ThreeDGravityModel &Model)
      {
        CheckModelConsistency(Model);
        return Imp->Calculate(Model, *this);
      }
  }
