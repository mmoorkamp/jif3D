//============================================================================
// Name        : ThreeDGravityImplementation.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "ThreeDGravMagImplementation.h"
#include "ThreeDGravMagCalculator.h"

namespace jif3D
  {

    ThreeDGravMagImplementation::ThreeDGravMagImplementation() :
      Transform(), XCoord(), YCoord(), ZCoord(), XSizes(), YSizes(), ZSizes()
      {

      }

    ThreeDGravMagImplementation::~ThreeDGravMagImplementation()
      {
      }

    void ThreeDGravMagImplementation::CacheGeometry(
        const ThreeDGravityModel &Model)
      {

        XCoord.resize(boost::extents[Model.GetXCoordinates().size()]);
        YCoord.resize(boost::extents[Model.GetYCoordinates().size()]);
        ZCoord.resize(boost::extents[Model.GetZCoordinates().size()]);
        XSizes.resize(boost::extents[Model.GetXCellSizes().size()]);
        YSizes.resize(boost::extents[Model.GetYCellSizes().size()]);
        ZSizes.resize(boost::extents[Model.GetZCellSizes().size()]);
        std::copy(Model.GetXCoordinates().begin(),
            Model.GetXCoordinates().end(), XCoord.begin());
        std::copy(Model.GetYCoordinates().begin(),
            Model.GetYCoordinates().end(), YCoord.begin());
        std::copy(Model.GetZCoordinates().begin(),
            Model.GetZCoordinates().end(), ZCoord.begin());
        std::copy(Model.GetXCellSizes().begin(), Model.GetXCellSizes().end(),
            XSizes.begin());
        std::copy(Model.GetYCellSizes().begin(), Model.GetYCellSizes().end(),
            YSizes.begin());
        std::copy(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
            ZSizes.begin());

      }

    /*! This function implements the grand structure of gravity forward calculation, i.e. processing
     * geometric information, looping over all measurements and combining the response
     * of the gridded part and the layered background. The details of the calculation,
     * the type of data (scalar or FTG) and platform specific code are implemented in derived classes
     * through the virtual functions CalcGridded and CalcBackground.
     *
     * This function also gets the calling calculator object as an argument. This way it can, if required, copy
     * the current row of the sensitivity matrix to the calculator object and
     * call the HandleSensitivity method that allows the calculator class to store, process
     * or discard the sensitivity information.
     *
     * If we override this function in a derived class, e.g. for CUDA specific issues, make sure to call CheckTransform
     * to guarantee a valid data transfrom object.
     *
     * @param Model The gravity model for which to calculate the data
     * @param Calculator A derived class of ThreeDGravityCalculator
     * @return A vector of measurements
     */
    rvec ThreeDGravMagImplementation::Calculate(
        const ThreeDGravityModel &Model, ThreeDGravMagCalculator &Calculator)
      {

        // get the number of measurements
        // this class is only called from a calculator object
        // which performs consistency check
        const size_t nmeas = Model.GetMeasPosX().size();
        //Cache the coordinate and size information
        //to avoid problems with thread safety
        CacheGeometry(Model);

        //allocate enough memory for all datapoints in the result vector
        rvec result(nmeas * GetDataPerMeasurement());

        // calculate the size of the modelling domain for the background adjustment
        // this way we do not have to recalculate within the loop
        const double modelxwidth = std::accumulate(
            Model.GetXCellSizes().begin(), Model.GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(
            Model.GetYCellSizes().begin(), Model.GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(
            Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(), 0.0);

        // for all measurement points add the responses of the discretized part and the 1D background
        for (size_t i = 0; i < nmeas; ++i)
          {
            // the vector to hold the result for the current measurement
            rvec currdata(RawDataPerMeasurement());

            currdata = CalcGridded(i, Model,
                Calculator.SetCurrentSensitivities());

            //adjust for the effects of finite extents of the grid
            currdata += CalcBackground(i, modelxwidth, modelywidth,
                modelzwidth, Model, Calculator.SetCurrentSensitivities());
            if (Transform)
              {
                //now apply the transformation to the data
                rvec currresult(Transform->Transform(currdata));
                std::copy(currresult.begin(), currresult.end(), result.begin()
                    + (i * GetDataPerMeasurement()));
              }
            else
              {
                std::copy(currdata.begin(), currdata.end(), result.begin() + (i
                    * GetDataPerMeasurement()));
              }
            //give the calculator object the chance to handle the sensitivities
            //for the current measurement
            Calculator.HandleSensitivities(i);

          }
        return result;
      }

    /*! This is the default implementation for the least squares derivatives. Derived classes can choose
     * to override this if the forward engine requires a different approach. However, the scheme is fairly
     * general and should work with most implementations.
     * @param Model The Model for which we want the derivative
     * @param Misfit The Misfit between observed and synthetic data
     * @return The derivative of a least-squares objective function with respect to the model parameters
     */
    rvec ThreeDGravMagImplementation::LQDerivative(
        const ThreeDGravityModel &Model, const rvec &Misfit)
      {
        // get the number of measurements
        // this class is only called from a calculator object
        // which performs consistency check
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t TransDataPerMeas = GetDataPerMeasurement();
        assert(Misfit.size() == nmeas * TransDataPerMeas);
        //Cache the coordinate and size information
        //to avoid problems with thread safety
        CacheGeometry(Model);

        // calculate the size of the modelling domain for the background adjustment
        // this way we do not have to recalculate within the loop
        const double modelxwidth = std::accumulate(
            Model.GetXCellSizes().begin(), Model.GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(
            Model.GetYCellSizes().begin(), Model.GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(
            Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(), 0.0);
        //the total number of model parameters, gridded domain + background
        const size_t nmod = Model.GetDensities().num_elements()
            + Model.GetBackgroundDensities().size();
        //allocate enough memory for all derivatives in the result vector
        rvec DerivMod(nmod);
        std::fill(DerivMod.begin(), DerivMod.end(), 0.0);

        //we need the Sensitivities
        rmat CurrentSensitivities(RawDataPerMeasurement(), nmod);
        rmat NewSensitivities;
        for (size_t i = 0; i < nmeas; ++i)
          {
            rvec currdata(RawDataPerMeasurement());
            //build up the full sensitivity matrix for the current measurement
            currdata = CalcGridded(i, Model, CurrentSensitivities);
            currdata += CalcBackground(i, modelxwidth, modelywidth,
                modelzwidth, Model, CurrentSensitivities);
            if (Transform)
              {
                rmat TransDeriv(Transform->Derivative(currdata));
                //treat the rows of the sensitivity matrix like columns
                //to implicitly perform the transpose for the adjoint
                //we might have more than one datum, e.g, for FTG data
                NewSensitivities
                    = ublas::prod(TransDeriv, CurrentSensitivities);
              }
            else
              {
                NewSensitivities = CurrentSensitivities;
              }
            for (size_t j = 0; j < TransDataPerMeas; ++j)
              {
                boost::numeric::ublas::matrix_row<rmat> CurrRow(
                    NewSensitivities, j);
                DerivMod += Misfit(i * TransDataPerMeas + j) * CurrRow;
              }
          }
        return 2.0 * DerivMod;
      }
  }
