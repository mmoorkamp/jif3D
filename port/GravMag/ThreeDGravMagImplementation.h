//============================================================================
// Name        : ThreeDGravityImplementation.h
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================

#ifndef THREEDGRAVITYIMPLEMENTATION_H_
#define THREEDGRAVITYIMPLEMENTATION_H_

#include <boost/shared_ptr.hpp>

#include "../Global/VecMat.h"
#include "../Global/VectorTransform.h"
#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "ThreeDGravMagCalculator.h"


namespace jif3D
  {
    /** \addtogroup gravity Gravity forward modeling, display and inversion */
    /* @{ */
    //this is just a forward declaration to avoid circular inclusions
    template<class ThreeDModelType> class ThreeDGravMagCalculator;
    //! The base class that provides the interface for the numerical implementation of the gravity forward calculations.
    /*! The calculation of the forward response is split into two class hierarchies that
     * have to be used in conjunction. The classes derived from ThreeDGravityImplementation
     * are responsible for the numerical implementation details of the calculation, i.e. using
     * parallelization with openmp or performing calculations on a graphics card with cuda, and
     * whether FTG or scalar data is calculated. These classes are not directly visible to the
     * user, but only through the calculator classes that determine whether sensitivity information
     * is stored, whether calculations are cached with stored sensitivity information and how this
     * information is processed and stored.
     *
     * This type of design resembles the bridge pattern and allows to freely combine optimized implementations
     * for different platforms with different sensitivity handlings.
     */


    template<class ThreeDModelType>
    class J3DEXPORT ThreeDGravMagImplementation
      {
    private:
      //! A data transform that we can apply to the calculated values, e.g. an invariant of the FTG data
      boost::shared_ptr<VectorTransform> Transform;
      //! Store the geometry information inside this object to avoid issues with thread safe access
      void CacheGeometry(const ThreeDModelType &Model);
      //! Calculate the response of the 1D background for a single measurement, this function has to be implemented in the derived class.
      virtual rvec CalcBackground(const size_t measindex, const double xwidth,
          const double ywidth, const double zwidth, const ThreeDModelType &Model,
          rmat &Sensitivities) = 0;
      //! Calculate the response of the gridded domain for a single measurement, this function has to be implemented in the derived class.
      virtual rvec CalcGridded(const size_t measindex, const ThreeDModelType &Model,
          rmat &Sensitivities) = 0;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & Transform;
          ar & XCoord;
          ar & YCoord;
          ar & ZCoord;
          ar & XSizes;
          ar & YSizes;
          ar & ZSizes;
        }
    protected:
      /*! The access functions for the coordinates are not thread safe
       * so we cache the values once in the default implementation of calculate.
       * Because CalcBackground and CalcGridded can only be called from within
       * there we are guaranteed correct values in the implementation of those functions
       */
      typedef typename ThreeDModelType::t3DModelDim ModelDimType;
      //! The cached values of the x-coordinates in m for each cell
      ModelDimType XCoord;
      //! The cached values of the y-coordinates in m for each cell
      ModelDimType YCoord;
      //! The cached values of the z-coordinates in m for each cell
      ModelDimType ZCoord;
      //! The cached values for the size of each cell in x-direction
      ModelDimType XSizes;
      //! The cached values for the size of each cell in y-direction
      ModelDimType YSizes;
      //! The cached values for the size of each cell in z-direction
      ModelDimType ZSizes;
    public:
      //! Returns the number of data before any transformation is applied
      virtual size_t RawDataPerMeasurement() = 0;
      //! Set a transformation class that should be applied to any data and gradient in the calculation
      void SetDataTransform(boost::shared_ptr<VectorTransform> DataTransform)
        {
          Transform = DataTransform;
        }
      //! We can implement tensor and scalar calculations in the derived classes, this function returns how many data values a single measurement yields and considers any transformation
      size_t GetDataPerMeasurement()
        {
          return Transform ? Transform->GetOutputSize() : RawDataPerMeasurement();
        }
      //! For a given Model calculate the forward response for all measurements and return it as a real vector, the calculator object is passed to process the sensitivity information
      virtual rvec Calculate(const ThreeDModelType &Model,
          jif3D::ThreeDGravMagCalculator<ThreeDModelType> &Calculator);
      //! Calculate the least-squres derivative vector for the given model and vector
      virtual rvec LQDerivative(const ThreeDModelType &Model, const rvec &Misfit);
      ThreeDGravMagImplementation();
      virtual ~ThreeDGravMagImplementation();
      };

    template<class ThreeDModelType>
    ThreeDGravMagImplementation<ThreeDModelType>::ThreeDGravMagImplementation() :
        Transform(), XCoord(), YCoord(), ZCoord(), XSizes(), YSizes(), ZSizes()
      {

      }

    template<class ThreeDModelType>
    ThreeDGravMagImplementation<ThreeDModelType>::~ThreeDGravMagImplementation()
      {
      }

    template<class ThreeDModelType>
    void ThreeDGravMagImplementation<ThreeDModelType>::CacheGeometry(
        const ThreeDModelType &Model)
      {

        XCoord.resize(Model.GetXCoordinates().size());
        YCoord.resize(Model.GetYCoordinates().size());
        ZCoord.resize(Model.GetZCoordinates().size());
        XSizes.resize(Model.GetXCellSizes().size());
        YSizes.resize(Model.GetYCellSizes().size());
        ZSizes.resize(Model.GetZCellSizes().size());
        std::copy(Model.GetXCoordinates().begin(), Model.GetXCoordinates().end(),
            XCoord.begin());
        std::copy(Model.GetYCoordinates().begin(), Model.GetYCoordinates().end(),
            YCoord.begin());
        std::copy(Model.GetZCoordinates().begin(), Model.GetZCoordinates().end(),
            ZCoord.begin());
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
    template<class ThreeDModelType>
    rvec ThreeDGravMagImplementation<ThreeDModelType>::Calculate(
        const ThreeDModelType &Model,
        jif3D::ThreeDGravMagCalculator<ThreeDModelType> &Calculator)
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
        const double modelxwidth = std::accumulate(Model.GetXCellSizes().begin(),
            Model.GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(Model.GetYCellSizes().begin(),
            Model.GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(Model.GetZCellSizes().begin(),
            Model.GetZCellSizes().end(), 0.0);

        // for all measurement points add the responses of the discretized part and the 1D background
        for (size_t i = 0; i < nmeas; ++i)
          {
            // the vector to hold the result for the current measurement
            rvec currdata(RawDataPerMeasurement());

            currdata = CalcGridded(i, Model, Calculator.SetCurrentSensitivities());

            //adjust for the effects of finite extents of the grid
            currdata += CalcBackground(i, modelxwidth, modelywidth, modelzwidth, Model,
                Calculator.SetCurrentSensitivities());
            if (Transform)
              {
                //now apply the transformation to the data
                rvec currresult(Transform->Transform(currdata));
                std::copy(currresult.begin(), currresult.end(),
                    result.begin() + (i * GetDataPerMeasurement()));
              }
            else
              {
                std::copy(currdata.begin(), currdata.end(),
                    result.begin() + (i * GetDataPerMeasurement()));
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
    template<class ThreeDModelType>
    rvec ThreeDGravMagImplementation<ThreeDModelType>::LQDerivative(
        const ThreeDModelType &Model, const rvec &Misfit)
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
        const double modelxwidth = std::accumulate(Model.GetXCellSizes().begin(),
            Model.GetXCellSizes().end(), 0.0);
        const double modelywidth = std::accumulate(Model.GetYCellSizes().begin(),
            Model.GetYCellSizes().end(), 0.0);
        const double modelzwidth = std::accumulate(Model.GetZCellSizes().begin(),
            Model.GetZCellSizes().end(), 0.0);
        //the total number of model parameters, gridded domain + background
        const size_t nmod = Model.GetNModelParm();
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
            currdata += CalcBackground(i, modelxwidth, modelywidth, modelzwidth, Model,
                CurrentSensitivities);
            if (Transform)
              {
                rmat TransDeriv(Transform->Derivative(currdata));
                //treat the rows of the sensitivity matrix like columns
                //to implicitly perform the transpose for the adjoint
                //we might have more than one datum, e.g, for FTG data
                NewSensitivities = ublas::prod(TransDeriv, CurrentSensitivities);
              }
            else
              {
                NewSensitivities = CurrentSensitivities;
              }
            for (size_t j = 0; j < TransDataPerMeas; ++j)
              {
                boost::numeric::ublas::matrix_row<rmat> CurrRow(NewSensitivities, j);
                DerivMod += Misfit(i * TransDataPerMeas + j) * CurrRow;
              }
          }
        return 2.0 * DerivMod;
      }
  /* @} */
  }
#endif /* THREEDGRAVITYIMPLEMENTATION_H_ */
