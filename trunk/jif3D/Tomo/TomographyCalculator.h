//============================================================================
// Name        : TomographyCalculator.h
// Author      : Apr 14, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef TOMOGRAPHYCALCULATOR_H_
#define TOMOGRAPHYCALCULATOR_H_

#include "../Tomo/tomo_types.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/TomographyData.h"
#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VectorTransform.h"
#include <boost/shared_ptr.hpp>

namespace jif3D
  {
    /** \addtogroup tomo Seismic tomography classes and functions */
    /* @{ */
    //! This class provides a convenient interface to the Podvin and Lecomte code and Bjoern's additions
    /*! This class calculates first arrival tomography data and the associated gradients for inversion.
     * It is mainly a driver for the eikonal solver by Podvin and Lecomte and the ray tracing routines
     * written by B. Heincke. Given a three-dimensional model of seismic slownesses in s/m it calculates
     * the time it takes the first arrival to travel for a set of specified source-receiver configurations.
     */
    class J3DEXPORT TomographyCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef ThreeDSeismicModel ModelType;
      typedef TomographyData DataType;

    private:
      //! Internally we only calculate traveltime in the area where we have sources and receivers, this is the shift between the model and the internal grid in x-direction
      int minxindex;
      //! Internally we only calculate traveltime in the area where we have sources and receivers, this is the shift between the model and the internal grid in y-direction
      int minyindex;
      //! The number of air layers on top of the model
      size_t nairlayers;
      //! Information about the source receiver geometry in the format of Bjoern's code
      jif3D::GEOMETRY geo;
      //! Information about the model grid in the format of Bjoern's code
      jif3D::GRID_STRUCT grid;
      //! Information about the data in the format of Bjoern's code
      jif3D::DATA_STRUCT data;
      //! The various segments of each raypath through the model cells
      std::vector<jif3D::RP_STRUCT> raypath;
      //! Set a transform for the observed data
      boost::shared_ptr<jif3D::VectorTransform> DataTransform;
      //! Perform the dynamic allocation for the c-structures above
      void Allocate(const size_t ngrid, const size_t ndata, const size_t npos);

    public:
      void  SetDataTransform(boost::shared_ptr<jif3D::VectorTransform> DT)
       {
         DataTransform = DT;
       }

      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & minxindex;
          ar & minyindex;
          ar & nairlayers;
          ar & geo;
          ar & grid;
          ar & data;
          ar & raypath;
          ar & DataTransform;
        }
      //! The constructor sets everything to default values
      TomographyCalculator();
      virtual ~TomographyCalculator();
      //! Return the raypath structure for the last forward modelling
      const std::vector<jif3D::RP_STRUCT> &GetRayPath() const
        {
          return raypath;
        }
      //! Write out ray distribution for the last forward modelling
      void WriteRays(const std::string &filename) const;
      //! Calculate the travel times for the given model
      /*! This is the core call to the forward modeling code. Has to be performed before LQDerivative can be called
       * @param Model The object containing the slowness distribution and measurement setup
       * @return The calculated travel times
       */
      rvec Calculate(const ModelType &Model, const jif3D::TomographyData &Data);
      //! Calculate the least-square derivative for the given model and data difference
      /*! For inversion we need the derivative of a least-squares objective function
       * with respect to the model parameters.
       * @param Model The object containing the slowness distribution and measurement setup, should be the same as the previous call to Calculate
       * @param Misfit The difference between observed and calculated data
       * @return The partial derivatives with respect to each model parameter
       */
      rvec LQDerivative(const ModelType &Model, const jif3D::TomographyData &Data, const rvec &Misfit);
      };
  /* @} */
  }



#endif /* TOMOGRAPHYCALCULATOR_H_ */
