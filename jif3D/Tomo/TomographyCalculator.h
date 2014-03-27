//============================================================================
// Name        : TomographyCalculator.h
// Author      : Apr 14, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef TOMOGRAPHYCALCULATOR_H_
#define TOMOGRAPHYCALCULATOR_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include "../Global/VecMat.h"
#include "ThreeDSeismicModel.h"
#include "tomo_types.h"

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
    class TomographyCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef ThreeDSeismicModel ModelType;
    private:
      //! Do we want to write out the file showing the ray distribution
      bool writerays;
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
      //! Perform the dynamic allocation for the c-structures above
      void Allocate(const size_t ngrid, const size_t ndata, const size_t npos);
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & writerays;
          ar & minxindex;
          ar & minyindex;
          ar & nairlayers;
          ar & geo;
          ar & grid;
          ar & data;
          ar & raypath;
        }
    public:
      //! We can tell the forward modelling object to write out the rays whenever it performs a calculation
      /*! Writing out rays is helpful to show coverage and identify problems, but the
       * files can be several Gigabytes, so we have this as an option.
       * @param saverays If True we write out a .vtk file showing the rays in the model
       */
      explicit TomographyCalculator(bool saverays = false);
      virtual ~TomographyCalculator();
      //! Return the raypath structure for the last forward modelling
      const std::vector<jif3D::RP_STRUCT> &GetRayPath() const
        {
          return raypath;
        }
      //! Calculate the travel times for the given model
      /*! This is the core call to the forward modeling code. Has to be performed before LQDerivative can be called
       * @param Model The object containing the slowness distribution and measurement setup
       * @return The calculated travel times
       */
      rvec Calculate(const ModelType &Model);
      //! Calculate the least-square derivative for the given model and data difference
      /*! For inversion we need the derivative of a least-squares objective function
       * with respect to the model parameters.
       * @param Model The object containing the slowness distribution and measurement setup, should be the same as the previous call to Calculate
       * @param Misfit The difference between observed and calculated data
       * @return The partial derivatives with respect to each model parameter
       */
      rvec LQDerivative(const ModelType &Model, const rvec &Misfit);
      };
  /* @} */
  }

BOOST_CLASS_EXPORT_KEY(jif3D::TomographyCalculator)

#endif /* TOMOGRAPHYCALCULATOR_H_ */
