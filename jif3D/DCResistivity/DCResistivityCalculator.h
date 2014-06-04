//============================================================================
// Name        : DCResistivityCalculator.h
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#ifndef DCRESISTIVITYCALCULATOR_H_
#define DCRESISTIVITYCALCULATOR_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/export.hpp>
#include "../Global/VecMat.h"
#include "ThreeDDCResistivityModel.h"
#include "DCResForwardBase.h"

namespace jif3D
  {
    /** \addtogroup DC Resistivity classes and functions */
    /* @{ */
    //! This class provides a convenient interface to the DCResForwarBase code
    /*! This class calculates DC Resistivity Forward Response data and the associated gradients for inversion.
     */
    class DCResistivityCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef ThreeDDCResistivityModel ModelType;
    private:
      //! Information about the source receiver geometry in the format of Bjoern's code
      jif3D::GEOMETRY_RES geo;
      //! Information about the model grid in the format of Bjoern's code
      jif3D::GRID_STRUCT_RES grid;
      //! Information about the data in the format of Bjoern's code
      jif3D::DATA_STRUCT_RES data;
      //! Information about the gradient in the format of Bjoern's code
      jif3D::GRADIENT_STRUCT_RES grad;
      //! Perform the dynamic allocation for the c-structures above
      void Allocate(const size_t ngrid, const size_t ndata, const size_t nshot, const size_t nmeaspoint);
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & geo;
          ar & grid;
          ar & data;
          ar & grad;
        }
    public:
      DCResistivityCalculator();
      virtual ~DCResistivityCalculator();
      //! Calculate the DC Resistivity Forward Response for the given model
      /*! This is the core call to the forward modelling code. Has to be performed before LQDerivative can be called
       * @param Model The object containing the resistivity distribution and measurement setup
       * @return The calculated DC Resistivity Forward Response
       */
      rvec Calculate(const ModelType &Model);
      //! Calculate the least-square derivative for the given model and data difference
      /*! For inversion we need the derivative of a least-squares objective function
       * with respect to the model parameters.
       * @param Model The object containing the resistivity distribution and measurement setup, should be the same as the previous call to Calculate
       * @param Misfit The difference between observed and calculated data
       * @return The partial derivatives with respect to each model parameter
       */
      rvec LQDerivative(const ModelType &Model, const rvec &Misfit);
      };
  /* @} */
  }

BOOST_CLASS_EXPORT_KEY(jif3D::DCResistivityCalculator)

#endif /* DCRESISTIVITYCALCULATOR_H_ */
