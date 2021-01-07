//============================================================================
// Name        : DCResistivityCalculator.h
// Author      : Zhanjie Shi ,Max Moorkamp and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi ,Max Moorkamp and Richard.W Hobbs
//============================================================================

#ifndef DCRESISTIVITYCALCULATOR_H_
#define DCRESISTIVITYCALCULATOR_H_

#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include "../DCResistivity/DCResistivityData.h"
#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VectorTransform.h"
#include <boost/shared_ptr.hpp>
#include "DCResForwardBase.h"

#include <boost/serialization/export.hpp>




namespace jif3D
  {
    /** \addtogroup DC Resistivity classes and functions */
    /* @{ */
    //! This class provides a convenient interface to the DCResForwarBase code
    /*! This class calculates DC Resistivity Forward Response data and the associated gradients for inversion.
     */
    class J3DEXPORT DCResistivityCalculator
      {
    public:
      //! This type definition is necessary so that ThreeDModelObjective can correctly deduce the native type for a model object for this class
      typedef ThreeDDCResistivityModel ModelType;
      typedef DCResistivityData DataType;

    private:
      //! Information about the source receiver geometry in the format of Bjoern's code
      jif3D::GEOMETRY_RES geo;
      //! Information about the model grid in the format of Bjoern's code
      jif3D::GRID_STRUCT_RES grid;
      //! Perform the dynamic allocation for the c-structures above
      void Allocate(const size_t ngrid, const size_t nshot, const size_t nmeaspoint, const size_t ndata);
      //! Translate a 3D model object to structures used for forward calculation
      void ModelToStruct(const ThreeDDCResistivityModel &Model, const DCResistivityData &Data, jif3D::GEOMETRY_RES &geo,
                jif3D::GRID_STRUCT_RES &grid);
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & geo;
          ar & grid;
        }
    public:
      DCResistivityCalculator();
      virtual ~DCResistivityCalculator();
      //! Calculate the DC Resistivity Forward Response for the given model
      /*! This is the core call to the forward modelling code. Has to be performed before LQDerivative can be called
       * @param Model The object containing the resistivity distribution and measurement setup
       * @return The calculated DC Resistivity Forward Response
       */
      rvec Calculate(const ModelType &Model, const DataType &Data);
      //! Calculate the least-square derivative for the given model and data difference
      /*! For inversion we need the derivative of a least-squares objective function
       * with respect to the model parameters.
       * @param Model The object containing the resistivity distribution and measurement setup, should be the same as the previous call to Calculate
       * @param Misfit The difference between observed and calculated data
       * @return The partial derivatives with respect to each model parameter
       */
      rvec LQDerivative(const ModelType &Model, const DataType &Data, const rvec &Misfit);
      };
  /* @} */
  }

BOOST_CLASS_EXPORT_KEY(jif3D::DCResistivityCalculator)

#endif /* DCRESISTIVITYCALCULATOR_H_ */
