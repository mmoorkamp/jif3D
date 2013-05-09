//============================================================================
// Name        : ModelRefiner.h
// Author      : Feb 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef MODELREFINER_H_
#define MODELREFINER_H_

#include <boost/serialization/serialization.hpp>
#include "ThreeDModelBase.h"
#include "../Global/VecMat.h"

namespace jiba
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */
    //! Refine a 3D model by inserting extra cells at specified coordinates
    /*! In the joint inversion every method has different gridding requirements.
     * In addition we generally might want to use a coarser inversion grid
     * than we use for forward modeling. This class performs the refinement
     * between the inversion grid and the modeling grid.
     *
     * The refinement procedure is relatively simple. We always keep the cell
     * boundaries of the input model (usually the inversion model) and add
     * new cell boundaries at coordinates specified for refinement in this class.
     * Therefore no interpolation is needed. Also we do not extend the modeling
     * domain beyond the original extent.
     */
    class ModelRefiner
      {
    private:
      //! The coordinates in x-direction in the refined model
      ThreeDModelBase::t3DModelDim RefiningXCoordinates;
      //! The coordinates in y-direction in the refined model
      ThreeDModelBase::t3DModelDim RefiningYCoordinates;
      //! The coordinates in z-direction in the refined model
      ThreeDModelBase::t3DModelDim RefiningZCoordinates;
      //! Refine the coordinate and size information for a single axis
      void RefineOneAxis(const ThreeDModelBase::t3DModelDim &OldCoordinates,
          const ThreeDModelBase::t3DModelDim &RefCoordinates,
          const ThreeDModelBase::t3DModelDim &OldSizes,
          ThreeDModelBase::t3DModelDim &NewSizes);
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void save(Archive & ar, const unsigned int version) const
        {
          //we have to split saving and loading because we cannot
          //directly serialize a multi_array
          //so we store the number of elements in each array
          ar & RefiningXCoordinates.shape()[0];
          ar & RefiningYCoordinates.shape()[0];
          ar & RefiningZCoordinates.shape()[0];
          //and then save the contents as an ordinary continuous array
          ar
              & boost::serialization::make_array(RefiningXCoordinates.origin(),
                  RefiningXCoordinates.num_elements());
          ar
              & boost::serialization::make_array(RefiningYCoordinates.origin(),
                  RefiningYCoordinates.num_elements());
          ar
              & boost::serialization::make_array(RefiningZCoordinates.origin(),
                  RefiningZCoordinates.num_elements());
        }
      template<class Archive>
      void load(Archive & ar, const unsigned int version)
        {
          //for loading we first need to get the number of elements
          // in each multi-array
          size_t nx, ny, nz;
          ar & nx;
          ar & ny;
          ar & nz;
          //allocate memory
          RefiningXCoordinates.resize(boost::extents[nx]);
          RefiningYCoordinates.resize(boost::extents[ny]);
          RefiningZCoordinates.resize(boost::extents[nz]);
          //and then do a raw read into the preallocated multi-arrays
          ar
              & boost::serialization::make_array(RefiningXCoordinates.origin(),
                  RefiningXCoordinates.num_elements());
          ar
              & boost::serialization::make_array(RefiningYCoordinates.origin(),
                  RefiningYCoordinates.num_elements());
          ar
              & boost::serialization::make_array(RefiningZCoordinates.origin(),
                  RefiningZCoordinates.num_elements());
        }
      BOOST_SERIALIZATION_SPLIT_MEMBER()
    public:
      //! Set the refinement coordinates for the x-Axis
      /*! This function takes a vector with the coordinates along the x-Axis
       * at which new interfaces should be inserted
       * @param XCoordinates The coordinates of the extra boundaries in m
       */
      void SetXCoordinates(const ThreeDModelBase::t3DModelDim &XCoordinates)
        {
          RefiningXCoordinates.resize(boost::extents[XCoordinates.size()]);
          RefiningXCoordinates = XCoordinates;
        }
      //! Set the refinement coordinates for the y-Axis
      /*! This function takes a vector with the coordinates along the y-Axis
       * at which new interfaces should be inserted
       * @param YCoordinates The coordinates of the extra boundaries in m
       */
      void SetYCoordinates(const ThreeDModelBase::t3DModelDim &YCoordinates)
        {
          RefiningYCoordinates.resize(boost::extents[YCoordinates.size()]);
          RefiningYCoordinates = YCoordinates;
        }
      //! Set the refinement coordinates for the Z-Axis
      /*! This function takes a vector with the coordinates along the z-Axis
       * at which new interfaces should be inserted
       * @param ZCoordinates The coordinates of the extra boundaries in m
       */
      void SetZCoordinates(const ThreeDModelBase::t3DModelDim &ZCoordinates)
        {
          RefiningZCoordinates.resize(boost::extents[ZCoordinates.size()]);
          RefiningZCoordinates = ZCoordinates;
        }
      //! Refine the axes of the InputModel and store the result in RefinedModel
      /*! This function uses the specified refinement coordinates and the
       * cell coordinates in the InputModel to create the output RefinedModel.
       * After this call the cell data structure will be properly allocated but
       * without meaningful values. Use ProjectValues to fill the grid.
       * @param InputModel The original coarse model
       * @param RefinedModel The refined model with extra cell boundaries
       */
      void RefineAxes(const ThreeDModelBase &InputModel, ThreeDModelBase &RefinedModel);
      //! Project the values from the coarse input model onto the refined model, assumes RefinedModel has allocated memory for the new grid
      /*! This function is usually called after RefineAxes. It assumes that the grid of cell values has been allocated
       * in RefinedModel to hold the refined values. It maps the cell values of the coarse model onto
       * one or more cells in the refined model.
       * @param InputModel The original coarse model
       * @param RefinedModel The refined model with projected parameters.
       */
      void ProjectValues(const ThreeDModelBase &InputModel,
          ThreeDModelBase &RefinedModel);
      //! A convenience function that combines the calls to RefineAxes and ProjectValues
      void RefineModel(const ThreeDModelBase &InputModel, ThreeDModelBase &RefinedModel)
        {
          RefineAxes(InputModel, RefinedModel);
          ProjectValues(InputModel, RefinedModel);
        }
      //! Combine the gradient information from the cells of the fine grid on the coarse grid
      /*! When we refine a coarse model for the forward calculation, we have to combine the gradient
       * that we get from this refined model for the coarse model. This function takes the refined gradient
       * and the coarse model as arguments and returns the gradient that corresponds to the cells of the coarse model.
       * @param FineGradient The gradient information for each cell of the refined model
       * @param CoarseModel The coarse model that was refined for the forward calculation
       * @param RefinedModel The grid information for the refined model is taken from this class
       * @return The gradient for each cell of the coarse model
       */
      jiba::rvec CombineGradient(const jiba::rvec &FineGradient,
          const ThreeDModelBase &CoarseModel, const ThreeDModelBase &RefinedModel);
      ModelRefiner();
      virtual ~ModelRefiner();
      };
  /* @} */
  }

#endif /* MODELREFINER_H_ */
