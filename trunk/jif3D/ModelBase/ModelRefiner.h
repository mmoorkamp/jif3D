//============================================================================
// Name        : ModelRefiner.h
// Author      : Feb 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef MODELREFINER_H_
#define MODELREFINER_H_

#include "../Global/Serialization.h"
#include "ThreeDModelBase.h"
#include "../Global/VecMat.h"

namespace jif3D
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
    public:
      //! Set the refinement coordinates for the x-Axis
      /*! This function takes a vector with the coordinates along the x-Axis
       * at which new interfaces should be inserted
       * @param XCoordinates The coordinates of the extra boundaries in m
       */
      void SetXCoordinates(const ThreeDModelBase::t3DModelDim &XCoordinates)
        {
          RefiningXCoordinates = XCoordinates;
        }
      //! Set the refinement coordinates for the y-Axis
      /*! This function takes a vector with the coordinates along the y-Axis
       * at which new interfaces should be inserted
       * @param YCoordinates The coordinates of the extra boundaries in m
       */
      void SetYCoordinates(const ThreeDModelBase::t3DModelDim &YCoordinates)
        {
          RefiningYCoordinates = YCoordinates;
        }
      //! Set the refinement coordinates for the Z-Axis
      /*! This function takes a vector with the coordinates along the z-Axis
       * at which new interfaces should be inserted
       * @param ZCoordinates The coordinates of the extra boundaries in m
       */
      void SetZCoordinates(const ThreeDModelBase::t3DModelDim &ZCoordinates)
        {
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
      jif3D::rvec CombineGradient(const jif3D::rvec &FineGradient,
          const ThreeDModelBase &CoarseModel, const ThreeDModelBase &RefinedModel);
      ModelRefiner();
      virtual ~ModelRefiner();
      };
  /* @} */
  }

#endif /* MODELREFINER_H_ */
