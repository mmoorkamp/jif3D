//============================================================================
// Name        : ModelRefiner.h
// Author      : Feb 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef MODELREFINER_H_
#define MODELREFINER_H_
#include "ThreeDModelBase.h"
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
      ThreeDModelBase::t3DModelDim RefiningXCoordinates;
      ThreeDModelBase::t3DModelDim RefiningYCoordinates;
      ThreeDModelBase::t3DModelDim RefiningZCoordinates;
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
      void SetXCoordinates(ThreeDModelBase::t3DModelDim &XCoordinates)
        {
          RefiningXCoordinates.resize(boost::extents[XCoordinates.size()]);
          RefiningXCoordinates = XCoordinates;
        }
      //! Set the refinement coordinates for the y-Axis
      /*! This function takes a vector with the coordinates along the y-Axis
       * at which new interfaces should be inserted
       * @param YCoordinates The coordinates of the extra boundaries in m
       */
      void SetYCoordinates(ThreeDModelBase::t3DModelDim &YCoordinates)
        {
          RefiningYCoordinates.resize(boost::extents[YCoordinates.size()]);
          RefiningYCoordinates = YCoordinates;
        }
      //! Set the refinement coordinates for the Z-Axis
      /*! This function takes a vector with the coordinates along the z-Axis
       * at which new interfaces should be inserted
       * @param ZCoordinates The coordinates of the extra boundaries in m
       */
      void SetZCoordinates(ThreeDModelBase::t3DModelDim &ZCoordinates)
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
      void RefineAxes(const ThreeDModelBase &InputModel,
          ThreeDModelBase &RefinedModel);
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
      void RefineModel(const ThreeDModelBase &InputModel,
          ThreeDModelBase &RefinedModel)
        {
          RefineAxes(InputModel, RefinedModel);
          ProjectValues(InputModel, RefinedModel);
        }
      ModelRefiner();
      virtual ~ModelRefiner();
      };
  /* @} */
  }

#endif /* MODELREFINER_H_ */
