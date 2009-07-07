//============================================================================
// Name        : X3DModel.h
// Author      : Jul 2, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef X3DMODEL_H_
#define X3DMODEL_H_

#include "ThreeDMTModel.h"

namespace jiba
  {

    class X3DModel: public ThreeDMTModel
      {
    private:
      std::vector<double> Frequencies;
      std::vector<double> bg_thicknesses;
      std::vector<double> bg_conductivities;
    public:
      enum ProblemType
        {
        MT, CSMT, EDIP, MDIP
        };
      const std::vector<double> &GetFrequencies() const
        {
          return Frequencies;
        }
      std::vector<double> &SetFrequencies()
        {
          return Frequencies;
        }
      //! Set the thicknesses of the background layers, the individual thicknesses are given in m
      void SetBackgroundThicknesses(const std::vector<double> &value)
        {
          bg_thicknesses.clear();
          copy(value.begin(), value.end(), back_inserter(bg_thicknesses));
        }
      //! Return the thicknesses of the background layers in m
      const std::vector<double> &GetBackgroundThicknesses() const
        {
          return bg_thicknesses;
        }
      //! Set the conductivities of the background layers, the individual thicknesses are given in S/m
      void SetBackgroundConductivities(const std::vector<double> &value)
        {
          bg_conductivities.clear();
          bg_conductivities.reserve(value.size());
          copy(value.begin(), value.end(), back_inserter(bg_conductivities));
        }
      //! Return the conductivities of the background layers in S/m
      const std::vector<double> &GetBackgroundConductivities() const
        {
          return bg_conductivities;
        }
      //! The MT model for X3D by Avdeev et al. has the same cell size for all cells in the two horizontal directions so we just have one function to set it
      /*! This function sets both the size of all cells as well as the number of cells in the horizontal (x and y) directions
       * @param Size The size of each cell in all directions in m
       * @param nx The number of cells in x-direction (North)
       * @param ny The number of cells in y-direction (East)
       */
      void SetHorizontalCellSize(const double Size, const size_t nx,
          const size_t ny)
        {
          ThreeDModelBase::SetXCellSizes().resize(boost::extents[nx]);
          std::fill_n(ThreeDModelBase::SetXCellSizes().begin(), nx, Size);
          ThreeDModelBase::SetYCellSizes().resize(boost::extents[ny]);
          std::fill_n(ThreeDModelBase::SetYCellSizes().begin(), ny, Size);
        }
      //! The vertical cells can all have different sizes so we allow direct access to the CellSize structure
      t3DModelDim &SetZCellSizes()
        {
          return ThreeDModelBase::SetZCellSizes();
        }
      X3DModel();
      virtual ~X3DModel();
      };

  }

#endif /* X3DMODEL_H_ */
