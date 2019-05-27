//============================================================================
// Name        : X3DModel.h
// Author      : Jul 2, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef X3DMODEL_H_
#define X3DMODEL_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "ThreeDMTModel.h"

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! This class stores MT models with the additional information required for calculations with x3d
    /*! The x3d forward calculation code requires the gridded "anomalous" domain to be embedded in a layered
     * background. This class provides the necessary functionality to store information about the background
     * so that it can seamlessly used with the x3d forward calculation classes.
     *
     */
    class J3DEXPORT X3DModel: public ThreeDMTModel
      {
    private:
      //! The thicknesses of the background layers in m
      std::vector<double> bg_thicknesses;
      //! The conductivities of the background layers in S/m
      std::vector<double> bg_conductivities;
    public:
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDMTModel>(*this);
          ar & bg_thicknesses;
          ar & bg_conductivities;
        }
      //! The problem type we want to perform the calculation for
      /*! We can use this enum to specify which type of forward calculation
       * we want to perform, magnetotelluric (MT), controlled source EM (CSEM),
       * electric dipole (EDIP), or magnetic dipole (MDIP).
       */
      enum ProblemType
        {
        MT, CSMT, EDIP, MDIP
        };

      //! Set the thicknesses of the background layers, the individual thicknesses are given in m
      void SetBackgroundThicknesses(const std::vector<double> &value)
        {
          bg_thicknesses = value;
        }
      //! Return the thicknesses of the background layers in m
      const std::vector<double> &GetBackgroundThicknesses() const
        {
          return bg_thicknesses;
        }
      //! Set the conductivities of the background layers, the individual thicknesses are given in S/m
      void SetBackgroundConductivities(const std::vector<double> &value)
        {
          bg_conductivities = value;
        }
      //! Return the conductivities of the background layers in S/m
      const std::vector<double> &GetBackgroundConductivities() const
        {
          return bg_conductivities;
        }
      //! The MT model for X3D by Avdeev et al. has the same cell size for all cells in each horizontal directions so we just have one function to set it
      /*! This function sets both the size of all cells as well as the number of cells in the horizontal (x and y) directions
       * @param XSize The size of each cell in x-direction (North) in m
       * @param YSize The size of each cell in y-direction (East) in m
       * @param nx The number of cells in x-direction (North)
       * @param ny The number of cells in y-direction (East)
       */
      void SetHorizontalCellSize(const double XSize, const double YSize, const size_t nx,
          const size_t ny)
        {
          ThreeDModelBase::t3DModelDim XCS(nx, XSize);
          ThreeDModelBase::t3DModelDim YCS(ny, YSize);
          SetXCellSizes(XCS);
          SetYCellSizes(YCS);
        }
      //! The vertical cells can all have different sizes so we allow direct access to the CellSize structure
      void SetZCellSizes(ThreeDModelBase::t3DModelDim &ZCS)
        {
          ThreeDModelBase::SetZCellSizes(ZCS);
        }
      //! Given three coordinates in m, find the indices of the model cell that corresponds to these coordinates, this is a more efficient implementation than the one in the base class
      virtual boost::array<ThreeDModelBase::t3DModelData::index, 3>
      FindAssociatedIndices(const double xcoord, const double ycoord,
          const double zcoord) const override;
      //! Write all model information to a netcdf file
      virtual void WriteNetCDF(const std::string &filename) const override;
      //! Read all model information from a netcdf file
      virtual void ReadNetCDF(const std::string &filename) override;

      //! Other models will be copied by the copy operator for the base class
      X3DModel& operator=(const ThreeDModelBase& source);
      bool operator ==(const X3DModel &b) const;
      X3DModel();

      virtual ~X3DModel();
      };
  /* @} */
  }

#endif /* X3DMODEL_H_ */
