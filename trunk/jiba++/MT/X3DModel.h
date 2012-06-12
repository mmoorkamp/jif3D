//============================================================================
// Name        : X3DModel.h
// Author      : Jul 2, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef X3DMODEL_H_
#define X3DMODEL_H_

#include <boost/serialization/serialization.hpp>
#include "ThreeDMTModel.h"

namespace jiba
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! This class stores MT models with the additional information required for calculations with x3d
    /*! The x3d forward calculation code requires the gridded "anomalous" domain to be embedded in a layered
     * background. This class provides the necessary functionality to store information about the background
     * so that it can seamlessly used with the x3d forward calculation classes.
     *
     */
    class X3DModel: public ThreeDMTModel
      {
    private:
      //! The thicknesses of the background layers in m
      std::vector<double> bg_thicknesses;
      //! The conductivities of the background layers in S/m
      std::vector<double> bg_conductivities;
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDMTModel>(*this);
          ar & bg_thicknesses;
          ar & bg_conductivities;
        }
    public:
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
      //! The MT model for X3D by Avdeev et al. has the same cell size for all cells in each horizontal directions so we just have one function to set it
      /*! This function sets both the size of all cells as well as the number of cells in the horizontal (x and y) directions
       * @param XSize The size of each cell in x-direction (North) in m
       * @param YSize The size of each cell in y-direction (East) in m
       * @param nx The number of cells in x-direction (North)
       * @param ny The number of cells in y-direction (East)
       */
      void SetHorizontalCellSize(const double XSize, const double YSize,
          const size_t nx, const size_t ny)
        {
          ThreeDModelBase::SetXCellSizes().resize(boost::extents[nx]);
          std::fill_n(ThreeDModelBase::SetXCellSizes().begin(), nx, XSize);
          ThreeDModelBase::SetYCellSizes().resize(boost::extents[ny]);
          std::fill_n(ThreeDModelBase::SetYCellSizes().begin(), ny, YSize);
        }
      //! The vertical cells can all have different sizes so we allow direct access to the CellSize structure
      t3DModelDim &SetZCellSizes()
        {
          return ThreeDModelBase::SetZCellSizes();
        }
      //! Copy the source and receiver positions and the indices for the source-receiver combinations from the source model, also copies frequencies for calculation
      void CopyMeasurementConfigurations(const X3DModel &Source)
        {
          ClearMeasurementPoints();
          const size_t nmeas = Source.GetMeasPosX().size();
          for (size_t i = 0; i < nmeas; ++i)
            {
              AddMeasurementPoint(Source.GetMeasPosX()[i],
                  Source.GetMeasPosY()[i], Source.GetMeasPosZ()[i]);
            }
          SetFrequencies() = Source.GetFrequencies();
        }
      //! Write all model information to a netcdf file
      void WriteNetCDF(const std::string filename) const;
      //! Read all model information from a netcdf file
      void ReadNetCDF(const std::string filename);
      //! The copy operator for X3DModels
      X3DModel& operator=(const X3DModel& source);
      //! Other models will be copied by the copy operator for the base class
      X3DModel& operator=(const ThreeDModelBase& source);
      X3DModel();
      //! We define our own copy constructor. This has to be updated if additional information is stored in this object.
      X3DModel(const X3DModel &source);
      virtual ~X3DModel();
      };
  /* @} */
  }

#endif /* X3DMODEL_H_ */
