//============================================================================
// Name        : ThreeDSeismicModel.h
// Author      : Apr 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef THREEDSEISMICMODEL_H_
#define THREEDSEISMICMODEL_H_

#include "../Global/Serialization.h"
#include "../Global/FatalException.h"
#include "../ModelBase/ThreeDModelBase.h"
#include "../Global/Jif3DGlobal.h"

namespace jif3D
  {
    /** \addtogroup tomo Seismic tomography classes and functions
     * This module contains the functionality for seismic refraction forward modeling and data handling.
     *
     * As for the other methods the 3D Models are stored in a class derived from ThreeDModelBase. The model
     * parameter for these models is slowness in s/m. At the moment we only have one forward modeling code based
     * on the algorithm of Podvin and Lecomte with ray tracing implemented by B. Heincke.
     *
     * We therefore only have a single forward modeling class TomographyCalculator that is used to calculate arrival
     * times from a given slowness model. In contrast to MT and Gravity models we also need to store the position
     * of the sources in the model class.
     * */
    /* @{ */
    class TomoExtraParameterSetter;

    class TomographyData;
    //! The seismic model class stores all the information necessary to calculate arrival times, i.e. slownesses as well as source and receiver positions
    /*! As this class is derived from ThreeDModelBase the overall handling is similar to the other model classes. There are a few notable exceptions.
     * The model has to have equal cell sizes in all three directions, therefore there is only a single function to set the size and
     * number of cells for all axes.
     *
     * We also need to store the position of the seismic sources. The receiver positions are the measurement positions of the base class and
     * we have additional functions for the receivers through the GetSourcePosX etc. functions. For efficient storage we only store each
     * receiver and source position once even though we usually have a multitude of source-receiver combinations that form our data.
     *
     * We therefore have an additional function AddMeasurementConfiguration that takes the index of the source and receiver for an individual measurement.
     *
     */
    class J3DEXPORT ThreeDSeismicModel: public jif3D::ThreeDModelBase
      {
    public:
      typedef TomoExtraParameterSetter ExtraParameterSetter;

    public:
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDModelBase>(*this);

        }
      //! Given three coordinates in m, find the indices of the model cell that correponds to these coordinates, this is a more efficient implementation than the one in the base class
      virtual boost::array<ThreeDModelBase::t3DModelData::index, 3>
      FindAssociatedIndices(const double xcoord, const double ycoord,
          const double zcoord) const override;
      //! The seismic model has the same cell size for all cells in all directions so we just have one function to set it
      /*! This function sets both the size of all cells as well as the number of cells in x,y and z-direction
       * @param Size The size of each cell in all directions in m
       * @param nx The number of cells in x-direction (North)
       * @param ny The number of cells in y-direction (East)
       * @param nz The number of cells in z-direction
       */
      void SetCellSize(const double Size, const size_t nx, const size_t ny,
          const size_t nz)
        {
          ThreeDModelBase::t3DModelDim XCS(nx,Size);
          ThreeDModelBase::t3DModelDim YCS(ny,Size);
          ThreeDModelBase::t3DModelDim ZCS(nz,Size);
          SetXCellSizes(XCS);
          SetYCellSizes(YCS);
          SetZCellSizes(ZCS);
          ThreeDModelBase::SetData().resize(boost::extents[nx][ny][nz]);
        }
      //! return read only access to the stored slowness values
      const t3DModelData &GetSlownesses() const
        {
          return ThreeDModelBase::GetData();
        }
      //! return a reference to stored slownesses
      t3DModelData &SetSlownesses()
        {
          return ThreeDModelBase::SetData();
        }

      //! Write the seismic model in VTK format, at the moment the best format for plotting, this only writes the slowness and not the source and receiver positions
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, "Slowness");
        }
      //! Write the seimic model and all associated information in a netcdf file
      virtual void WriteNetCDF(const std::string &filename) const override;
      //! Read the seismic model and all associated information from a netcdf file, this version does not check the validity of the grid
      virtual void ReadNetCDF(const std::string &filename) override;
      //! Read the seismic model and all associated information from a netcdf file, optionally check if grid is valid for tomography
      /*! A valid grid for tomography has equal cell sizes in all directions. However,
       * we reuse grid files for other purposes to make it easier to handle, so this
       * routine can optionally check if all cells have the same size in all directions
       * @param filename The name of the netcdf file
       * @param checkgrid Perform check that all cells have the same size
       */
      virtual void ReadNetCDF(const std::string &filename, bool checkgrid);
      ThreeDSeismicModel();
      ThreeDSeismicModel& operator=(const ThreeDModelBase& source);
      virtual ~ThreeDSeismicModel();
      };

    class TomoExtraParameterSetter
      {
    public:
      void operator()(ThreeDSeismicModel &Model, TomographyData &Data, const std::vector<double> &Dens)
        {
          throw jif3D::FatalException(
              "Tomography class does not support extra parameters !", __FILE__, __LINE__);
        }
      };
  /* @} */
  }

#endif /* THREEDSEISMICMODEL_H_ */
