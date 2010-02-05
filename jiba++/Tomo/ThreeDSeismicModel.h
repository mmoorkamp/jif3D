//============================================================================
// Name        : ThreeDSeismicModel.h
// Author      : Apr 7, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef THREEDSEISMICMODEL_H_
#define THREEDSEISMICMODEL_H_

#include "../ModelBase/ThreeDModelBase.h"

namespace jiba
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
    class ThreeDSeismicModel: public jiba::ThreeDModelBase
      {
    public:
      //! We need vectors of indices to associate a source with a receiver for a given shot
      typedef std::vector<int> tIndexVec;
    private:
      //The x-position of the sources
      ThreeDModelBase::tMeasPosVec SourcePosX;
      //The y-position of the sources
      ThreeDModelBase::tMeasPosVec SourcePosY;
      //The z-position of the sources
      ThreeDModelBase::tMeasPosVec SourcePosZ;
      //Each source can can correspond to a number of measurements with different receivers
      //here we record for each datum the index of the source position in the above arrays
      tIndexVec SourceIndices;
      //and similarly for the receivers in the MeasPos* arrays
      tIndexVec ReceiverIndices;
    public:
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
          ThreeDModelBase::SetXCellSizes().resize(boost::extents[nx]);
          std::fill_n(ThreeDModelBase::SetXCellSizes().begin(), nx, Size);
          ThreeDModelBase::SetYCellSizes().resize(boost::extents[ny]);
          std::fill_n(ThreeDModelBase::SetYCellSizes().begin(), ny, Size);
          ThreeDModelBase::SetZCellSizes().resize(boost::extents[nz]);
          std::fill_n(ThreeDModelBase::SetZCellSizes().begin(), nz, Size);
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
      //! Add a source to the model
      void AddSource(const double xcoord, const double ycoord,
          const double zcoord)
        {
          SourcePosX.push_back(xcoord + XOrigin);
          SourcePosY.push_back(ycoord + YOrigin);
          SourcePosZ.push_back(zcoord + ZOrigin);
        }
      //! read only access to the x-positions of the sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourcePosX() const
        {
          return SourcePosX;
        }
      //! read only access to the y-positions of the sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourcePosY() const
        {
          return SourcePosY;
        }
      //! read only access to the z-positions of the sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourcePosZ() const
        {
          return SourcePosZ;
        }
      //! read only access to the indices of the sources for each data point
      const tIndexVec &GetSourceIndices() const
        {
          return SourceIndices;
        }
      //! read only access to the indices of the receivers for each data point
      const tIndexVec &GetReceiverIndices() const
        {
          return ReceiverIndices;
        }
      //! add a source-receiver combination for which a travel time will be calculated
      /*! We only store the positions of the sources and receivers once. With this
       * function we can add a configuration for which we want to calculate a travel time.
       * @param SourceIndex The index of the source in the vector of positions
       * @param ReceiverIndex The index of the receiver in the vector of positions
       */
      void AddMeasurementConfiguration(const size_t SourceIndex,
          const size_t ReceiverIndex)
        {
          assert(SourceIndex < SourcePosX.size());
          assert(ReceiverIndex < GetMeasPosX().size());
          SourceIndices.push_back(SourceIndex);
          ReceiverIndices.push_back(ReceiverIndex);
        }
      //! Remove all the index values that determine which source-receiver combinations are active
      void ClearMeasurementConfigurations()
        {
          SourceIndices.clear();
          ReceiverIndices.clear();
        }
      //! Clear all the source position definitions including the source-receiver combinations
      void ClearSourcePos()
        {
          SourcePosX.clear();
          SourcePosY.clear();
          SourcePosZ.clear();
          ClearMeasurementConfigurations();
        }
      //! Copy the source and receiver positions and the indices for the source-receiver combinations from the source model
      void CopyMeasurementConfigurations(const ThreeDSeismicModel &Source)
        {

          SourcePosX.resize(Source.SourcePosX.size());
          std::copy(Source.SourcePosX.begin(), Source.SourcePosX.end(),
              SourcePosX.begin());
          SourcePosY.resize(Source.SourcePosY.size());
          std::copy(Source.SourcePosY.begin(), Source.SourcePosY.end(),
              SourcePosY.begin());
          SourcePosZ.resize(Source.SourcePosZ.size());
          std::copy(Source.SourcePosZ.begin(), Source.SourcePosZ.end(),
              SourcePosZ.begin());
          SourceIndices.resize(Source.SourceIndices.size());
          std::copy(Source.SourceIndices.begin(), Source.SourceIndices.end(),
              SourceIndices.begin());
          ReceiverIndices.resize(Source.ReceiverIndices.size());
          std::copy(Source.ReceiverIndices.begin(),
              Source.ReceiverIndices.end(), ReceiverIndices.begin());
        }
      //! Set the origin of the coordinate system, this is a reimplementation from the base class to also change the source positions
      virtual void SetOrigin(const double x, const double y, const double z);
      //! Write the seismic model in VTK format, at the moment the best format for plotting, this only writes the slowness and not the source and receiver positions
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, "Slowness");
        }
      //! Write the seimic model and all associated information in a netcdf file
      void WriteNetCDF(const std::string filename) const;
      //! Read the seismic model and all associated information from a netcdf file
      void ReadNetCDF(const std::string filename);
      ThreeDSeismicModel();
      virtual ~ThreeDSeismicModel();
      };
  /* @} */
  }

#endif /* THREEDSEISMICMODEL_H_ */
