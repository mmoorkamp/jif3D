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

    class ThreeDSeismicModel: public jiba::ThreeDModelBase
      {
    private:
      //The x-position of the sources
      ThreeDModelBase::tMeasPosVec SourcePosX;
      //The y-position of the sources
      ThreeDModelBase::tMeasPosVec SourcePosY;
      //The z-position of the sources
      ThreeDModelBase::tMeasPosVec SourcePosZ;
      //Each source can can correspond to a number of measurements with different receivers
      //here we record for each datum the index of the source position in the above arrays
      std::vector<size_t> SourceIndices;
      //and similarly for the receivers in the MeasPos* arrays
      std::vector<size_t> ReceiverIndices;
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
      //! Set the origin of the coordinate system, this is a reimplementation from the base class to also change the source positions
      virtual void SetOrigin(const double x, const double y, const double z);
      //! Write the seismic model in VTK format, at the moment the best format for plotting, this only writes the slowness and not the source and receiver positions
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, "Slowness");
        }
      //! Write the seimic model and all associated information in a netcdf file
      void WriteNetCDF(const std::string filename) const;
      //! Read the density model and all associated information from a netcdf file
      void ReadNetCDF(const std::string filename);
      ThreeDSeismicModel();
      virtual ~ThreeDSeismicModel();
      };

  }

#endif /* THREEDSEISMICMODEL_H_ */
