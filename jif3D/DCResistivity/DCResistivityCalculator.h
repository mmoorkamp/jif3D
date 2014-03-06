//============================================================================
// Name        : DCResistivityCalculator.h
// Author      : 6 Mar 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#ifndef DCRESISTIVITYCALCULATOR_H_
#define DCRESISTIVITYCALCULATOR_H_

#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
    /** \addtogroup dcresistivity DC Resistivity classes and functions */
    /* @{ */

    //! This class provides
    /*! This class calculates
     */
    class DCResistivityCalculator: public ThreeDModelBase
      {
    public:
      typedef std::vector<int> tIndexVec;
    private:
      //! The x-positions of the sources
      ThreeDModelBase::tMeasPosVec SourcePosX;
      //! The y-positions of the sources
      ThreeDModelBase::tMeasPosVec SourcePosY;
      //! The z-positions of the sources
      ThreeDModelBase::tMeasPosVec SourcePosZ;
      //! Each source can can correspond to a number of measurements with different receiver
      //! configurations here we record for each datum the index of the first injection point
      tIndexVec SourceIndices1;
      //! Each source can can correspond to a number of measurements with different receiver
      //! configurations here we record for each datum the index of the second injection point
      tIndexVec SourceIndices2;
      //! The index in the position vectors for the first measurement electrode
      tIndexVec ReceiverIndices1;
      //! The index in the position vectors for the second measurement electrode
      tIndexVec ReceiverIndices2;
    public:
      //! Add a source to the model
      void AddSource(const double xcoord, const double ycoord, const double zcoord)
        {
          SourcePosX.push_back(xcoord - XOrigin);
          SourcePosY.push_back(ycoord - YOrigin);
          SourcePosZ.push_back(zcoord - ZOrigin);
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
      //! read only access to the indices of the first source for each data point
      const tIndexVec &GetSourceIndices1() const
        {
          return SourceIndices1;
        }
      //! read only access to the indices of the first source for each data point
      const tIndexVec &GetSourceIndices2() const
        {
          return SourceIndices2;
        }
      //! read only access to the indices of the receivers for each data point
      const tIndexVec &GetReceiverIndices1() const
        {
          return ReceiverIndices1;
        }
      //! read only access to the indices of the receivers for each data point
      const tIndexVec &GetReceiverIndices2() const
        {
          return ReceiverIndices2;
        }
      //! add a source-receiver combination for which a travel time will be calculated
      /*! We only store the positions of the sources and receivers once. With this
       * function we can add a configuration for which we want to calculate a travel time.
       * @param SourceIndex The index of the source in the vector of positions
       * @param ReceiverIndex The index of the receiver in the vector of positions
       */
      void AddMeasurementConfiguration(size_t SourceIndex1, size_t SourceIndex2,
          size_t ReceiverIndex1, size_t ReceiverIndex2)
        {
          assert(SourceIndex1 < SourcePosX.size());
          assert(SourceIndex2 < SourcePosX.size());
          assert(ReceiverIndex1 < GetMeasPosX().size());
          assert(ReceiverIndex2 < GetMeasPosX().size());
          SourceIndices1.push_back(SourceIndex1);
          SourceIndices2.push_back(SourceIndex2);
          ReceiverIndices1.push_back(ReceiverIndex1);
          ReceiverIndices2.push_back(ReceiverIndex2);
        }
      //! Remove all the index values that determine which source-receiver combinations are active
      void ClearMeasurementConfigurations()
        {
          SourceIndices1.clear();
          SourceIndices2.clear();
          ReceiverIndices1.clear();
          ReceiverIndices2.clear();
        }
      //! Clear all the source position definitions including the source-receiver combinations
      void ClearSourcePos()
        {
          SourcePosX.clear();
          SourcePosY.clear();
          SourcePosZ.clear();
          ClearMeasurementConfigurations();
        }
      void SetHorizontalCellSize(const double Size, const size_t nx, const size_t ny)
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
      //! Write the seismic model in VTK format, at the moment the best format for plotting, this only writes the slowness and not the source and receiver positions
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, "Slowness");
        }
      //! Write the seimic model and all associated information in a netcdf file
      void WriteNetCDF(const std::string filename) const;
      //! Read the seismic model and all associated information from a netcdf file
      void ReadNetCDF(const std::string filename, bool checkgrid = true);

      DCResistivityCalculator();
      virtual ~DCResistivityCalculator();
      };
  /* @} */
  } /* namespace jif3D */

#endif /* DCRESISTIVITYCALCULATOR_H_ */
