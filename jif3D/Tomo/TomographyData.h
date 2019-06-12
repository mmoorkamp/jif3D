/*
 * TomographyData.h
 *
 *  Created on: May 21, 2019
 *      Author: max
 */

#ifndef TOMO_TOMOGRAPHYDATA_H_
#define TOMO_TOMOGRAPHYDATA_H_

#include "../Tomo/ThreeDSeismicModel.h"
#include "../Global/Serialization.h"
#include "../DataBase/GeneralData.h"
namespace jif3D
  {

    class TomographyData: public GeneralData
      {
    public:
      typedef ThreeDSeismicModel ModelType;

      //! We need vectors of indices to associate a source with a receiver for a given shot
      typedef std::vector<int> tIndexVec;
    private:
      //! The x-position of the sources
      std::vector<double> SourcePosX;
      //! The y-position of the sources
      std::vector<double> SourcePosY;
      //! The z-position of the sources
      std::vector<double> SourcePosZ;
      //! Each source can can correspond to a number of measurements with different receivers
      //! here we record for each datum the index of the source position in the above arrays
      tIndexVec SourceIndices;
      //! The index in the position vectors for the receivers
      tIndexVec ReceiverIndices;
    public:
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & base_object<GeneralData>(*this);
          ar & SourcePosX;
          ar & SourcePosY;
          ar & SourcePosZ;
          ar & SourceIndices;
          ar & ReceiverIndices;
        }
      //! Add a source to the model
      void AddSource(const double xcoord, const double ycoord, const double zcoord)
        {
          SourcePosX.push_back(xcoord);
          SourcePosY.push_back(ycoord);
          SourcePosZ.push_back(zcoord);
        }
      //! read only access to the x-positions of the sources in m
      const std::vector<double> &GetSourcePosX() const
        {
          return SourcePosX;
        }
      //! read only access to the y-positions of the sources in m
      const std::vector<double> &GetSourcePosY() const
        {
          return SourcePosY;
        }
      //! read only access to the z-positions of the sources in m
      const std::vector<double> &GetSourcePosZ() const
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
      void CopyMeasurementConfigurations(const TomographyData &Source)
        {

          SourcePosX = Source.SourcePosX;
          SourcePosY = Source.SourcePosY;
          SourcePosZ = Source.SourcePosZ;
          SourceIndices = Source.SourceIndices;

          ClearMeasurementPoints();
          const size_t nmeas = Source.GetMeasPosX().size();
          for (size_t i = 0; i < nmeas; ++i)
            {
              AddMeasurementPoint(Source.GetMeasPosX()[i], Source.GetMeasPosY()[i],
                  Source.GetMeasPosZ()[i]);
            }

          ReceiverIndices = Source.ReceiverIndices;
        }
      void SetSourcePoints(const std::vector<double> &xcoords,
          const std::vector<double> &ycoords, const std::vector<double> &zcoords)
        {
          const size_t nsource = xcoords.size();
          if (ycoords.size() != nsource)
            {
              throw jif3D::FatalException(
                  "Trying to set coordinates with different number of values", __FILE__,
                  __LINE__);
            }
          if (zcoords.size() != nsource)
            {
              throw jif3D::FatalException(
                  "Trying to set coordinates with different number of values", __FILE__,
                  __LINE__);
            }
          SourcePosX = xcoords;
          SourcePosY = ycoords;
          SourcePosZ = zcoords;
        }
      void WriteSourcePoints(const std::string &filename);
      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) override;
      TomographyData();
      virtual ~TomographyData();
      };

  } /* namespace jif3D */

#endif /* TOMO_TOMOGRAPHYDATA_H_ */
