/*
 * DCResistivityData.h
 *
 *  Created on: Nov 13, 2020
 *      Author: zhanjie
 */

#ifndef DC_DCRESISTIVITYDATA_H_
#define DC_DCRESISTIVITYDATA_H_

#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include "../Global/Serialization.h"
#include "../DataBase/GeneralData.h"

namespace jif3D
  {

    class DCResistivityData: public GeneralData
      {
    public:
      typedef ThreeDDCResistivityModel ModelType;

      //! We need vectors of indices to associate a source with a receiver for a given shot
      typedef std::vector<int> tIndexVec;
    private:
      //! The x-position of the positive-sources
      std::vector<double> SourcePosPosX;
      //! The y-position of the positive-sources
      std::vector<double> SourcePosPosY;
      //! The z-position of the positive-sources
      std::vector<double> SourcePosPosZ;
      //! The x-position of the negative-sources
      std::vector<double> SourceNegPosX;
      //! The y-position of the negative-sources
      std::vector<double> SourceNegPosY;
      //! The z-position of the negative-sources
      std::vector<double> SourceNegPosZ;
      //! The x-position of the second receivers
      std::vector<double> MeasSecPosX;
      //! The y-position of the second receivers
      std::vector<double> MeasSecPosY;
      //! The z-position of the second receivers
      std::vector<double> MeasSecPosZ;
      //! Each source can correspond to a number of measurements with different receivers
      //! here we record for each datum the index of the source position in the above arrays
      tIndexVec SourceIndices;

    public:
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
    	  ar & base_object<GeneralData>(*this);
          ar & SourcePosPosX;
          ar & SourcePosPosY;
          ar & SourcePosPosZ;
          ar & SourceNegPosX;
          ar & SourceNegPosY;
          ar & SourceNegPosZ;
          ar & MeasSecPosX;
          ar & MeasSecPosY;
          ar & MeasSecPosZ;
          ar & SourceIndices;
        }

      //! Add a pair of source to the model
      void AddSource(double xcoord, double ycoord, double zcoord, double negxcoord,
          double negycoord, double negzcoord)
        {
          SourcePosPosX.push_back(xcoord);
          SourcePosPosY.push_back(ycoord);
          SourcePosPosZ.push_back(zcoord);
          SourceNegPosX.push_back(negxcoord);
          SourceNegPosY.push_back(negycoord);
          SourceNegPosZ.push_back(negzcoord);
        }

      //! Add a the second receiver electrode to the model
      void AddMeasurementPoint(const double xcoord, const double ycoord,
          const double zcoord, double negxcoord, double negycoord, double negzcoord,
          int sourceindex)
        {
    	  GeneralData::AddMeasurementPoint(xcoord, ycoord, zcoord);
          MeasSecPosX.push_back(negxcoord);
          MeasSecPosY.push_back(negycoord);
          MeasSecPosZ.push_back(negzcoord);
          SourceIndices.push_back(sourceindex);
        }
      //! read only access to the x-positions of the Positive sources in m
      const std::vector<double> &GetSourcePosPosX() const
        {
          return SourcePosPosX;
        }
      //! read only access to the y-positions of the Positive sources in m
      const std::vector<double> &GetSourcePosPosY() const
        {
          return SourcePosPosY;
        }
      //! read only access to the z-positions of the Positive sources in m
      const std::vector<double> &GetSourcePosPosZ() const
        {
          return SourcePosPosZ;
        }
      //! read only access to the x-positions of the Negative sources in m
      const std::vector<double> &GetSourceNegPosX() const
        {
          return SourceNegPosX;
        }
      //! read only access to the y-positions of the Negative sources in m
      const std::vector<double> &GetSourceNegPosY() const
        {
          return SourceNegPosY;
        }
      //! read only access to the z-positions of the Negative sources in m
      const std::vector<double> &GetSourceNegPosZ() const
        {
          return SourceNegPosZ;
        }
      //! read only access to the indices of the sources for each data point
      const tIndexVec &GetSourceIndices() const
        {
          return SourceIndices;
        }
      //! read only access to the x-positions of the second receiver electrode in m
      const std::vector<double> &GetMeasSecPosX() const
        {
          return MeasSecPosX;
        }
      //! read only access to the y-positions of the second receiver electrode in m
      const std::vector<double> &GetMeasSecPosY() const
        {
          return MeasSecPosY;
        }
      //! read only access to the z-positions of the second receiver electrode in m
      const std::vector<double> &GetMeasSecPosZ() const
        {
          return MeasSecPosZ;
        }
      //! Clear all the source position definitions including the source-receiver combinations
      void ClearSourcePosPos()
        {
          SourcePosPosX.clear();
          SourcePosPosY.clear();
          SourcePosPosZ.clear();
        }
      void ClearSourceNegPos()
        {
          SourceNegPosX.clear();
          SourceNegPosY.clear();
          SourceNegPosZ.clear();
        }
      //! Clear the second receiver electrode position definitions including the source-receiver combinations
      void ClearSecMeasurementPoint()
        {
          MeasSecPosX.clear();
          MeasSecPosY.clear();
          MeasSecPosZ.clear();
        }
      //! Clear the sourceindics
      void ClearSourceIndices()
        {
          SourceIndices.clear();
        }

      //! Copy the source and receiver positions and the indices for the source-receiver combinations from the source model
      void CopyMeasurementConfigurations(const DCResistivityData &Source)
        {
    	  SourcePosPosX = Source.SourcePosPosX;
    	  SourcePosPosY = Source.SourcePosPosY;
    	  SourcePosPosZ = Source.SourcePosPosZ;
    	  SourceNegPosX = Source.SourceNegPosX;
    	  SourceNegPosY = Source.SourceNegPosY;
    	  SourceNegPosZ = Source.SourceNegPosZ;

    	  MeasSecPosX = Source.MeasSecPosX;
    	  MeasSecPosY = Source.MeasSecPosY;
    	  MeasSecPosZ = Source.MeasSecPosZ;
    	  SourceIndices = Source.SourceIndices;

    	  GeneralData::ClearMeasurementPoints();
          const size_t nmeas = Source.GetMeasPosX().size();
          for (size_t i = 0; i < nmeas; ++i)
            {
        	  GeneralData::AddMeasurementPoint(Source.GetMeasPosX()[i], Source.GetMeasPosY()[i],
                  Source.GetMeasPosZ()[i]);
            }
        }

      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) const override;
      DCResistivityData();
      virtual ~DCResistivityData();
      };

  /* @} */
  } /* namespace jif3D */

#endif /* DC_DCRESISTIVITYDATA_H_ */
