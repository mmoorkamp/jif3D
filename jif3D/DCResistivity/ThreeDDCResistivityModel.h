//============================================================================
// Name        : ThreeDDCResistivityModel.h
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#ifndef THREEDDCRESISTIVITYMODEL_H_
#define THREEDDCRESISTIVITYMODEL_H_

#include "../Global/Serialization.h"
#include "../Global/FatalException.h"
#include "../Global/Jif3DGlobal.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
    /** \addtogroup DC Resistivity classes and functions
     * This module contains the functionality for DC resistivity forward modelling and data handling.
     * As for the other methods the 3D Models are stored in a class derived from ThreeDModelBase. The model parameter for these models is resistivity in ohm.m.
     * The DC resistivity model class stores all the information necessary to calculate DCresistivity response, i.e. resistivity as well as source and receiver positions
     * As this class is derived from ThreeDModelBase the overall handling is similar to the other model classes. There are a few notable exceptions.
     * The model has to have variable cell sizes in z direction, which is thinner for shallow layers but thicker for deeper layers.
     * Moreover, there are often two sources and two receivers for each datum, so we store the Positive and Negative source at the same time. And the other receiver
     * electrode need to be stored, too. We store the source position only once, but store all receiver positions and source index for each datum.
     */
    /* @{ */
    class DCExtraParameterSetter;

    class J3DEXPORT ThreeDDCResistivityModel: public jif3D::ThreeDModelBase
      {
    public:
      typedef DCExtraParameterSetter ExtraParameterSetter;
      //! We need vectors of indices to associate a source with a receiver for a given shot
      typedef std::vector<int> tIndexVec;
    private:
      //! The x-position of the positive-sources
      ThreeDModelBase::tMeasPosVec SourcePosPosX;
      //! The y-position of the positive-sources
      ThreeDModelBase::tMeasPosVec SourcePosPosY;
      //! The z-position of the positive-sources
      ThreeDModelBase::tMeasPosVec SourcePosPosZ;
      //! The x-position of the negative-sources
      ThreeDModelBase::tMeasPosVec SourceNegPosX;
      //! The y-position of the negative-sources
      ThreeDModelBase::tMeasPosVec SourceNegPosY;
      //! The z-position of the negative-sources
      ThreeDModelBase::tMeasPosVec SourceNegPosZ;
      //! The x-position of the second receivers
      ThreeDModelBase::tMeasPosVec MeasSecPosX;
      //! The y-position of the second receivers
      ThreeDModelBase::tMeasPosVec MeasSecPosY;
      //! The z-position of the second receivers
      ThreeDModelBase::tMeasPosVec MeasSecPosZ;
      //! Each source can correspond to a number of measurements with different receivers
      //! here we record for each datum the index of the source position in the above arrays
      tIndexVec SourceIndices;
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDModelBase>(*this);
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
    public:
      //! The DC Resistivity model for DCResistivityCalculator has the same cell size for all cells in each horizontal directions so we just have one function
      //! to set it
      /*! This function sets both the size of all cells as well as the number of cells in the horizontal (x and y) directions
       * @param XSize The size of each cell in x-direction (North) in m
       * @param YSize The size of each cell in y-direction (East) in m
       * @param nx The number of cells in x-direction (North)
       * @param ny The number of cells in y-direction (East)
       */
      void SetHorizontalCellSize(const double XSize, const double YSize, const size_t nx,
          const size_t ny)
        {
          ThreeDModelBase::t3DModelDim XCS(nx,XSize);
          ThreeDModelBase::t3DModelDim YCS(ny,YSize);
          ThreeDModelBase::SetXCellSizes(XCS);
          ThreeDModelBase::SetYCellSizes(YCS);
        }
      //! The vertical cells can all have different sizes so we allow direct access to the CellSize structure
      void SetZCellSizes(ThreeDModelBase::t3DModelDim &ZCS)
        {
          ThreeDModelBase::SetZCellSizes(ZCS);
        }
      //! return read only access to the stored resistivity values
      const t3DModelData &GetResistivities() const
        {
          return ThreeDModelBase::GetData();
        }
      //! return a reference to stored resistivities
      t3DModelData &SetResistivities()
        {
          return ThreeDModelBase::SetData();
        }
      //! Add a Positive source to the model
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
          ThreeDModelBase::AddMeasurementPoint(xcoord, ycoord, zcoord);
          MeasSecPosX.push_back(negxcoord);
          MeasSecPosY.push_back(negycoord);
          MeasSecPosZ.push_back(negzcoord);
          SourceIndices.push_back(sourceindex);
        }
      //! read only access to the x-positions of the Positive sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourcePosPosX() const
        {
          return SourcePosPosX;
        }
      //! read only access to the y-positions of the Positive sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourcePosPosY() const
        {
          return SourcePosPosY;
        }
      //! read only access to the z-positions of the Positive sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourcePosPosZ() const
        {
          return SourcePosPosZ;
        }
      //! read only access to the x-positions of the Negative sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourceNegPosX() const
        {
          return SourceNegPosX;
        }
      //! read only access to the y-positions of the Negative sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourceNegPosY() const
        {
          return SourceNegPosY;
        }
      //! read only access to the z-positions of the Negative sources in m
      const ThreeDModelBase::tMeasPosVec &GetSourceNegPosZ() const
        {
          return SourceNegPosZ;
        }
      //! read only access to the indices of the sources for each data point
      const tIndexVec &GetSourceIndices() const
        {
          return SourceIndices;
        }
      //! read only access to the x-positions of the second receiver electrode in m
      const ThreeDModelBase::tMeasPosVec &GetMeasSecPosX() const
        {
          return MeasSecPosX;
        }
      //! read only access to the y-positions of the second receiver electrode in m
      const ThreeDModelBase::tMeasPosVec &GetMeasSecPosY() const
        {
          return MeasSecPosY;
        }
      //! read only access to the z-positions of the second receiver electrode in m
      const ThreeDModelBase::tMeasPosVec &GetMeasSecPosZ() const
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
      void CopyMeasurementConfigurations(const ThreeDDCResistivityModel &Source)
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

    	  jif3D::ThreeDModelBase::ClearMeasurementPoints();
          const size_t nmeas = Source.GetMeasPosX().size();
          for (size_t i = 0; i < nmeas; ++i)
            {
              jif3D::ThreeDModelBase::AddMeasurementPoint(Source.GetMeasPosX()[i], Source.GetMeasPosY()[i],
                  Source.GetMeasPosZ()[i]);
            }
        }
      //! Write the resistivity model in VTK format, at the moment the best format for plotting, this only writes the resistivity and not the source and receiver positions
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, "Resistivity");
        }
      //! Write the resistivity model and all associated information in a netcdf file
      virtual void WriteNetCDF(const std::string &filename) const override;
      //! Read the resistivity model and all associated information from a netcdf file
      virtual void ReadNetCDF(const std::string &filename) override;
      ThreeDDCResistivityModel();
      //! We define our own copy constructor
      ThreeDDCResistivityModel(const ThreeDDCResistivityModel &source);
      //! We define our own copy operator
      ThreeDDCResistivityModel& operator=(const ThreeDDCResistivityModel& source);
      //! Other models will be copied by the copy operator for the base class
      ThreeDDCResistivityModel& operator=(const ThreeDModelBase& source);
      virtual ~ThreeDDCResistivityModel();
      };

    class J3DEXPORT DCExtraParameterSetter
      {
    public:
      void operator()(ThreeDDCResistivityModel &Model, const std::vector<double> &Dens)
        {
          throw jif3D::FatalException(
              "Tomography class does not support extra parameters !");
        }
      };
  /* @} */
  }

#endif /* THREEDDCRESISTIVITYMODEL_H_ */
