//============================================================================
// Name        : ThreeDMagneticModel.h
// Author      : 11 Oct 2013
// Version     :
// Copyright   : 2013, mm489
//============================================================================

#ifndef THREEDMAGNETICMODEL_H_
#define THREEDMAGNETICMODEL_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include "../Global/FatalException.h"
#include "../ModelBase/ThreeDModelBase.h"
#include "../ModelBase/GridOnlyModelCache.h"

namespace jif3D
  {

    class MagneticsExtraParameterSetter;
    class MagneticData;

    static const std::string SusceptibilityName = "Susceptibility";
    static const std::string SusceptibilityUnit = " ";
    class J3DEXPORT ThreeDMagneticModel: public jif3D::ThreeDModelBase
      {
    public:
      typedef jif3D::GridOnlyModelCache<ThreeDMagneticModel> ModelCacheType;
      typedef MagneticsExtraParameterSetter ExtraParameterSetter;

      jif3D::rvec GetModelParameters() const
        {
          jif3D::rvec parms(GetData().num_elements());
          std::copy(GetData().origin(), GetData().origin() + GetData().num_elements(),
              parms.begin());
          return parms;
        }
      //! The implementation class for magnetic and gravity data needs to know the number of model parameters
      size_t GetNModelParm() const
        {
          return GetData().num_elements();
        }
      ThreeDMagneticModel();
      //! Other models will be copied by the copy operator for the base class
      ThreeDMagneticModel& operator=(const ThreeDModelBase& source);
      virtual ~ThreeDMagneticModel();
      //! Set the sizes of the grid cells in x-direction in m
      void SetXCellSizes(ThreeDModelBase::t3DModelDim &XCS)
        {
          ThreeDModelBase::SetXCellSizes(XCS);
        }
      //! Set the sizes of the grid cells in y-direction in m
      void SetYCellSizes(ThreeDModelBase::t3DModelDim &YCS)
        {
          ThreeDModelBase::SetYCellSizes(YCS);
        }
      //! Set the sizes of the grid cells in z-direction in m
      void SetZCellSizes(ThreeDModelBase::t3DModelDim &ZCS)
        {
          ThreeDModelBase::SetZCellSizes(ZCS);
        }
      //! return read only access to the stored susceptibility values
      const t3DModelData &GetSusceptibilities() const
        {
          return ThreeDModelBase::GetData();
        }
      //! return a reference to stored susceptibilities
      t3DModelData &SetSusceptibilities()
        {
          return ThreeDModelBase::SetData();
        }

      //! Write the density model and all associated information in a netcdf file
      virtual void WriteNetCDF(const std::string &filename) const override
        {
          netCDF::NcFile DataFile(filename, netCDF::NcFile::replace);
          //first write the 3D discretized part
          WriteDataToNetCDF(DataFile, SusceptibilityName, SusceptibilityUnit);
        }

      //! Write the density model in VTK format, at the moment the best format for plotting
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, SusceptibilityName);
        }

      //! Read the density model and all associated information from a netcdf file
      virtual void ReadNetCDF(const std::string &filename) override
        {
          //create the netcdf file object
          netCDF::NcFile DataFile(filename, netCDF::NcFile::read);
          //read in the 3D gridded data
          ReadDataFromNetCDF(DataFile, SusceptibilityName, SusceptibilityUnit);
        }

    private:
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDModelBase>(*this);
        }
      };

    class J3DEXPORT MagneticsExtraParameterSetter
      {
    public:
      void operator()(ThreeDMagneticModel &Model, MagneticData &Data, const std::vector<double> &Dens)
        {
          throw jif3D::FatalException(
              "Magnetics class does not support extra parameters !", __FILE__, __LINE__);
        }
      };
  } /* namespace jif3D */
#endif /* THREEDMAGNETICMODEL_H_ */
