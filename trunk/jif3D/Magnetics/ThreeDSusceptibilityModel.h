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
#include "../Global/NetCDFPortHelper.h"

#include "../ModelBase/ThreeDModelBase.h"
#include "../ModelBase/GridOnlyModelCache.h"
#include "../ModelBase/GridBackgroundModelCache.h"

namespace jif3D
  {

    class MagneticsExtraParameterSetter;
    class TotalFieldMagneticData;

    static const std::string SusceptibilityName = "Susceptibility";
    static const std::string SusceptibilityUnit = " ";
    class J3DEXPORT ThreeDSusceptibilityModel: public jif3D::ThreeDModelBase
      {
    public:
      //! The type of the background thickness and susceptibility vector, this is a std::vector because we want to easily append elements
      typedef std::vector<double> tBackgroundVec;
      typedef jif3D::GridBackgroundModelCache<ThreeDSusceptibilityModel> ModelCacheType;
    private:
      //! The susceptibilities of the background layers
      tBackgroundVec bg_susceptibilities;
      //! The thicknesses of the background layers
      tBackgroundVec bg_thicknesses;
    public:
      typedef MagneticsExtraParameterSetter ExtraParameterSetter;
      jif3D::rvec GetModelParameters() const
        {
          jif3D::rvec parms(GetData().num_elements() + bg_susceptibilities.size());
          std::copy(GetData().origin(), GetData().origin() + GetData().num_elements(),
              parms.begin());
          std::copy(bg_susceptibilities.begin(), bg_susceptibilities.end(),
              parms.begin() + GetData().num_elements());
          return parms;
        }
      //! The implementation class for magnetic and gravity data needs to know the number of model parameters
      size_t GetNModelParm() const
        {
          return GetData().num_elements() + bg_susceptibilities.size();
        }
      ThreeDSusceptibilityModel();
      //! Other models will be copied by the copy operator for the base class
      ThreeDSusceptibilityModel& operator=(const ThreeDModelBase &source);
      virtual ~ThreeDSusceptibilityModel();
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
      const t3DModelData& GetSusceptibilities() const
        {
          return ThreeDModelBase::GetData();
        }
      //! return a reference to stored susceptibilities
      t3DModelData& SetSusceptibilities()
        {
          return ThreeDModelBase::SetData();
        }
      //! Set the susceptibility of the background, it extends to infinity in horizontal directions and to the depth specified by the thicknesses in vertical direction
      void SetBackgroundSusceptibilities(const tBackgroundVec &value)
        {
          bg_susceptibilities.clear();
          copy(value.begin(), value.end(), back_inserter(bg_susceptibilities));
        }
      //! Return the susceptibility of the background layers
      const tBackgroundVec& GetBackgroundSusceptibilities() const
        {
          return bg_susceptibilities;
        }
      const tBackgroundVec& GetBackgroundValues() const
        {
          return bg_susceptibilities;
        }
      //! Set the thicknesses of the background layers, the individual thicknesses are given in m the total thickness of the background layers does not need to coincide with the gridded domain
      void SetBackgroundThicknesses(const tBackgroundVec &value)
        {
          bg_thicknesses.clear();
          copy(value.begin(), value.end(), back_inserter(bg_thicknesses));
        }
      //! Return the thicknesses of the background layers in m
      const tBackgroundVec& GetBackgroundThicknesses() const
        {
          return bg_thicknesses;
        }
      //! Write the susceptibility model and all associated information in a netcdf file
      virtual void WriteNetCDF(const std::string &filename) const override
        {
          netCDF::NcFile DataFile(filename, netCDF::NcFile::replace);
          //first write the 3D discretized part
          WriteDataToNetCDF(DataFile, SusceptibilityName, SusceptibilityUnit);

          if (!bg_susceptibilities.empty())
            {
              assert(bg_susceptibilities.size() == bg_thicknesses.size());
              //Create the matching dimension for the background layers
              netCDF::NcDim BackgroundDim = DataFile.addDim("bg_layers",
                  bg_thicknesses.size());
              //now we can write the actual parameters for the layers
              netCDF::NcVar bgSusVar = DataFile.addVar("bg_susceptibilities",
                  netCDF::ncDouble, BackgroundDim);
              bgSusVar.putAtt("long_name", "Background Susceptibilities");
              bgSusVar.putAtt("units", SusceptibilityUnit);
              jif3D::cxxport::put_legacy_ncvar(bgSusVar, bg_susceptibilities.data(),
                  bg_susceptibilities.size());

              netCDF::NcVar bgThickVar = DataFile.addVar("bg_thicknesses",
                  netCDF::ncDouble, BackgroundDim);
              bgThickVar.putAtt("long_name", "Background Thicknesses");
              bgThickVar.putAtt("units", "m");
              jif3D::cxxport::put_legacy_ncvar(bgThickVar, bg_thicknesses.data(),
                  bg_thicknesses.size());
              //            bgThickVar.put(&bg_thicknesses[0], bg_thicknesses.size()); // old
            }

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

          try
            {
              //try to read the background data
              netCDF::NcDim BackgroundDim = DataFile.getDim("bg_layers");
              //if we succeed
              if (!BackgroundDim.isNull())
                {
                  //find out how many layers the background has
                  const size_t nbglayers = BackgroundDim.getSize();
                  //allocate memory
                  bg_susceptibilities.assign(nbglayers, 0.0);
                  bg_thicknesses.assign(nbglayers, 0.0);
                  //and read in from the file
                  netCDF::NcVar bgSusVar = DataFile.getVar("bg_susceptibilities");
                  netCDF::NcVar bgThickVar = DataFile.getVar("bg_thicknesses");

                  jif3D::cxxport::get_legacy_ncvar(bgSusVar, bg_susceptibilities.data(),
                      nbglayers);
                  jif3D::cxxport::get_legacy_ncvar(bgThickVar, bg_thicknesses.data(),
                      nbglayers);
                }
            } catch (const netCDF::exceptions::NcException &ex)
            {
              /*
               * Ignore exception. Before porting, the legacy NcError class was used to temporarily
               * set error handling to "silent_nonfatal".
               */
            }
        }

    private:
      friend class access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & base_object<ThreeDModelBase>(*this);
          ar & bg_thicknesses;
          ar & bg_susceptibilities;
        }
      };

    class J3DEXPORT MagneticsExtraParameterSetter
      {
    public:
      void operator()(ThreeDSusceptibilityModel &Model, TotalFieldMagneticData &Data,
          const std::vector<double> &Sus)
        {
          Model.SetBackgroundSusceptibilities(Sus);
        }
      };
  } /* namespace jif3D */
#endif /* THREEDMAGNETICMODEL_H_ */
