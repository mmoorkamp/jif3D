//============================================================================
// Name        : ThreeDMagneticModel.h
// Author      : 11 Oct 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#ifndef THREEDMAGNETICMODEL_H_
#define THREEDMAGNETICMODEL_H_

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {

    static const std::string SusceptibilityName = "Susceptibility";
    static const std::string SusceptibilityUnit = " ";
    class ThreeDMagneticModel: public jif3D::ThreeDModelBase
      {
    private:
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & boost::serialization::base_object<ThreeDModelBase>(*this);
        }
    public:
      //! The implementation class for magnetic and gravity data needs to know the number of model parameters
      size_t GetNModelParm() const
      {
    	  return GetData().num_elements();
      }
      ThreeDMagneticModel();
      virtual ~ThreeDMagneticModel();
      //! Set the sizes of the grid cells in x-direction in m
      t3DModelDim &SetXCellSizes()
        {
          return ThreeDModelBase::SetXCellSizes();
        }
      //! Set the sizes of the grid cells in y-direction in m
      t3DModelDim &SetYCellSizes()
        {
          return ThreeDModelBase::SetYCellSizes();
        }
      //! Set the sizes of the grid cells in z-direction in m
      t3DModelDim &SetZCellSizes()
        {
          return ThreeDModelBase::SetZCellSizes();
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
      void WriteNetCDF(const std::string filename) const
        {
          NcFile DataFile(filename.c_str(), NcFile::Replace);
          //write the 3D discretized part
          WriteDataToNetCDF(DataFile, SusceptibilityName, SusceptibilityUnit);
        }
      //! Write the density model in VTK format, at the moment the best format for plotting
      void WriteVTK(const std::string filename) const
        {
          ThreeDModelBase::WriteVTK(filename, SusceptibilityName);
        }
      //! Read the density model and all associated information from a netcdf file
      void ReadNetCDF(const std::string filename)
        {
          //create the netcdf file object
          NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
          //read in the 3D gridded data
          ReadDataFromNetCDF(DataFile, SusceptibilityName, SusceptibilityUnit);
        }
      };

  } /* namespace jif3D */
#endif /* THREEDMAGNETICMODEL_H_ */
