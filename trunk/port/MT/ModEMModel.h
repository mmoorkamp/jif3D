//============================================================================
// Name        : ModEMModel.h
// Author      : 14 Jan 2016
// Version     :
// Copyright   : 2016, mm489
//============================================================================

#ifndef MT_MODEMMODEL_H_
#define MT_MODEMMODEL_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "ThreeDMTModel.h"

namespace jif3D
  {

    class J3DEXPORT ModEMModel : public ThreeDMTModel
      {
    public:
      ModEMModel();
      virtual ~ModEMModel();

      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<ThreeDMTModel>(*this);
        }
      //! Write all model information to a netcdf file
      virtual void WriteNetCDF(const std::string &filename) const;
      //! Read all model information from a netcdf file
      virtual void ReadNetCDF(const std::string &filename);
      //! Other models will be copied by the copy operator for the base class
      ModEMModel& operator=(const ThreeDModelBase& source);
      };

  } /* namespace jif3D */

#endif /* MT_MODEMMODEL_H_ */
