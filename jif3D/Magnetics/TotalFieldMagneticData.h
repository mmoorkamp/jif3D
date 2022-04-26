/*
 * MagneticData.h
 *
 *  Created on: May 21, 2019
 *      Author: max
 */

#ifndef MAGNETICS_TOTALFIELDMAGNETICDATA_H_
#define MAGNETICS_TOTALFIELDMAGNETICDATA_H_

#include "../DataBase/GeneralData.h"
#include "ThreeDSusceptibilityModel.h"

namespace jif3D
  {

    class TotalFieldMagneticData: public GeneralData
      {
    public:
      typedef ThreeDSusceptibilityModel ModelType;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralData>(*this);
        }
      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) const override;
      void WriteVTK(const std::string &filename) const;
      TotalFieldMagneticData();
      virtual ~TotalFieldMagneticData();
      };

  } /* namespace jif3D */

#endif /* MAGNETICS_TOTALFIELDMAGNETICDATA_H_ */
