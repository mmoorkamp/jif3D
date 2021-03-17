/*
 * MagneticData.h
 *
 *  Created on: May 21, 2019
 *      Author: max
 */

#ifndef MAGNETICS_MAGNETICDATA_H_
#define MAGNETICS_MAGNETICDATA_H_

#include "ThreeDMagneticModel.h"
#include "../DataBase/GeneralData.h"

namespace jif3D
  {

    class MagneticData: public GeneralData
      {
    public:
      typedef ThreeDMagneticModel ModelType;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralData>(*this);
        }
      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) const override;
      void WriteVTK(const std::string &filename) const;
      MagneticData();
      virtual ~MagneticData();
      };

  } /* namespace jif3D */

#endif /* MAGNETICS_MAGNETICDATA_H_ */
