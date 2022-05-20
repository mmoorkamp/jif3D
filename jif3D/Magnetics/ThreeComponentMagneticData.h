/*
 * ThreeComponentMagneticData.h
 *
 *  Created on: Apr 26, 2022
 *      Author: max
 */

#ifndef MAGNETICS_THREECOMPONENTMAGNETICDATA_H_
#define MAGNETICS_THREECOMPONENTMAGNETICDATA_H_

#include "../DataBase/GeneralData.h"
#include "ThreeDMagnetizationModel.h"

namespace jif3D
  {

    class ThreeComponentMagneticData: public GeneralData
      {
    public:

      typedef ThreeDMagnetizationModel ModelType;
      template<class Archive>
      void serialize(Archive &ar, const unsigned int version)
        {
          ar & base_object < GeneralData > (*this);
        }
      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) const override;
      void WriteVTK(const std::string &filename) const;
      ThreeComponentMagneticData();
      virtual ~ThreeComponentMagneticData();
      };

  } /* namespace jif3D */

#endif /* MAGNETICS_THREECOMPONENTMAGNETICDATA_H_ */
