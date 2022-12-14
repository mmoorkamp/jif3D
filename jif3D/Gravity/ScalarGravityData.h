/*
 * ScalarGravityData.h
 *
 *  Created on: May 20, 2019
 *      Author: max
 */

#ifndef GRAVITY_SCALARGRAVITYDATA_H_
#define GRAVITY_SCALARGRAVITYDATA_H_

#include "../Global/Serialization.h"
#include "../DataBase/GeneralData.h"
#include "ThreeDGravityModel.h"

namespace jif3D
  {

    class ScalarGravityData: public GeneralData
      {
    public:
      typedef ThreeDGravityModel ModelType;
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralData>(*this);
        }
      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) const override;
      void WriteVTK(const std::string &filename) const;
      ScalarGravityData();
      virtual ~ScalarGravityData();

      };

  } /* namespace jif3D */

#endif /* GRAVITY_SCALARGRAVITYDATA_H_ */
