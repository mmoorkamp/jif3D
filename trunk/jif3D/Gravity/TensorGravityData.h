/*
 * TensorGravityData.h
 *
 *  Created on: May 20, 2019
 *      Author: max
 */

#ifndef GRAVITY_TENSORGRAVITYDATA_H_
#define GRAVITY_TENSORGRAVITYDATA_H_

#include "../Global/Serialization.h"
#include "../DataBase/GeneralData.h"
#include "ThreeDGravityModel.h"


namespace jif3D
  {

    class TensorGravityData: public GeneralData
      {
    public:
      typedef ThreeDGravityModel ModelType;

      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & base_object<GeneralData>(*this);
        }
      virtual void ReadNetCDF(const std::string &filename) override;
      virtual void WriteNetCDF(const std::string &filename) override;
      TensorGravityData();
      virtual ~TensorGravityData();
      };

  } /* namespace jif3D */

#endif /* GRAVITY_TENSORGRAVITYDATA_H_ */
