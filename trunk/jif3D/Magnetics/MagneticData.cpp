/*
 * MagneticData.cpp
 *
 *  Created on: May 21, 2019
 *      Author: max
 */

#include "MagneticData.h"
#include "ReadWriteMagneticData.h"

namespace jif3D
  {

    void MagneticData::ReadNetCDF(const std::string &filename)
      {
        std::vector<double> Data, Error, MeasX, MeasY, MeasZ;
        ReadTotalFieldMagneticMeasurements(filename, Data, MeasX, MeasY, MeasZ, Error);
        SetMeasurementPoints(MeasX, MeasY, MeasZ);
        SetDataAndErrors(Data, Error);

      }
    void MagneticData::WriteNetCDF(const std::string &filename)
      {
        SaveTotalFieldMagneticMeasurements(filename, GetData(), GetMeasPosX(),
            GetMeasPosY(), GetMeasPosZ(), GetErrors());
      }

    MagneticData::MagneticData()
      {
        // TODO Auto-generated constructor stub

      }

    MagneticData::~MagneticData()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
