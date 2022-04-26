/*
 * MagneticData.cpp
 *
 *  Created on: May 21, 2019
 *      Author: max
 */

#include "ReadWriteMagneticData.h"
#include "../ModelBase/VTKTools.h"
#include "TotalFieldMagneticData.h"


namespace jif3D
  {

    void TotalFieldMagneticData::WriteVTK(const std::string &filename) const
      {
        jif3D::Write3DDataToVTK(filename + ".vtk", "T", GetData(), GetMeasPosX(),
            GetMeasPosY(), GetMeasPosZ());
      }
    void TotalFieldMagneticData::ReadNetCDF(const std::string &filename)
      {
        std::vector<double> Data, Error, MeasX, MeasY, MeasZ;
        ReadTotalFieldMagneticMeasurements(filename, Data, MeasX, MeasY, MeasZ, Error);
        SetMeasurementPoints(MeasX, MeasY, MeasZ);
        SetDataAndErrors(Data, Error);

      }
    void TotalFieldMagneticData::WriteNetCDF(const std::string &filename) const
      {
        SaveTotalFieldMagneticMeasurements(filename, GetData(), GetMeasPosX(),
            GetMeasPosY(), GetMeasPosZ(), GetErrors());
      }

    TotalFieldMagneticData::TotalFieldMagneticData()
      {
        // TODO Auto-generated constructor stub

      }

    TotalFieldMagneticData::~TotalFieldMagneticData()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
