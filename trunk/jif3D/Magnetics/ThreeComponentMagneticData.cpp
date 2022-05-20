/*
 * ThreeComponentMagneticData.cpp
 *
 *  Created on: Apr 26, 2022
 *      Author: max
 */

#include "ThreeComponentMagneticData.h"
#include "ReadWriteMagneticData.h"
#include "../ModelBase/VTKTools.h"

namespace jif3D
  {

    void ThreeComponentMagneticData::ReadNetCDF(const std::string &filename)
      {
        std::vector<double> PosX, PosY, PosZ;
        std::vector<double> Data, Error;

        ReadMagneticComponentMeasurements(filename, Data, PosX, PosY, PosZ, Error);
        const size_t nmeas = PosX.size();
        assert(nmeas == PosY.size());
        assert(nmeas == PosZ.size());
        assert(nmeas * 3 == Data.size());
        assert(nmeas * 3 == Error.size());
        for (size_t i = 0; i < nmeas; ++i)
          {
            AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
          }
        SetDataAndErrors(Data, Error);
      }
    void ThreeComponentMagneticData::WriteNetCDF(const std::string &filename) const
      {
        SaveMagneticComponentMeasurements(filename, GetData(), GetMeasPosX(), GetMeasPosY(),
            GetMeasPosZ(), GetErrors());
      }
    void ThreeComponentMagneticData::WriteVTK(const std::string &filename) const
      {
        jif3D::Write3DVectorDataToVTK(filename + ".vtk", "B", GetData(), GetMeasPosX(),
            GetMeasPosY(), GetMeasPosZ());
      }

    ThreeComponentMagneticData::ThreeComponentMagneticData()
      {
        // TODO Auto-generated constructor stub

      }

    ThreeComponentMagneticData::~ThreeComponentMagneticData()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
