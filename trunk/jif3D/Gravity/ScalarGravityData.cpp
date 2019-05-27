/*
 * ScalarGravityData.cpp
 *
 *  Created on: May 20, 2019
 *      Author: max
 */

#include "ScalarGravityData.h"
#include "ReadWriteGravityData.h"
#include <cassert>
namespace jif3D
  {

    void ScalarGravityData::ReadNetCDF(const std::string &filename)
      {
        std::vector<double> PosX, PosY, PosZ, Data, Error;

        ReadScalarGravityMeasurements(filename, Data, PosX, PosY, PosZ, Error);
        const size_t nmeas = PosX.size();
        assert(nmeas == PosY.size());
        assert(nmeas == PosZ.size());
        assert(nmeas == Data.size());
        assert(nmeas == Error.size());
        for (size_t i = 0; i < nmeas; ++i)
          {
            AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
          }
        SetDataAndErrors(Data, Error);
      }

    void ScalarGravityData::WriteNetCDF(const std::string &filename)
      {
        SaveScalarGravityMeasurements(filename, GetData(), GetMeasPosX(), GetMeasPosY(),
            GetMeasPosZ(), GetErrors());

      }

    ScalarGravityData::ScalarGravityData()
      {
        // TODO Auto-generated constructor stub

      }

    ScalarGravityData::~ScalarGravityData()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
