/*
 * TensorGravityData.cpp
 *
 *  Created on: May 20, 2019
 *      Author: max
 */

#include "TensorGravityData.h"
#include "ReadWriteGravityData.h"


namespace jif3D
  {

    void TensorGravityData::ReadNetCDF(const std::string &filename)
      {
        std::vector<double> PosX, PosY, PosZ, Data, Error;

        ReadTensorGravityMeasurements(filename, Data, PosX, PosY, PosZ, Error);
        const size_t nmeas = PosX.size();
        assert(nmeas == PosY.size());
        assert(nmeas == PosZ.size());
        assert(nmeas * 9 == Data.size());
        assert(nmeas * 9 == Error.size());
        for (size_t i = 0; i < nmeas; ++i)
          {
            AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
          }
        SetDataAndErrors(Data, Error);
      }

    void TensorGravityData::WriteNetCDF(const std::string &filename)
      {
        SaveTensorGravityMeasurements(filename, GetData(), GetMeasPosX(), GetMeasPosY(),
            GetMeasPosZ(), GetErrors());

      }
    TensorGravityData::TensorGravityData()
      {
        // TODO Auto-generated constructor stub

      }

    TensorGravityData::~TensorGravityData()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
