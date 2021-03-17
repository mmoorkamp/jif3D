/*
 * GeneralData.cpp
 *
 *  Created on: May 20, 2019
 *      Author: max
 */

#include "GeneralData.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/VTKTools.h"
#include <fstream>
#include <cassert>

namespace jif3D
  {

    GeneralData::GeneralData()
      {
        // TODO Auto-generated constructor stub

      }

    GeneralData::~GeneralData()
      {
        // TODO Auto-generated destructor stub
      }

    void GeneralData::WriteMeasurementPoints(const std::string &filename) const
      {
        std::vector<double> RecNum(GetMeasPosX().size());
        std::iota(RecNum.begin(), RecNum.end(), 1);
        jif3D::Write3DDataToVTK(filename, "Receiver", RecNum, GetMeasPosX(),
            GetMeasPosY(), GetMeasPosZ());
      }

//    void GeneralData::ReadMeasPosNetCDF(const std::string filename)
//      {
//        jif3D::ReadMeasPosNetCDF(filename, MeasPosX, MeasPosY, MeasPosZ);
//      }
//
//    void GeneralData::ReadMeasPosAscii(const std::string filename)
//      {
//        //we assume that the measurement positions are simply
//        //in the format x,y,z
//        std::ifstream infile(filename.c_str());
//        double posx, posy, posz;
//        while (infile.good())
//          {
//            infile >> posx >> posy >> posz;
//            if (infile.good())
//              {
//                MeasPosX.push_back(posx);
//                MeasPosY.push_back(posy);
//                MeasPosZ.push_back(posz);
//              }
//          }
//        assert(MeasPosX.size() == MeasPosY.size());
//        assert(MeasPosX.size() == MeasPosZ.size());
//      }
  }
/* namespace jif3D */
