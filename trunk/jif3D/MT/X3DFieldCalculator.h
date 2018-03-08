/*
 * X3DFieldCalculator.h
 *
 *  Created on: 26 Feb 2018
 *      Author: mm489
 */

#ifndef MT_X3DFIELDCALCULATOR_H_
#define MT_X3DFIELDCALCULATOR_H_

#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/VectorTransform.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/convert.h"
#include "../Global/Jif3DPlatformHelper.h"

#include <limits>
#include <vector>
#include <utility>

#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include "X3DModel.h"
#include "X3DTypes.h"

namespace jif3D
  {

    class X3DFieldCalculator
      {
    private:
      jif3D::cmat Ex1_obs, Ex2_obs, Ey1_obs, Ey2_obs, Hx1_obs, Hx2_obs, Hy1_obs, Hy2_obs;
      jif3D::cmat Ex1_all, Ex2_all, Ey1_all, Ey2_all, Ez1_all, Ez2_all;
      bool HaveForwardFields;
    public:
      void CalculateFields(const X3DModel &Model, size_t minfreqindex,
        size_t maxfreqindex);
      X3DFieldCalculator(boost::filesystem::path TDir = boost::filesystem::current_path(),
          std::string x3d = "x3d", bool Clean = true);
      virtual ~X3DFieldCalculator();
      };

  } /* namespace jif3D */

#endif /* MT_X3DFIELDCALCULATOR_H_ */
