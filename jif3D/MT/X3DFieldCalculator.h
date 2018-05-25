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
//#include "X3DTypes.h"
#include "ReadWriteX3D.h"

namespace jif3D
  {

    class X3DFieldCalculator
      {
    private:
      std::vector<std::complex<double> > Ex1, Ex2, Ey1, Ey2, Hx1, Hx2, Hy1, Hy2, Hz1, Hz2;
      boost::filesystem::path TempDir;
      std::string X3DName;
      std::string NameRoot;
      bool CleanFiles;
      jif3D::GreenCalcType GreenStage1;
      jif3D::GreenCalcType GreenStage4;
      jif3D::X3DModel OldModel;
      std::string ForwardDirName;
      std::string ObjectID()
        {
          //a unique ID created on construction
          boost::uuids::uuid tag = boost::uuids::random_generator()();
          //make a unique filename for the sensitivity file created by this object
          //we use boost uuid to generate a unique identifier tag
          //and translate it to a string to generate the filename
          return "mt" + jif3D::stringify(jif3D::platform::get_process_id()) + "x"
              + jif3D::stringify(this) + "t" + jif3D::stringify(tag);
        }
    public:
      const std::vector<std::complex<double> > &GetEx1() const
        {
          return Ex1;
        }
      const std::vector<std::complex<double> > &GetEx2() const
        {
          return Ex2;
        }
      const std::vector<std::complex<double> > &GetEy1() const
        {
          return Ey1;
        }
      const std::vector<std::complex<double> > &GetEy2() const
        {
          return Ey2;
        }
      const std::vector<std::complex<double> > &GetHx1() const
        {
          return Hx1;
        }
      const std::vector<std::complex<double> > &GetHx2() const
        {
          return Hx2;
        }
      const std::vector<std::complex<double> > &GetHy1() const
        {
          return Hy1;
        }
      const std::vector<std::complex<double> > &GetHy2() const
        {
          return Hy2;
        }
      const std::vector<std::complex<double> > &GetHz1() const
        {
          return Hz1;
        }
      const std::vector<std::complex<double> > &GetHz2() const
        {
          return Hz2;
        }
      std::string GetForwardDirName() const
        {
          return ForwardDirName;
        }
      void CalculateFields(const X3DModel &Model, size_t freqindex);
      X3DFieldCalculator(boost::filesystem::path TDir = boost::filesystem::current_path(),
          std::string x3d = "x3d", bool Clean = true, jif3D::GreenCalcType GS1 = hst,
          jif3D::GreenCalcType GS4 = hst);
      virtual ~X3DFieldCalculator();
      };

  } /* namespace jif3D */

#endif /* MT_X3DFIELDCALCULATOR_H_ */
