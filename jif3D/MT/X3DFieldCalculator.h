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
#include "../Global/Jif3DGlobal.h"
#include "../Global/convert.h"
#include "../Global/Jif3DPlatformHelper.h"

#include <vector>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/shared_ptr.hpp>
#include "X3DModel.h"
#include "MTData.h"
#include "ReadWriteX3D.h"
//#include "X3DTypes.h"

namespace jif3D
  {

    class X3DFieldCalculator
      {
    private:
      void CleanUp();
      std::vector<std::vector<std::complex<double>>> Ex1, Ex2, Ey1, Ey2, Hx1, Hx2, Hy1,
          Hy2, Hz1, Hz2;
      std::vector<std::vector<std::complex<double>>> Ex1_all, Ex2_all, Ey1_all, Ey2_all,
          Ez1_all, Ez2_all;
      std::vector<bool> HaveCurrentFields;
      std::map<double, int> FrequencyMap;
      boost::filesystem::path TempDir;
      std::string X3DName;
      std::string NameRoot;
      bool CleanFiles;
      double zshift;
      jif3D::GreenCalcType GreenStage1;
      jif3D::GreenCalcType GreenStage4;
      jif3D::X3DModel OldModel;
      std::vector<std::pair<size_t, size_t>> ForwardExecTime;
      //! A file to store statistics of execution time for the forward, can help to find problematic frequencies
      std::ofstream ForwardTimesFile;
      std::vector<std::string> ForwardDirName;
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
      const std::vector<std::complex<double> >& ReturnField(double Freq,
          const std::vector<std::vector<std::complex<double>>> &Field) const;
      void CalculateFields(const X3DModel &Model, const std::vector<double> &Frequencies,
          const std::vector<double> MeasPosZ, size_t freqindex);
    public:
      //! Provide serialization to be able to store objects and, more importantly for hpx parallelization
      template<class Archive>
      void save(Archive &ar, unsigned int version) const
        {
          ar & Ex1;
          ar & Ex2;
          ar & Ey1;
          ar & Ey2;
          ar & Hx1;
          ar & Hx2;
          ar & Hy1;
          ar & Hy2;
          ar & Hz1;
          ar & Hz2;
          ar & Ex1_all;
          ar & Ex2_all;
          ar & Ey1_all;
          ar & Ey2_all;
          ar & Ez1_all;
          ar & Ez2_all;
          std::string DirName(TempDir.string());
          ar & DirName;
          ar & CleanFiles;
          ar & zshift;
          ar & X3DName;
          ar & GreenStage1;
          ar & GreenStage4;
          ar & OldModel;
        }
      template<class Archive>
      void load(Archive &ar, unsigned int version)
        {
          ar & Ex1;
          ar & Ex2;
          ar & Ey1;
          ar & Ey2;
          ar & Hx1;
          ar & Hx2;
          ar & Hy1;
          ar & Hy2;
          ar & Hz1;
          ar & Hz2;
          ar & Ex1_all;
          ar & Ex2_all;
          ar & Ey1_all;
          ar & Ey2_all;
          ar & Ez1_all;
          ar & Ez2_all;
          std::string DirName;
          ar & DirName;
          TempDir = DirName;
          ar & CleanFiles;
          ar & zshift;
          ar & X3DName;
          NameRoot = ObjectID();
          ar & GreenStage1;
          ar & GreenStage4;
          ar & OldModel;
        }
#ifdef HAVEHPX
      HPX_SERIALIZATION_SPLIT_MEMBER()
#else
      BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
      const std::vector<std::complex<double> >& GetEx1(double Freq) const
        {
          return ReturnField(Freq, Ex1);
        }
      const std::vector<std::complex<double> >& GetEx2(double Freq) const
        {
          return ReturnField(Freq, Ex2);
        }
      const std::vector<std::complex<double> >& GetEy1(double Freq) const
        {
          return ReturnField(Freq, Ey1);
        }
      const std::vector<std::complex<double> >& GetEy2(double Freq) const
        {
          return ReturnField(Freq, Ey2);
        }
      const std::vector<std::complex<double> >& GetHx1(double Freq) const
        {
          return ReturnField(Freq, Hx1);
        }
      const std::vector<std::complex<double> >& GetHx2(double Freq) const
        {
          return ReturnField(Freq, Hx2);
        }
      const std::vector<std::complex<double> >& GetHy1(double Freq) const
        {
          return ReturnField(Freq, Hy1);
        }
      const std::vector<std::complex<double> >& GetHy2(double Freq) const
        {
          return ReturnField(Freq, Hy2);
        }
      const std::vector<std::complex<double> >& GetHz1(double Freq) const
        {
          return ReturnField(Freq, Hz1);
        }
      const std::vector<std::complex<double> >& GetHz2(double Freq) const
        {
          return ReturnField(Freq, Hz2);
        }
      const std::vector<std::complex<double> >& GetEx1_all(double Freq) const
        {
          return ReturnField(Freq, Ex1_all);
        }
      const std::vector<std::complex<double> >& GetEx2_all(double Freq) const
        {
          return ReturnField(Freq, Ex2_all);
        }
      const std::vector<std::complex<double> >& GetEy1_all(double Freq) const
        {
          return ReturnField(Freq, Ey1_all);
        }
      const std::vector<std::complex<double> >& GetEy2_all(double Freq) const
        {
          return ReturnField(Freq, Ey2_all);
        }
      const std::vector<std::complex<double> >& GetEz1_all(double Freq) const
        {
          return ReturnField(Freq, Ez1_all);
        }
      const std::vector<std::complex<double> >& GetEz2_all(double Freq) const
        {
          return ReturnField(Freq, Ez2_all);
        }
      std::string GetForwardDirName(double freq) const
        {
          size_t freqindex = FrequencyMap.at(freq);
          return ForwardDirName.at(freqindex);
        }
      void CalculateFields(const X3DModel &Model, const std::vector<double> &Frequencies,
          const std::vector<double> MeasPosZ);
      X3DFieldCalculator(boost::filesystem::path TDir = boost::filesystem::current_path(),
          std::string x3d = "x3d", double zs = 0.0, bool Clean = true,
          jif3D::GreenCalcType GS1 = hst, jif3D::GreenCalcType GS4 = hst);
      X3DFieldCalculator(const jif3D::X3DFieldCalculator &source);
      virtual ~X3DFieldCalculator();
      };

  }
/* namespace jif3D */

#endif /* MT_X3DFIELDCALCULATOR_H_ */
