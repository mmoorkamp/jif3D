/*
 * X3DFreqFunctions.h
 *
 *  Created on: Mar 17, 2014
 *      Author: mmoorkamp
 */

#ifndef X3DFREQFUNCTIONS_H_
#define X3DFREQFUNCTIONS_H_

#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#endif

#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/VectorTransform.h"
#include "ReadWriteX3D.h"
#include "X3DModel.h"
#include <string>
#include <vector>

struct ForwardResult
  {
  jif3D::rvec DistImpedance;
  jif3D::rvec RawImpedance;
  //! Provide serialization to be able to store objects
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
    {
      ar & std::vector<double>(DistImpedance.begin(), DistImpedance.end());
      ar & std::vector<double>(RawImpedance.begin(), RawImpedance.end());
    }
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
    {
      std::vector<double> v;
      ar & v;
      DistImpedance.resize(v.size());
      std::copy(v.begin(), v.end(), DistImpedance.begin());
      ar & v;
      RawImpedance.resize(v.size());
      std::copy(v.begin(), v.end(), RawImpedance.begin());
    }
#ifdef HAVEHPX
  HPX_SERIALIZATION_SPLIT_MEMBER()
#else
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
}  ;

struct GradInfo
  {
  jif3D::rvec Misfit;
  jif3D::rvec RawImpedance;
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
    {
      ar & std::vector<double>(Misfit.begin(), Misfit.end());
      ar & std::vector<double>(RawImpedance.begin(), RawImpedance.end());
    }
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
    {
      std::vector<double> v;
      ar & v;
      Misfit.resize(v.size());
      std::copy(v.begin(), v.end(), Misfit.begin());
      ar & v;
      RawImpedance.resize(v.size());
      std::copy(v.begin(), v.end(), RawImpedance.begin());
    }
#ifdef HAVEHPX
  HPX_SERIALIZATION_SPLIT_MEMBER()
#else
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
GradInfo  () :
  Misfit(), RawImpedance()
    {

    }
  GradInfo(jif3D::rvec Mf, jif3D::rvec RI) :
  Misfit(Mf), RawImpedance(RI)
    {

    }
};

struct GradResult
  {
  jif3D::rvec Gradient;
  template<class Archive>
  void save(Archive & ar, const unsigned int version) const
    {
      ar & std::vector<double>(Gradient.begin(), Gradient.end());
    }
  template<class Archive>
  void load(Archive & ar, const unsigned int version)
    {
      std::vector<double> v;
      ar & v;
      Gradient.resize(v.size());
      std::copy(v.begin(), v.end(), Gradient.begin());

    }
#ifdef HAVEHPX
  HPX_SERIALIZATION_SPLIT_MEMBER()
#else
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
GradResult  () :
  Gradient()
    {
    }
  GradResult(jif3D::rvec G) :
  Gradient(G)
    {
    }
};

struct ForwardInfo
  {
  jif3D::X3DModel Model;
  std::vector<double> C;
  size_t freqindex;
  std::string TempDirName;
  std::string X3DName;
  std::string NameRoot;
  jif3D::GreenCalcType GreenStage1;
  jif3D::GreenCalcType GreenStage4;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
    {
      ar & Model;
      ar & C;
      ar & freqindex;
      ar & TempDirName;
      ar & X3DName;
      ar & NameRoot;
      ar & GreenStage1;
      ar & GreenStage4;
    }
  ForwardInfo(jif3D::X3DModel M, std::vector<double> Dist, size_t f, std::string TD,
      std::string XN, std::string NR, jif3D::GreenCalcType G1 = jif3D::GreenCalcType::hst,
      jif3D::GreenCalcType G4 = jif3D::GreenCalcType::hst) :
      Model(M), C(Dist), freqindex(f), TempDirName(TD), X3DName(XN), NameRoot(NR), GreenStage1(
          G1), GreenStage4(G4)
    {

    }
  ForwardInfo() :
      Model(), C(), freqindex(), TempDirName(), X3DName(), NameRoot(), GreenStage1(
          jif3D::GreenCalcType::hst), GreenStage4(jif3D::GreenCalcType::hst)
    {

    }
  };

ForwardResult CalculateFrequency(const ForwardInfo &Info);

GradResult LQDerivativeFreq(const ForwardInfo &Info, const GradInfo &GI, const jif3D::VectorTransform &DataTransform);

#ifdef HAVEHPX
HPX_DEFINE_PLAIN_ACTION(CalculateFrequency, CalculateFrequency_action);
HPX_DEFINE_PLAIN_ACTION(LQDerivativeFreq, LQDerivativeFreq_action);
#endif

#endif /* X3DFREQFUNCTIONS_H_ */
