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

#include <string>
#include <vector>

#include "../Global/Jif3DGlobal.h"
#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/VectorTransform.h"
#include "ReadWriteX3D.h"
#include "X3DModel.h"
#include "X3DFieldCalculator.h"

struct J3DEXPORT ForwardResult
  {
  std::vector<double> DistImpedance;
  std::vector<double> RawImpedance;
  //! Provide serialization to be able to store objects
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
    {
      ar & DistImpedance;
      ar & RawImpedance;
    }
  };

struct GradInfo
  {
  std::vector<double> Misfit;
  std::vector<double> RawImpedance;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
    {
      ar & Misfit;
      ar & RawImpedance;
    }
  GradInfo() :
      Misfit(), RawImpedance()
    {

    }
  GradInfo(jif3D::rvec Mf, jif3D::rvec RI)
    {
      Misfit.resize(Mf.size());
      std::copy(Mf.begin(), Mf.end(), Misfit.begin());
      RawImpedance.resize(RI.size());
      std::copy(RI.begin(), RI.end(), RawImpedance.begin());
    }
  };

struct GradResult
  {
  std::vector<double> Gradient;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
    {
      ar & Gradient;
    }
  GradResult() :
      Gradient()
    {
    }
  GradResult(jif3D::rvec G)
    {
      Gradient.resize(G.size());
      std::copy(G.begin(), G.end(), Gradient.begin());
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

ForwardResult CalculateFrequency(const ForwardInfo &Info,
    boost::shared_ptr<jif3D::X3DFieldCalculator> Calc);

GradResult LQDerivativeFreq(const ForwardInfo &Info, const GradInfo &GI,
    boost::shared_ptr<jif3D::X3DFieldCalculator> Calc);

#ifdef HAVEHPX
HPX_DEFINE_PLAIN_ACTION(CalculateFrequency, CalculateFrequency_action);
HPX_DEFINE_PLAIN_ACTION(LQDerivativeFreq, LQDerivativeFreq_action);
#endif

#endif /* X3DFREQFUNCTIONS_H_ */
