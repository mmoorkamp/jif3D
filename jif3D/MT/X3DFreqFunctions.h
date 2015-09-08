/*
 * X3DFreqFunctions.h
 *
 *  Created on: Mar 17, 2014
 *      Author: mmoorkamp
 */

#ifndef X3DFREQFUNCTIONS_H_
#define X3DFREQFUNCTIONS_H_

#ifdef HAVEHPX
#include <hpx/config.hpp>
#include <hpx/include/actions.hpp>
#include <hpx/include/serialization.hpp>
#else
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#endif

#include <string>
#include <vector>
#include "../Global/VecMat.h"
#include "ReadWriteX3D.h"
#include "X3DModel.h"

struct ForwardResult
  {
  jif3D::rvec DistImpedance;
  jif3D::rvec RawImpedance;
  //! Provide serialization to be able to store objects
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
    {
      ar & DistImpedance;
      ar & RawImpedance;
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

jif3D::rvec LQDerivativeFreq(const ForwardInfo &Info, const jif3D::rvec &Misfit,
    const jif3D::rvec &RawImpedance);

#ifdef HAVEHPX
HPX_PLAIN_ACTION(CalculateFrequency, CalculateFrequency_action)
HPX_PLAIN_ACTION(LQDerivativeFreq, LQDerivativeFreq_action)
#endif

#endif /* X3DFREQFUNCTIONS_H_ */
