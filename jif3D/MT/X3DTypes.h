/*
 * X3DTypes.h
 *
 *  Created on: 26 Feb 2018
 *      Author: mm489
 */

#ifndef MT_X3DTYPES_H_
#define MT_X3DTYPES_H_


#include "../Global/Jif3DGlobal.h"
#include "../Global/Serialization.h"
#include "../Global/VecMat.h"
#include "../Global/VectorTransform.h"
#include "X3DModel.h"

#include <string>
#include <vector>


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
      std::copy(G.begin(),G.end(),Gradient.begin());
    }
  };

struct ForwardInfo
  {
  jif3D::X3DModel Model;
  std::vector<double> C;
  std::vector<double> Angles;
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
      ar & Angles;
      ar & freqindex;
      ar & TempDirName;
      ar & X3DName;
      ar & NameRoot;
      ar & GreenStage1;
      ar & GreenStage4;
    }
  ForwardInfo(jif3D::X3DModel M, std::vector<double> Dist, td::vector<double> Rot, size_t f, std::string TD,
      std::string XN, std::string NR, jif3D::GreenCalcType G1 = jif3D::GreenCalcType::hst,
      jif3D::GreenCalcType G4 = jif3D::GreenCalcType::hst) :
      Model(M), C(Dist), Angles(Rot), freqindex(f), TempDirName(TD), X3DName(XN), NameRoot(NR), GreenStage1(
          G1), GreenStage4(G4)
    {

    }
  ForwardInfo() :
      Model(), C(), Angles(), freqindex(), TempDirName(), X3DName(), NameRoot(), GreenStage1(
          jif3D::GreenCalcType::hst), GreenStage4(jif3D::GreenCalcType::hst)
    {

    }
  };



#endif /* MT_X3DTYPES_H_ */
