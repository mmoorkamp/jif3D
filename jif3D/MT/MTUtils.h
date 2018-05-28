//============================================================================
// Name        : MTUtils.h
// Author      : 10 Mar 2014
// Version     :
// Copyright   : 2014, mm489
//============================================================================

#ifndef MTUTILS_H_
#define MTUTILS_H_

#include "../Global/Jif3DGlobal.h"
#include "../ModelBase/ThreeDModelBase.h"
#include "../Global/VecMat.h"
#include "X3DModel.h"

#include <vector>
#include <string>
#include <complex>

#include <boost/filesystem/operations.hpp>

namespace jif3D
  {

//define some names that are always the same
//either because x3d uses this convention
//or because we use them in their own directory
//and want to keep them simple to make sure x3d can handle them
    const std::string modelfilename("x3d.model");
    const std::string resultfilename("x3d.result");
    const std::string emaname = resultfilename + "0.ema";
    const std::string sourceafilename = modelfilename + "0a.source";
    const std::string sourcebfilename = modelfilename + "0b.source";
    const std::string emoAname = resultfilename + "0a.emo";
    const std::string emoBname = resultfilename + "0b.emo";
    const std::string emaAname = resultfilename + "0a.ema";
    const std::string emaBname = resultfilename + "0b.ema";
    const std::string dirext = "_dir";

    J3DEXPORT std::string MakeUniqueName(const std::string &NameRoot,
        X3DModel::ProblemType Type, const size_t FreqIndex);
    J3DEXPORT void MakeRunFile(const std::string &NameRoot, const std::string &DirName,
        const std::string &X3DName);
    J3DEXPORT void RunX3D(const std::string &NameRoot);

    J3DEXPORT void CopyHNK(const boost::filesystem::path &SourceDir,
        const boost::filesystem::path &TargetDir);
    J3DEXPORT bool CheckHNK(const boost::filesystem::path &TargetDir);

    J3DEXPORT jif3D::rvec AdaptDist(const std::vector<double> &C,
        const jif3D::rvec &RawImpedance, const jif3D::rvec &Misfit);
    //! Calculate the strength of the electric dipole necessary for MT adjoint calculations
    J3DEXPORT cmat CalcEExt(const rvec &Misfit, const std::vector<double> &C,
        const size_t startindex, const size_t freq_start_index,
        const std::complex<double> &Hx1, const std::complex<double> &Hx2,
        const std::complex<double> &Hy1, const std::complex<double> &Hy2);

    J3DEXPORT void CalcU(const std::string &RootName,
        const std::vector<std::complex<double> > &XPolMoments,
        const std::vector<std::complex<double> > &YPolMoments,
        const std::vector<std::complex<double> > &ZPolMoments,
        std::vector<std::complex<double> > &Ux, std::vector<std::complex<double> > &Uy,
        std::vector<std::complex<double> > &Uz, const std::vector<size_t> &XSourceXIndex,
        const std::vector<size_t> &XSourceYIndex,
        const std::vector<double> &XSourceDepths,
        const std::vector<size_t> &YSourceXIndex,
        const std::vector<size_t> &YSourceYIndex,
        const std::vector<double> &YSourceDepths,
        const std::vector<size_t> &ZSourceXIndex,
        const std::vector<size_t> &ZSourceYIndex,
        const std::vector<double> &ZSourceDepths,
        const jif3D::ThreeDModelBase::t3DModelDim &ZCellBoundaries,
        const jif3D::ThreeDModelBase::t3DModelDim &ZCellSizes, const size_t ncellsx,
        const size_t ncellsy, const size_t ncellsz);

    J3DEXPORT void CheckField(const std::vector<std::complex<double> > &Field,
        size_t nelem);
    J3DEXPORT void CompareDepths(const std::vector<double> &BGDepths,
        const jif3D::ThreeDModelBase::t3DModelDim &ModelDepths);
  }
#endif /* MTUTILS_H_ */
