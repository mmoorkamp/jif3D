//============================================================================
// Name        : MTUtils.h
// Author      : 10 Mar 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#ifndef MTUTILS_H_
#define MTUTILS_H_

#include <vector>
#include <string>
#include <complex>
#include <boost/filesystem/operations.hpp>
#include "../ModelBase/ThreeDModelBase.h"
#include "../Global/VecMat.h"

namespace jif3D
  {
    void MakeRunFile(const std::string &NameRoot, const std::string &DirName,
        const std::string &X3DName);
    void RunX3D(const std::string &NameRoot);

    void CopyHNK(const boost::filesystem::path &SourceDir,
        const boost::filesystem::path &TargetDir);
    bool CheckHNK(const boost::filesystem::path &TargetDir);

    jif3D::rvec AdaptDist(const std::vector<double> &C, const jif3D::rvec &RawImpedance,
        const jif3D::rvec &Misfit);
    cmat CalcEExt(const rvec &Misfit, const std::vector<double> &C,
        const size_t startindex, const size_t freq_start_index,
        const std::complex<double> &Hx1, const std::complex<double> &Hx2,
        const std::complex<double> &Hy1, const std::complex<double> &Hy2);
    void CalcHext(const std::complex<double> &omega_mu, std::complex<double> &Xp1,
        std::complex<double> &Xp2, std::complex<double> &Yp1, std::complex<double> &Yp2,
        const std::complex<double> &Zxx, const std::complex<double> &Zxy,
        const std::complex<double> &Zyx, const std::complex<double> &Zyy);

    void CheckField(const std::vector<std::complex<double> > &Field, size_t nelem);
    void CompareDepths(const std::vector<double> &BGDepths,
        const jif3D::ThreeDModelBase::t3DModelDim &ModelDepths);
  }
#endif /* MTUTILS_H_ */
