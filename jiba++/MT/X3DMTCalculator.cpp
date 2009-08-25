//============================================================================
// Name        : X3DMTCalculator.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <unistd.h>
#include <cassert>
#include <complex>
#include <fstream>
#include <map>
#include <omp.h>
#include <boost/assign/list_of.hpp> // for 'map_list_of()'
#include <boost/multi_array.hpp>
#include <boost/filesystem/operations.hpp>
#include "../Global/FatalException.h"
#include "../Global/convert.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"

namespace jiba
  {

    const std::string modelfilename("x3d.model");
    const std::string resultfilename("x3d.result");
    const std::string emaname = resultfilename + "0.ema";
    const std::string sourcefilename = modelfilename + "0a.source";
    const std::string emoAname = resultfilename + "0a.emo";
    const std::string emoBname = resultfilename + "0b.emo";
    const std::string emaAname = resultfilename + "0a.ema";
    const std::string emaBname = resultfilename + "0b.ema";
    const std::map<X3DModel::ProblemType, std::string> Extension =
        boost::assign::map_list_of(X3DModel::MT, "MT")(X3DModel::EDIP, "EDIP")(
            X3DModel::MDIP, "MDIP");
    X3DMTCalculator::X3DMTCalculator()
      {

      }

    X3DMTCalculator::~X3DMTCalculator()
      {

      }

    bool CheckHNK(const boost::filesystem::path &TargetDir)
      {
        return boost::filesystem::exists(TargetDir / "ndec15.hnk")
            && boost::filesystem::exists(TargetDir / "ndec20.hnk")
            && boost::filesystem::exists(TargetDir / "ndec30.hnk")
            && boost::filesystem::exists(TargetDir / "ndec40.hnk");
      }

    void CopyHNK(const boost::filesystem::path &SourceDir,
        const boost::filesystem::path &TargetDir)
      {
        if (!CheckHNK(TargetDir))
          {
            boost::filesystem::copy_file(SourceDir / "ndec15.hnk", TargetDir
                / "ndec15.hnk");
            boost::filesystem::copy_file(SourceDir / "ndec20.hnk", TargetDir
                / "ndec20.hnk");
            boost::filesystem::copy_file(SourceDir / "ndec30.hnk", TargetDir
                / "ndec30.hnk");
            boost::filesystem::copy_file(SourceDir / "ndec40.hnk", TargetDir
                / "ndec40.hnk");
          }
      }

    std::string X3DMTCalculator::MakeUniqueName(X3DModel::ProblemType Type,
        const size_t FreqIndex)
      {
        std::string result("p" + jiba::stringify(getpid()) + jiba::stringify(
            this));

        result += Extension.at(Type);

        result += jiba::stringify(FreqIndex);
        return result;
      }

    void MakeRunFile(const std::string &NameRoot)
      {
        std::string DirName = NameRoot + "_dir";
        std::string RunFileName = NameRoot + "_run";
        boost::filesystem::create_directory(DirName);
        std::ofstream runfile;
        runfile.open(RunFileName.c_str());
        runfile << "#!/bin/bash\n";
        runfile << "cd " << DirName << "\n";
        runfile << "x3d > /dev/null\n";
        runfile << "cd ..\n";
        runfile.close();
        system((std::string("chmod u+x ./") + RunFileName).c_str());
        CopyHNK(boost::filesystem::current_path(), DirName);
      }

    void RunX3D(const std::string &NameRoot)
      {
        system(("./" + NameRoot + "_run").c_str());
      }

    void CleanFiles(const std::string &NameRoot)
      {
        boost::filesystem::remove_all(NameRoot + "_dir");
        boost::filesystem::remove_all(NameRoot + "_run");
      }

    rvec X3DMTCalculator::Calculate(const X3DModel &Model)
      {

        const size_t nfreq = Model.GetFrequencies().size();
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nmodx = Model.GetXCoordinates().size();
        const size_t nmody = Model.GetYCoordinates().size();
        const size_t nmodz = Model.GetZCoordinates().size();

        jiba::rvec result(nmeas * nfreq * 8);
        result.clear();
        omp_lock_t lck;
        omp_init_lock(&lck);
#pragma omp parallel default(shared)
          {
#pragma omp for
            for (int i = 0; i < nfreq; ++i)
              {
                std::string RootName = MakeUniqueName(X3DModel::MT, i);
                std::string DirName = RootName + "_dir";
                std::vector<double> CurrFreq(1, Model.GetFrequencies()[i]);
                //writing out files causes problems in parallel
                // so we make sure it is done one at a time
#pragma omp critical
                  {
                    MakeRunFile(RootName);
                    WriteProjectFile(DirName, CurrFreq, X3DModel::MT,
                        resultfilename, modelfilename);
                    Write3DModelForX3D(DirName + "/" + modelfilename,
                        Model.GetXCellSizes(), Model.GetYCellSizes(),
                        Model.GetZCellSizes(), Model.GetConductivities(),
                        Model.GetBackgroundConductivities(),
                        Model.GetBackgroundThicknesses());
                  }
                RunX3D(RootName);
                std::vector<std::complex<double> > Ex1, Ex2, Ey1, Ey2, Hx1,
                    Hx2, Hy1, Hy2;
                std::complex<double> Zxx, Zxy, Zyx, Zyy;
                ReadEMO(DirName + "/" + emoAname, Ex1, Ey1, Hx1, Hy1);
                ReadEMO(DirName + "/" + emoBname, Ex2, Ey2, Hx2, Hy2);
                const size_t freq_index = nmeas * i * 8;
                for (size_t j = 0; j < nmeas; ++j)
                  {
                    boost::array<ThreeDModelBase::t3DModelData::index, 3>
                        StationIndex = Model.FindAssociatedIndices(
                            Model.GetMeasPosX()[j], Model.GetMeasPosY()[j],
                            Model.GetMeasPosZ()[j]);
                    //at the moment we ignore the depth/elevation of the site
                    const size_t offset = StationIndex[1] * nmodx
                        + StationIndex[0];
                    const size_t meas_index = freq_index + j * 8;
                    FieldsToImpedance(Ex1[offset], Ex2[offset], Ey1[offset],
                        Ey2[offset], Hx1[offset], Hx2[offset], Hy1[offset],
                        Hy2[offset], Zxx, Zxy, Zyx, Zyy);
                    omp_set_lock(&lck);
                    result(meas_index) = Zxx.real();
                    result(meas_index + 1) = Zxx.imag();
                    result(meas_index + 2) = Zxy.real();
                    result(meas_index + 3) = Zxy.imag();
                    result(meas_index + 4) = Zyx.real();
                    result(meas_index + 5) = Zyx.imag();
                    result(meas_index + 6) = Zyy.real();
                    result(meas_index + 7) = Zyy.imag();
                    omp_unset_lock(&lck);
                  }
                //finished with one frequency
              }
            //end of parallel section
          }
        omp_destroy_lock(&lck);
        return result;
      }

    cmat CalcATimesH(const cmat &A, const cmat &H)
      {
        cmat result(2, 2);
        const std::complex<double> magdet = 1. / (H(0, 0) * H(1, 1) - H(0, 1)
            * H(1, 0));
        result(0, 0) = magdet * (conj(A(0, 0)) * H(1, 1) - conj(A(0, 1)) * H(0,
            1));
        result(0, 1) = magdet * (-conj(A(0, 0)) * H(1, 0) + conj(A(0, 1)) * H(
            0, 0));
        result(1, 0) = magdet * (conj(A(1, 0)) * H(1, 1) - conj(A(1, 1)) * H(0,
            1));
        result(1, 1) = magdet * (-conj(A(1, 0)) * H(1, 0) + conj(A(1, 1)) * H(
            0, 0));
        return result;
      }

    cmat MisfitToA(const rvec &Misfit, const size_t startindex)
      {
        assert(startindex <= Misfit.size()-8);
        cmat result(2, 2);
        result(0, 0) = std::complex<double>(Misfit(startindex), Misfit(
            startindex + 1));
        result(0, 1) = std::complex<double>(Misfit(startindex + 2), Misfit(
            startindex + 3));
        result(1, 0) = std::complex<double>(Misfit(startindex + 4), Misfit(
            startindex + 5));
        result(1, 1) = std::complex<double>(Misfit(startindex + 6), Misfit(
            startindex + 7));
        return result;
      }

    cmat MakeH(const std::complex<double> &Hx1,
        const std::complex<double> &Hx2, const std::complex<double> &Hy1,
        const std::complex<double> &Hy2)
      {
        cmat result(2, 2);
        result(0, 0) = Hx1;
        result(0, 1) = Hx2;
        result(1, 0) = Hy1;
        result(1, 1) = Hy2;
        return result;
      }

    void CalcHext(const std::complex<double> &omega_mu,
        std::complex<double> &Xp1, std::complex<double> &Xp2, std::complex<
            double> &Yp1, std::complex<double> &Yp2,
        const std::complex<double> &Zxx, const std::complex<double> &Zxy,
        const std::complex<double> &Zyx, const std::complex<double> &Zyy)
      {
        std::complex<double> temp1 = conj(Zxx) * Xp1 + conj(Zyx) * Yp1;
        std::complex<double> temp2 = conj(Zxy) * Xp1 + conj(Zyy) * Yp1;
        Xp1 = temp1 * omega_mu;
        Yp1 = temp2 * omega_mu;
        temp1 = conj(Zxx) * Xp2 + conj(Zyx) * Yp2;
        temp2 = conj(Zxy) * Xp2 + conj(Zyy) * Yp2;
        Xp2 = temp1 * omega_mu;
        Yp2 = temp2 * omega_mu;
      }

    void CalcU(const std::string &RootName, const boost::multi_array<
        std::complex<double>, 2> &XPolMoments, const boost::multi_array<
        std::complex<double>, 2> &YPolMoments,
        std::vector<std::complex<double> > &Ux, std::vector<
            std::complex<double> > &Uy, std::vector<std::complex<double> > &Uz,
        const size_t ncellsx, const size_t ncellsy, const size_t ncellsz)
      {
        std::string DirName = RootName + "_dir/";
#pragma omp critical
          {
            WriteSourceFile(DirName + sourcefilename, 0.0, XPolMoments,
                YPolMoments);
          }
        RunX3D(RootName);
        ReadEMA(DirName + emaname, Ux, Uy, Uz, ncellsx, ncellsy, ncellsz);
        //boost::filesystem::remove_all(emaname);
      }

    rvec X3DMTCalculator::LQDerivative(const X3DModel &Model,
        const rvec &Misfit)
      {
        const size_t nfreq = Model.GetFrequencies().size();
        const size_t ncellsx = Model.GetConductivities().shape()[0];
        const size_t ncellsy = Model.GetConductivities().shape()[1];
        const size_t ncellsz = Model.GetConductivities().shape()[2];
        const size_t nobs = ncellsx * ncellsy;
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t ndata = nmeas * nfreq * 8;
        const size_t nmod = ncellsx * ncellsy * ncellsz;
        assert(Misfit.size() == ndata);
        //std::cout << "Misfit: " << Misfit << std::endl;
        jiba::rvec Gradient(nmod);
        Gradient.clear();;
#pragma omp parallel default(shared)
          {
#pragma omp for
            for (int i = 0; i < nfreq; ++i)
              {
                std::string ForwardName = MakeUniqueName(X3DModel::MT, i);
                std::string ForwardDirName = ForwardName + "_dir/";

                //read the fields from the forward calculation
                std::vector<std::complex<double> > Ex1_obs, Ex2_obs, Ey1_obs,
                    Ey2_obs, Hx1_obs, Hx2_obs, Hy1_obs, Hy2_obs;
                ReadEMO(ForwardDirName + emoAname, Ex1_obs, Ey1_obs, Hx1_obs,
                    Hy1_obs);
                ReadEMO(ForwardDirName + emoBname, Ex2_obs, Ey2_obs, Hx2_obs,
                    Hy2_obs);
                assert(Ex1_obs.size()==nobs);
                assert(Ex2_obs.size()==nobs);
                std::vector<std::complex<double> > Ex1_all, Ex2_all, Ey1_all,
                    Ey2_all, Ez1_all, Ez2_all;
                ReadEMA(ForwardDirName + emaAname, Ex1_all, Ey1_all, Ez1_all,
                    ncellsx, ncellsy, ncellsz);
                ReadEMA(ForwardDirName + emaBname, Ex2_all, Ey2_all, Ez2_all,
                    ncellsx, ncellsy, ncellsz);

                //create variables for the adjoint field calculation
                const size_t freq_index = nmeas * i * 8;
                boost::multi_array<std::complex<double>, 2> XPolMoments1(
                    boost::extents[ncellsx][ncellsy]), YPolMoments1(
                    boost::extents[ncellsx][ncellsy]);
                boost::multi_array<std::complex<double>, 2> XPolMoments2(
                    boost::extents[ncellsx][ncellsy]), YPolMoments2(
                    boost::extents[ncellsx][ncellsy]);
                boost::multi_array<std::complex<double>, 2> Zeros(
                    boost::extents[ncellsx][ncellsy]);
                std::fill(Zeros.origin(),
                    Zeros.origin() + Zeros.num_elements(), 0.0);
                std::fill(XPolMoments1.origin(), XPolMoments1.origin() + nobs,
                    0.0);
                std::fill(XPolMoments2.origin(), XPolMoments2.origin() + nobs,
                    0.0);
                std::fill(YPolMoments1.origin(), YPolMoments1.origin() + nobs,
                    0.0);
                std::fill(YPolMoments2.origin(), YPolMoments2.origin() + nobs,
                    0.0);
                //make the sources for the electric dipoles
                for (size_t j = 0; j < nmeas; ++j)
                  {
                    boost::array<ThreeDModelBase::t3DModelData::index, 3>
                        StationIndex = Model.FindAssociatedIndices(
                            Model.GetMeasPosX()[j], Model.GetMeasPosY()[j],
                            Model.GetMeasPosZ()[j]);
                    //at the moment we ignore the depth/elevation of the site
                    const size_t offset = StationIndex[1] * ncellsx
                        + StationIndex[0];

                    const size_t siteindex = freq_index + j * 8;
                    cmat j_ext = CalcATimesH(MisfitToA(Misfit, siteindex),
                        MakeH(Hx1_obs[offset], Hx2_obs[offset],
                            Hy1_obs[offset], Hy2_obs[offset]));
                    //std::cout << result << std::endl;
                    XPolMoments1.data()[offset] = conj(j_ext(0, 0));
                    YPolMoments1.data()[offset] = conj(j_ext(1, 0));
                    XPolMoments2.data()[offset] = conj(j_ext(0, 1));
                    YPolMoments2.data()[offset] = conj(j_ext(1, 1));
                  }
                std::vector<double> CurrFreq(1, 0.0);
                CurrFreq[0] = Model.GetFrequencies()[i];
                std::string EdipName = MakeUniqueName(X3DModel::EDIP, i);
                std::string EdipDirName = EdipName + "_dir/";
#pragma omp critical
                  {
                    MakeRunFile(EdipName);
                    WriteProjectFile(EdipDirName, CurrFreq, X3DModel::EDIP,
                        resultfilename, modelfilename);
                    Write3DModelForX3D(EdipDirName + modelfilename,
                        Model.GetXCellSizes(), Model.GetYCellSizes(),
                        Model.GetZCellSizes(), Model.GetConductivities(),
                        Model.GetBackgroundConductivities(),
                        Model.GetBackgroundThicknesses());
                    //write an empty source file for the second source polarization
                    WriteSourceFile(EdipDirName + modelfilename + "0b.source",
                        0.0, Zeros, Zeros);
                  }
                std::vector<std::complex<double> > Ux1_el, Ux2_el, Uy1_el,
                    Uy2_el, Uz1_el, Uz2_el;
                CalcU(EdipName, XPolMoments1, YPolMoments1, Ux1_el, Uy1_el,
                    Uz1_el, ncellsx, ncellsy, ncellsz);
                // boost::filesystem::rename(sourcefilename, "model.edipa.source");
                //calculate the second polarization
                CalcU(EdipName, XPolMoments2, YPolMoments2, Ux2_el, Uy2_el,
                    Uz2_el, ncellsx, ncellsy, ncellsz);
                //boost::filesystem::rename(sourcefilename, "model.edipb.source");


                //first polarization of the magnetic dipole
                const std::complex<double> omega_mu = -1.0 / (std::complex<
                    double>(0.0, jiba::mag_mu) * 2.0 * M_PI
                    * Model.GetFrequencies()[i]);

                std::string MdipName = MakeUniqueName(X3DModel::MDIP, i);
                std::string MdipDirName = MdipName + "_dir/";
#pragma omp critical
                  {
                    MakeRunFile(MdipName);
                    WriteProjectFile(MdipDirName, CurrFreq, X3DModel::MDIP,
                        resultfilename, modelfilename);
                    Write3DModelForX3D(MdipDirName + modelfilename,
                        Model.GetXCellSizes(), Model.GetYCellSizes(),
                        Model.GetZCellSizes(), Model.GetConductivities(),
                        Model.GetBackgroundConductivities(),
                        Model.GetBackgroundThicknesses());
                  }
                for (size_t j = 0; j < nobs; ++j)
                  {
                    cmat Z(2, 2);
                    FieldsToImpedance(Ex1_obs[j], Ex2_obs[j], Ey1_obs[j],
                        Ey2_obs[j], Hx1_obs[j], Hx2_obs[j], Hy1_obs[j],
                        Hy2_obs[j], Z(0, 0), Z(0, 1), Z(1, 0), Z(1, 1));
                    CalcHext(omega_mu, XPolMoments1.data()[j],
                        XPolMoments2.data()[j], YPolMoments1.data()[j],
                        YPolMoments2.data()[j], Z(0, 0), Z(0, 1), Z(1, 0), Z(1,
                            1));
                  }

                std::vector<std::complex<double> > Ux1_mag, Ux2_mag, Uy1_mag,
                    Uy2_mag, Uz1_mag, Uz2_mag;
                CalcU(MdipName, XPolMoments1, YPolMoments1, Ux1_mag, Uy1_mag,
                    Uz1_mag, ncellsx, ncellsy, ncellsz);
                //boost::filesystem::rename(sourcefilename, "model.mdipa.source");
                CalcU(MdipName, XPolMoments2, YPolMoments2, Ux2_mag, Uy2_mag,
                    Uz2_mag, ncellsx, ncellsy, ncellsz);
                // boost::filesystem::rename(sourcefilename, "model.mdipb.source");
                const double cell_sizex = Model.GetXCellSizes()[0];
                const double cell_sizey = Model.GetYCellSizes()[0];
                for (size_t j = 0; j < nmod; ++j)
                  {
                    const double Volume = cell_sizex * cell_sizey
                        * Model.GetZCellSizes()[j % ncellsz];
#pragma omp atomic
                    Gradient(j) += std::real((Ux1_el[j] + Ux1_mag[j])
                        * Ex1_all[j] + (Uy1_el[j] + Uy1_mag[j]) * Ey1_all[j]
                        + (Uz1_el[j] + Uz1_mag[j]) * Ez1_all[j] + (Ux2_el[j]
                        + Ux2_mag[j]) * Ex2_all[j] + (Uy2_el[j] + Uy2_mag[j])
                        * Ey2_all[j] + (Uz2_el[j] + Uz2_mag[j]) * Ez2_all[j])
                        * Volume;
                  }
                CleanFiles(ForwardDirName);
                CleanFiles(EdipName);
                CleanFiles(MdipName);
                //finished with one frequency
              }
            //end of parallel region
          }
        X3DModel GradMod(Model);
        std::copy(Gradient.begin(), Gradient.end(),
            GradMod.SetConductivities().origin());
        GradMod.WriteVTK("grad.vtk");
        return Gradient;
      }
  }
