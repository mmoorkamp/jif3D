//============================================================================
// Name        : X3DMTCalculator.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <cassert>
#include <complex>
#include <boost/multi_array.hpp>
#include <boost/filesystem/operations.hpp>
#include "../Global/convert.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"

namespace jiba
  {

    const std::string modelfilename("x3d.model");
    const std::string resultfilename("x3d.result");
    X3DMTCalculator::X3DMTCalculator()
      {

      }

    X3DMTCalculator::~X3DMTCalculator()
      {

      }

    rvec X3DMTCalculator::Calculate(const X3DModel &Model)
      {

        WriteProjectFile(Model.GetFrequencies(), X3DModel::MT, resultfilename,
            modelfilename);
        Write3DModelForX3D(modelfilename, Model.GetXCellSizes(),
            Model.GetYCellSizes(), Model.GetZCellSizes(),
            Model.GetConductivities(), Model.GetBackgroundConductivities(),
            Model.GetBackgroundThicknesses());
        system("x3d> /dev/null");
        jiba::rvec result(Model.GetXCellSizes().size()
            * Model.GetYCellSizes().size() * Model.GetFrequencies().size() * 8);
        result.clear();
        const size_t nfreq = Model.GetFrequencies().size();
        for (size_t i = 0; i < nfreq; ++i)
          {
            std::string emoAname = resultfilename + jiba::stringify(i)
                + "a.emo";
            std::string emoBname = resultfilename + jiba::stringify(i)
                + "b.emo";
            std::vector<std::complex<double> > Ex1, Ex2, Ey1, Ey2, Hx1, Hx2,
                Hy1, Hy2;
            std::complex<double> Zxx, Zxy, Zyx, Zyy;
            ReadEMO(emoAname, Ex1, Ey1, Hx1, Hy1);
            ReadEMO(emoBname, Ex2, Ey2, Hx2, Hy2);
            const size_t nobs = Ex1.size();
            const size_t freq_index = nobs * i * 8;
            for (size_t j = 0; j < nobs; ++j)
              {
                const size_t obs_index = freq_index + j * 8;
                FieldsToImpedance(Ex1[j], Ex2[j], Ey1[j], Ey2[j], Hx1[j],
                    Hx2[j], Hy1[j], Hy2[j], Zxx, Zxy, Zyx, Zyy);
                result(obs_index) = Zxx.real();
                result(obs_index + 1) = Zxx.imag();
                result(obs_index + 2) = Zxy.real();
                result(obs_index + 3) = Zxy.imag();
                result(obs_index + 4) = Zyx.real();
                result(obs_index + 5) = Zyx.imag();
                result(obs_index + 6) = Zyy.real();
                result(obs_index + 7) = Zyy.imag();
              }
            //finished with one frequency
          }
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

    void CalcU(const std::string &sourcefilename, const std::string &emaname,
        const boost::multi_array<std::complex<double>, 2> &XPolMoments,
        const boost::multi_array<std::complex<double>, 2> &YPolMoments,
        std::vector<std::complex<double> > &Ux, std::vector<
            std::complex<double> > &Uy, std::vector<std::complex<double> > &Uz,
        const size_t ncellsx, const size_t ncellsy, const size_t ncellsz)
      {
        WriteSourceFile(sourcefilename, 0.0, XPolMoments, YPolMoments);
        system("x3d > /dev/null");
        ReadEMA(emaname, Ux, Uy, Uz);
        Ux = ResortFields(Ux, ncellsx, ncellsy, ncellsz);
        Uy = ResortFields(Uy, ncellsx, ncellsy, ncellsz);
        Uz = ResortFields(Uz, ncellsx, ncellsy, ncellsz);
      }

    rvec X3DMTCalculator::LQDerivative(const X3DModel &Model,
        const rvec &Misfit)
      {
        const size_t nfreq = Model.GetFrequencies().size();
        const size_t ncellsx = Model.GetConductivities().shape()[0];
        const size_t ncellsy = Model.GetConductivities().shape()[1];
        const size_t ncellsz = Model.GetConductivities().shape()[2];
        const size_t nobs = ncellsx * ncellsy;
        const size_t ndata = nobs * nfreq * 8;
        const size_t nmod = ncellsx * ncellsy * ncellsz;
        assert(Misfit.size() == ndata);
        //std::cout << "Misfit: " << Misfit << std::endl;
        jiba::rvec Gradient(Model.GetConductivities().num_elements());
        Gradient.clear();
        for (size_t i = 0; i < nfreq; ++i)
          {
            std::string emoAname = resultfilename + jiba::stringify(i)
                + "a.emo";
            std::string emoBname = resultfilename + jiba::stringify(i)
                + "b.emo";
            std::string emaname = resultfilename + jiba::stringify(i) + ".ema";

            std::vector<std::complex<double> > Ex1_obs, Ex2_obs, Ey1_obs,
                Ey2_obs, Hx1_obs, Hx2_obs, Hy1_obs, Hy2_obs;

            ReadEMO(emoAname, Ex1_obs, Ey1_obs, Hx1_obs, Hy1_obs);
            ReadEMO(emoBname, Ex2_obs, Ey2_obs, Hx2_obs, Hy2_obs);
            assert(Ex1_obs.size()==nobs);
            assert(Ex2_obs.size()==nobs);
            std::vector<std::complex<double> > Ex1_all, Ex2_all, Ey1_all,
                Ey2_all, Ez1_all, Ez2_all;
            ReadEMA(resultfilename + jiba::stringify(i) + "a.ema", Ex1_all,
                Ey1_all, Ez1_all);
            ReadEMA(resultfilename + jiba::stringify(i) + "b.ema", Ex2_all,
                Ey2_all, Ez2_all);
            assert(Ex1_all.size() == nmod);
            assert(Ex2_all.size() == nmod);
            Ex1_all = ResortFields(Ex1_all, ncellsx, ncellsy, ncellsz);
            Ex2_all = ResortFields(Ex2_all, ncellsx, ncellsy, ncellsz);
            Ey1_all = ResortFields(Ey1_all, ncellsx, ncellsy, ncellsz);
            Ey2_all = ResortFields(Ey2_all, ncellsx, ncellsy, ncellsz);
            Ez1_all = ResortFields(Ez1_all, ncellsx, ncellsy, ncellsz);
            Ez2_all = ResortFields(Ez2_all, ncellsx, ncellsy, ncellsz);
            const size_t freq_index = nobs * i * 8;
            boost::multi_array<std::complex<double>, 2> XPolMoments1(
                boost::extents[ncellsx][ncellsy]), YPolMoments1(
                boost::extents[ncellsx][ncellsy]);
            boost::multi_array<std::complex<double>, 2> XPolMoments2(
                boost::extents[ncellsx][ncellsy]), YPolMoments2(
                boost::extents[ncellsx][ncellsy]);
            boost::multi_array_ref<std::complex<double>, 1> XPolVec1(
                XPolMoments1.data(),
                boost::extents[XPolMoments1.num_elements()]);
            boost::multi_array_ref<std::complex<double>, 1> YPolVec1(
                YPolMoments1.data(),
                boost::extents[YPolMoments1.num_elements()]);
            boost::multi_array_ref<std::complex<double>, 1> XPolVec2(
                XPolMoments2.data(),
                boost::extents[XPolMoments2.num_elements()]);
            boost::multi_array_ref<std::complex<double>, 1> YPolVec2(
                YPolMoments2.data(),
                boost::extents[YPolMoments2.num_elements()]);
            for (size_t j = 0; j < nobs; ++j)
              {
                const size_t siteindex = freq_index + j * 8;
                cmat j_ext = CalcATimesH(MisfitToA(Misfit, siteindex), MakeH(
                    Hx1_obs[j], Hx2_obs[j], Hy1_obs[j], Hy2_obs[j]));
                //std::cout << result << std::endl;
                XPolVec1[j] = j_ext(0, 0);
                YPolVec1[j] = j_ext(1, 0);
                XPolVec2[j] = j_ext(0, 1);
                YPolVec2[j] = j_ext(1, 1);
              }
            WriteProjectFile(Model.GetFrequencies(), X3DModel::EDIP,
                resultfilename, modelfilename);
            std::string sourcefilename = modelfilename + stringify(i)
                + "a.source";
            //write an empty source file for the second source polarization
            boost::multi_array<std::complex<double>, 2> Zeros(
                boost::extents[ncellsx][ncellsy]);
            std::fill(Zeros.origin(), Zeros.origin() + Zeros.num_elements(),
                0.0);

            WriteSourceFile(modelfilename + stringify(i) + "b.source", 0.0,
                Zeros, Zeros);
            std::vector<std::complex<double> > Ux1_el, Ux2_el, Uy1_el, Uy2_el,
                Uz1_el, Uz2_el;
            CalcU(sourcefilename, emaname, XPolMoments1, YPolMoments1, Ux1_el,
                Uy1_el, Uz1_el, ncellsx, ncellsy, ncellsz);

            //calculate the second polarization
            CalcU(sourcefilename, emaname, XPolMoments2, YPolMoments2, Ux2_el,
                Uy2_el, Uz2_el, ncellsx, ncellsy, ncellsz);

            //first polarization of the magnetic dipole
            const std::complex<double> omega_mu = -1.0 / (std::complex<double>(
                0.0, jiba::mag_mu) * 2.0 * M_PI * Model.GetFrequencies()[i]);
            WriteProjectFile(Model.GetFrequencies(), X3DModel::MDIP,
                resultfilename, modelfilename);
            for (size_t j = 0; j < nobs; ++j)
              {
                cmat Z(2, 2);
                const size_t siteindex = freq_index + j * 8;
                FieldsToImpedance(Ex1_obs[i], Ex2_obs[i], Ey1_obs[i],
                    Ey2_obs[i], Hx1_obs[i], Hx2_obs[i], Hy1_obs[i], Hy2_obs[i],
                    Z(0, 0), Z(0, 1), Z(1, 0), Z(1, 1));
                cmat AH = CalcATimesH(MisfitToA(Misfit, siteindex), MakeH(
                    Hx1_obs[j], Hx2_obs[j], Hy1_obs[j], Hy2_obs[j]));
                cmat h_ext = ublas::prod(trans(Z), AH);
                XPolVec1[j] = h_ext(0, 0) * omega_mu;
                YPolVec1[j] = h_ext(1, 0) * omega_mu;
                XPolVec2[j] = h_ext(0, 1) * omega_mu;
                YPolVec2[j] = h_ext(1, 1) * omega_mu;
              }

            std::vector<std::complex<double> > Ux1_mag, Ux2_mag, Uy1_mag,
                Uy2_mag, Uz1_mag, Uz2_mag;
            CalcU(sourcefilename, emaname, XPolMoments1, YPolMoments1, Ux1_mag,
                Uy1_mag, Uz1_mag, ncellsx, ncellsy, ncellsz);

            CalcU(sourcefilename, emaname, XPolMoments2, YPolMoments2, Ux2_mag,
                Uy2_mag, Uz2_mag, ncellsx, ncellsy, ncellsz);
            const double cell_sizex = Model.GetXCellSizes()[0];
            const double cell_sizey = Model.GetYCellSizes()[0];
            for (size_t i = 0; i < nmod; ++i)
              {
                const double Volume = cell_sizex * cell_sizey
                    * Model.GetZCellSizes()[i % ncellsz];
                Gradient(i) += std::real((Ux1_el[i] + Ux1_mag[i]) * Ex1_all[i]
                    + (Uy1_el[i] + Uy1_mag[i]) * Ey1_all[i] + (Uz1_el[i]
                    + Uz1_mag[i]) * Ez1_all[i] + (Ux2_el[i] + Ux2_mag[i])
                    * Ex2_all[i] + (Uy2_el[i] + Uy2_mag[i]) * Ey2_all[i]
                    + (Uz2_el[i] + Uz2_mag[i]) * Ez2_all[i]) * Volume;
              }
            //finished with one frequency
          }
        return Gradient;
      }
  }
