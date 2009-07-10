//============================================================================
// Name        : X3DMTCalculator.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

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
        system("x3d");
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
            const size_t freq_index = nobs * i;
            for (size_t j = 0; j < nobs; ++j)
              {
                const size_t obs_index = freq_index + j * 8;
                FieldsToImpedance(Ex1[i], Ex2[i], Ey1[i], Ey2[i], Hx1[i],
                    Hx2[i], Hy1[i], Hy2[i], Zxx, Zxy, Zyx, Zyy);
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
        result(0, 0) = magdet * (conj(A(0, 0)) * H(1, 1) - conj(A(0, 1) * H(0,
            1)));
        result(0, 1) = magdet * (-conj(A(0, 0)) * H(1, 0) + conj(A(0, 1) * H(0,
            0)));
        result(1, 0) = magdet * (conj(A(1, 0)) * H(1, 1) - conj(A(1, 1) * H(0,
            1)));
        result(1, 1) = magdet * (-conj(A(1, 0)) * H(1, 0) + conj(A(1, 1) * H(0,
            0)));
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

    rvec X3DMTCalculator::LQDerivative(const X3DModel &Model,
        const rvec &Misfit)
      {
        const size_t nfreq = Model.GetFrequencies().size();
        const size_t ncellsx = Model.GetConductivities().shape()[0];
        const size_t ncellsy = Model.GetConductivities().shape()[1];
        const size_t nobs = ncellsx * ncellsy;
        const size_t ndata = nobs * nfreq * 8;
        const size_t nmod = ncellsx * ncellsy
            * Model.GetConductivities().shape()[2];
        assert(Misfit.size() == ndata);
        std::cout << "Misfit: " << Misfit << std::endl;
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

            std::vector<std::complex<double> > Ex1_all, Ex2_all, Ey1_all,
                Ey2_all, Ez1_all, Ez2_all;
            ReadEMA(resultfilename + jiba::stringify(i) + "a.ema", Ex1_all,
                Ey1_all, Ez1_all);
            ReadEMA(resultfilename + jiba::stringify(i) + "b.ema", Ex2_all,
                Ey2_all, Ez2_all);
            const size_t freq_index = nobs * i;
            boost::multi_array<std::complex<double>, 2> XPolMoments(
                boost::extents[ncellsx][ncellsy]), YPolMoments(
                boost::extents[ncellsx][ncellsy]);
            boost::multi_array_ref<std::complex<double>, 1> XPolVec(
                XPolMoments.data(), boost::extents[XPolMoments.num_elements()]);
            boost::multi_array_ref<std::complex<double>, 1> YPolVec(
                YPolMoments.data(), boost::extents[YPolMoments.num_elements()]);
            for (size_t j = 0; j < nobs; ++j)
              {
                const size_t siteindex = freq_index + j * 8;
                cmat H(2, 2);
                H(0, 0) = Hx1_obs[j];
                H(0, 1) = Hx2_obs[j];
                H(1, 0) = Hy1_obs[j];
                H(1, 1) = Hy2_obs[j];
                cmat result = CalcATimesH(MisfitToA(Misfit, siteindex), H);
                std::cout << result << std::endl;
                XPolVec[j] = result(0, 0);
                YPolVec[j] = result(0, 1);
              }
            WriteProjectFile(Model.GetFrequencies(), X3DModel::EDIP,
                resultfilename, modelfilename);
            std::string sourcefilename = modelfilename + stringify(i)
                + "a.source";
            WriteSourceFile(sourcefilename, 0.0, XPolMoments, YPolMoments);
            //write an empty source file for the second source polarization
            std::fill(XPolVec.begin(), XPolVec.end(), 0.0);
            std::fill(YPolVec.begin(), YPolVec.end(), 0.0);
            WriteSourceFile(modelfilename + stringify(i) + "b.source", 0.0,
                XPolMoments, YPolMoments);
            boost::filesystem::remove_all(resultfilename+"*");
            system("x3d");
            std::vector<std::complex<double> > Ux1_el, Ux2_el, Uy1_el, Uy2_el,
                Uz1_el, Uz2_el;
            ReadEMA(emaname, Ux1_el, Uy1_el, Uz1_el);

            for (size_t j = 0; j < nobs; ++j)
              {
                const size_t siteindex = freq_index + j * 8;
                cmat H(2, 2);
                H(0, 0) = Hx1_obs[j];
                H(0, 1) = Hx2_obs[j];
                H(1, 0) = Hy1_obs[j];
                H(1, 1) = Hy2_obs[j];
                cmat result = CalcATimesH(MisfitToA(Misfit, siteindex), H);
                XPolVec[j] = result(1, 0);
                YPolVec[j] = result(1, 1);
              }
            WriteSourceFile(sourcefilename, 0.0, XPolMoments, YPolMoments);
            boost::filesystem::remove_all(resultfilename+"*");
            system("x3d");
            ReadEMA(emaname, Ux2_el, Uy2_el, Uz2_el);

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
                cmat H(2, 2);
                H(0, 0) = Hx1_obs[j];
                H(0, 1) = Hx2_obs[j];
                H(1, 0) = Hy1_obs[j];
                H(1, 1) = Hy2_obs[j];
                cmat AH = CalcATimesH(MisfitToA(Misfit, siteindex), H);
                cmat result = ublas::prod(trans(Z), AH);
                XPolVec[j] = result(0, 0) * omega_mu;
                YPolVec[j] = result(0, 1) * omega_mu;
              }
            WriteSourceFile(sourcefilename, 0.0, XPolMoments, YPolMoments);
            boost::filesystem::remove_all(resultfilename+"*");
            system("x3d");
            std::vector<std::complex<double> > Ux1_mag, Ux2_mag, Uy1_mag,
                Uy2_mag, Uz1_mag, Uz2_mag;
            ReadEMA(emaname, Ux1_mag, Uy1_mag, Uz1_mag);
            for (size_t j = 0; j < nobs; ++j)
              {
                cmat Z(2, 2);
                const size_t siteindex = freq_index + j * 8;
                FieldsToImpedance(Ex1_obs[i], Ex2_obs[i], Ey1_obs[i],
                    Ey2_obs[i], Hx1_obs[i], Hx2_obs[i], Hy1_obs[i], Hy2_obs[i],
                    Z(0, 0), Z(0, 1), Z(1, 0), Z(1, 1));
                cmat H(2, 2);
                H(0, 0) = Hx1_obs[j];
                H(0, 1) = Hx2_obs[j];
                H(1, 0) = Hy1_obs[j];
                H(1, 1) = Hy2_obs[j];
                cmat AH = CalcATimesH(MisfitToA(Misfit, siteindex), H);
                cmat result = ublas::prod(trans(Z), AH);
                XPolVec[j] = result(1, 0) * omega_mu;
                YPolVec[j] = result(1, 1) * omega_mu;
              }
            WriteSourceFile(sourcefilename, 0.0, XPolMoments, YPolMoments);
            boost::filesystem::remove_all(resultfilename+"*");
            system("x3d");
            ReadEMA(emaname, Ux2_mag, Uy2_mag, Uz2_mag);
            const double cell_sizex = Model.GetXCellSizes()[0];
            const double cell_sizey = Model.GetYCellSizes()[0];
            for (size_t i = 0; i < nmod; ++i)
              {
                int xindex, yindex, zindex;
                Model.OffsetToIndex(i, xindex, yindex, zindex);
                const double Volume = cell_sizex * cell_sizey * Model.GetZCellSizes()[zindex];
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
