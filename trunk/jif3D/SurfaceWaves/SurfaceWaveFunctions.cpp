/*
 * SurfaceWaveFunctions.cpp
 *
 *  Created on: 20 Sep 2019
 *      Author: bweise
 */

#include <vector>
#include <complex>
#include <tuple>
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Rhumb.hpp>
#include <GeographicLib/Constants.hpp>
#include "../SurfaceWaves/SurfaceWaveFunctions.h"
#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
    typedef std::complex<double> dcomp;
    //const std::complex<double> i(0, 1.0);

    std::vector<double> T2w(const std::vector<double> &periods)
      {
        // Conversion of periods to angular frequencies
        std::vector<double> w(periods.size());
        int nperiods = periods.size();
        for (int n = 0; n < nperiods; n++)
          {
            w[n] = 2.0 * M_PI / periods[n];
          }
        return w;
      }

    std::vector<double> array2vector(const ThreeDModelBase::t3DModelData &array,
        const int &NX, const int &NY, const int &NZ)
      {
        std::vector<double> tmp(NX * NY * NZ);
        for (int i = 0; i < NZ; ++i)
          {
            for (int j = 0; j < NY; ++j)
              {
                for (int k = 0; k < NX; ++k)
                  {
                    //DataVp[k][j][i] = tmp_data[k + j * NX + i * (NX * NY)];
                    tmp[i + j * NZ + k * (NZ * NY)] = array[k][j][i];
                  }
              }
          }
        return tmp;
      }

    double newton_vr(const double &vp, const double &vs, const double &tolerance)
      {
        // variables to store rayleigh velocity
        double vrlast = vs * 0.99;
        double vr = vrlast;
        double diff = tolerance + 99999.0;
        // calculate Rayleigh velocity for homogenous half space
        while (diff > tolerance)
          {
            double fvr = 4.0 - 4.0 * (pow(vrlast, 2) / pow(vs, 2))
                + pow(vrlast, 4) / pow(vs, 4)
                - 4.0 * sqrt(1 - pow(vrlast, 2) / pow(vp, 2))
                    * sqrt(1.0 - pow(vrlast, 2) / pow(vs, 2));
            double dfvr = -8.0 * vrlast / pow(vs, 2) + 4.0 * pow(vrlast, 3) / pow(vs, 4)
                + (4.0 * vrlast * (pow(vp, 2) + pow(vs, 2) - 2.0 * pow(vrlast, 2)))
                    / (vp * vs * sqrt((vp - vrlast) * (vp + vrlast))
                        * sqrt((vs - vrlast) * (vs + vrlast)));
            vr = vrlast - fvr / dfvr;
            diff = sqrt(pow(vr - vrlast, 2));
            vrlast = vr;
          }
        return vr;
      }

    std::tuple<dcomp, dcomp, double, double, double, double, double, dcomp, dcomp, double,
        double, double> compute_util(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const bool &botlay)
      {
        // computes some constants for each layer (vertical wave numbers and some derived properties)
        dcomp hv, kv;
        double SH, CH, SK, CK, mh, mk;
        const double k = w / c;
        if (c < vp)
          {
            mh = sqrt(pow(w / c, 2) - pow(w / vp, 2));
            hv = mh;
            if (botlay == 0)
              {
                SH = (k / mh) * sinh(mh * dn);
                CH = cosh(mh * dn);
              }
          }
        else
          {
            mh = sqrt(pow(w / vp, 2) - pow(w / c, 2));
            hv = mh;
            if (botlay == 0)
              {
                hv = hv * i;
                SH = (k / mh) * sin(mh * dn);
                CH = cos(mh * dn);
              }
          }
        if (c < vs)
          {
            mk = sqrt(pow(w / c, 2) - pow(w / vs, 2));
            kv = mk;
            if (botlay == 0)
              {
                SK = (k / mk) * sinh(mk * dn);
                CK = cosh(mk * dn);
              }
          }
        else
          {
            mk = sqrt(pow(w / vs, 2) - pow(w / c, 2));
            kv = mk;
            if (botlay == 0)
              {
                kv = kv * i;
                SK = (k / mk) * sin(mk * dn);
                CK = cos(mk * dn);
              }
          }
        const double gam = 2.0 * pow(vs, 2) / pow(c, 2);
        const dcomp hvnorm = hv / k;
        const dcomp kvnorm = kv / k;
        const double l = 2.0 * pow(k, 2) - pow(w / vs, 2);

        return std::make_tuple(hv, kv, SH, CH, SK, CK, gam, hvnorm, kvnorm, l, mh, mk);
      }

    std::tuple<dcomp, dcomp, double, double, double, double, double, double, double,
        double, dcomp, dcomp, double, double, double, double, double, double, double,
        double> compute_util_grads(const double &w, const double &vs, const double &vp,
        const double &c, const double &thck, const bool &botlay)
      {
        const auto util = compute_util(w, c, vp, vs, thck, botlay);

        const double SH = std::get<2>(util);
        const double CH = std::get<3>(util);
        const double SK = std::get<4>(util);
        const double CK = std::get<5>(util);
        const double mh = std::get<10>(util);
        const double mk = std::get<11>(util);

        dcomp hv_h, kv_k, hv_c, kv_c;
        double mkk, mhh, mhc, mkc;
        const double k = w / c;
        if (c < vp)
          {
            mhh = pow(w, 2) / (mh * pow(vp, 3));
            hv_h = mhh;
            mhc = (-1.0) * (pow(w, 2) / (mh * pow(c, 3)));
            hv_c = mhc;
          }
        else
          {
            mhh = ((-1.0) * pow(w, 2)) / (mh * pow(vp, 3));
            hv_h = mhh;
            mhc = pow(w, 2) / (mh * pow(c, 3));
            hv_c = mhc;
            if (botlay == 0)
              {
                hv_h = i * mhh;
                hv_c = i * mhc;
              }
          }

        if (c < vs)
          {
            mkk = pow(w, 2) / (mk * pow(vs, 3));
            kv_k = mkk;
            mkc = (-1.0) * (pow(w, 2) / (mk * pow(c, 3)));
            kv_c = mkc;
          }
        else
          {
            mkk = ((-1.0) * pow(w, 2)) / (mk * pow(vs, 3));
            kv_k = mkk;
            mkc = pow(w, 2) / (mk * pow(c, 3));
            kv_c = mkc;
            if (botlay == 0)
              {
                kv_k = i * mkk;
                kv_c = i * mkc;
              }
          }

        const double hvnorm_c = (-1.0) * (2.0 * c) / pow(vp, 2);
        const double hvnorm_h = (2.0 * pow(c, 2)) / pow(vp, 3);
        const double kvnorm_c = (-1.0) * (2.0 * c) / pow(vs, 2);
        const double kvnorm_k = (2.0 * pow(c, 2)) / pow(vs, 3);
        const double gam_k = (4.0 * vs) / pow(c, 2);
        const double gam_c = (-4.0 * pow(vs, 2)) / pow(c, 3);
        const double l_k = 2.0 * pow(w, 2) / pow(vs, 3);
        const double l_c = -4.0 * pow(w, 2) / pow(c, 3);
        const double CHH = (w * c * thck * SH) / pow(vp, 3);
        const double CKK = (w * c * thck * SK) / pow(vs, 3);
        const double CHC = ((-1.0) * k * thck * SH) / c;
        const double CKC = ((-1.0) * k * thck * SK) / c;
        const double SHC = (-1.0) * ((mhc / mh) + (1 / c)) * SH
            + (k * thck * mhc * CH / mh);
        const double SKC = (-1.0) * ((mkc / mk) + (1 / c)) * SK
            + (k * thck * mkc * CK / mk);
        const double SHH = (mhh / mh) * (k * thck * CH - SH);
        const double SKK = (mkk / mk) * (k * thck * CK - SK);

        return std::make_tuple(hv_h, kv_k, hvnorm_h, kvnorm_k, gam_k, l_k, CHH, CKK, SHH,
            SKK, hv_c, kv_c, hvnorm_c, kvnorm_c, gam_c, l_c, CHC, CKC, SHC, SKC);
      }

    std::tuple<double, double, double, double, double> compute_T(const double &w,
        const double &c, const double &vp, const double &vs, const double &mu)
      {
        // computes layer matrix for bottom layer (an i denotes an imaginary subdeterminant e.g. iT1214)
        const double k = w / c;
        const auto util = compute_util(w, c, vp, vs, 99999.0, 1);
        const dcomp hv = std::get<0>(util);
        const dcomp kv = std::get<1>(util);
        const double l = std::get<9>(util);

        const dcomp fact = pow((-1.0) * pow(vs, 2) / (2 * mu * hv * kv * pow(w, 2)), 2);

        const double T1212 = std::real(
            pow(mu, 2) * kv * hv * (pow(l, 2) - 4.0 * pow(k, 2) * kv * hv) * fact);
        const double T1213 = std::real(
            mu * pow(hv, 2) * kv * (l - 2.0 * pow(k, 2)) * fact);
        const double iT1214 = std::real(k * mu * hv * kv * (l - 2.0 * hv * kv) * fact);
        const double T1224 = std::real(
            mu * hv * pow(kv, 2) * (2.0 * pow(k, 2) - l) * fact);
        const double T1234 = std::real(hv * kv * (pow(k, 2) - hv * kv) * fact);

        return std::make_tuple(T1212, T1213, iT1214, T1224, T1234);
      }

    std::tuple<double, double, double, double, double> compute_T_vs(const double &w,
        const double &c, const double &vp, const double &vs, const double &mu)
      {
        const auto T = compute_T(w, c, vp, vs, mu);
        const double T1212 = std::get<0>(T);
        const double iT1214 = std::get<2>(T);

        const auto util = compute_util(w, c, vp, vs, 99999.0, 1);
        const dcomp hv = std::get<0>(util);
        const dcomp kv = std::get<1>(util);
        const double l = std::get<9>(util);

        const auto util_grads = compute_util_grads(w, vs, vp, c, 99999.0, 1);
        const dcomp kv_k = std::get<1>(util_grads);
        const double l_k = std::get<5>(util_grads);

        const double gT1212 = std::real(
            (4.0 * T1212 / vs)
                + pow(vs, 4) * (2.0 * l * l_k * kv - pow(l, 2) * kv_k)
                    / (4.0 * pow(w, 4) * hv * pow(kv, 2)));
        const double gT1213 = std::real(
            kv_k * pow(vs, 2) / (4.0 * mu * pow(w, 2) * pow(kv, 2)));
        const double igT1214 = std::real(
            (2.0 * iT1214 / vs)
                + pow(vs, 4) * (l_k * kv - l * kv_k)
                    / (4.0 * mu * pow(w, 3) * c * hv * pow(kv, 2)));
        const double gT1224 = 0.0;
        const double gT1234 = std::real(
            (-1.0) * kv_k * pow(vs, 4)
                / (4.0 * pow(mu, 2) * pow(w, 2) * pow(c, 2) * hv * pow(kv, 2)));

        return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
      }

    std::tuple<double, double, double, double, double> compute_T_vp(const double &w,
        const double &c, const double &vp, const double &vs, const double &mu)
      {
        const auto util = compute_util(w, c, vp, vs, 99999.0, 1);
        const dcomp hv = std::get<0>(util);
        const dcomp kv = std::get<1>(util);
        const double l = std::get<9>(util);

        const auto util_grads = compute_util_grads(w, vs, vp, c, 99999.0, 1);
        const dcomp hv_h = std::get<0>(util_grads);

        const double gT1212 = std::real(
            (-1.0) * pow(vs, 4) * pow(l, 2) * hv_h / (4.0 * pow(w, 4) * pow(hv, 2) * kv));
        const double gT1213 = 0.0;
        const double igT1214 = std::real(
            (-1.0) * pow(vs, 4) * l * hv_h
                / (4.0 * mu * pow(w, 3) * c * pow(hv, 2) * kv));
        const double gT1224 = std::real(
            (-1.0) * hv_h * pow(vs, 2) / (4.0 * mu * pow(hv, 2)));
        const double gT1234 = std::real(
            (-1.0) * hv_h * pow(vs, 4)
                / (4.0 * pow(mu, 2) * pow(w, 2) * pow(c, 2) * pow(hv, 2) * kv));

        return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
      }

    std::tuple<double, double, double, double, double> compute_T_rho(const double &w,
        const double &c, const double &vp, const double &vs, const double &mu)
      {
        const auto T = compute_T(w, c, vp, vs, mu);
        const double T1213 = std::get<1>(T);
        const double iT1214 = std::get<2>(T);
        const double T1224 = std::get<3>(T);
        const double T1234 = std::get<4>(T);

        const double gT1212 = 0.0;
        const double gT1213 = (-1.0) * T1213 * pow(vs, 2) / mu;
        const double igT1214 = (-1.0) * iT1214 * pow(vs, 2) / mu;
        const double gT1224 = (-1.0) * T1224 * pow(vs, 2) / mu;
        const double gT1234 = (-2.0) * T1234 * pow(vs, 2) / mu;

        return std::make_tuple(gT1212, gT1213, igT1214, gT1224, gT1234);
      }

    std::vector<double> compute_T_c(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu)
      {
        auto util = compute_util(w, c, vp, vs, 99999.0, 1);
        const dcomp hv = std::get<0>(util);
        const dcomp kv = std::get<1>(util);
        const double l = std::get<9>(util);

        auto util_g = compute_util_grads(w, vs, vp, c, 99999.0, 1);
        const dcomp hv_c = std::get<10>(util_g);
        const dcomp kv_c = std::get<11>(util_g);
        const double l_c = std::get<15>(util_g);

        const auto T = compute_T(w, c, vp, vs, mu);
        const double iT1214 = std::get<2>(T);

        const double dens = mu / pow(vs, 2);

        const double gT1212 = std::real(
            (pow(vs, 4) / 4.0 * pow(w, 4))
                * (((2.0 * l * l_c * hv * kv - pow(l, 2) * (hv_c * kv + hv * kv_c))
                    / (pow(hv * kv, 2))) + (8.0 * pow(w, 2) / pow(c, 3))));
        const double gT1213 = std::real(kv_c / (4.0 * dens * pow(w * kv, 2)));
        const double igT1214 = std::real(
            ((-1.0) * iT1214 / c)
                + ((pow(vs, 2) * (l_c * kv * hv - l * (hv_c * kv + hv * kv_c)))
                    / (4.0 * dens * pow(w, 3) * c * pow(kv * hv, 2))));
        const double gT1224 = std::real(((-1.0) * hv_c) / (4.0 * dens * pow(hv * w, 2)));
        const double gT1234 = std::real(
            (-2.0 * hv * kv + c * (hv_c * kv + hv * kv_c))
                / (4.0 * pow(dens * w * hv * kv, 2) * pow(c, 3)));

        std::vector<double> gT(5);
        gT[0] = gT1212;
        gT[1] = gT1213;
        gT[2] = igT1214;
        gT[3] = gT1224;
        gT[4] = gT1234;
        return gT;
      }

    std::vector<double> compute_G(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens)
      {
        // computes subdeterminants of G matrix (an i denotes an imaginary subdeterminant e.g. iG1214)
        const auto kvert = compute_util(w, c, vp, vs, dn, 0);
        const double SH = std::get<2>(kvert);
        const double CH = std::get<3>(kvert);
        const double SK = std::get<4>(kvert);
        const double CK = std::get<5>(kvert);
        const double gam = std::get<6>(kvert);
        const dcomp hvnorm = std::get<7>(kvert);
        const dcomp kvnorm = std::get<8>(kvert);

        const double G1212 = std::real(
            2.0 * gam * (1.0 - gam) + (2.0 * pow(gam, 2) - 2.0 * gam + 1.0) * CH * CK
                - (pow(1.0 - gam, 2) + pow(gam, 2) * pow(hvnorm, 2) * pow(kvnorm, 2)) * SH
                    * SK);
        const double G1213 = std::real(
            (1.0 / (dens * w * c)) * (CH * SK - SH * CK * pow(hvnorm, 2)));
        const double iG1214 = std::real(
            (1.0 / (dens * w * c))
                * ((1.0 - 2.0 * gam) * (1.0 - CK * CH)
                    + (1.0 - gam - gam * pow(hvnorm, 2) * pow(kvnorm, 2)) * SH * SK));
        const double G1224 = std::real(
            (1.0 / (dens * w * c)) * (pow(kvnorm, 2) * CH * SK - SH * CK));
        const double G1234 = std::real(
            (-1.0 / (pow(dens, 2) * pow(w, 2) * pow(c, 2)))
                * (2.0 * (1.0 - CH * CK)
                    + (1.0 + pow(kvnorm, 2) * pow(hvnorm, 2)) * SK * SH));
        const double G1312 = std::real(
            dens * w * c
                * (pow(gam, 2) * pow(kvnorm, 2) * CH * SK - pow(1.0 - gam, 2) * SH * CK));
        const double G1313 = std::real(CH * CK);
        const double iG1314 = std::real(
            (1.0 - gam) * SH * CK + gam * pow(kvnorm, 2) * CH * SK);
        const double G1324 = std::real((-1.0) * pow(kvnorm, 2) * SH * SK);
        const double iG1412 = std::real(
            dens * w * c
                * ((3.0 * pow(gam, 2) - 2.0 * pow(gam, 3) - gam) * (1.0 - CH * CK)
                    + (pow(1.0 - gam, 3) - pow(gam, 3) * pow(hvnorm, 2) * pow(kvnorm, 2))
                        * SH * SK));
        const double iG1413 = std::real(
            (-1.0) * ((1.0 - gam) * CH * SK + gam * pow(hvnorm, 2) * SH * CK));
        const double G1414 = std::real(
            1.0 - 2.0 * gam * (1.0 - gam) * (1.0 - CH * CK)
                + (pow(1.0 - gam, 2) + pow(gam, 2) * pow(kvnorm, 2) * pow(hvnorm, 2)) * SH
                    * SK);
        const double G2412 = std::real(
            dens * w * c
                * (pow(1.0 - gam, 2) * CH * SK - pow(gam, 2) * SH * CK * pow(hvnorm, 2)));
        const double G2413 = std::real((-1.0) * pow(hvnorm, 2) * SH * SK);
        const double G3412 = std::real(
            (-1.0) * pow(dens, 2) * pow(w, 2) * pow(c, 2)
                * (2.0 * pow(gam, 2) * pow(1.0 - gam, 2) * (1.0 - CH * CK)
                    + (pow(1.0 - gam, 4) + pow(gam, 4) * pow(hvnorm, 2) * pow(kvnorm, 2))
                        * SH * SK));

        std::vector<double> G(15);
        G[0] = G1212;
        G[1] = G1213;
        G[2] = iG1214;
        G[3] = G1224;
        G[4] = G1234;
        G[5] = G1312;
        G[6] = G1313;
        G[7] = iG1314;
        G[8] = G1324;
        G[9] = iG1412;
        G[10] = iG1413;
        G[11] = G1414;
        G[12] = G2412;
        G[13] = G2413;
        G[14] = G3412;
        return G;
      }

    std::vector<double> compute_G_vs(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens)
      {
        const auto kvert = compute_util(w, c, vp, vs, dn, 0);
        const double SH = std::get<2>(kvert);
        const double CH = std::get<3>(kvert);
        const double SK = std::get<4>(kvert);
        const double CK = std::get<5>(kvert);
        const double gam = std::get<6>(kvert);
        const dcomp hvnorm = std::get<7>(kvert);
        const dcomp kvnorm = std::get<8>(kvert);

        const auto util_grads = compute_util_grads(w, vs, vp, c, dn, 0);
        const double kvnorm_k = std::get<3>(util_grads);
        const double gam_k = std::get<4>(util_grads);
        const double CKK = std::get<7>(util_grads);
        const double SKK = std::get<9>(util_grads);

        const double gG1212 = std::real(
            2.0 * gam_k * (1.0 - 2.0 * gam) * (1.0 - CH * CK)
                + (2.0 * pow(gam, 2) - 2.0 * gam + 1.0) * CH * CKK
                - (2.0 * gam_k * (gam - 1.0)
                    + 2.0 * gam * gam_k * pow(hvnorm, 2) * pow(kvnorm, 2)
                    + pow(gam, 2) * pow(hvnorm, 2) * kvnorm_k) * SH * SK
                - (pow((1.0 - gam), 2) + pow(gam, 2) * pow(hvnorm, 2) * pow(kvnorm, 2))
                    * SH * SKK);
        const double gG1213 = std::real(
            (1.0 / (dens * w * c)) * (CH * SKK - pow(hvnorm, 2) * SH * CKK));
        const double igG1214 =
            std::real(
                (1.0 / (dens * w * c))
                    * ((-2.0) * gam_k * (1.0 - CH * CK) - (1.0 - 2.0 * gam) * CH * CKK)
                    + (1.0 / (dens * w * c))
                        * (((-1.0) * gam_k - gam_k * pow(hvnorm, 2) * pow(kvnorm, 2)
                            - gam * pow(hvnorm, 2) * kvnorm_k) * SH * SK
                            + (1.0 - gam - gam * pow(hvnorm, 2) * pow(kvnorm, 2)) * SH
                                * SKK));
        const double gG1224 = std::real(
            (1.0 / (dens * w * c))
                * (kvnorm_k * CH * SK + pow(kvnorm, 2) * CH * SKK - SH * CKK));
        const double gG1234 = std::real(
            (1.0 / pow((dens * w * c), 2))
                * (-2.0 * CH * CKK + (1.0 + pow(hvnorm, 2) * pow(kvnorm, 2)) * SH * SKK
                    + pow(hvnorm, 2) * kvnorm_k * SH * SK));
        const double gG1312 = std::real(
            dens * w * c
                * (2.0 * gam * gam_k * pow(kvnorm, 2) * CH * SK
                    + pow(gam, 2) * kvnorm_k * CH * SK
                    + pow(gam, 2) * pow(kvnorm, 2) * CH * SKK)
                + dens * w * c
                    * (2.0 * gam_k * (1.0 - gam) * SH * CK
                        - pow((1.0 - gam), 2) * SH * CKK));
        const double gG1313 = std::real(CH * CKK);
        const double igG1314 = std::real(
            (-1.0) * gam_k * SH * CK + (1.0 - gam) * SH * CKK
                + gam_k * pow(kvnorm, 2) * CH * SK + gam * kvnorm_k * CH * SK
                + gam * pow(kvnorm, 2) * CH * SKK);
        const double gG1324 = std::real(
            (-1.0) * kvnorm_k * SH * SK - pow(kvnorm, 2) * SH * SKK);
        const double igG1412 =
            std::real(
                dens * w * c
                    * (gam_k * (-6.0 * pow(gam, 2) + 6.0 * gam - 1.0) * (1.0 - CH * CK)
                        - (pow(gam, 2) - gam) * (1.0 - 2.0 * gam) * CH * CKK)
                    + dens * w * c
                        * (-3.0 * pow((1.0 - gam), 2) * gam_k
                            - 3.0 * pow(gam, 2) * gam_k * pow(hvnorm, 2) * pow(kvnorm, 2)
                            - pow(gam, 3) * pow(hvnorm, 2) * kvnorm_k) * SH * SK
                    + dens * w * c
                        * (pow((1.0 - gam), 3)
                            - pow(gam, 3) * pow(hvnorm, 2) * pow(kvnorm, 2)) * SH * SKK);
        const double igG1413 =
            std::real(
                (-1.0)
                    * ((1.0 - gam) * CH * SKK - gam_k * CH * SK
                        + gam_k * pow(hvnorm, 2) * SH * CK
                        + gam * pow(hvnorm, 2) * SH * CKK));
        const double gG1414 = std::real(
            2.0 * gam_k * (2.0 * gam - 1.0) * (1.0 - CH * CK)
                - 2.0 * (pow(gam, 2) - gam) * CH * CKK
                + (2.0 * (gam - 1) * gam_k
                    + 2.0 * gam * gam_k * pow(hvnorm, 2) * pow(kvnorm, 2)
                    + pow(gam, 2) * pow(hvnorm, 2) * kvnorm_k) * SH * SK
                + (pow((1.0 - gam), 2) + pow(gam, 2) * pow(hvnorm, 2) * pow(kvnorm, 2))
                    * SH * SKK);
        const double gG2412 = std::real(
            dens * w * c
                * (2.0 * (gam - 1.0) * gam_k * CH * SK + pow((1.0 - gam), 2) * CH * SKK
                    - 2.0 * gam * gam_k * pow(hvnorm, 2) * SH * CK
                    - pow(gam, 2) * pow(hvnorm, 2) * SH * CKK));
        const double gG2413 = std::real((-1.0) * pow(hvnorm, 2) * SH * SKK);
        const double gG3412 =
            std::real(
                (-1.0) * pow((dens * c * w), 2)
                    * (4.0 * gam * gam_k * (1.0 + 2.0 * pow(gam, 2) - 3.0 * gam)
                        * (1.0 - CH * CK)
                        - 2.0 * pow(gam, 2) * pow((1.0 - gam), 2) * CH * CKK)
                    - pow((dens * c * w), 2)
                        * ((pow((1.0 - gam), 4)
                            + pow(gam, 4) * pow(hvnorm, 2) * pow(kvnorm, 2)) * SH * SKK)
                    - pow((dens * c * w), 2)
                        * ((-4.0 * gam_k * pow((1.0 - gam), 3)
                            + 4.0 * pow(gam, 3) * gam_k * pow(hvnorm, 2) * pow(kvnorm, 2)
                            + pow(gam, 4) * pow(hvnorm, 2) * kvnorm_k) * SH * SK));

        std::vector<double> G(15);
        G[0] = gG1212;
        G[1] = gG1213;
        G[2] = igG1214;
        G[3] = gG1224;
        G[4] = gG1234;
        G[5] = gG1312;
        G[6] = gG1313;
        G[7] = igG1314;
        G[8] = gG1324;
        G[9] = igG1412;
        G[10] = igG1413;
        G[11] = gG1414;
        G[12] = gG2412;
        G[13] = gG2413;
        G[14] = gG3412;
        return G;
      }

    std::vector<double> compute_G_vp(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens)
      {
        const auto kvert = compute_util(w, c, vp, vs, dn, 0);
        const double SH = std::get<2>(kvert);
        const double SK = std::get<4>(kvert);
        const double CK = std::get<5>(kvert);
        const double gam = std::get<6>(kvert);
        const dcomp hvnorm = std::get<7>(kvert);
        const dcomp kvnorm = std::get<8>(kvert);

        const auto util_grads = compute_util_grads(w, vs, vp, c, dn, 0);
        const double hvnorm_h = std::get<2>(util_grads);
        const double CHH = std::get<6>(util_grads);
        const double SHH = std::get<8>(util_grads);

        const double gG1212 = std::real(
            (2.0 * pow(gam, 2) - 2.0 * gam + 1.0) * CHH * CK
                - pow(gam, 2) * hvnorm_h * pow(kvnorm, 2) * SH * SK
                - (pow((1.0 - gam), 2) + pow(gam, 2) * pow(hvnorm, 2) * pow(kvnorm, 2))
                    * SHH * SK);
        const double gG1213 = std::real(
            (1.0 / (dens * w * c))
                * (CHH * SK - hvnorm_h * SH * CK - pow(hvnorm, 2) * SHH * CK));
        const double igG1214 = std::real(
            (1.0 / (dens * w * c))
                * ((2.0 * gam - 1.0) * CHH * CK
                    + (1.0 - gam - gam * pow(hvnorm, 2) * pow(kvnorm, 2)) * SHH * SK
                    - gam * hvnorm_h * pow(kvnorm, 2) * SH * SK));
        const double gG1224 = std::real(
            (1.0 / (dens * w * c)) * (pow(kvnorm, 2) * CHH * SK - SHH * CK));
        const double gG1234 = std::real(
            ((-1.0) / pow((w * dens * c), 2))
                * (-2.0 * CHH * CK + (1.0 + pow(kvnorm, 2) * pow(hvnorm, 2)) * SHH * SK
                    + hvnorm_h * pow(kvnorm, 2) * SH * SK));
        const double gG1312 = std::real(
            dens * w * c
                * (pow(gam, 2) * pow(kvnorm, 2) * CHH * SK
                    - pow((1.0 - gam), 2) * SHH * CK));
        const double gG1313 = std::real(CHH * CK);
        const double igG1314 = std::real(
            (1.0 - gam) * SHH * CK + gam * pow(kvnorm, 2) * CHH * SK);
        const double gG1324 = std::real((-1.0) * pow(kvnorm, 2) * SHH * SK);
        const double igG1412 = std::real(
            dens * w * c
                * (-1.0 * gam * (gam - 1.0) * (1.0 - 2.0 * gam) * CHH * CK
                    - pow(gam, 3) * hvnorm_h * pow(kvnorm, 2) * SH * SK)
                + dens * w * c
                    * ((pow((1.0 - gam), 3)
                        - pow(gam, 3) * pow(hvnorm, 2) * pow(kvnorm, 2)) * SHH * SK));
        const double igG1413 = std::real(
            (gam - 1.0) * CHH * SK - gam * hvnorm_h * SH * CK
                - gam * pow(hvnorm, 2) * SHH * CK);
        const double gG1414 = std::real(
            2.0 * gam * (1.0 - gam) * CHH * CK
                + (pow((1.0 - gam), 2) + pow(gam, 2) * pow(hvnorm, 2) * pow(kvnorm, 2))
                    * SHH * SK + pow(gam, 2) * hvnorm_h * pow(kvnorm, 2) * SH * SK);
        const double gG2412 = std::real(
            dens * w * c
                * (pow((1.0 - gam), 2) * CHH * SK - pow(gam, 2) * hvnorm_h * SH * CK
                    - pow(gam, 2) * pow(hvnorm, 2) * SHH * CK));
        const double gG2413 = std::real(
            (-1.0) * (hvnorm_h * SH + pow(hvnorm, 2) * SHH) * SK);
        const double gG3412 = std::real(
            (-1.0) * pow((dens * c * w), 2)
                * (-2.0 * pow(gam, 2) * pow((1 - gam), 2) * CHH * CK)
                - pow((dens * c * w), 2)
                    * ((pow((1 - gam), 4) + pow(gam, 4) * pow(hvnorm, 2) * pow(kvnorm, 2))
                        * SHH * SK + pow(gam, 4) * hvnorm_h * pow(kvnorm, 2) * SH * SK));

        std::vector<double> G(15);
        G[0] = gG1212;
        G[1] = gG1213;
        G[2] = igG1214;
        G[3] = gG1224;
        G[4] = gG1234;
        G[5] = gG1312;
        G[6] = gG1313;
        G[7] = igG1314;
        G[8] = gG1324;
        G[9] = igG1412;
        G[10] = igG1413;
        G[11] = gG1414;
        G[12] = gG2412;
        G[13] = gG2413;
        G[14] = gG3412;
        return G;
      }

    std::vector<double> compute_G_rho(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens)
      {
        const auto G = compute_G(c, dn, w, vp, vs, dens);
        const double G1213 = G[1];
        const double iG1214 = G[2];
        const double G1224 = G[3];
        const double G1234 = G[4];
        const double G1312 = G[5];
        const double iG1412 = G[9];
        const double G2412 = G[12];
        const double G3412 = G[14];

        const double gG1212 = 0.0;
        const double gG1213 = std::real((-1.0) * G1213 / dens);
        const double igG1214 = std::real((-1.0) * iG1214 / dens);
        const double gG1224 = std::real((-1.0) * G1224 / dens);
        const double gG1234 = std::real((-2.0) * G1234 / dens);
        const double gG1312 = std::real(G1312 / dens);
        const double gG1313 = 0.0;
        const double igG1314 = 0.0;
        const double gG1324 = 0.0;
        const double igG1412 = std::real(iG1412 / dens);
        const double igG1413 = 0.0;
        const double gG1414 = 0.0;
        const double gG2412 = std::real(G2412 / dens);
        const double gG2413 = 0.0;
        const double gG3412 = std::real(2.0 * G3412 / dens);

        std::vector<double> Gout(15);
        Gout[0] = gG1212;
        Gout[1] = gG1213;
        Gout[2] = igG1214;
        Gout[3] = gG1224;
        Gout[4] = gG1234;
        Gout[5] = gG1312;
        Gout[6] = gG1313;
        Gout[7] = igG1314;
        Gout[8] = gG1324;
        Gout[9] = igG1412;
        Gout[10] = igG1413;
        Gout[11] = gG1414;
        Gout[12] = gG2412;
        Gout[13] = gG2413;
        Gout[14] = gG3412;
        return Gout;
      }

    std::vector<double> compute_G_c(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens)
      {
        const auto utils = compute_util(w, c, vp, vs, dn, 0);
        const double SH = std::get<2>(utils);
        const double CH = std::get<3>(utils);
        const double SK = std::get<4>(utils);
        const double CK = std::get<5>(utils);
        const double gam = std::get<6>(utils);
        const dcomp hvnorm = std::get<7>(utils);
        const dcomp kvnorm = std::get<8>(utils);

        const auto utils_c = compute_util_grads(w, vs, vp, c, dn, 0);
        const double SHC = std::get<18>(utils_c);
        const double CHC = std::get<16>(utils_c);
        const double SKC = std::get<19>(utils_c);
        const double CKC = std::get<17>(utils_c);
        const double gam_c = std::get<14>(utils_c);
        const double hvnorm_c = std::get<12>(utils_c);
        const double kvnorm_c = std::get<13>(utils_c);

        const double dCC = CHC * CK + CH * CKC;
        const double dSS = SHC * SK + SH * SKC;
        const double dCS = CHC * SK + CH * SKC;
        const double dSC = SHC * CK + SH * CKC;
        const double dk = std::real(
            hvnorm_c * pow(kvnorm, 2) + pow(hvnorm, 2) * kvnorm_c);

        const auto G = compute_G(c, dn, w, vp, vs, dens);
        const double G1213 = G[1];
        const double iG1214 = G[2];
        const double G1224 = G[3];
        const double G1234 = G[4];
        const double G1312 = G[5];
        const double iG1412 = G[9];
        const double G2412 = G[12];
        const double G3412 = G[14];

        const double gG1212 = std::real(
            2.0 * gam_c * (1.0 - 2.0 * gam) * (1.0 - CH * CK)
                + (2.0 * pow(gam, 2) - 2.0 * gam + 1.0) * dCC
                - (2.0 * gam_c * (gam - 1.0 + gam * pow(hvnorm * kvnorm, 2))
                    + pow(gam, 2) * dk) * SH * SK
                - (pow(1.0 - gam, 2) + pow(gam * hvnorm * kvnorm, 2)) * dSS);
        const double gG1213 = std::real(
            ((-1.0) * G1213 / c)
                + (1.0 / (dens * w * c))
                    * (dCS - hvnorm_c * SH * CK - pow(hvnorm, 2) * dSC));
        const double igG1214 =
            std::real(
                ((-1.0) * iG1214 / c)
                    + (1.0 / (dens * w * c))
                        * ((2.0 * gam - 1.0) * dCC - 2.0 * gam_c * (1.0 - CH * CK))
                    + (1.0 / (dens * w * c))
                        * ((1.0 - gam - gam * pow(hvnorm * kvnorm, 2)) * dSS
                            - (gam_c + gam_c * pow(hvnorm * kvnorm, 2) + gam * dk) * SH
                                * SK));
        const double gG1224 = std::real(
            ((-1.0) * G1224 / c)
                + (1.0 / (dens * w * c))
                    * (kvnorm_c * CH * SK + pow(kvnorm, 2) * dCS - dSC));
        const double gG1234 = std::real(
            ((-2.0) * G1234 / c)
                + pow((1.0 / (dens * w * c)), 2)
                    * (2.0 * dCC - dk * SH * SK - (1.0 + pow(hvnorm * kvnorm, 2)) * dSS));
        const double gG1312 = std::real(
            (G1312 / c)
                + dens * w * c
                    * (2.0 * gam * gam_c * pow(kvnorm, 2) * CH * SK
                        + pow(gam, 2) * kvnorm_c * CH * SK + pow(gam * kvnorm, 2) * dCS)
                + dens * w * c
                    * (2.0 * gam_c * (1.0 - gam) * SH * CK - pow(1.0 - gam, 2) * dSC));
        const double gG1313 = std::real(dCC);
        const double igG1314 = std::real(
            (-1.0) * gam_c * SH * CK + (1.0 - gam) * dSC
                + (gam_c * pow(kvnorm, 2) + gam * kvnorm_c) * CH * SK
                + gam * pow(kvnorm, 2) * dCS);
        const double gG1324 = std::real(
            (-1.0) * (kvnorm_c * SH * SK + pow(kvnorm, 2) * dSS));
        const double igG1412 = std::real(
            (iG1412 / c)
                + dens * w * c
                    * (gam_c * (-6.0 * pow(gam, 2) + 6.0 * gam - 1.0) * (1.0 - CH * CK)
                        - (pow(gam, 2) - gam) * (1.0 - 2.0 * gam) * dCC)
                + dens * c * w
                    * ((pow(1.0 - gam, 3) - pow(gam, 3) * pow(hvnorm * kvnorm, 2)) * dSS
                        + (-3.0 * gam_c * pow(1.0 - gam, 2)
                            - 3.0 * pow(gam, 2) * gam_c * pow(hvnorm * kvnorm, 2)
                            - pow(gam, 3) * dk) * SH * SK));
        const double igG1413 = std::real(
            gam_c * CH * SK - (1.0 - gam) * dCS
                - (gam_c * pow(hvnorm, 2) + gam * hvnorm_c) * SH * CK
                - gam * pow(hvnorm, 2) * dSC);
        const double gG1414 = std::real(
            2.0 * gam * (1.0 - gam) * dCC
                + 2.0 * gam_c * (2.0 * gam - 1.0) * (1.0 - CH * CK)
                + (pow(1 - gam, 2) + pow(gam * hvnorm * kvnorm, 2)) * dSS
                + (2 * gam_c * (gam - 1.0 + gam * pow(hvnorm * kvnorm, 2))
                    + pow(gam, 2) * dk) * SH * SK);
        const double gG2412 = std::real(
            (G2412 / c)
                + dens * w * c
                    * (-2.0 * (1.0 - gam) * gam_c * CH * SK + pow(1.0 - gam, 2) * dCS)
                + dens * w * c
                    * (-1.0
                        * (2.0 * gam * gam_c * pow(hvnorm, 2) + pow(gam, 2) * hvnorm_c)
                        * SH * CK - pow(gam * hvnorm, 2) * dSC));
        const double gG2413 = std::real(
            (-1.0) * (hvnorm_c * SH * SK + pow(hvnorm, 2) * dSS));
        const double gG3412 = std::real(
            (2.0 * G3412 / c)
                - pow(dens * c * w, 2)
                    * (4.0 * gam * gam_c * (gam - 1.0) * (2.0 * gam - 1.0)
                        * (1.0 - CH * CK) - 2.0 * pow(gam, 2) * pow(1.0 - gam, 2) * dCC)
                - pow(dens * c * w, 2)
                    * ((pow(1.0 - gam, 4) + pow(gam, 4) * pow(hvnorm * kvnorm, 2)) * dSS
                        + (-4.0 * gam_c * pow(1.0 - gam, 3)
                            + 4.0 * pow(gam, 3) * gam_c * pow(hvnorm * kvnorm, 2)
                            + pow(gam, 4) * dk) * SH * SK));

        std::vector<double> Gout(15);
        Gout[0] = gG1212;
        Gout[1] = gG1213;
        Gout[2] = igG1214;
        Gout[3] = gG1224;
        Gout[4] = gG1234;
        Gout[5] = gG1312;
        Gout[6] = gG1313;
        Gout[7] = igG1314;
        Gout[8] = gG1324;
        Gout[9] = igG1412;
        Gout[10] = igG1413;
        Gout[11] = gG1414;
        Gout[12] = gG2412;
        Gout[13] = gG2413;
        Gout[14] = gG3412;
        return Gout;
      }

    std::vector<double> compute_R_c(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const double &dens,
        const std::tuple<double, double, double, double, double> &T,
        const std::vector<double> &T_c)
      {

        const double T1212 = std::get<0>(T);
        const double T1213 = std::get<1>(T);
        const double iT1214 = std::get<2>(T);
        const double T1224 = std::get<3>(T);
        const double T1234 = std::get<4>(T);

        const double T1212_c = T_c[0];
        const double T1213_c = T_c[1];
        const double iT1214_c = T_c[2];
        const double T1224_c = T_c[3];
        const double T1234_c = T_c[4];

        std::vector<double> G = compute_G(c, dn, w, vp, vs, dens);
        const double G1212 = G[0];
        const double G1213 = G[1];
        const double iG1214 = G[2];
        const double G1224 = G[3];
        const double G1234 = G[4];
        const double G1312 = G[5];
        const double G1313 = G[6];
        const double iG1314 = G[7];
        const double G1324 = G[8];
        const double iG1412 = G[9];
        const double iG1413 = G[10];
        const double G1414 = G[11];
        const double G2412 = G[12];
        const double G2413 = G[13];
        const double G3412 = G[14];

        std::vector<double> G_c = compute_G_c(c, dn, w, vp, vs, dens);
        const double G1212_c = G_c[0];
        const double G1213_c = G_c[1];
        const double iG1214_c = G_c[2];
        const double G1224_c = G_c[3];
        const double G1234_c = G_c[4];
        const double G1312_c = G_c[5];
        const double G1313_c = G_c[6];
        const double iG1314_c = G_c[7];
        const double G1324_c = G_c[8];
        const double iG1412_c = G_c[9];
        const double iG1413_c = G_c[10];
        const double G1414_c = G_c[11];
        const double G2412_c = G_c[12];
        const double G2413_c = G_c[13];
        const double G3412_c = G_c[14];

        double R1212 = T1212_c * G1212 + T1212 * G1212_c + T1213_c * G1312
            + T1213 * G1312_c - 2.0 * iT1214_c * iG1412 - 2.0 * iT1214 * iG1412_c
            + T1224_c * G2412 + T1224 * G2412_c + T1234_c * G3412 + T1234 * G3412_c;
        double R1213 = T1212_c * G1213 + T1212 * G1213_c + T1213_c * G1313
            + T1213 * G1313_c - 2.0 * iT1214_c * iG1413 - 2.0 * iT1214 * iG1413_c
            + T1224_c * G2413 + T1224 * G2413_c + T1234_c * G2412 + T1234 * G2412_c;
        double iR1214 = T1212_c * iG1214 + T1212 * iG1214_c + T1213_c * iG1314
            + T1213 * iG1314_c + iT1214_c * (2.0 * G1414 - 1.0) + iT1214 * (2.0 * G1414_c)
            + T1224_c * iG1413 + T1224 * iG1413_c + T1234_c * iG1412 + T1234 * iG1412_c;
        double R1224 = T1212_c * G1224 + T1212 * G1224_c + T1213_c * G1324
            + T1213 * G1324_c - 2.0 * iT1214_c * iG1314 - 2.0 * iT1214 * iG1314_c
            + T1224_c * G1313 + T1224 * G1313_c + T1234_c * G1312 + T1234 * G1312_c;
        double R1234 = T1212_c * G1234 + T1212 * G1234_c + T1213_c * G1224
            + T1213 * G1224_c - 2.0 * iT1214_c * iG1214 - 2.0 * iT1214 * iG1214_c
            + T1224_c * G1213 + T1224 * G1213_c + T1234_c * G1212 + T1234 * G1212_c;

        std::vector<double> R(5);
        R[0] = R1212;
        R[1] = R1213;
        R[2] = iR1214;
        R[3] = R1224;
        R[4] = R1234;
        return R;
      }

    std::tuple<double, double, double, double, double> compute_R(const double &w,
        const double &c, const double &vp, const double &vs, const double &dn,
        const double &dens, const std::tuple<double, double, double, double, double> &T,
        const int &param)
      {
        // Recursive layer stacking from bottom to top layer (an i denotes an imaginary subdeterminant e.g. iR1214)
        const double T1212 = std::get<0>(T);
        const double T1213 = std::get<1>(T);
        const double iT1214 = std::get<2>(T);
        const double T1224 = std::get<3>(T);
        const double T1234 = std::get<4>(T);

        std::vector<double> G;
        if (param == 0)
          {
            G = compute_G(c, dn, w, vp, vs, dens);
          }
        else if (param == 1)
          {
            G = compute_G_vs(c, dn, w, vp, vs, dens);
          }
        else if (param == 2)
          {
            G = compute_G_vp(c, dn, w, vp, vs, dens);
          }
        else if (param == 3)
          {
            G = compute_G_rho(c, dn, w, vp, vs, dens);
          }

        const double G1212 = G[0];
        const double G1213 = G[1];
        const double iG1214 = G[2];
        const double G1224 = G[3];
        const double G1234 = G[4];
        const double G1312 = G[5];
        const double G1313 = G[6];
        const double iG1314 = G[7];
        const double G1324 = G[8];
        const double iG1412 = G[9];
        const double iG1413 = G[10];
        const double G1414 = G[11];
        const double G2412 = G[12];
        const double G2413 = G[13];
        const double G3412 = G[14];

        double R1212 = T1212 * G1212 + T1213 * G1312 - 2.0 * iT1214 * iG1412
            + T1224 * G2412 + T1234 * G3412;
        double R1213 = T1212 * G1213 + T1213 * G1313 - 2.0 * iT1214 * iG1413
            + T1224 * G2413 + T1234 * G2412;
        double iR1214 = T1212 * iG1214 + T1213 * iG1314 + iT1214 * (2.0 * G1414 - 1.0)
            + T1224 * iG1413 + T1234 * iG1412;
        double R1224 = T1212 * G1224 + T1213 * G1324 - 2.0 * iT1214 * iG1314
            + T1224 * G1313 + T1234 * G1312;
        double R1234 = T1212 * G1234 + T1213 * G1224 - 2.0 * iT1214 * iG1214
            + T1224 * G1213 + T1234 * G1212;

        if (param == 0)
          {
            // Normalize R matrix components to +-100000
            double tmpmax, maxR = R1212;
            if (maxR < 0)
              maxR = (-1.0) * maxR;
            tmpmax = R1213;
            if (tmpmax < 0)
              tmpmax = (-1.0) * tmpmax;
            if (tmpmax > maxR)
              maxR = tmpmax;
            tmpmax = iR1214;
            if (tmpmax < 0)
              tmpmax = (-1.0) * tmpmax;
            if (tmpmax > maxR)
              maxR = tmpmax;
            tmpmax = R1224;
            if (tmpmax < 0)
              tmpmax = (-1.0) * tmpmax;
            if (tmpmax > maxR)
              maxR = tmpmax;
            tmpmax = R1234;
            if (tmpmax < 0)
              tmpmax = (-1.0) * tmpmax;
            if (tmpmax > maxR)
              maxR = tmpmax;
            if (maxR > 1.0e5)
              {
                maxR = 1.0e5 / maxR;
                R1212 = maxR * R1212;
                R1213 = maxR * R1213;
                iR1214 = maxR * iR1214;
                R1224 = maxR * R1224;
                R1234 = maxR * R1234;
              }
          }
        return std::make_tuple(R1212, R1213, iR1214, R1224, R1234);
      }

    double compute_R1212(const double &w, const double &c, const std::vector<double> &vp,
        const std::vector<double> &vs, const double &mu, const std::vector<double> &depth,
        const std::vector<double> &dens, const int &nlay)
      {
        // Recursive layer stacking from bottom to top to get R1212
        std::tuple<double, double, double, double, double> R;
        R = compute_T(w, c, vp[nlay - 1], vs[nlay - 1], mu);
        for (int n = nlay - 2; n >= 0; n--)
          {
            double dn = depth[n + 1] - depth[n];
            R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 0);
          }
        return std::get<0>(R);
      }

    double compute_R1212_vs(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs, const double &mu,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &nlay, const int &gradlay)
      {
        // Recursive layer stacking from bottom to top to get R1212
        std::tuple<double, double, double, double, double> R;
        if (nlay - 1 == gradlay)
          {
            R = compute_T_vs(w, c, vp[nlay - 1], vs[nlay - 1], mu);
          }
        else
          {
            R = compute_T(w, c, vp[nlay - 1], vs[nlay - 1], mu);
          }

        for (int n = nlay - 2; n >= 0; n--)
          {
            double dn = depth[n + 1] - depth[n];
            if (n == gradlay)
              R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 1);
            else
              R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 0);
          }
        return std::get<0>(R);
      }

    double compute_R1212_vp(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs, const double &mu,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &nlay, const int &gradlay)
      {
        // Recursive layer stacking from bottom to top to get R1212
        std::tuple<double, double, double, double, double> R;
        if (nlay - 1 == gradlay)
          {
            R = compute_T_vp(w, c, vp[nlay - 1], vs[nlay - 1], mu);
          }
        else
          {
            R = compute_T(w, c, vp[nlay - 1], vs[nlay - 1], mu);
          }

        for (int n = nlay - 2; n >= 0; n--)
          {
            double dn = depth[n + 1] - depth[n];
            if (n == gradlay)
              R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 2);
            else
              R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 0);
          }
        return std::get<0>(R);
      }

    double compute_R1212_dens(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs, const double &mu,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &nlay, const int &gradlay)
      {
        // Recursive layer stacking from bottom to top to get R1212
        std::tuple<double, double, double, double, double> R;
        if (nlay - 1 == gradlay)
          {
            R = compute_T_rho(w, c, vp[nlay - 1], vs[nlay - 1], mu);
          }
        else
          {
            R = compute_T(w, c, vp[nlay - 1], vs[nlay - 1], mu);
          }

        for (int n = nlay - 2; n >= 0; n--)
          {
            double dn = depth[n + 1] - depth[n];
            if (n == gradlay)
              R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 3);
            else
              R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 0);
          }
        return std::get<0>(R);
      }

    double compute_R1212_c(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs, const double &mu,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &nlay)
      {
        // Recursive layer stacking from bottom to top to get R1212
        std::tuple<double, double, double, double, double> R;
        std::vector<double> R_c(5, 0);
        R = compute_T(w, c, vp[nlay - 1], vs[nlay - 1], mu);
        R_c = compute_T_c(w, c, vp[nlay - 1], vs[nlay - 1], mu);
        for (int n = nlay - 2; n >= 0; n--)
          {
            double dn = depth[n + 1] - depth[n];
            R_c = compute_R_c(w, c, vp[n], vs[n], dn, dens[n], R, R_c);
            R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 0);
          }
        return R_c[0];
      }

    std::vector<std::vector<double>> get_gc_segments(const double &east0,
        const double &north0, const double &east1, const double &north1,
        const double &lon_centr, const double &false_east, const double &length_tolerance)
      {
        // Approximates great circle path with loxodrome segments.

        // Set up of projections
        GeographicLib::TransverseMercatorExact proj(GeographicLib::Constants::WGS84_a(),
            GeographicLib::Constants::WGS84_f(), GeographicLib::Constants::UTM_k0());
        GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(),
            GeographicLib::Constants::WGS84_f());
        GeographicLib::Rhumb rhumb(GeographicLib::Constants::WGS84_a(),
            GeographicLib::Constants::WGS84_f());

        // Geographic coordinates for start & end point
        double lon0, lat0, lon1, lat1;
        proj.Reverse(lon_centr, east0 - false_east, north0, lat0, lon0);
        proj.Reverse(lon_centr, east1 - false_east, north1, lat1, lon1);
        std::vector<double> lon;
        std::vector<double> lat;
        std::vector<double> easting;
        std::vector<double> northing;
        lon.push_back(lon0);
        lon.push_back(lon1);
        lat.push_back(lat0);
        lat.push_back(lat1);
        easting.push_back(east0);
        easting.push_back(east1);
        northing.push_back(north0);
        northing.push_back(north1);

        // Great circle between start & end point
        const GeographicLib::GeodesicLine line = geod.InverseLine(lat[0], lon[0], lat[1],
            lon[1]);

        // Loxodrome distance between start & end point
        double azi_tmp, dist_lox;
        rhumb.Inverse(lat[0], lon[0], lat[1], lon[1], dist_lox, azi_tmp);

        // start with two points, check if distance difference is already smaller than tolerance
        int npts = 2;
        while (dist_lox - line.Distance() > length_tolerance)
          {
            // In each iteration keep only the starting point
            lon.erase(lon.begin() + 1, lon.end());
            lat.erase(lat.begin() + 1, lat.end());
            easting.erase(easting.begin() + 1, easting.end());
            northing.erase(northing.begin() + 1, northing.end());

            // calculate segment length
            const double segment_length = line.Distance() / (npts + 1);
            for (int pt = 1; pt <= npts; pt++)
              {
                // calculate point along great circle
                double lon_tmp, lat_tmp;
                line.Position(pt * segment_length, lat_tmp, lon_tmp);
                lat.push_back(lat_tmp);
                lon.push_back(lon_tmp);

                // transform to utm
                double x, y;
                proj.Forward(lon_centr, lat[pt], lon[pt], x, y);
                easting.push_back(x + false_east);
                northing.push_back(y);
              }
            // add end point to vectors
            lon.push_back(lon1);
            lat.push_back(lat1);
            easting.push_back(east1);
            northing.push_back(north1);

            //caculate loxodrome distance
            dist_lox = 0.0;
            double dist_tmp;
            int easting_size = easting.size();
            for (int segment = 0; segment < easting_size - 1; segment++)
              {
                rhumb.Inverse(lat[segment], lon[segment], lat[segment + 1],
                    lon[segment + 1], dist_tmp, azi_tmp);
                dist_lox = dist_lox + dist_tmp;
              }

            // increase umber of points/segments
            npts = npts * 2;
          }

        // write eastings/northings to pts vector
        std::vector<std::vector<double>> pts(2);
        pts[0] = easting;
        pts[1] = northing;
        return pts;
      }

    std::vector<std::vector<double>> get_t_segments(double &east0, double &north0,
        const double &east1, const double &north1, const double &event_e,
        const double &event_n, const double &lon_centr, const std::vector<double> &origin,
        const double &deast, const double &dnorth, const std::vector<double> &c,
        const int &ncells_east, const std::vector<double> &dcdvs,
        const std::vector<double> &dcdvp, const std::vector<double> &dcdrho,
        const int &nlay, const double &false_east)
      {
        // Computes relative phase delay for a station pair and a given earthquake

        // Set up coordinate transformations
        GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(),
            GeographicLib::Constants::WGS84_f());
        GeographicLib::TransverseMercatorExact proj(GeographicLib::Constants::WGS84_a(),
            GeographicLib::Constants::WGS84_f(), GeographicLib::Constants::UTM_k0());
        double event_lat, event_lon;
        proj.Reverse(lon_centr, event_e - false_east, event_n, event_lat, event_lon);

        // Path segment is interpreted as line segment (not loxodrome)
        // because we don't care about its length (only orientation).
        // Check which model cells are crossed by inter-station path segment.
        const double slope = (north1 - north0) / (east1 - east0);
        const double intercept = north0 - slope * east0;
        double ecell0 = floor((east0 - origin[0]) / deast);
        double ncell0 = floor((north0 - origin[1]) / dnorth);
        const double ecell1 = floor((east1 - origin[0]) / deast);
        const double ncell1 = floor((north1 - origin[1]) / dnorth);
        std::vector<std::vector<double>> times_grads(4);
        std::vector<double> times(1);
        double time = 0.0;
        double mid_e, mid_n, mid_lon, mid_lat, s12, az1, az2, dist_segment_e,
            dist_segment_n;
        std::vector<double> rhograd(dcdrho.size(), 0.0), vsgrad(dcdvs.size(), 0.0),
            vpgrad(dcdvp.size(), 0.0);

        while ((abs(ecell0 - ecell1) > 0.0) || (abs(ncell0 - ncell1) > 0.0))
          {
            double north_intercept, east_intercept, estep, nstep;
            if (east0 < east1)
              {
                east_intercept = origin[0] + (ecell0 + 1) * deast;
                estep = 1.0;
                if (slope < 0)
                  {
                    north_intercept = ((origin[1] + (ncell0 + 1) * dnorth) - intercept)
                        / slope;
                    if (north_intercept != east_intercept)
                      north_intercept = ((origin[1] + ncell0 * dnorth) - intercept)
                          / slope;
                    nstep = -1.0;
                  }
                else if (slope > 0)
                  {
                    north_intercept = ((origin[1] + (ncell0 + 1) * dnorth) - intercept)
                        / slope;
                    nstep = 1.0;
                  }
                else
                  {
                    north_intercept = east_intercept + 9999999.9;
                    nstep = 0.0;
                  }
              }
            else if (east0 > east1)
              {
                east_intercept = origin[0] + ecell0 * deast;
                estep = -1.0;
                if (slope < 0)
                  {
                    north_intercept = ((origin[1] + ncell0 * dnorth) - intercept) / slope;
                    if (north_intercept != east_intercept)
                      north_intercept = ((origin[1] + (ncell0 + 1) * dnorth) - intercept)
                          / slope;
                    nstep = 1.0;
                  }
                else if (slope > 0)
                  {
                    north_intercept = ((origin[1] + ncell0 * dnorth) - intercept) / slope;
                    nstep = -1.0;
                  }
                else
                  {
                    north_intercept = east_intercept + 9999999.9;
                    nstep = 0.0;
                  }
              }
            else
              {
                estep = 0.0;
                north_intercept = east0;
                east_intercept = north_intercept + 9999999.9;
                if (north0 < north1)
                  nstep = 1.0;
                else
                  nstep = -1.0;

              }
            double east_dist = abs(east0 - east_intercept);
            double north_dist = abs(east0 - north_intercept);
            if (north_dist < east_dist)
              {
                dist_segment_e = north_intercept - east0;
                dist_segment_n = (slope * north_intercept + intercept) - north0;
                if (east0 == east1)
                  {
                    if (north0 < north1)
                      {
                        dist_segment_n = (origin[1] + (ncell0 + 1) * dnorth) - north0;
                      }
                    else
                      {
                        dist_segment_n = (origin[1] + ncell0 * dnorth) - north0;
                      }
                  }
                mid_e = east0 + dist_segment_e / 2.0;
                mid_n = north0 + dist_segment_n / 2.0;
                north0 = north0 + dist_segment_n;
                east0 = north_intercept;
                ncell0 = ncell0 + nstep;
              }
            else if (east_dist < north_dist)
              {
                dist_segment_e = east_intercept - east0;
                dist_segment_n = (slope * east_intercept + intercept) - north0;
                mid_e = east0 + dist_segment_e / 2.0;
                mid_n = north0 + dist_segment_n / 2.0;
                east0 = east_intercept;
                north0 = north0 + dist_segment_n;
                ecell0 = ecell0 + estep;
              }
            else
              {
                dist_segment_e = east_intercept - east0;
                dist_segment_n = (slope * east_intercept + intercept) - north0;
                mid_e = east0 + dist_segment_e / 2.0;
                mid_n = north0 + dist_segment_n / 2.0;
                east0 = east_intercept;
                north0 = north0 + dist_segment_n;
                ecell0 = ecell0 + estep;
                ncell0 = ncell0 + nstep;
              }
            proj.Reverse(lon_centr, mid_e - false_east, mid_n, mid_lat, mid_lon);
            geod.Inverse(event_lat, event_lon, mid_lat, mid_lon, s12, az1, az2);
            double a = 450 - az2;
            if (a > 360)
              {
                a = a - 360;
              }
            double ti = (-1.0)
                * (cos(a * M_PI / 180.0) * dist_segment_e
                    + sin(a * M_PI / 180.0) * dist_segment_n)
                / c[ncell0 * ncells_east + ecell0];
            time = time + ti;

            for (int n = 0; n < nlay; n++)
              {
                vsgrad[ncell0 * ncells_east * nlay + ecell0 * nlay + n] = vsgrad[ncell0
                    * ncells_east * nlay + ecell0 * nlay + n]
                    - dcdvs[ncell0 * ncells_east * nlay + ecell0 * nlay + n]
                        * (ti / c[ncell0 * ncells_east + ecell0]);
                vpgrad[ncell0 * ncells_east * nlay + ecell0 * nlay + n] = vpgrad[ncell0
                    * ncells_east * nlay + ecell0 * nlay + n]
                    - dcdvp[ncell0 * ncells_east * nlay + ecell0 * nlay + n]
                        * (ti / c[ncell0 * ncells_east + ecell0]);
                rhograd[ncell0 * ncells_east * nlay + ecell0 * nlay + n] = rhograd[ncell0
                    * ncells_east * nlay + ecell0 * nlay + n]
                    - dcdrho[ncell0 * ncells_east * nlay + ecell0 * nlay + n]
                        * (ti / c[ncell0 * ncells_east + ecell0]);
                ;
              }
          }
        dist_segment_e = east1 - east0;
        dist_segment_n = north1 - north0;
        mid_e = east0 + dist_segment_e / 2.0;
        mid_n = north0 + dist_segment_n / 2.0;
        proj.Reverse(lon_centr, mid_e - false_east, mid_n, mid_lat, mid_lon);
        geod.Inverse(event_lat, event_lon, mid_lat, mid_lon, s12, az1, az2);
        double a = 450 - az2;
        if (a > 360)
          {
            a = a - 360;
          }
        double ti = (-1.0)
            * (cos(a * M_PI / 180.0) * dist_segment_e
                + sin(a * M_PI / 180.0) * dist_segment_n)
            / c[ncell1 * ncells_east + ecell1];
        time = time + ti;

        for (int n = 0; n < nlay; n++)
          {
            vsgrad[ncell1 * ncells_east * nlay + ecell1 * nlay + n] = vsgrad[ncell1
                * ncells_east * nlay + ecell1 * nlay + n]
                - dcdvs[ncell1 * ncells_east * nlay + ecell1 * nlay + n]
                    * (ti / c[ncell1 * ncells_east + ecell1]);
            vpgrad[ncell1 * ncells_east * nlay + ecell1 * nlay + n] = vpgrad[ncell1
                * ncells_east * nlay + ecell1 * nlay + n]
                - dcdvp[ncell1 * ncells_east * nlay + ecell1 * nlay + n]
                    * (ti / c[ncell1 * ncells_east + ecell1]);
            rhograd[ncell1 * ncells_east * nlay + ecell1 * nlay + n] = rhograd[ncell1
                * ncells_east * nlay + ecell1 * nlay + n]
                - dcdrho[ncell1 * ncells_east * nlay + ecell1 * nlay + n]
                    * (ti / c[ncell1 * ncells_east + ecell1]);
          }

        times[0] = time;
        times_grads[0] = times;
        times_grads[1] = vsgrad;
        times_grads[2] = vpgrad;
        times_grads[3] = rhograd;
        return times_grads;
      }
  }
