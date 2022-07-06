/*
 * SurfaceWaveFunctions.cpp
 *
 *  Created on: 20 Sep 2019
 *      Author: bweise
 */

#include "../SurfaceWaves/SurfaceWaveFunctions.h"
#include "../ModelBase/ThreeDModelBase.h"
#include "../Global/NumUtil.h"

#include <vector>
#include <complex>
#include <tuple>
#include <string>
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Rhumb.hpp>
#include <GeographicLib/Constants.hpp>
#include <boost/math/special_functions/sign.hpp>

namespace jif3D
  {
    std::vector<double> T2w(const std::vector<double> &periods)
      {
        // Conversion of periods to angular frequencies
        std::vector<double> w(periods.size());
        int nperiods = periods.size();
        for (int n = 0; n < nperiods; n++)
          {
            w[n] = 2.0 * M_PI / periods[n];
            //w[n] = 1.0 / periods[n];
          }
        return w;
      }

    std::vector<double> array2vector(const ThreeDModelBase::t3DModelData &array)
      {
        const size_t NX = array.shape()[0];
        const size_t NY = array.shape()[1];
        const size_t NZ = array.shape()[2];
        std::vector<double> tmp(NX * NY * NZ);
        for (size_t i = 0; i < NZ; ++i)
          {
            for (size_t j = 0; j < NY; ++j)
              {
                for (size_t k = 0; k < NX; ++k)
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

    SWUtilities compute_util(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const bool &botlay)
      {
        // computes some constants for each layer (vertical wave numbers and some derived properties)
        SWUtilities u;

        const double k = w / c;
        if (c < vp)
          {
            u.ma = sqrt(jif3D::pow2(w / c) - jif3D::pow2(w / vp));
            u.nu_a = u.ma;
            if (botlay == 0)
              {
                u.SA = (k / u.ma) * sinh(u.ma * dn);
                u.CA = cosh(u.ma * dn);
              }
          }
        else
          {
            u.ma = sqrt(jif3D::pow2(w / vp) - jif3D::pow2(w / c));
            u.nu_a = u.ma;
            if (botlay == 0)
              {
                u.nu_a = u.nu_a * i;
                u.SA = (k / u.ma) * sin(u.ma * dn);
                u.CA = cos(u.ma * dn);
              }
          }
        if (c < vs)
          {
            u.mb = sqrt(jif3D::pow2(w / c) - jif3D::pow2(w / vs));
            u.nu_b = u.mb;
            if (botlay == 0)
              {
                u.SB = (k / u.mb) * sinh(u.mb * dn);
                u.CB = cosh(u.mb * dn);
              }
          }
        else
          {
            u.mb = sqrt(jif3D::pow2(w / vs) - jif3D::pow2(w / c));
            u.nu_b = u.mb;
            if (botlay == 0)
              {
                u.nu_b = u.nu_b * i;
                u.SB = (k / u.mb) * sin(u.mb * dn);
                u.CB = cos(u.mb * dn);
              }
          }
        u.gam = 2.0 * jif3D::pow2(vs) / jif3D::pow2(c);
        u.nu_a_nrm = u.nu_a / k;
        u.nu_b_nrm = u.nu_b / k;
        u.l = 2.0 * jif3D::pow2(k) - jif3D::pow2(w / vs);

        /*if (botlay == 1)
         {
         u.PrintSWUtilities("BottomLayerUtilities_" + std::to_string(w));
         }
         else
         {
         u.PrintSWUtilities(
         "Utilities_" + std::to_string(w) + "_" + std::to_string(vs) + "_"
         + std::to_string(vp));
         }*/

        return u;
      }

    GradientUtilities compute_util_grads(const double &w, const double &vs,
        const double &vp, const double &c, const double &thck, const bool &botlay)
      {
        const SWUtilities u = compute_util(w, c, vp, vs, thck, botlay);

        GradientUtilities gu;

        const double k = w / c;

        if (c < vp)
          {
            gu.maa = pow(w, 2) / (u.ma * pow(vp, 3));
            gu.mac = (-1.0) * (pow(w, 2) / (u.ma * pow(c, 3)));
          }
        else
          {
            gu.maa = (-1.0 * pow(w, 2)) / (u.ma * pow(vp, 3));
            gu.mac = pow(w, 2) / (u.ma * pow(c, 3));
          }

        if (c < vs)
          {
            gu.mbb = pow(w, 2) / (u.mb * pow(vs, 3));
            gu.mbc = (-1.0) * (pow(w, 2) / (u.mb * pow(c, 3)));
          }
        else
          {
            gu.mbb = (-1.0 * pow(w, 2)) / (u.mb * pow(vs, 3));
            gu.mbc = pow(w, 2) / (u.mb * pow(c, 3));
          }

        gu.nu_a_nrm2_c = (-2.0 * c) / pow(vp, 2);
        gu.nu_a_nrm2_a = (2.0 * pow(c, 2)) / pow(vp, 3);
        gu.nu_b_nrm2_c = (-2.0 * c) / pow(vs, 2);
        gu.nu_b_nrm2_b = (2.0 * pow(c, 2)) / pow(vs, 3);
        gu.gam_b = (4.0 * vs) / pow(c, 2);
        gu.gam_c = (-4.0 * pow(vs, 2)) / pow(c, 3);
        gu.l_b = 2.0 * pow(w, 2) / pow(vs, 3);
        gu.l_c = -4.0 * pow(w, 2) / pow(c, 3);
        gu.CAA = (w * c * thck * u.SA) / pow(vp, 3);
        gu.CBB = (w * c * thck * u.SB) / pow(vs, 3);
        gu.CAC = ((-1.0) * k * thck * u.SA) / c;
        gu.CBC = ((-1.0) * k * thck * u.SB) / c;
        gu.SAC = (-1.0) * ((gu.mac / u.ma) + (1 / c)) * u.SA
            + (k * thck * gu.mac * u.CA / u.ma);
        gu.SBC = (-1.0) * ((gu.mbc / u.mb) + (1 / c)) * u.SB
            + (k * thck * gu.mbc * u.CB / u.mb);
        gu.SAA = (gu.maa / u.ma) * (k * thck * u.CA - u.SA);
        gu.SBB = (gu.mbb / u.mb) * (k * thck * u.CB - u.SB);

        /*if (botlay == 1)
         {
         gu.PrintGradientUtilities(
         "BottomLayerGradientUtilities_" + std::to_string(w));
         }
         else
         {
         gu.PrintGradientUtilities(
         "GradientUtilities_" + std::to_string(w) + "_" + std::to_string(vs) + "_"
         + std::to_string(vp));
         }*/

        return gu;
      }

    LayerSubdeterminants compute_T(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu)
      {
        // computes layer matrix for bottom layer (an i denotes an imaginary subdeterminant e.g. iT1214)
        const double k = w / c;
        const SWUtilities u = compute_util(w, c, vp, vs, 99999.0, 1);

        const dcomp fact = pow(
            (-1.0) * pow(vs, 2) / (2 * mu * u.nu_a * u.nu_b * pow(w, 2)), 2);

        LayerSubdeterminants T;

        T.R1212 = std::real(
            pow(mu, 2) * u.nu_b * u.nu_a
                * (pow(u.l, 2) - 4.0 * pow(k, 2) * u.nu_b * u.nu_a) * fact);
        T.R1213 = std::real(
            mu * pow(u.nu_a, 2) * u.nu_b * (u.l - 2.0 * pow(k, 2)) * fact);
        T.iR1214 = std::real(
            k * mu * u.nu_a * u.nu_b * (u.l - 2.0 * u.nu_a * u.nu_b) * fact);
        T.R1224 = std::real(
            mu * u.nu_a * pow(u.nu_b, 2) * (2.0 * pow(k, 2) - u.l) * fact);
        T.R1234 = std::real(u.nu_a * u.nu_b * (pow(k, 2) - u.nu_a * u.nu_b) * fact);

        //T.PrintLayerSubdeterminants("BottomLayerSubdeterminants_" + std::to_string(w));

        return T;
      }

    LayerSubdeterminants compute_T_vs(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu)
      {
        const LayerSubdeterminants T = compute_T(w, c, vp, vs, mu);

        LayerSubdeterminants Tvs;

        const SWUtilities u = compute_util(w, c, vp, vs, 99999.0, 1);

        const GradientUtilities ug = compute_util_grads(w, vs, vp, c, 99999.0, 1);

        const double dens = mu / pow(vs, 2);

        Tvs.R1212 = std::real(
            (4.0 * T.R1212 / vs)
                + pow(vs, 4) * (2.0 * u.l * ug.l_b * u.nu_b - pow(u.l, 2) * ug.mbb)
                    / (4.0 * pow(w, 4) * u.nu_a * pow(u.nu_b, 2)));
        Tvs.R1213 = std::real(ug.mbb / (4.0 * dens * pow(w, 2) * pow(u.nu_b, 2)));
        Tvs.iR1214 = std::real(
            (2.0 * T.iR1214 / vs)
                + pow(vs, 2) * (ug.l_b * u.nu_b - u.l * ug.mbb)
                    / (4.0 * dens * pow(w, 3) * c * u.nu_a * pow(u.nu_b, 2)));
        Tvs.R1224 = 0.0;
        Tvs.R1234 = std::real(
            (-1.0) * ug.mbb
                / (4.0 * pow(dens, 2) * pow(w, 2) * pow(c, 2) * u.nu_a * pow(u.nu_b, 2)));

        /*Tvs.PrintLayerSubdeterminants(
         "BottomLayerSubdeterminants_vs_" + std::to_string(w));*/

        return Tvs;
      }

    LayerSubdeterminants compute_T_vp(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu)
      {
        const SWUtilities u = compute_util(w, c, vp, vs, 99999.0, 1);

        const GradientUtilities ug = compute_util_grads(w, vs, vp, c, 99999.0, 1);

        LayerSubdeterminants Tvp;

        Tvp.R1212 = std::real(
            (-1.0) * pow(vs, 4) * pow(u.l, 2) * ug.maa
                / (4.0 * pow(w, 4) * pow(u.nu_a, 2) * u.nu_b));
        Tvp.R1213 = 0.0;
        Tvp.iR1214 = std::real(
            (-1.0) * pow(vs, 4) * u.l * ug.maa
                / (4.0 * mu * pow(w, 3) * c * pow(u.nu_a, 2) * u.nu_b));
        Tvp.R1224 = std::real(
            (-1.0) * ug.maa * pow(vs, 2) / (4.0 * mu * pow(w, 2) * pow(u.nu_a, 2)));
        Tvp.R1234 = std::real(
            (-1.0) * ug.maa * pow(vs, 4)
                / (4.0 * pow(mu, 2) * pow(w, 2) * pow(c, 2) * pow(u.nu_a, 2) * u.nu_b));

        /*Tvp.PrintLayerSubdeterminants(
         "BottomLayerSubdeterminants_vp_" + std::to_string(w));*/

        return Tvp;
      }

    LayerSubdeterminants compute_T_rho(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu)
      {
        const LayerSubdeterminants T = compute_T(w, c, vp, vs, mu);

        LayerSubdeterminants Tdens;

        const double dens = mu / pow(vs, 2);

        Tdens.R1212 = 0.0;
        Tdens.R1213 = (-1.0) * T.R1213 / dens;
        Tdens.iR1214 = (-1.0) * T.iR1214 / dens;
        Tdens.R1224 = (-1.0) * T.R1224 / dens;
        Tdens.R1234 = (-2.0) * T.R1234 / dens;

        /*Tdens.PrintLayerSubdeterminants(
         "BottomLayerSubdeterminants_dens_" + std::to_string(w));*/

        return Tdens;
      }

    LayerSubdeterminants compute_T_c(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu)
      {
        const SWUtilities u = compute_util(w, c, vp, vs, 99999.0, 1);

        const GradientUtilities ug = compute_util_grads(w, vs, vp, c, 99999.0, 1);

        const LayerSubdeterminants T = compute_T(w, c, vp, vs, mu);

        LayerSubdeterminants Tc;

        const double dens = mu / pow(vs, 2);

        Tc.R1212 = std::real(
            (pow(vs, 4) / (4.0 * pow(w, 4)))
                * (((2.0 * u.l * ug.l_c * u.nu_a * u.nu_b
                    - pow(u.l, 2) * (ug.mac * u.nu_b + u.nu_a * ug.mbc))
                    / (pow(u.nu_a * u.nu_b, 2))) + (8.0 * pow(w, 2) / pow(c, 3))));
        Tc.R1213 = std::real(ug.mbc / (4.0 * dens * pow(w * u.nu_b, 2)));
        Tc.iR1214 = std::real(
            ((-1.0) * T.iR1214 / c)
                + ((pow(vs, 2)
                    * (ug.l_c * u.nu_b * u.nu_a
                        - u.l * (ug.mac * u.nu_b + u.nu_a * ug.mbc)))
                    / (4.0 * dens * pow(w, 3) * c * pow(u.nu_b * u.nu_a, 2))));
        Tc.R1224 = std::real(((-1.0) * ug.mac) / (4.0 * dens * pow(u.nu_a * w, 2)));
        Tc.R1234 = std::real(
            (-1.0) * (2.0 * u.nu_a * u.nu_b + c * (ug.mac * u.nu_b + u.nu_a * ug.mbc))
                / (4.0 * pow(dens * w * u.nu_a * u.nu_b, 2) * pow(c, 3)));

        //Tc.PrintLayerSubdeterminants("BottomLayerSubdeterminants_c_" + std::to_string(w));

        return Tc;
      }

    PropagatorSubdeterminants compute_G(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens)
      {
        // computes subdeterminants of G matrix (an i denotes an imaginary subdeterminant e.g. iG1214)
        const SWUtilities u = compute_util(w, c, vp, vs, dn, 0);

        PropagatorSubdeterminants G;

        G.G1212 = std::real(
            2.0 * u.gam * (1.0 - u.gam)
                + (2.0 * pow(u.gam, 2) - 2.0 * u.gam + 1.0) * u.CA * u.CB
                - (pow(1.0 - u.gam, 2)
                    + pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * u.SA
                    * u.SB);
        G.G1213 = std::real(
            (1.0 / (dens * w * c)) * (u.CA * u.SB - u.SA * u.CB * pow(u.nu_a_nrm, 2)));
        G.iG1214 = std::real(
            (1.0 / (dens * w * c))
                * ((1.0 - 2.0 * u.gam) * (1.0 - u.CB * u.CA)
                    + (1.0 - u.gam - u.gam * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2))
                        * u.SA * u.SB));
        G.G1224 = std::real(
            (1.0 / (dens * w * c)) * (pow(u.nu_b_nrm, 2) * u.CA * u.SB - u.SA * u.CB));
        G.G1234 = std::real(
            (-1.0 / (pow(dens, 2) * pow(w, 2) * pow(c, 2)))
                * (2.0 * (1.0 - u.CA * u.CB)
                    + (1.0 + pow(u.nu_b_nrm, 2) * pow(u.nu_a_nrm, 2)) * u.SB * u.SA));
        G.G1312 = std::real(
            dens * w * c
                * (pow(u.gam, 2) * pow(u.nu_b_nrm, 2) * u.CA * u.SB
                    - pow(1.0 - u.gam, 2) * u.SA * u.CB));
        G.G1313 = std::real(u.CA * u.CB);
        G.iG1314 = std::real(
            (1.0 - u.gam) * u.SA * u.CB + u.gam * pow(u.nu_b_nrm, 2) * u.CA * u.SB);
        G.G1324 = std::real((-1.0) * pow(u.nu_b_nrm, 2) * u.SA * u.SB);
        G.iG1412 = std::real(
            dens * w * c
                * ((3.0 * pow(u.gam, 2) - 2.0 * pow(u.gam, 3) - u.gam)
                    * (1.0 - u.CA * u.CB)
                    + (pow(1.0 - u.gam, 3)
                        - pow(u.gam, 3) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * u.SA
                        * u.SB));
        G.iG1413 =
            std::real(
                (-1.0)
                    * ((1.0 - u.gam) * u.CA * u.SB
                        + u.gam * pow(u.nu_a_nrm, 2) * u.SA * u.CB));
        G.G1414 = std::real(
            1.0 - 2.0 * u.gam * (1.0 - u.gam) * (1.0 - u.CA * u.CB)
                + (pow(1.0 - u.gam, 2)
                    + pow(u.gam, 2) * pow(u.nu_b_nrm, 2) * pow(u.nu_a_nrm, 2)) * u.SA
                    * u.SB);
        G.G2412 = std::real(
            dens * w * c
                * (pow(1.0 - u.gam, 2) * u.CA * u.SB
                    - pow(u.gam, 2) * u.SA * u.CB * pow(u.nu_a_nrm, 2)));
        G.G2413 = std::real((-1.0) * pow(u.nu_a_nrm, 2) * u.SA * u.SB);
        G.G3412 = std::real(
            (-1.0) * pow(dens, 2) * pow(w, 2) * pow(c, 2)
                * (2.0 * pow(u.gam, 2) * pow(1.0 - u.gam, 2) * (1.0 - u.CA * u.CB)
                    + (pow(1.0 - u.gam, 4)
                        + pow(u.gam, 4) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * u.SA
                        * u.SB));

        /*G.PrintPropagatorSubdeterminants(
         "PropagatorSubdeterminants_" + std::to_string(w) + "_" + std::to_string(vs)
         + "_" + std::to_string(vp) + "_" + std::to_string(dens));*/

        return G;
      }

    PropagatorSubdeterminants compute_G_vs(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens)
      {
        const SWUtilities u = compute_util(w, c, vp, vs, dn, 0);

        const GradientUtilities ug = compute_util_grads(w, vs, vp, c, dn, 0);

        PropagatorSubdeterminants G;

        G.G1212 = std::real(
            2.0 * ug.gam_b * (1.0 - 2.0 * u.gam) * (1.0 - u.CA * u.CB)
                + (2.0 * pow(u.gam, 2) - 2.0 * u.gam + 1.0) * u.CA * ug.CBB
                - (2.0 * ug.gam_b * (u.gam - 1.0)
                    + 2.0 * u.gam * ug.gam_b * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)
                    + pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * ug.nu_b_nrm2_b) * u.SA * u.SB
                - (pow((1.0 - u.gam), 2)
                    + pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * u.SA
                    * ug.SBB);
        G.G1213 = std::real(
            (1.0 / (dens * w * c))
                * (u.CA * ug.SBB - pow(u.nu_a_nrm, 2) * u.SA * ug.CBB));
        G.iG1214 = std::real(
            (1.0 / (dens * w * c))
                * ((-2.0) * ug.gam_b * (1.0 - u.CA * u.CB)
                    - (1.0 - 2.0 * u.gam) * u.CA * ug.CBB)
                + (1.0 / (dens * w * c))
                    * (((-1.0) * ug.gam_b
                        - ug.gam_b * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)
                        - u.gam * pow(u.nu_a_nrm, 2) * ug.nu_b_nrm2_b) * u.SA * u.SB
                        + (1.0 - u.gam - u.gam * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2))
                            * u.SA * ug.SBB));
        G.G1224 = std::real(
            (1.0 / (dens * w * c))
                * (ug.nu_b_nrm2_b * u.CA * u.SB + pow(u.nu_b_nrm, 2) * u.CA * ug.SBB
                    - u.SA * ug.CBB));
        G.G1234 = std::real(
            (-1.0 / pow((dens * w * c), 2))
                * (-2.0 * u.CA * ug.CBB
                    + (1.0 + pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * u.SA * ug.SBB
                    + pow(u.nu_a_nrm, 2) * ug.nu_b_nrm2_b * u.SA * u.SB));
        G.G1312 = std::real(
            dens * w * c
                * (2.0 * u.gam * ug.gam_b * pow(u.nu_b_nrm, 2) * u.CA * u.SB
                    + pow(u.gam, 2) * ug.nu_b_nrm2_b * u.CA * u.SB
                    + pow(u.gam, 2) * pow(u.nu_b_nrm, 2) * u.CA * ug.SBB)
                + dens * w * c
                    * (2.0 * ug.gam_b * (1.0 - u.gam) * u.SA * u.CB
                        - pow((1.0 - u.gam), 2) * u.SA * ug.CBB));
        G.G1313 = std::real(u.CA * ug.CBB);
        G.iG1314 = std::real(
            (-1.0) * ug.gam_b * u.SA * u.CB + (1.0 - u.gam) * u.SA * ug.CBB
                + ug.gam_b * pow(u.nu_b_nrm, 2) * u.CA * u.SB
                + u.gam * ug.nu_b_nrm2_b * u.CA * u.SB
                + u.gam * pow(u.nu_b_nrm, 2) * u.CA * ug.SBB);
        G.G1324 = std::real(
            (-1.0) * ug.nu_b_nrm2_b * u.SA * u.SB - pow(u.nu_b_nrm, 2) * u.SA * ug.SBB);
        G.iG1412 = std::real(
            dens * w * c
                * (ug.gam_b * (-6.0 * pow(u.gam, 2) + 6.0 * u.gam - 1.0)
                    * (1.0 - u.CA * u.CB)
                    - (pow(u.gam, 2) - u.gam) * (1.0 - 2.0 * u.gam) * u.CA * ug.CBB)
                + dens * w * c
                    * (-3.0 * pow((1.0 - u.gam), 2) * ug.gam_b
                        - 3.0 * pow(u.gam, 2) * ug.gam_b * pow(u.nu_a_nrm, 2)
                            * pow(u.nu_b_nrm, 2)
                        - pow(u.gam, 3) * pow(u.nu_a_nrm, 2) * ug.nu_b_nrm2_b) * u.SA
                    * u.SB
                + dens * w * c
                    * (pow((1.0 - u.gam), 3)
                        - pow(u.gam, 3) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * u.SA
                    * ug.SBB);
        G.iG1413 = std::real(
            (-1.0)
                * ((1.0 - u.gam) * u.CA * ug.SBB - ug.gam_b * u.CA * u.SB
                    + ug.gam_b * pow(u.nu_a_nrm, 2) * u.SA * u.CB
                    + u.gam * pow(u.nu_a_nrm, 2) * u.SA * ug.CBB));
        G.G1414 = std::real(
            2.0 * ug.gam_b * (2.0 * u.gam - 1.0) * (1.0 - u.CA * u.CB)
                - 2.0 * (pow(u.gam, 2) - u.gam) * u.CA * ug.CBB
                + (2.0 * (u.gam - 1.0) * ug.gam_b
                    + 2.0 * u.gam * ug.gam_b * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)
                    + pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * ug.nu_b_nrm2_b) * u.SA * u.SB
                + (pow(1.0 - u.gam, 2)
                    + pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * u.SA
                    * ug.SBB);
        G.G2412 = std::real(
            dens * w * c
                * (2.0 * (u.gam - 1.0) * ug.gam_b * u.CA * u.SB
                    + pow((1.0 - u.gam), 2) * u.CA * ug.SBB
                    - 2.0 * u.gam * ug.gam_b * pow(u.nu_a_nrm, 2) * u.SA * u.CB
                    - pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * u.SA * ug.CBB));
        G.G2413 = std::real((-1.0) * pow(u.nu_a_nrm, 2) * u.SA * ug.SBB);
        G.G3412 = std::real(
            (-1.0) * pow((dens * c * w), 2)
                * (4.0 * u.gam * ug.gam_b * (1.0 + 2.0 * pow(u.gam, 2) - 3.0 * u.gam)
                    * (1.0 - u.CA * u.CB)
                    - 2.0 * pow(u.gam, 2) * pow((1.0 - u.gam), 2) * u.CA * ug.CBB)
                - pow((dens * c * w), 2)
                    * ((pow((1.0 - u.gam), 4)
                        + pow(u.gam, 4) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * u.SA
                        * ug.SBB)
                - pow((dens * c * w), 2)
                    * ((-4.0 * ug.gam_b * pow((1.0 - u.gam), 3)
                        + 4.0 * pow(u.gam, 3) * ug.gam_b * pow(u.nu_a_nrm, 2)
                            * pow(u.nu_b_nrm, 2)
                        + pow(u.gam, 4) * pow(u.nu_a_nrm, 2) * ug.nu_b_nrm2_b) * u.SA
                        * u.SB));

        /*G.PrintPropagatorSubdeterminants(
         "PropagatorSubdeterminants_vs_" + std::to_string(w) + "_" + std::to_string(vs)
         + "_" + std::to_string(vp) + "_" + std::to_string(dens));*/

        return G;
      }

    PropagatorSubdeterminants compute_G_vp(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens)
      {
        const SWUtilities u = compute_util(w, c, vp, vs, dn, 0);

        const GradientUtilities ug = compute_util_grads(w, vs, vp, c, dn, 0);

        PropagatorSubdeterminants G;

        G.G1212 = std::real(
            (2.0 * pow(u.gam, 2) - 2.0 * u.gam + 1.0) * ug.CAA * u.CB
                - pow(u.gam, 2) * ug.nu_a_nrm2_a * pow(u.nu_b_nrm, 2) * u.SA * u.SB
                - (pow((1.0 - u.gam), 2)
                    + pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * ug.SAA
                    * u.SB);
        G.G1213 = std::real(
            (1.0 / (dens * w * c))
                * (ug.CAA * u.SB - ug.nu_a_nrm2_a * u.SA * u.CB
                    - pow(u.nu_a_nrm, 2) * ug.SAA * u.CB));
        G.iG1214 = std::real(
            (1.0 / (dens * w * c))
                * ((2.0 * u.gam - 1.0) * ug.CAA * u.CB
                    + (1.0 - u.gam - u.gam * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2))
                        * ug.SAA * u.SB
                    - u.gam * ug.nu_a_nrm2_a * pow(u.nu_b_nrm, 2) * u.SA * u.SB));
        G.G1224 = std::real(
            (1.0 / (dens * w * c))
                * (pow(u.nu_b_nrm, 2) * ug.CAA * u.SB - ug.SAA * u.CB));
        G.G1234 = std::real(
            ((-1.0) / pow(w * dens * c, 2))
                * (-2.0 * ug.CAA * u.CB
                    + (1.0 + pow(u.nu_b_nrm, 2) * pow(u.nu_a_nrm, 2)) * ug.SAA * u.SB
                    + ug.nu_a_nrm2_a * pow(u.nu_b_nrm, 2) * u.SA * u.SB));
        G.G1312 = std::real(
            dens * w * c
                * (pow(u.gam, 2) * pow(u.nu_b_nrm, 2) * ug.CAA * u.SB
                    - pow((1.0 - u.gam), 2) * ug.SAA * u.CB));
        G.G1313 = std::real(ug.CAA * u.CB);
        G.iG1314 = std::real(
            (1.0 - u.gam) * ug.SAA * u.CB + u.gam * pow(u.nu_b_nrm, 2) * ug.CAA * u.SB);
        G.G1324 = std::real((-1.0) * pow(u.nu_b_nrm, 2) * ug.SAA * u.SB);
        G.iG1412 = std::real(
            dens * w * c
                * (-1.0 * u.gam * (u.gam - 1.0) * (1.0 - 2.0 * u.gam) * ug.CAA * u.CB
                    - pow(u.gam, 3) * ug.nu_a_nrm2_a * pow(u.nu_b_nrm, 2) * u.SA * u.SB)
                + dens * w * c
                    * ((pow((1.0 - u.gam), 3)
                        - pow(u.gam, 3) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2))
                        * ug.SAA * u.SB));
        G.iG1413 = std::real(
            (u.gam - 1.0) * ug.CAA * u.SB - u.gam * ug.nu_a_nrm2_a * u.SA * u.CB
                - u.gam * pow(u.nu_a_nrm, 2) * ug.SAA * u.CB);
        G.G1414 = std::real(
            2.0 * u.gam * (1.0 - u.gam) * ug.CAA * u.CB
                + (pow((1.0 - u.gam), 2)
                    + pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2)) * ug.SAA
                    * u.SB
                + pow(u.gam, 2) * ug.nu_a_nrm2_a * pow(u.nu_b_nrm, 2) * u.SA * u.SB);
        G.G2412 = std::real(
            dens * w * c
                * (pow((1.0 - u.gam), 2) * ug.CAA * u.SB
                    - pow(u.gam, 2) * ug.nu_a_nrm2_a * u.SA * u.CB
                    - pow(u.gam, 2) * pow(u.nu_a_nrm, 2) * ug.SAA * u.CB));
        G.G2413 = std::real(
            (-1.0) * (ug.nu_a_nrm2_a * u.SA + pow(u.nu_a_nrm, 2) * ug.SAA) * u.SB);
        G.G3412 = std::real(
            (-1.0) * pow((dens * c * w), 2)
                * (-2.0 * pow(u.gam, 2) * pow((1 - u.gam), 2) * ug.CAA * u.CB)
                - pow((dens * c * w), 2)
                    * ((pow((1 - u.gam), 4)
                        + pow(u.gam, 4) * pow(u.nu_a_nrm, 2) * pow(u.nu_b_nrm, 2))
                        * ug.SAA * u.SB
                        + pow(u.gam, 4) * ug.nu_a_nrm2_a * pow(u.nu_b_nrm, 2) * u.SA
                            * u.SB));

        /*G.PrintPropagatorSubdeterminants(
         "PropagatorSubdeterminants_vp_" + std::to_string(w) + "_" + std::to_string(vs)
         + "_" + std::to_string(vp) + "_" + std::to_string(dens));*/

        return G;
      }

    PropagatorSubdeterminants compute_G_rho(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens)
      {
        const PropagatorSubdeterminants G = compute_G(c, dn, w, vp, vs, dens);

        PropagatorSubdeterminants Gout;

        Gout.G1212 = 0.0;
        Gout.G1213 = (-1.0) * G.G1213 / dens;
        Gout.iG1214 = (-1.0) * G.iG1214 / dens;
        Gout.G1224 = (-1.0) * G.G1224 / dens;
        Gout.G1234 = (-2.0) * G.G1234 / dens;
        Gout.G1312 = G.G1312 / dens; // Check this out!!
        Gout.G1313 = 0.0;
        Gout.iG1314 = 0.0;
        Gout.G1324 = 0.0;
        Gout.iG1412 = G.iG1412 / dens;
        Gout.iG1413 = 0.0;
        Gout.G1414 = 0.0;
        Gout.G2412 = G.G2412 / dens;
        Gout.G2413 = 0.0;
        Gout.G3412 = 2.0 * G.G3412 / dens;

        /*Gout.PrintPropagatorSubdeterminants(
         "PropagatorSubdeterminants_dens_" + std::to_string(w) + "_"
         + std::to_string(vs) + "_" + std::to_string(vp) + "_"
         + std::to_string(dens));*/

        return Gout;
      }

    PropagatorSubdeterminants compute_G_c(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens)
      {
        const SWUtilities u = compute_util(w, c, vp, vs, dn, 0);

        const GradientUtilities ug = compute_util_grads(w, vs, vp, c, dn, 0);

        const double dCC = ug.CAC * u.CB + u.CA * ug.CBC;
        const double dSS = ug.SAC * u.SB + u.SA * ug.SBC;
        const double dCS = ug.CAC * u.SB + u.CA * ug.SBC;
        const double dSC = ug.SAC * u.CB + u.SA * ug.CBC;
        const double dk = std::real(
            ug.nu_a_nrm2_c * pow(u.nu_b_nrm, 2) + pow(u.nu_a_nrm, 2) * ug.nu_b_nrm2_c);

        const PropagatorSubdeterminants G = compute_G(c, dn, w, vp, vs, dens);

        PropagatorSubdeterminants Gout;

        Gout.G1212 = std::real(
            2.0 * ug.gam_c * (1.0 - 2.0 * u.gam) * (1.0 - u.CA * u.CB)
                + (2.0 * pow(u.gam, 2) - 2.0 * u.gam + 1.0) * dCC
                - (2.0 * ug.gam_c
                    * (u.gam - 1.0 + u.gam * pow(u.nu_a_nrm * u.nu_b_nrm, 2))
                    + pow(u.gam, 2) * dk) * u.SA * u.SB
                - (pow(1.0 - u.gam, 2) + pow(u.gam * u.nu_a_nrm * u.nu_b_nrm, 2)) * dSS);
        Gout.G1213 = std::real(
            ((-1.0) * G.G1213 / c)
                + (1.0 / (dens * w * c))
                    * (dCS - ug.nu_a_nrm2_c * u.SA * u.CB - pow(u.nu_a_nrm, 2) * dSC));
        Gout.iG1214 = std::real(
            ((-1.0) * G.iG1214 / c)
                + (1.0 / (dens * w * c))
                    * ((2.0 * u.gam - 1.0) * dCC - 2.0 * ug.gam_c * (1.0 - u.CA * u.CB))
                + (1.0 / (dens * w * c))
                    * ((1.0 - u.gam - u.gam * pow(u.nu_a_nrm * u.nu_b_nrm, 2)) * dSS
                        - (ug.gam_c + ug.gam_c * pow(u.nu_a_nrm * u.nu_b_nrm, 2)
                            + u.gam * dk) * u.SA * u.SB));
        Gout.G1224 = std::real(
            ((-1.0) * G.G1224 / c)
                + (1.0 / (dens * w * c))
                    * (ug.nu_b_nrm2_c * u.CA * u.SB + pow(u.nu_b_nrm, 2) * dCS - dSC));
        Gout.G1234 = std::real(
            ((-2.0) * G.G1234 / c)
                + pow((1.0 / (dens * w * c)), 2)
                    * (2.0 * dCC - dk * u.SA * u.SB
                        - (1.0 + pow(u.nu_a_nrm * u.nu_b_nrm, 2)) * dSS));
        Gout.G1312 = std::real(
            (G.G1312 / c)
                + dens * w * c
                    * (2.0 * u.gam * ug.gam_c * pow(u.nu_b_nrm, 2) * u.CA * u.SB
                        + pow(u.gam, 2) * ug.nu_b_nrm2_c * u.CA * u.SB
                        + pow(u.gam * u.nu_b_nrm, 2) * dCS)
                + dens * w * c
                    * (2.0 * ug.gam_c * (1.0 - u.gam) * u.SA * u.CB
                        - pow(1.0 - u.gam, 2) * dSC));
        Gout.G1313 = std::real(dCC);
        Gout.iG1314 = std::real(
            (-1.0) * ug.gam_c * u.SA * u.CB + (1.0 - u.gam) * dSC
                + (ug.gam_c * pow(u.nu_b_nrm, 2) + u.gam * ug.nu_b_nrm2_c) * u.CA * u.SB
                + u.gam * pow(u.nu_b_nrm, 2) * dCS);
        Gout.G1324 = std::real(
            (-1.0) * (ug.nu_b_nrm2_c * u.SA * u.SB + pow(u.nu_b_nrm, 2) * dSS));
        Gout.iG1412 = std::real(
            (G.iG1412 / c)
                + dens * w * c
                    * (ug.gam_c * (-6.0 * pow(u.gam, 2) + 6.0 * u.gam - 1.0)
                        * (1.0 - u.CA * u.CB)
                        - (pow(u.gam, 2) - u.gam) * (1.0 - 2.0 * u.gam) * dCC)
                + dens * c * w
                    * ((pow(1.0 - u.gam, 3)
                        - pow(u.gam, 3) * pow(u.nu_a_nrm * u.nu_b_nrm, 2)) * dSS
                        + (-3.0 * ug.gam_c * pow(1.0 - u.gam, 2)
                            - 3.0 * pow(u.gam, 2) * ug.gam_c
                                * pow(u.nu_a_nrm * u.nu_b_nrm, 2) - pow(u.gam, 3) * dk)
                            * u.SA * u.SB));
        Gout.iG1413 = std::real(
            ug.gam_c * u.CA * u.SB - (1.0 - u.gam) * dCS
                - (ug.gam_c * pow(u.nu_a_nrm, 2) + u.gam * ug.nu_a_nrm2_c) * u.SA * u.CB
                - u.gam * pow(u.nu_a_nrm, 2) * dSC);
        Gout.G1414 = std::real(
            2.0 * u.gam * (1.0 - u.gam) * dCC
                + 2.0 * ug.gam_c * (2.0 * u.gam - 1.0) * (1.0 - u.CA * u.CB)
                + (pow(1.0 - u.gam, 2) + pow(u.gam * u.nu_a_nrm * u.nu_b_nrm, 2)) * dSS
                + (2.0 * ug.gam_c
                    * (u.gam - 1.0 + u.gam * pow(u.nu_a_nrm * u.nu_b_nrm, 2))
                    + pow(u.gam, 2) * dk) * u.SA * u.SB);
        Gout.G2412 = std::real(
            (G.G2412 / c)
                + dens * w * c
                    * (-2.0 * (1.0 - u.gam) * ug.gam_c * u.CA * u.SB
                        + pow(1.0 - u.gam, 2) * dCS)
                + dens * w * c
                    * (-1.0
                        * (2.0 * u.gam * ug.gam_c * pow(u.nu_a_nrm, 2)
                            + pow(u.gam, 2) * ug.nu_a_nrm2_c) * u.SA * u.CB
                        - pow(u.gam * u.nu_a_nrm, 2) * dSC));
        Gout.G2413 = std::real(
            (-1.0) * (ug.nu_a_nrm2_c * u.SA * u.SB + pow(u.nu_a_nrm, 2) * dSS));
        Gout.G3412 = std::real(
            (2.0 * G.G3412 / c)
                - pow(dens * c * w, 2)
                    * (4.0 * u.gam * ug.gam_c * (u.gam - 1.0) * (2.0 * u.gam - 1.0)
                        * (1.0 - u.CA * u.CB)
                        - 2.0 * pow(u.gam, 2) * pow(1.0 - u.gam, 2) * dCC)
                - pow(dens * c * w, 2)
                    * ((pow(1.0 - u.gam, 4)
                        + pow(u.gam, 4) * pow(u.nu_a_nrm * u.nu_b_nrm, 2)) * dSS
                        + (-4.0 * ug.gam_c * pow(1.0 - u.gam, 3)
                            + 4.0 * pow(u.gam, 3) * ug.gam_c
                                * pow(u.nu_a_nrm * u.nu_b_nrm, 2) + pow(u.gam, 4) * dk)
                            * u.SA * u.SB));

        /*Gout.PrintPropagatorSubdeterminants(
         "PropagatorSubdeterminants_c_" + std::to_string(w) + "_" + std::to_string(vs)
         + "_" + std::to_string(vp) + "_" + std::to_string(dens));*/

        return Gout;
      }

    LayerSubdeterminants compute_R_c(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const double &dens,
        const LayerSubdeterminants &T, const LayerSubdeterminants &Tc)
      {
        PropagatorSubdeterminants G = compute_G(c, dn, w, vp, vs, dens);
        PropagatorSubdeterminants G_c = compute_G_c(c, dn, w, vp, vs, dens);

        LayerSubdeterminants Rc;

        Rc.R1212 = Tc.R1212 * G.G1212 + T.R1212 * G_c.G1212 + Tc.R1213 * G.G1312
            + T.R1213 * G_c.G1312 - 2.0 * Tc.iR1214 * G.iG1412
            - 2.0 * T.iR1214 * G_c.iG1412 + Tc.R1224 * G.G2412 + T.R1224 * G_c.G2412
            + Tc.R1234 * G.G3412 + T.R1234 * G_c.G3412;
        //      std::cout << Tc.R1212 * G.G1212 << " " << T.R1212 * G_c.G1212 << " " << Tc.R1213 * G.G1312 << " " <<
        //          T.R1213 * G_c.G1312 << " " << - 2.0 * Tc.iR1214 * G.iG1412 << " " << - 2.0 * T.iR1214 * G_c.iG1412
        //          << " " << Tc.R1224 * G.G2412 << " " << T.R1224 * G_c.G2412 << " " << Tc.R1234 * G.G3412 << " " << T.R1234 * G_c.G3412 << std::endl;
        Rc.R1213 = Tc.R1212 * G.G1213 + T.R1212 * G_c.G1213 + Tc.R1213 * G.G1313
            + T.R1213 * G_c.G1313 - 2.0 * Tc.iR1214 * G.iG1413
            - 2.0 * T.iR1214 * G_c.iG1413 + Tc.R1224 * G.G2413 + T.R1224 * G_c.G2413
            + Tc.R1234 * G.G2412 + T.R1234 * G_c.G2412;
        Rc.iR1214 = Tc.R1212 * G.iG1214 + T.R1212 * G_c.iG1214 + Tc.R1213 * G.iG1314
            + T.R1213 * G_c.iG1314 + Tc.iR1214 * (2.0 * G.G1414 - 1.0)
            + 2.0 * T.iR1214 * G_c.G1414 + Tc.R1224 * G.iG1413 + T.R1224 * G_c.iG1413
            + Tc.R1234 * G.iG1412 + T.R1234 * G_c.iG1412;
        Rc.R1224 = Tc.R1212 * G.G1224 + T.R1212 * G_c.G1224 + Tc.R1213 * G.G1324
            + T.R1213 * G_c.G1324 - 2.0 * Tc.iR1214 * G.iG1314
            - 2.0 * T.iR1214 * G_c.iG1314 + Tc.R1224 * G.G1313 + T.R1224 * G_c.G1313
            + Tc.R1234 * G.G1312 + T.R1234 * G_c.G1312;
        Rc.R1234 = Tc.R1212 * G.G1234 + T.R1212 * G_c.G1234 + Tc.R1213 * G.G1224
            + T.R1213 * G_c.G1224 - 2.0 * Tc.iR1214 * G.iG1214
            - 2.0 * T.iR1214 * G_c.iG1214 + Tc.R1224 * G.G1213 + T.R1224 * G_c.G1213
            + Tc.R1234 * G.G1212 + T.R1234 * G_c.G1212;

        /*Rc.PrintLayerSubdeterminants(
         "LayerSubdeterminants_c_" + std::to_string(w) + "_" + std::to_string(vs) + "_"
         + std::to_string(vp) + "_" + std::to_string(dens));*/

        return Rc;
      }

    LayerSubdeterminants compute_R(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const double &dens,
        const LayerSubdeterminants &T, const int &param)
      {
        // Recursive layer stacking from bottom to top layer (an i denotes an imaginary subdeterminant e.g. iR1214)

        PropagatorSubdeterminants G;
        if (param == 0 || param == 4)
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

        LayerSubdeterminants R;

        if (param == 0 || param == 4)
          {
            R.R1212 = T.R1212 * G.G1212 + T.R1213 * G.G1312 - 2.0 * T.iR1214 * G.iG1412
                + T.R1224 * G.G2412 + T.R1234 * G.G3412;
            R.R1213 = T.R1212 * G.G1213 + T.R1213 * G.G1313 - 2.0 * T.iR1214 * G.iG1413
                + T.R1224 * G.G2413 + T.R1234 * G.G2412;
            R.iR1214 = T.R1212 * G.iG1214 + T.R1213 * G.iG1314
                + T.iR1214 * (2.0 * G.G1414 - 1.0) + T.R1224 * G.iG1413
                + T.R1234 * G.iG1412;
            R.R1224 = T.R1212 * G.G1224 + T.R1213 * G.G1324 - 2.0 * T.iR1214 * G.iG1314
                + T.R1224 * G.G1313 + T.R1234 * G.G1312;
            R.R1234 = T.R1212 * G.G1234 + T.R1213 * G.G1224 - 2.0 * T.iR1214 * G.iG1214
                + T.R1224 * G.G1213 + T.R1234 * G.G1212;
          }
        else
          {
            R.R1212 = T.R1212 * G.G1212 + T.R1213 * G.G1312 - 2.0 * T.iR1214 * G.iG1412
                + T.R1224 * G.G2412 + T.R1234 * G.G3412;
            R.R1213 = T.R1212 * G.G1213 + T.R1213 * G.G1313 - 2.0 * T.iR1214 * G.iG1413
                + T.R1224 * G.G2413 + T.R1234 * G.G2412;
            R.iR1214 = T.R1212 * G.iG1214 + T.R1213 * G.iG1314 + 2.0 * T.iR1214 * G.G1414
                + T.R1224 * G.iG1413 + T.R1234 * G.iG1412;
            R.R1224 = T.R1212 * G.G1224 + T.R1213 * G.G1324 - 2.0 * T.iR1214 * G.iG1314
                + T.R1224 * G.G1313 + T.R1234 * G.G1312;
            R.R1234 = T.R1212 * G.G1234 + T.R1213 * G.G1224 - 2.0 * T.iR1214 * G.iG1214
                + T.R1224 * G.G1213 + T.R1234 * G.G1212;
          }

        if (param == 0)
          {
            // Normalize R matrix components to +-100000
            double tmpmax, maxR = R.R1212;
            if (maxR < 0)
              maxR = (-1.0) * maxR;
            tmpmax = R.R1213;
            if (tmpmax < 0)
              tmpmax = (-1.0) * tmpmax;
            if (tmpmax > maxR)
              maxR = tmpmax;
            tmpmax = R.iR1214;
            if (tmpmax < 0)
              tmpmax = (-1.0) * tmpmax;
            if (tmpmax > maxR)
              maxR = tmpmax;
            tmpmax = R.R1224;
            if (tmpmax < 0)
              tmpmax = (-1.0) * tmpmax;
            if (tmpmax > maxR)
              maxR = tmpmax;
            tmpmax = R.R1234;
            if (tmpmax < 0)
              tmpmax = (-1.0) * tmpmax;
            if (tmpmax > maxR)
              maxR = tmpmax;
            if (maxR > 1.0e5)
              {
                maxR = 1.0e5 / maxR;
                R.R1212 = maxR * R.R1212;
                R.R1213 = maxR * R.R1213;
                R.iR1214 = maxR * R.iR1214;
                R.R1224 = maxR * R.R1224;
                R.R1234 = maxR * R.R1234;
              }
          }

        /*if (param == 0)
         {
         R.PrintLayerSubdeterminants(
         "LayerSubdeterminants_normalized_" + std::to_string(w) + "_"
         + std::to_string(vs) + "_" + std::to_string(vp) + "_"
         + std::to_string(dens));
         }
         else if (param == 1)
         {
         R.PrintLayerSubdeterminants(
         "LayerSubdeterminants_vs_" + std::to_string(w) + "_" + std::to_string(vs)
         + "_" + std::to_string(vp) + "_" + std::to_string(dens));
         }
         else if (param == 2)
         {
         R.PrintLayerSubdeterminants(
         "LayerSubdeterminants_vp_" + std::to_string(w) + "_" + std::to_string(vs)
         + "_" + std::to_string(vp) + "_" + std::to_string(dens));
         }
         else if (param == 3)
         {
         R.PrintLayerSubdeterminants(
         "LayerSubdeterminants_dens_" + std::to_string(w) + "_"
         + std::to_string(vs) + "_" + std::to_string(vp) + "_"
         + std::to_string(dens));
         }
         else if (param == 4)
         {
         R.PrintLayerSubdeterminants(
         "LayerSubdeterminants_" + std::to_string(w) + "_" + std::to_string(vs)
         + "_" + std::to_string(vp) + "_" + std::to_string(dens));
         }*/

        return R;
      }

    double compute_R1212(const double &w, const double &c, const std::vector<double> &vp,
        const std::vector<double> &vs, const std::vector<double> &depth,
        const std::vector<double> &dens)
      {
        const int nlay = vs.size();
        const double mu = dens[nlay - 1] * pow(vs[nlay - 1], 2);
        // Recursive layer stacking from bottom to top to get R1212
        LayerSubdeterminants R;
        R = compute_T(w, c, vp[nlay - 1], vs[nlay - 1], mu);
        for (int n = nlay - 2; n >= 0; n--)
          {
            double dn = 0.0;
            if (n > 0)
              {
                dn = depth[n] - depth[n - 1];
              }
            else
              {
                dn = depth[n];
              }
            R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 0);
          }
        return R.R1212;
      }

    double compute_R1212_vs(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &gradlay)
      {
        const int nlay = vs.size();
        const double mu = dens[nlay - 1] * pow(vs[nlay - 1], 2);
        // Recursive layer stacking from bottom to top to get R1212
        LayerSubdeterminants R;
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
            double dn = 0.0;
            if (n > 0)
              {
                dn = depth[n] - depth[n - 1];
              }
            else
              {
                dn = depth[n];
              }
            if (n == gradlay)
              {
                R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 1);
              }
            else
              {
                R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 4);
              }
          }
        return R.R1212;
      }

    double compute_R1212_vp(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &gradlay)
      {
        const int nlay = vs.size();
        const double mu = dens[nlay - 1] * pow(vs[nlay - 1], 2);
        // Recursive layer stacking from bottom to top to get R1212
        LayerSubdeterminants R;
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
            double dn = 0.0;
            if (n > 0)
              {
                dn = depth[n] - depth[n - 1];
              }
            else
              {
                dn = depth[n];
              }
            if (n == gradlay)
              {
                R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 2);
              }
            else
              {
                R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 4);
              }
          }
        return R.R1212;
      }

    double compute_R1212_dens(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &gradlay)
      {
        const int nlay = vs.size();
        const double mu = dens[nlay - 1] * pow(vs[nlay - 1], 2);
        // Recursive layer stacking from bottom to top to get R1212
        LayerSubdeterminants R;
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
            double dn = 0.0;
            if (n > 0)
              {
                dn = depth[n] - depth[n - 1];
              }
            else
              {
                dn = depth[n];
              }
            if (n == gradlay)
              {
                R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 3);
              }
            else
              {
                R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 4);
              }
          }
        return R.R1212;
      }

    double compute_R1212_c(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs,
        const std::vector<double> &depth, const std::vector<double> &dens)
      {
        const int nlay = vs.size();
        const double mu = dens[nlay - 1] * pow(vs[nlay - 1], 2);
        // Recursive layer stacking from bottom to top to get R1212
        LayerSubdeterminants R;
        LayerSubdeterminants R_c;
        R = compute_T(w, c, vp[nlay - 1], vs[nlay - 1], mu);
        R_c = compute_T_c(w, c, vp[nlay - 1], vs[nlay - 1], mu);
        for (int n = nlay - 2; n >= 0; n--)
          {
            double dn = 0.0;
            if (n > 0)
              {
                dn = depth[n] - depth[n - 1];
              }
            else
              {
                dn = depth[n];
              }
            R_c = compute_R_c(w, c, vp[n], vs[n], dn, dens[n], R, R_c);
            R = compute_R(w, c, vp[n], vs[n], dn, dens[n], R, 4);
          }
        return R_c.R1212;
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

    t_segments get_t_segments(const double &east0_in, const double &north0_in,
        const double &east1, const double &north1, const double &event_lat,
        const double &event_lon, const double &lon_centr,
        const std::vector<double> &origin, const double &deast, const double &dnorth,
        const int &ncells_east, const double &false_east)
      {
        // Computes relative phase delay for a station pair and a given earthquake

        // Set up coordinate transformations
        GeographicLib::Geodesic geod(GeographicLib::Constants::WGS84_a(),
            GeographicLib::Constants::WGS84_f());
        GeographicLib::TransverseMercatorExact proj(GeographicLib::Constants::WGS84_a(),
            GeographicLib::Constants::WGS84_f(), GeographicLib::Constants::UTM_k0());
        //double event_lat, event_lon;
        //proj.Reverse(lon_centr, event_e - false_east, event_n, event_lat, event_lon);

        // Path segment is interpreted as line segment (not loxodrome)
        // because we don't care about its length (only orientation).
        // Check which model cells are crossed by inter-station path segment.
        double north0 = north0_in, east0 = east0_in;
        const double slope = (north1 - north0) / (east1 - east0);
        const double intercept = north0 - slope * east0;
        double ecell0 = floor((east0 - origin[0]) / deast);
        double ncell0 = floor((north0 - origin[1]) / dnorth);
        const double ecell1 = floor((east1 - origin[0]) / deast);
        const double ncell1 = floor((north1 - origin[1]) / dnorth);

        double mid_e, mid_n, mid_lon, mid_lat, s12, az1, az2, dist_segment_e,
            dist_segment_n;

        t_segments result;
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
            double ncell = ncell0;
            double ecell = ecell0;
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

            //  std::cout << dist_segment_n << " " << dist_segment_e << " "
            //      << std::sqrt(
            //          dist_segment_n * dist_segment_n + dist_segment_e * dist_segment_e)
            //      << std::endl;
            double ti = (cos(a * M_PI / 180.0) * dist_segment_e
                + sin(a * M_PI / 180.0) * dist_segment_n);

            //std::cout << "Mid lat: " << mid_lat << " Mid lon: " << mid_lon << " Az1: "
            //                << az1 << " Az2: " << az2 << " a: " << a << " ti: " << ti << std::endl;
            //            std::cout << "Ev lat: " << event_lat << " Ev lon: " << event_lon << " N Cell: " << ncell  << " ECell: " << ecell<< std::endl;

            result.seglength.push_back(ti);
            result.cellindices.push_back(ncell * ncells_east + ecell);
            /*const double tic = ti / c[ncell * ncells_east + ecell];
             for (int n = 0; n < nlay; n++)
             {
             const size_t index = ncell * ncells_east * nlay + ecell * nlay + n;
             vsgrad[index] += 2.0 * dcdvs[index] * tic;
             vpgrad[index] += 2.0 * dcdvp[index] * tic;
             rhograd[index] += 2.0 * dcdrho[index] * tic;
             }*/
          }
        dist_segment_e = east1 - east0;
        dist_segment_n = north1 - north0;
        mid_e = east0 + dist_segment_e / 2.0;
        mid_n = north0 + dist_segment_n / 2.0;
        proj.Reverse(lon_centr, mid_e - false_east, mid_n, mid_lat, mid_lon);
        geod.Inverse(event_lat, event_lon, mid_lat, mid_lon, s12, az1, az2);
        //std::cout << "Mid lat: " << mid_lat << " Mid lon: " << mid_lon << " Az1: " << az1
        //    << " Az2: " << az2 << std::endl;
        double a = 450 - az2;
        if (a > 360)
          {
            a = a - 360;
          }
        // std::cout << dist_segment_n << " " << dist_segment_e << " "
//            << std::sqrt(
        //       dist_segment_n * dist_segment_n + dist_segment_e * dist_segment_e)
        //  << std::endl;

        double ti = (cos(a * M_PI / 180.0) * dist_segment_e
            + sin(a * M_PI / 180.0) * dist_segment_n);
        result.seglength.push_back(ti);
        result.cellindices.push_back(ncell1 * ncells_east + ecell1);
        /*double ti = (cos(a * M_PI / 180.0) * dist_segment_e
         + sin(a * M_PI / 180.0) * dist_segment_n) / c[ncell1 * ncells_east + ecell1];
         time = time - ti;

         const double tic = (ti / c[ncell1 * ncells_east + ecell1]);
         for (int n = 0; n < nlay; n++)
         {
         const size_t index = ncell1 * ncells_east * nlay + ecell1 * nlay + n;
         vsgrad[index] += 2.0 * dcdvs[index] * tic;
         vpgrad[index] += 2.0 * dcdvp[index] * tic;
         rhograd[index] += 2.0 * dcdrho[index] * tic;
         }*/

        return result;
      }
  }
