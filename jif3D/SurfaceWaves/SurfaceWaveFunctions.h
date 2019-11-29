/*
 * SurfaceWaveFunctions.h
 *
 *  Created on: 20 Sep 2019
 *      Author: bweise
 */

#ifndef SURFACEWAVES_SURFACEWAVEFUNCTIONS_H_
#define SURFACEWAVES_SURFACEWAVEFUNCTIONS_H_
#include<vector>
#include <complex>
#include <tuple>
#include "../ModelBase/ThreeDModelBase.h"
#include <GeographicLib/TransverseMercatorExact.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Rhumb.hpp>
#include <GeographicLib/Constants.hpp>

namespace jif3D
  {
    typedef std::complex<double> dcomp;
    const std::complex<double> i(0, 1.0);
    std::vector<double> T2w(const std::vector<double> &periods);
    std::vector<double> array2vector(const ThreeDModelBase::t3DModelData &array,
        const int &NX, const int &NY, const int &NZ);
    double newton_vr(const double &vp, const double &vs, const double &tolerance);
    std::tuple<dcomp, dcomp, double, double, double, double, double, dcomp, dcomp, double,
        double, double> compute_util(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const bool &botlay);
    std::tuple<dcomp, dcomp, double, double, double, double, double, double, double,
        double, dcomp, dcomp, double, double, double, double, double, double, double,
        double> compute_util_grads(const double &w, const double &vs, const double &vp,
        const double &c, const double &thck, const bool &botlay);
    std::tuple<double, double, double, double, double> compute_T(const double &w,
        const double &c, const double &vp, const double &vs, const double &mu);
    std::tuple<double, double, double, double, double> compute_T_vs(const double &w,
        const double &c, const double &vp, const double &vs, const double &mu);
    std::tuple<double, double, double, double, double> compute_T_vp(const double &w,
        const double &c, const double &vp, const double &vs, const double &mu);
    std::tuple<double, double, double, double, double> compute_T_rho(const double &w,
        const double &c, const double &vp, const double &vs, const double &mu);
    std::vector<double> compute_T_c(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu);
    std::vector<double> compute_G(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens);
    std::vector<double> compute_G_vs(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens);
    std::vector<double> compute_G_vp(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens);
    std::vector<double> compute_G_rho(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens);
    std::vector<double> compute_G_c(const double &c, const double &dn, const double &w,
        const double &vp, const double &vs, const double &dens);
    std::vector<double> compute_R_c(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const double &dens,
        const std::tuple<double, double, double, double, double> &T,
        const std::vector<double> &T_c);
    std::tuple<double, double, double, double, double> compute_R(const double &w,
        const double &c, const double &vp, const double &vs, const double &dn,
        const double &dens, const std::tuple<double, double, double, double, double> &T,
        const int &param);
    std::vector<double> compute_R1212(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs, const double &mu,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &nlay, const int &param, const int &gradlay);
    std::vector<std::vector<double>> get_gc_segments(const double &east0,
        const double &north0, const double &east1, const double &north1,
        const double &lon_centr, const double &false_east,
        const double &length_tolerance);
    std::vector<std::vector<double>> get_t_segments(double &east0, double &north0,
        const double &east1, const double &north1, const double &event_e,
        const double &event_n, const double &lon_centr, const std::vector<double> &origin,
        const double &deast, const double &dnorth, const std::vector<double> &c,
        const int &ncells_east, const std::vector<double> &dsdvs,
        const std::vector<double> &dsdvp, const std::vector<double> &dsdrho,
        const int &nlay, const double &false_east);

    class R1212_root
      {
      // Root finding functor for R1212. Keeps all other varibles constants and only changes c (phase velocity) to finds roots. Uses the TOMS algorithm from the boost libraries
    public:
      R1212_root(const double &w_, const std::vector<double> &vp_,
          const std::vector<double> &vs_, const double &mu_,
          const std::vector<double> &depth_, const std::vector<double> &dens_,
          const int &nlay_) :
          w(w_), vp(vp_), vs(vs_), mu(mu_), depth(depth_), dens(dens_), nlay(nlay_)
        {
        }
      ;
      double operator()(const double c)
        {
          std::vector<double> R1212 = compute_R1212(w, c, vp, vs, mu, depth, dens, nlay,
              0, -999);
          ;
          return R1212[0];
        }
    private:
      const double &w;
      const std::vector<double> &vp;
      const std::vector<double> &vs;
      const double &mu;
      const std::vector<double> &depth;
      const std::vector<double> &dens;
      const int &nlay;
      };

    struct TerminationCondition
      {
      // checks whether root bracketing has sufficiently converged
      TerminationCondition(const double &tolerance_) :
          tolerance(tolerance_)
        {
        }
      ;
      const double tolerance;
      bool operator()(const double &min, const double &max)
        {
          return abs(min - max) < tolerance;
        }
      };

    struct weighted_add
      {
      const double weight;
      double operator()(const double &aa, const double &bb)
        {
          return aa + bb * weight;
        }
      weighted_add(const double weight_) :
          weight(weight_)
        {
        }
      };
  }

#endif /* SURFACEWAVES_SURFACEWAVEFUNCTIONS_H_ */
