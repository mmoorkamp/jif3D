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

    struct PropagatorSubdeterminants
      {
      double G1212;
      double G1213;
      double iG1214;
      double G1224;
      double G1234;
      double G1312;
      double G1313;
      double iG1314;
      double G1324;
      double iG1412;
      double iG1413;
      double G1414;
      double G2412;
      double G2413;
      double G3412;
      };

    struct LayerSubdeterminants
      {
      double R1212;
      double R1213;
      double iR1214;
      double R1224;
      double R1234;
      };

    struct SWUtilities
      {
      dcomp nu_a;
      dcomp nu_b;
      double SA;
      double CA;
      double SB;
      double CB;
      double gam;
      dcomp nu_a_nrm;
      dcomp nu_b_nrm;
      double l;
      double ma;
      double mb;
      };

    struct GradientUtilities
      {
      double maa;
      double mbb;
      double nu_a_nrm_a;
      double nu_b_nrm_b;
      double gam_b;
      double l_b;
      double CAA;
      double CBB;
      double SAA;
      double SBB;
      double mac;
      double mbc;
      double nu_a_nrm_c;
      double nu_b_nrm_c;
      double gam_c;
      double l_c;
      double CAC;
      double CBC;
      double SAC;
      double SBC;
      };

    std::vector<double> T2w(const std::vector<double> &periods);
    std::vector<double> array2vector(const ThreeDModelBase::t3DModelData &array,
        const int &NX, const int &NY, const int &NZ);
    double newton_vr(const double &vp, const double &vs, const double &tolerance);
    SWUtilities compute_util(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const bool &botlay);
    GradientUtilities compute_util_grads(const double &w, const double &vs,
        const double &vp, const double &c, const double &thck, const bool &botlay);
    LayerSubdeterminants compute_T(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu);
    LayerSubdeterminants compute_T_vs(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu);
    LayerSubdeterminants compute_T_vp(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu);
    LayerSubdeterminants compute_T_rho(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu);
    LayerSubdeterminants compute_T_c(const double &w, const double &c, const double &vp,
        const double &vs, const double &mu);
    PropagatorSubdeterminants compute_G(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens);
    PropagatorSubdeterminants compute_G_vs(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens);
    PropagatorSubdeterminants compute_G_vp(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens);
    PropagatorSubdeterminants compute_G_rho(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens);
    PropagatorSubdeterminants compute_G_c(const double &c, const double &dn,
        const double &w, const double &vp, const double &vs, const double &dens);
    LayerSubdeterminants compute_R_c(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const double &dens,
        const LayerSubdeterminants &T, const LayerSubdeterminants &T_c);
    LayerSubdeterminants compute_R(const double &w, const double &c, const double &vp,
        const double &vs, const double &dn, const double &dens,
        const LayerSubdeterminants &T, const int &param);
    double compute_R1212(const double &w, const double &c, const std::vector<double> &vp,
        const std::vector<double> &vs, const std::vector<double> &depth,
        const std::vector<double> &dens);
    double compute_R1212_vs(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &gradlay);
    double compute_R1212_vp(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &gradlay);
    double compute_R1212_dens(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs,
        const std::vector<double> &depth, const std::vector<double> &dens,
        const int &gradlay);
    double compute_R1212_c(const double &w, const double &c,
        const std::vector<double> &vp, const std::vector<double> &vs,
        const std::vector<double> &depth, const std::vector<double> &dens);
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
      // Root finding functor for R1212. Keeps all other variables constants and only changes c (phase velocity) to finds roots. Uses the TOMS algorithm from the boost libraries
    public:
      R1212_root(const double &w_, const std::vector<double> &vp_,
          const std::vector<double> &vs_, const std::vector<double> &depth_,
          const std::vector<double> &dens_) :
          w(w_), vp(vp_), vs(vs_), depth(depth_), dens(dens_)
        {
        }
      ;
      double operator()(const double c)
        {
          double R1212 = compute_R1212(w, c, vp, vs, depth, dens);
          return R1212;
        }
    private:
      const double &w;
      const std::vector<double> &vp;
      const std::vector<double> &vs;
      const std::vector<double> &depth;
      const std::vector<double> &dens;
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
