/*
 * SurfaceWaveFunctions.h
 *
 *  Created on: 20 Sep 2019
 *      Author: bweise
 */

#ifndef SURFACEWAVES_SURFACEWAVEFUNCTIONS_H_
#define SURFACEWAVES_SURFACEWAVEFUNCTIONS_H_
#include <vector>
#include <complex>
#include <tuple>
#include <fstream>
#include <iomanip>
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
      PropagatorSubdeterminants() :
          G1212(0.0), G1213(0.0), iG1214(0.0), G1224(0.0), G1234(0.0), G1312(0.0), G1313(
              0.0), iG1314(0.0), G1324(0.0), iG1412(0.0), iG1413(0.0), G1414(0.0), G2412(
              0.0), G2413(0.0), G3412(0.0)
        {
        }
      void PrintPropagatorSubdeterminants(const std::string &filename)
        {
          std::ofstream PropSubDetFile;
          PropSubDetFile.open(filename + ".out");
          PropSubDetFile << std::fixed;
          PropSubDetFile << std::setprecision(15);
          PropSubDetFile << "G1212: \t" << G1212 << "\n";
          PropSubDetFile << "G1213: \t" << G1213 << "\n";
          PropSubDetFile << "iG1214: \t" << iG1214 << "\n";
          PropSubDetFile << "G1224: \t" << G1224 << "\n";
          PropSubDetFile << "G1234: \t" << G1234 << "\n";
          PropSubDetFile << "G1312: \t" << G1312 << "\n";
          PropSubDetFile << "G1313: \t" << G1313 << "\n";
          PropSubDetFile << "iG1314: \t" << iG1314 << "\n";
          PropSubDetFile << "G1324: \t" << G1324 << "\n";
          PropSubDetFile << "iG1412: \t" << iG1412 << "\n";
          PropSubDetFile << "iG1413: \t" << iG1413 << "\n";
          PropSubDetFile << "G1414: \t" << G1414 << "\n";
          PropSubDetFile << "G2412: \t" << G2412 << "\n";
          PropSubDetFile << "G2413: \t" << G2413 << "\n";
          PropSubDetFile << "G3412: \t" << G3412 << "\n";
        }
      };

    struct LayerSubdeterminants
      {
      double R1212;
      double R1213;
      double iR1214;
      double R1224;
      double R1234;
      LayerSubdeterminants() :
          R1212(0.0), R1213(0.0), iR1214(0.0), R1224(0.0), R1234(0.0)
        {
        }
      void PrintLayerSubdeterminants(const std::string &filename)
        {
          std::ofstream PropLayDetFile;
          PropLayDetFile.open(filename + ".out");
          PropLayDetFile << std::fixed;
          PropLayDetFile << std::setprecision(15);
          PropLayDetFile << "R1212: \t" << R1212 << "\n";
          PropLayDetFile << "R1213: \t" << R1213 << "\n";
          PropLayDetFile << "iR1214: \t" << iR1214 << "\n";
          PropLayDetFile << "R1224: \t" << R1224 << "\n";
          PropLayDetFile << "R1234: \t" << R1234 << "\n";
        }
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
      SWUtilities() :
          nu_a(0.0, 0.0), nu_b(0.0, 0.0), SA(0.0), CA(0.0), SB(0.0), CB(0.0), gam(0.0), nu_a_nrm(
              0.0, 0.0), nu_b_nrm(0.0, 0.0), l(0.0), ma(0.0), mb(0.0)
        {
        }
      void PrintSWUtilities(const std::string &filename)
        {
          std::ofstream SWUtilFile;
          SWUtilFile.open(filename + ".out");
          SWUtilFile << std::fixed;
          SWUtilFile << std::setprecision(15);
          SWUtilFile << "nu_a: \t" << nu_a << "\n";
          SWUtilFile << "nu_b: \t" << nu_b << "\n";
          SWUtilFile << "nu_a_nrm: \t" << nu_a_nrm << "\n";
          SWUtilFile << "nu_b_nrm: \t" << nu_b_nrm << "\n";
          SWUtilFile << "Sa: \t" << SA << "\n";
          SWUtilFile << "Ca: \t" << CA << "\n";
          SWUtilFile << "Sb: \t" << SB << "\n";
          SWUtilFile << "Cb: \t" << CB << "\n";
          SWUtilFile << "gam: \t" << gam << "\n";
          SWUtilFile << "l: \t" << l << "\n";
          SWUtilFile << "ma: \t" << ma << "\n";
          SWUtilFile << "mb: \t" << mb << "\n";
        }
      };

    struct GradientUtilities
      {
      double maa;
      double mbb;
      double nu_a_nrm2_a;
      double nu_b_nrm2_b;
      double gam_b;
      double l_b;
      double CAA;
      double CBB;
      double SAA;
      double SBB;
      double mac;
      double mbc;
      double nu_a_nrm2_c;
      double nu_b_nrm2_c;
      double gam_c;
      double l_c;
      double CAC;
      double CBC;
      double SAC;
      double SBC;
      GradientUtilities() :
          maa(0.0), mbb(0.0), nu_a_nrm2_a(0.0), nu_b_nrm2_b(0.0), gam_b(0.0), l_b(0.0), CAA(
              0.0), CBB(0.0), SAA(0.0), SBB(0.0), mac(0.0), mbc(0.0), nu_a_nrm2_c(0.0), nu_b_nrm2_c(
              0.0), gam_c(0.0), l_c(0.0), CAC(0.0), CBC(0.0), SAC(0.0), SBC(0.0)
        {
        }
      void PrintGradientUtilities(const std::string &filename)
        {
          std::ofstream GradUtilFile;
          GradUtilFile.open(filename + ".out");
          GradUtilFile << std::fixed;
          GradUtilFile << std::setprecision(15);
          GradUtilFile << "ma_a: \t" << maa << "\n";
          GradUtilFile << "mb_b: \t" << mbb << "\n";
          GradUtilFile << "ma_c: \t" << mac << "\n";
          GradUtilFile << "mb_c: \t" << mbc << "\n";

          GradUtilFile << "nu_a_nrm_a: \t" << nu_a_nrm2_a << "\n";
          GradUtilFile << "nu_b_nrm_b: \t" << nu_b_nrm2_b << "\n";
          GradUtilFile << "nu_a_nrm2_c: \t" << nu_a_nrm2_c << "\n";
          GradUtilFile << "nu_b_nrm2_c: \t" << nu_b_nrm2_c << "\n";

          GradUtilFile << "gam_b: \t" << gam_b << "\n";
          GradUtilFile << "gam_c: \t" << gam_c << "\n";

          GradUtilFile << "l_b: \t" << l_b << "\n";
          GradUtilFile << "l_c: \t" << l_c << "\n";

          GradUtilFile << "Ca_a: \t" << CAA << "\n";
          GradUtilFile << "Cb_b: \t" << CBB << "\n";
          GradUtilFile << "Sa_a: \t" << SAA << "\n";
          GradUtilFile << "Sb_b: \t" << SBB << "\n";

          GradUtilFile << "Ca_c: \t" << CAC << "\n";
          GradUtilFile << "Cb_c: \t" << CBC << "\n";
          GradUtilFile << "Sa_c: \t" << SAC << "\n";
          GradUtilFile << "Sb_c: \t" << SBC << "\n";
        }
      };

    std::vector<double> T2w(const std::vector<double> &periods);
    std::vector<double> array2vector(const ThreeDModelBase::t3DModelData &array);
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
    std::vector<std::vector<double>> get_t_segments(const double &east0, const double &north0,
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
