/*
 * SurfaceWaveCalculator.cpp
 *
 *  Created on: 20 Sep 2019
 *      Author: bweise
 */

#include <string>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <functional>
#include <fstream>
#include <stdlib.h>
#include <netcdf>
#include <boost/math/tools/roots.hpp>
#include "../Global/NumUtil.h"
#include "../SurfaceWaves/SurfaceWaveCalculator.h"
#include "../SurfaceWaves/SurfaceWaveFunctions.h"

namespace jif3D
  {
    SurfaceWaveCalculator::SurfaceWaveCalculator() :
        false_east(500000.0), tolerance(0.01), length_tolerance(1.0), mode_skip_it(2), toms_max_iter(
            50)
      {
      }

    void SurfaceWaveCalculator::forward(const SurfaceWaveModel &Model,
        const SurfaceWaveData &Data)
      {
        const std::vector<double> eventx = Data.GetEventPosX();
        const std::vector<double> eventy = Data.GetEventPosY();
        const std::vector<double> periods = Data.GetPeriods();
        const std::vector<double> mpn = Data.GetMeasPosX();
        const std::vector<double> mpe = Data.GetMeasPosY();
        const std::vector<double> mpz = Data.GetMeasPosZ();
        const std::vector<double> src_rcvr_cmb = Data.GetStatComb();
        const std::vector<double> event_stat_cmb = Data.GetEventStatComb();
        const std::vector<double> dtp = Data.GetData();
        const double lon_centr = Data.GetCentrLon();
        const double dtp_dummy = Data.GetDummy();

        std::vector<double> depth(Model.GetZCoordinates().size() - 1);
        std::vector<double> easting(Model.GetYCoordinates().size() - 1);
        std::vector<double> northing(Model.GetXCoordinates().size() - 1);
        std::copy(Model.GetXCoordinates().begin() + 1, Model.GetXCoordinates().end(),
            northing.begin());
        std::copy(Model.GetYCoordinates().begin() + 1, Model.GetYCoordinates().end(),
            easting.begin());
        std::copy(Model.GetZCoordinates().begin() + 1, Model.GetZCoordinates().end(),
            depth.begin());

        ThreeDModelBase::t3DModelData vp_all = Model.GetVp();
        ThreeDModelBase::t3DModelData vs_all = Model.GetData();
        ThreeDModelBase::t3DModelData dens_all = Model.GetDens();

        const int nperiods = periods.size();
        const int nsrcs = src_rcvr_cmb.size() / 2;
        const int nevents_per_src = dtp.size() / (nsrcs * nperiods);
        const int NX = northing.size();
        const int NY = easting.size();
        const int NZ = depth.size();
        const int nmod = NX * NY * NZ;
        if (dens_grad.size() != nmod)
          {
            dens_grad.resize(nmod);
            vs_grad.resize(nmod);
            vp_grad.resize(nmod);
          }
        if (dcdrho.size() != nmod)
          {
            dcdrho.resize(nmod);
            dcdvs.resize(nmod);
            dcdvp.resize(nmod);
          }
        if (dtp_mod.size() != dtp.size())
          {
            dtp_mod.resize(dtp.size());
          }
        std::fill(dens_grad.begin(), dens_grad.end(), 0.0);
        std::fill(vs_grad.begin(), vs_grad.end(), 0.0);
        std::fill(vp_grad.begin(), vp_grad.end(), 0.0);
        std::fill(dcdrho.begin(), dcdrho.end(), 0.0);
        std::fill(dcdvs.begin(), dcdvs.end(), 0.0);
        std::fill(dcdvp.begin(), dcdvp.end(), 0.0);
        std::fill(dtp_mod.begin(), dtp_mod.end(), 0.0);
        const std::vector<double> vs = array2vector(vs_all, NX, NY, NZ);
        const std::vector<double> vp = array2vector(vp_all, NX, NY, NZ);
        const std::vector<double> dens = array2vector(dens_all, NX, NY, NZ);

        const double deast = easting[1] - easting[0];
        const double dnorth = northing[1] - northing[0];
        const std::vector<double> model_origin =
          { easting[0] - deast, northing[0] - dnorth };

        const std::vector<double> w = T2w(periods);

        for (int freq = 0; freq < nperiods; freq++)
          {
            //cout << "Period: " << periods[freq] << " s.";
            //cout << "\n";
            //Vectors to store gradients, dispersion curves
            std::vector<double> vph_map(NX * NY);
            for (int nstep = 0; nstep < NX; nstep++)
              {
                for (int estep = 0; estep < NY; estep++)
                  {
                    std::vector<double> dens_1D(NZ);
                    std::vector<double> vs_1D(NZ);
                    std::vector<double> vp_1D(NZ);
                    bool lvz = 0;
                    for (int n = 0; n < NZ; n++)
                      {
                        // sort velocities, densities into 1D models
                        dens_1D[n] = dens[n + NZ * estep + NY * NZ * nstep];
                        vp_1D[n] = vp[n + NZ * estep + NY * NZ * nstep];
                        vs_1D[n] = vs[n + NZ * estep + NY * NZ * nstep];
                        // check if there's a low velocity zone
                        if (n > 0 && vs_1D[n] < vs_1D[n - 1])
                          {
                            lvz = 1;
                          }
                      }

                    if (vs_1D[0] <= 0)
                      {
                        // This still needs some work. Computation of dispersion curves does not work if top layer has vs=0
                        vph_map[estep + NY * nstep] = 0.0;
                        //resultfile << "\n" << easting[estep] << "\t" << northing[nstep] << "\t" << (2.0*M_PI)/w[freq] << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0;
                        continue;
                      }
                    else
                      {
                        // Calculation of velocity limits
                        const double vsmin = *std::min_element(vs_1D.begin(),
                            vs_1D.end());
                        const double vsmax = *std::max_element(vs_1D.begin(),
                            vs_1D.end());
                        const double vpmin = *std::min_element(vp_1D.begin(),
                            vp_1D.end());
                        const double vpmax = *std::max_element(vp_1D.begin(),
                            vp_1D.end());
                        const std::vector<double> c_lim
                          { newton_vr(vpmin, vsmin, tolerance) / 1.05, newton_vr(vpmax,
                              vsmax, tolerance) * 1.05 };

                        // step ratio for root bracketing
                        const double stepratio = (vsmin - c_lim[0]) / (2.0 * vsmin);

                        // Shear modulus bottom layer
                        const double mu = pow(vs_1D[NZ - 1], 2) * dens_1D[NZ - 1];

                        // Compute initial R1212 polarization for large period below fundamental mode
                        double R1212 = compute_R1212(w[nperiods - 1] / 10.0, c_lim[0],
                            vp_1D, vs_1D, mu, depth, dens_1D, NZ);
                        const bool pol0 = signbit(R1212);

                        double c_last = c_lim[0]; //initial value for c to start search

                        // Loop over all periods/frequencies
                        //for(int freq=0; freq<nperiods; freq++){
                        double c0, c1 = c_last; // stores brackets
                        bool pol1 = pol0; // initial polarization of R1212
                        double precision = 1;

                        if (lvz == 1) // If there is a LVZ, we have to start at the lowest possible velocity
                          c1 = c_lim[0];

                        // Loop to find root brackets, breaks when sign of R1212 changes
                        while (pol0 == pol1)
                          {
                            cnt: ;
                            c0 = c1; // set lower bracket to value of last iteration's upper bracket
                            c1 = c0 + c0 * (stepratio / precision); // increase upper bracket by step ratio

                            // Check polarization of R1212 for the upper bracket
                            R1212 = compute_R1212(w[freq], c1, vp_1D, vs_1D, mu, depth,
                                dens_1D, NZ);
                            pol1 = signbit(R1212);

                            // If a sign change is found check for mode skipping
                            if (pol0 != pol1 && (c1 - c0) > (2.0 * tolerance))
                              {
                                double c2 = (c1 + c0) / 2.0; // set new speed between brackets
                                double delta = (c2 - c0) / mode_skip_it;
                                if (delta < (mode_skip_it * tolerance))
                                  delta = tolerance;
                                // check sign of R1212 between brackets
                                while (tolerance < (c2 - c0))
                                  {
                                    R1212 = compute_R1212(w[freq], c2, vp_1D, vs_1D, mu,
                                        depth, dens_1D, NZ);
                                    const bool pol2 = signbit(R1212);
                                    // if mode skipping detected increase precision (-> decrease step ratio) and return to bracket search
                                    if (pol2 == pol1)
                                      {
                                        precision = precision * mode_skip_it;
                                        c1 = c0;
                                        goto cnt;
                                      }
                                    // "Downward" search along c-axis for mode skipping (2 runs per default)
                                    c2 = c2 - delta;
                                  }
                              }
                          }
                        // If a sign change is found, brackets c0 & c2 are passed to an instance of R1212_root (-> Boost TOMS root finding algorithm)
                        std::pair<double, double> brackets; // stores refinded root brackets
                        boost::uintmax_t max_iter = toms_max_iter; // Maximum number of TOMS iterations
                        R1212_root root(w[freq], vp_1D, vs_1D, mu, depth, dens_1D, NZ);
                        TerminationCondition tc(tolerance);
                        brackets = boost::math::tools::toms748_solve(root, c0, c1, tc,
                            max_iter);
                        if (lvz == 0 && (brackets.first + brackets.second) / 2.0 < c_last)
                          {
                            c1 = c_lim[0];
                            precision = precision * mode_skip_it;
                            goto cnt;
                          }
                        c_last = (brackets.first + brackets.second) / 2.0;

                        vph_map[estep + NY * nstep] = c_last;

                        // Write output to file
                        //resultfile << "\n" << easting[estep] << "\t" << northing[nstep] << "\t" << (2.0*M_PI)/w[freq] << "\t" << c_last << "\t" << brackets.second-brackets.first << "\t" << max_iter;

                        const double R_c = compute_R1212_c(w[freq], c_last, vp_1D, vs_1D,
                            mu, depth, dens_1D, NZ);
                        for (int n = 0; n < NZ; n++)
                          {
                            //Computation of Gradients
                            double R_tmp = compute_R1212_vs(w[freq], c_last, vp_1D, vs_1D,
                                mu, depth, dens_1D, NZ, n);
                            dcdvs[n + NZ * estep + NY * NZ * nstep] = ((-1.0) * R_tmp
                                / R_c);
                            R_tmp = compute_R1212_vp(w[freq], c_last, vp_1D, vs_1D, mu,
                                depth, dens_1D, NZ, n);
                            dcdvp[n + NZ * estep + NY * NZ * nstep] = ((-1.0) * R_tmp
                                / R_c);
                            R_tmp = compute_R1212_dens(w[freq], c_last, vp_1D, vs_1D, mu,
                                depth, dens_1D, NZ, n);
                            dcdrho[n + NZ * estep + NY * NZ * nstep] = ((-1.0) * R_tmp
                                / R_c);
                          }
                      }
                  } //end loop over northing
              } //end loop over easting

            // loop over all rays, computes phase delays
            for (int src = 0; src < nsrcs; src++)
              {
                std::vector<std::vector<double>> segments = get_gc_segments(
                    mpe[src_rcvr_cmb[src]], mpn[src_rcvr_cmb[src]],
                    mpe[src_rcvr_cmb[src + nsrcs]], mpn[src_rcvr_cmb[src + nsrcs]],
                    lon_centr, false_east, length_tolerance);
                std::vector<double> seg_east = segments[0];
                std::vector<double> seg_north = segments[1];
                int seg_east_size = seg_east.size();
                std::vector<std::vector<double>>().swap(segments);
                for (int event = 0; event < nevents_per_src; event++)
                  { //loop over events for each station pair
                    if (dtp[freq * nsrcs * nevents_per_src + event * nsrcs + src]
                        == dtp_dummy || event_stat_cmb[event * nsrcs + src] == dtp_dummy)
                      {
                        // if there is only a dummy value we can skip this period
                        dtp_mod[src + nsrcs * event + freq * nsrcs * nevents_per_src] =
                            dtp_dummy;
                        //delayfile << "\n" << event_stat_cmb[event*nsrcs+src] << "\t" << src_rcvr_cmb[src] << "\t" << src_rcvr_cmb[src+nsrcs] << "\t" << (2.0*M_PI)/w[freq] << "\t" << dtp_dummy;
                        continue;
                      }
                    else
                      {
                        std::vector<double> vs_tmpgrd(NX * NY * NZ, 0.0), vp_tmpgrd(
                            NX * NY * NZ, 0.0), dens_tmpgrd(NX * NY * NZ, 0.0);
                        double time_total = 0.0;
                        for (int seg = 0; seg < seg_east_size - 1; seg++)
                          { //loop over great circle segments
                            //cout << "Period: " << periods[freq] << " s\t Station combination: " << src << " of " << nsrcs << ",\t Event: " << event << " of " << nevents_per_src << ",\t Segment: " << seg << "\n";
                            std::vector<std::vector<double>> time_segment =
                                get_t_segments(seg_east[seg], seg_north[seg],
                                    seg_east[seg + 1], seg_north[seg + 1],
                                    eventy[event_stat_cmb[event * nsrcs + src]],
                                    eventx[event_stat_cmb[event * nsrcs + src]],
                                    lon_centr, model_origin, deast, dnorth, vph_map, NY,
                                    dcdvs, dcdvp, dcdrho, NZ, false_east);
                            std::vector<double> tmp = time_segment[0];
                            time_total = time_total + tmp[0];

                            tmp = time_segment[1];
                            std::transform(vs_tmpgrd.begin(), vs_tmpgrd.end(),
                                tmp.begin(), vs_tmpgrd.begin(), std::plus<double>());
                            tmp = time_segment[2];
                            std::transform(vp_tmpgrd.begin(), vp_tmpgrd.end(),
                                tmp.begin(), vp_tmpgrd.begin(), std::plus<double>());
                            tmp = time_segment[3];
                            std::transform(dens_tmpgrd.begin(), dens_tmpgrd.end(),
                                tmp.begin(), dens_tmpgrd.begin(), std::plus<double>());
                          } //end loop ever path segments
                        dtp_mod[src + nsrcs * event + freq * nsrcs * nevents_per_src] =
                            time_total;
                        //delayfile << "\n" << event_stat_cmb[event*nsrcs+src] << "\t" << src_rcvr_cmb[src] << "\t" << src_rcvr_cmb[src+nsrcs] << "\t" << (2.0*M_PI)/w[freq] << "\t" << time_total;

                        const double residual = (time_total
                            - dtp[freq * nsrcs * nevents_per_src + event * nsrcs + src])
                            / jif3D::pow2(
                                dtp_err[freq * nsrcs * nevents_per_src + event * nsrcs
                                    + src]);
                        std::transform(vs_grad.begin(), vs_grad.end(), vs_tmpgrd.begin(),
                            vs_grad.begin(), weighted_add(residual));
                        std::transform(vp_grad.begin(), vp_grad.end(), vp_tmpgrd.begin(),
                            vp_grad.begin(), weighted_add(residual));
                        std::transform(dens_grad.begin(), dens_grad.end(),
                            dens_tmpgrd.begin(), dens_grad.begin(),
                            weighted_add(residual));
                      }
                  } //end loop over events
              } // end loop rays
          } //end loop over periods
        /*resultfile << "\n";
         resultfile.close();

         delayfile << "\n";
         delayfile.close();*/
      }
  }

