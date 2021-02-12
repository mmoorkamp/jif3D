/*
 * SurfaceWaveCalculator.cpp
 *
 *  Created on: 20 Sep 2019
 *      Author: bweise
 */

#include <string>
#include <iostream>
#include <vector>
#include <iterator>
#include <complex>
#include <cmath>
#include <tuple>
#include <utility>
#include <algorithm>
#include <functional>
#include <fstream>
#include <map>
#include <stdlib.h>
#include <netcdf>
#include <omp.h>
#include <boost/math/tools/roots.hpp>
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "../Global/NumUtil.h"
#include "../Global/FatalException.h"
#include "../SurfaceWaves/SurfaceWaveCalculator.h"
#include "../SurfaceWaves/SurfaceWaveFunctions.h"
#include "../ModelBase/VTKTools.h"

namespace jif3D
  {
    SurfaceWaveCalculator::SurfaceWaveCalculator() :
        false_east(500000.0), tolerance(0.001), length_tolerance(1.0), mode_skip_it(2), toms_max_iter(
            50)
      {
      }

    void SurfaceWaveCalculator::WriteGradient(const int &NX, const int &NY, const int &NZ,
        const std::vector<double> &vs_grad,
        const ThreeDModelBase::t3DModelDim &XCellCoords,
        const ThreeDModelBase::t3DModelDim &YCellCoords,
        const ThreeDModelBase::t3DModelDim &ZCellCoords, const std::string &filename,
        const std::string &dataname)
      {
        ThreeDModelBase::t3DModelData dtdvs;
        dtdvs.resize(boost::extents[NX][NY][NZ]);
        for (int i = 0; i < NZ; ++i)
          {
            for (int j = 0; j < NY; ++j)
              {
                for (int k = 0; k < NX; ++k)
                  {
                    const size_t offset = i + NZ * j + NY * NZ * k;
                    dtdvs[k][j][i] = vs_grad[offset];
                  }
              }
          }

        Write3DModelToVTK(filename, dataname, XCellCoords, YCellCoords, ZCellCoords,
            dtdvs);
      }

    SurfaceWaveCalculator::Surf1DResult SurfaceWaveCalculator::CalcSurf1D(
        const std::vector<double> &w, size_t freqindex,
        const std::vector<double> &dens_1D, const std::vector<double> &vs_1D,
        const std::vector<double> &vp_1D, const std::vector<double> &depth)
      {
        bool lvz = 0;
        const size_t NZ = dens_1D.size();
        if (vs_1D.size() != NZ)
          {
            jif3D::FatalException(
                "1D models in surface wave calculation do not have the same size",
                __FILE__,
                __LINE__);
          }
        if (vp_1D.size() != NZ)
          {
            jif3D::FatalException(
                "1D models in surface wave calculation do not have the same size",
                __FILE__,
                __LINE__);
          }
        for (size_t n = 0; n < NZ; n++)
          {
            if (n > 0 && vs_1D[n] < vs_1D[n - 1])
              {
                lvz = 1;
              }
          }

        Surf1DResult Result;

        // Calculation of velocity limits
        const double vsmin = *std::min_element(vs_1D.begin(), vs_1D.end());
        const double vsmax = *std::max_element(vs_1D.begin(), vs_1D.end());
        const double vpmin = *std::min_element(vp_1D.begin(), vp_1D.end());
        const double vpmax = *std::max_element(vp_1D.begin(), vp_1D.end());
        const std::vector<double> c_lim
          { newton_vr(vpmin, vsmin, tolerance) / 1.05, newton_vr(vpmax, vsmax, tolerance)
              * 1.05 };

        // step ratio for root bracketing
        const double stepratio = (vsmin - c_lim[0]) / (2.0 * vsmin);

        // Compute initial R1212 polarization for large period below fundamental mode
        double R1212 = compute_R1212(w.back() / 10.0, c_lim[0], vp_1D, vs_1D, depth,
            dens_1D);
        const bool pol0 = signbit(R1212);

        double c_last = c_lim[0]; //initial value for c to start search

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
            R1212 = compute_R1212(w[freqindex], c1, vp_1D, vs_1D, depth, dens_1D);
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
                    R1212 = compute_R1212(w[freqindex], c2, vp_1D, vs_1D, depth, dens_1D);
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
        R1212_root root(w[freqindex], vp_1D, vs_1D, depth, dens_1D);
        TerminationCondition tc(tolerance);
        brackets = boost::math::tools::toms748_solve(root, c0, c1, tc, max_iter);
        if (lvz == 0 && (brackets.first + brackets.second) / 2.0 < c_last)
          {
            c1 = c_lim[0];
            precision = precision * mode_skip_it;
            goto cnt;
          }
        c_last = (brackets.first + brackets.second) / 2.0;
        Result.c = c_last;
        //R1212 = compute_R1212(w[freq], c_last, vp_1D, vs_1D, depth, dens_1D);

        // Write output to file
        //resultfile << "\n" << easting[estep] << "\t" << northing[nstep] << "\t" << (2.0*M_PI)/w[freq] << "\t" << c_last << "\t" << brackets.second-brackets.first << "\t" << max_iter;

        //double test = compute_R1212(w[freq], c_last, vp_1D, vs_1D, depth, dens_1D);
        const double R_c = compute_R1212_c(w[freqindex], c_last, vp_1D, vs_1D, depth,
            dens_1D);
        Result.dcdvp.resize(NZ);
        Result.dcdvs.resize(NZ);
        Result.dcdrho.resize(NZ);
        Result.rc = R_c;
        //std::cout << 2.0 * M_PI/w[freqindex] << " " << Result.c <<  " " <<  R_c << " ";
        double R_vs, R_vp, R_rho;
        for (int n = NZ - 1; n >= 0; n--)
          {
            //Computation of Gradients
            R_vs = compute_R1212_vs(w[freqindex], c_last, vp_1D, vs_1D, depth, dens_1D,
                n);
            Result.dcdvs[n] = ((-1.0) * R_vs / R_c);
            R_vp = compute_R1212_vp(w[freqindex], c_last, vp_1D, vs_1D, depth, dens_1D,
                n);
            Result.dcdvp[n] = ((-1.0) * R_vp / R_c);
            R_rho = compute_R1212_dens(w[freqindex], c_last, vp_1D, vs_1D, depth, dens_1D,
                n);
            Result.dcdrho[n] = ((-1.0) * R_rho / R_c);
          }
        //std::cout << R_vs << " " << R_vp << " " << R_rho << std::endl;
        return Result;
      }

    void SurfaceWaveCalculator::forward(const SurfaceWaveModel &Model,
        const SurfaceWaveData &Data)
      {
        const std::vector<double> eventlat = Data.GetEventPosLat();
        const std::vector<double> eventlon = Data.GetEventPosLon();
        const std::vector<double> periods = Data.GetPeriods();
        const std::vector<double> mpn = Data.GetMeasPosX();
        const std::vector<double> mpe = Data.GetMeasPosY();
        //const std::vector<double> mpz = Data.GetMeasPosZ();
        //const std::vector<int> NDataPerT = Data.GetDataPerT();
        const std::vector<int> StatPairs = Data.GetStatPairs();
        const std::vector<double> dtp_obs = Data.GetData();
        const std::vector<double> err_obs = Data.GetErrors();
        auto indexmap = Data.GetIndexMap();
        const double lon_centr = Data.GetCentrLon();

        const std::vector<double> depth(Model.GetZCoordinates().begin() + 1,
            Model.GetZCoordinates().end());
        const std::vector<double> easting(Model.GetYCoordinates().begin() + 1,
            Model.GetYCoordinates().end());
        const std::vector<double> northing(Model.GetXCoordinates().begin() + 1,
            Model.GetXCoordinates().end());

        const size_t nperiods = periods.size();
        const size_t npairs = StatPairs.size() / 2;
        const size_t ndata = dtp_obs.size();
        const size_t NX = northing.size();
        const size_t NY = easting.size();
        const size_t NZ = depth.size();
        const size_t nmod = NX * NY * NZ;
        if (dens_grad.size() != nmod)
          {
            dens_grad.resize(nmod);
            vs_grad.resize(nmod);
            vp_grad.resize(nmod);
          }
        std::fill(dens_grad.begin(), dens_grad.end(), 0.0);
        std::fill(vs_grad.begin(), vs_grad.end(), 0.0);
        std::fill(vp_grad.begin(), vp_grad.end(), 0.0);
        if (dtp_mod.size() != ndata)
          {
            dtp_mod.resize(ndata);
          }
        std::fill(dtp_mod.begin(), dtp_mod.end(), 0.0);
        const std::vector<double> vs = array2vector(Model.GetData());
        const std::vector<double> vp = array2vector(Model.GetVp());
        const std::vector<double> dens = array2vector(Model.GetDens());

        const double deast = easting[1] - easting[0];
        const double dnorth = northing[1] - northing[0];
        const std::vector<double> model_origin =
          { easting[0] - deast, northing[0] - dnorth };
        /*boost::posix_time::ptime starttime =
         boost::posix_time::microsec_clock::local_time();
         std::cout << "Starting calculation " << starttime << std::endl;*/
        const std::vector<double> w = T2w(periods);

        //Vectors to store gradients, dispersion curves
        std::vector<double> vph_map(NX * NY * nperiods, 0.0);
        std::vector<double> dcdrho(nmod * nperiods, 0.0), dcdvs(nmod * nperiods, 0.0),
            dcdvp(nmod * nperiods, 0.0);

        omp_lock_t lck;
        omp_init_lock(&lck);

#pragma omp parallel for default(shared)
        for (int freq = 0; freq < nperiods; freq++)
          {
            /*std::cout << "Period: " << periods[freq] << " s.";
             std::cout << "\n";*/

            std::vector<double> dens_1D(NZ);
            std::vector<double> vs_1D(NZ);
            std::vector<double> vp_1D(NZ);

            for (size_t nstep = 0; nstep < NX; nstep++)
              {
                for (size_t estep = 0; estep < NY; estep++)
                  {

                    for (size_t n = 0; n < NZ; n++)
                      {
                        const size_t offset = n + NZ * estep + NY * NZ * nstep;
                        // sort velocities, densities into 1D models
                        dens_1D[n] = dens[offset];
                        vp_1D[n] = vp[offset];
                        vs_1D[n] = vs[offset];
                        // check if there's a low velocity zone
                      }
                    if (vs_1D[0] > 0)
                      {
                        SurfaceWaveCalculator::Surf1DResult Result = CalcSurf1D(w, freq,
                            dens_1D, vs_1D, vp_1D, depth);
                        vph_map[estep + NY * nstep + NY * NX * freq] = Result.c;
                        for (size_t n = 0; n < NZ; n++)
                          {
                            const size_t offset = n + NZ * estep + NY * NZ * nstep
                                + nmod * freq;

                            dcdvs[offset] = Result.dcdvs[n];
                            dcdvp[offset] = Result.dcdvp[n];
                            dcdrho[offset] = Result.dcdrho[n];
                          }
                      }
                  } //end loop over easting
                /*if (nstep % 10 == 0)
                 {

                 std::cout << " Finished " << nstep << " cells, took "
                 << (boost::posix_time::microsec_clock::local_time() - starttime).total_seconds()
                 << " s" << std::endl;
                 }*/
              } //end loop over northing
          } // end frequency loop
        /*boost::posix_time::ptime currtime =
         boost::posix_time::microsec_clock::local_time();
         std::cout << " Finished 1D calculation took "
         << (currtime - starttime).total_seconds() << " s" << std::endl;*/
        // loop over all station pairs, computes phase delays
        int LastPair = -1;
        //std::multimap<int, std::vector<std::vector<double>>> pathmap;

        bool NewGeo = !Data.CompareGeometry(PreviousData);

        if (NewGeo)
          {
            //std::cout << "Creating geometry" << std::endl;
            pathmap_n.resize(npairs);
            pathmap_e.resize(npairs);
            SegLengths.resize(ndata);
            CellIndices.resize(ndata);
//#pragma omp parallel for default(shared), firstprivate(LastPair), schedule(static)
            for (size_t datacounter = 0; datacounter < ndata; datacounter++)
              {
                auto data_it = indexmap.begin();
                std::advance(data_it, datacounter);
                int pairindex = (*data_it).first;
                if (LastPair != pairindex)
                  {
                    std::vector<std::vector<double>> segments = get_gc_segments(
                        mpe[StatPairs[2 * pairindex]], mpn[StatPairs[2 * pairindex]],
                        mpe[StatPairs[2 * pairindex + 1]],
                        mpn[StatPairs[2 * pairindex + 1]], lon_centr, false_east,
                        length_tolerance);
                    std::vector<double> seg_east = segments[0];
                    std::vector<double> seg_north = segments[1];
                    LastPair = pairindex;
                    /*pathmap.insert(
                     std::pair<int, std::vector<std::vector<double>>>(pairindex,
                     segments));*/
                    pathmap_n[pairindex] = seg_north;
                    pathmap_e[pairindex] = seg_east;
                  }
                auto IndicesData = (*data_it).second;
                int eventid = std::get<0>(IndicesData);

                int seg_east_size = pathmap_e[pairindex].size();
                //std::cout << "Data: " << datacounter << std::endl;
                std::vector<double> CurrLengths;
                std::vector<int> CurrIndices;
                for (int seg = 0; seg < seg_east_size - 1; seg++)
                  { //loop over great circle segments
                    //std::cout << "Seg: " << seg << std::endl;

                    std::vector<double> seg_east = pathmap_e[pairindex];
                    std::vector<double> seg_north = pathmap_n[pairindex];

                    t_segments time_segment = get_t_segments(seg_east[seg],
                        seg_north[seg], seg_east[seg + 1], seg_north[seg + 1],
                        eventlat[eventid], eventlon[eventid], lon_centr, model_origin,
                        deast, dnorth, NY, false_east);
                    std::copy(time_segment.cellindices.begin(),
                        time_segment.cellindices.end(), std::back_inserter(CurrIndices));
                    std::copy(time_segment.seglength.begin(),
                        time_segment.seglength.end(), std::back_inserter(CurrLengths));
                  }
                SegLengths.at(datacounter) = CurrLengths;
                CellIndices.at(datacounter) = CurrIndices;
              }
            PreviousData = Data;
            Data.WriteRayPaths(pathmap_e, pathmap_n);
          }

#pragma omp parallel for default(shared)
        for (int datacounter = 0; datacounter < ndata; datacounter++)
          {
            std::vector<double> vs_tmpgrd(NX * NY * NZ, 0.0), vp_tmpgrd(NX * NY * NZ,
                0.0), dens_tmpgrd(NX * NY * NZ, 0.0);

            auto data_it = indexmap.begin();
            std::advance(data_it, datacounter);
            auto IndicesData = (*data_it).second;
            int periodid = std::get<1>(IndicesData);

            double time = 0;
            std::vector<double> rhograd(nmod, 0.0), vsgrad(nmod, 0.0), vpgrad(nmod, 0.0);
            for (size_t i = 0; i < SegLengths.at(datacounter).size(); ++i)
              {
                double ti = SegLengths.at(datacounter)[i]
                    / vph_map[CellIndices.at(datacounter)[i] + periodid * NX * NY];
                time = time - ti;

                const double tic = (ti
                    / vph_map[CellIndices.at(datacounter)[i] + periodid * NX * NY]);
                for (size_t n = 0; n < NZ; n++)
                  {
                    const size_t index = CellIndices.at(datacounter)[i] * NZ + n;
                    const size_t index2 = index + periodid * nmod;
                    vs_tmpgrd[index] += 2.0 * dcdvs[index2] * tic;
                    vp_tmpgrd[index] += 2.0 * dcdvp[index2] * tic;
                    dens_tmpgrd[index] += 2.0 * dcdrho[index2] * tic;
                  }
              }

            dtp_mod[datacounter] = time;
            //delayfile << "\n" << event_stat_cmb[event*nsrcs+src] << "\t" << src_rcvr_cmb[src] << "\t" << src_rcvr_cmb[src+nsrcs] << "\t" << (2.0*M_PI)/w[freq] << "\t" << time_total;

            const double residual = (time - dtp_obs[datacounter])
                / jif3D::pow2(err_obs[datacounter]);

            /*WriteGradient(NX, NY, NZ, vs_tmpgrd, Model.GetXCoordinates(),
             Model.GetYCoordinates(), Model.GetZCoordinates(),
             "dtdvs" + std::to_string(datacounter) + ".vtk", "dtdvs");*/
            omp_set_lock(&lck);
            std::transform(vs_grad.begin(), vs_grad.end(), vs_tmpgrd.begin(),
                vs_grad.begin(), weighted_add(residual));
            std::transform(vp_grad.begin(), vp_grad.end(), vp_tmpgrd.begin(),
                vp_grad.begin(), weighted_add(residual));
            std::transform(dens_grad.begin(), dens_grad.end(), dens_tmpgrd.begin(),
                dens_grad.begin(), weighted_add(residual));
            omp_unset_lock(&lck);

            /*std::vector<double> vs_weighted(nmod);
             std::transform(vs_tmpgrd.begin(), vs_tmpgrd.end(), vs_weighted.begin(),
             [residual](double val)
             { return residual *val;});*/

            /*WriteGradient(NX, NY, NZ, vs_weighted, Model.GetXCoordinates(),
             Model.GetYCoordinates(), Model.GetZCoordinates(),
             "dtdvs_total" + std::to_string(datacounter) + ".vtk", "dtdvs_weighted");*/

            /*if (std::distance(datamap.begin(), data_it) % 10 == 0)
             {

             std::cout << " Finished " << std::distance(datamap.begin(), data_it)
             << " out of " << datamap.size() << " station pairs, took "
             << (boost::posix_time::microsec_clock::local_time() - starttime).total_seconds()
             << " s" << std::endl;
             }*/
          } // end loop over data

      }
  }

