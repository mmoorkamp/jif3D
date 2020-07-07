/*
 * SurfaceWaveCalculator.h
 *
 *  Created on: 20 Sep 2019
 *      Author: bweise
 */

#ifndef SURFACEWAVES_SURFACEWAVECALCULATOR_H_
#define SURFACEWAVES_SURFACEWAVECALCULATOR_H_

#include "../SurfaceWaves/SurfaceWaveModel.h"
#include "../SurfaceWaves/SurfaceWaveData.h"
#include "../SurfaceWaves/SurfaceWaveFunctions.h"
#include <boost/math/tools/roots.hpp>

namespace jif3D
  {
    class J3DEXPORT SurfaceWaveCalculator
      {
    public:
      typedef SurfaceWaveModel ModelType;
      typedef SurfaceWaveData DataType;
      SurfaceWaveCalculator();
      jif3D::rvec Calculate(const ModelType &Model, const DataType &Data)
        {
          forward(Model, Data);
          std::vector<double> dtp_mod;
          std::multimap<int, std::tuple<int, int, double, double>>::iterator it;
          for (it = datamap_mod.begin(); it != datamap_mod.end(); ++it)
            {
              auto tmptple = (*it).second;
              double tmpdat = std::get<2>(tmptple);
              dtp_mod.push_back(tmpdat);
            }
          jif3D::rvec result(dtp_mod.size());
          std::copy(dtp_mod.begin(), dtp_mod.end(), result.begin());
          return result;
        }
      jif3D::rvec LQDerivative(const ModelType &Model, const DataType &Data,
          const rvec &Misfit) const
        {
          const size_t nparm = vs_grad.size();
          jif3D::rvec result(3 * nparm);
          std::copy(vs_grad.begin(), vs_grad.end(), result.begin());
          std::copy(vp_grad.begin(), vp_grad.end(), result.begin() + nparm);
          std::copy(dens_grad.begin(), dens_grad.end(), result.begin() + 2 * nparm);
          return result;
        }
      void set_false_east(const double &fe)
        {
          false_east = fe;
        }
      void set_vel_tolerance(const double &vtol)
        {
          tolerance = vtol;
        }
      void set_distance_tolerance(const double &dtol)
        {
          length_tolerance = dtol;
        }
      void set_root_search_iterations(const int &root_it)
        {
          toms_max_iter = root_it;
        }
      void set_mode_skip_iterations(const int &mskip_it)
        {
          mode_skip_it = mskip_it;
        }
      void forward(const ModelType &Model, const DataType &Data);
      struct Surf1DResult
        {
        double c;
        double rc;
        std::vector<double> dcdvs;
        std::vector<double> dcdvp;
        std::vector<double> dcdrho;
        };
      Surf1DResult CalcSurf1D(const std::vector<double> &w, size_t freqindex,
          const std::vector<double> &dens_1D, const std::vector<double> &vs_1D,
          const std::vector<double> &vp_1D, const std::vector<double> &depth);
    private:
      double false_east, tolerance, length_tolerance;
      int mode_skip_it, toms_max_iter;
      std::vector<double> dens_grad, vs_grad, vp_grad;
      std::multimap<int, std::tuple<int, int, double, double>> datamap_mod;
      };
  }

#endif /* SURFACEWAVES_SURFACEWAVECALCULATOR_H_ */
