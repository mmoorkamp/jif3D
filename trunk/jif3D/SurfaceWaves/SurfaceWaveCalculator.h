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
      void set_data_err(const std::vector<double> &err)
        {
          dtp_err.resize(err.size());
          dtp_err = err;
        }
      std::vector<double> GetDcdrho() const
        {
          return dcdrho;
        }
      std::vector<double> GetDcdvs() const
        {
          return dcdvs;
        }
      std::vector<double> GetDcdvp() const
        {
          return dcdvp;
        }
    private:
      double false_east, tolerance, length_tolerance;
      int mode_skip_it, toms_max_iter;
      std::vector<double> dcdrho, dcdvs, dcdvp, vph_map;
      std::vector<double> dens_grad, vs_grad, vp_grad, dtp_mod, dtp_err;
      void forward(const ModelType &Model, const DataType &Data);
      void WritePhaseVelocityMaps(const std::string filename, const std::vector<double> &vph_map);
      };
  }

#endif /* SURFACEWAVES_SURFACEWAVECALCULATOR_H_ */
