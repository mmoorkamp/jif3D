//============================================================================
// Name        : MT2DForward.cpp
// Author      : Sep 24, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "MT2DForward.h"
#include "Tarits2DMT.h"
namespace jiba
  {

    MT2DForward::MT2DForward()
      {
      }

    MT2DForward::~MT2DForward()
      {
      }
    void MT2DForward::CalcEpol(const std::vector<double> &Periods)
      {
        const long nionos = 20;
        const long natmos = 20;

        const int nperiods = Periods.size();
        const long nx = XSizes.size();
        const long nz = ZSizes.size();
        const int nzearth = nz;
        const int modelsize = nx * nzearth;
        const int nelements = modelsize * nperiods;
        const double rionos = 1.0;

        Hx_real.resize(boost::extents[nperiods][nx][nzearth]);
        Hx_imag.resize(boost::extents[nperiods][nx][nzearth]);
        Ey_real.resize(boost::extents[nperiods][nx][nzearth]);
        Ey_imag.resize(boost::extents[nperiods][nx][nzearth]);
        std::fill_n(Hx_real.origin(), nelements, 0.0);
        std::fill_n(Hx_imag.origin(), nelements, 0.0);
        std::fill_n(Ey_real.origin(), nelements, 0.0);
        std::fill_n(Ey_imag.origin(), nelements, 0.0);
#pragma omp parallel default(shared)
          {
            const t2DModelDim XS(XSizes);
            t2DModelDim ZS(boost::extents[nz + nionos + natmos]);
            std::fill_n(ZS.begin(), nionos + natmos, 100.0);
            std::copy(ZSizes.begin(), ZSizes.end(), ZS.begin() + nionos
                + natmos);
            const t2DModelData Res(Resistivities);
#pragma omp for
            for (int i = 0; i < nperiods; ++i)
              {
                int startingindex = i * (nx * nzearth);
                epol_(&Periods[i], &nx, &nz, &nionos, &natmos, XS.origin(),
                    ZS.origin(), Res.origin(), &rionos, Hx_real.origin()
                        + startingindex, Hx_imag.origin() + startingindex,
                    Ey_real.origin() + startingindex, Ey_imag.origin()
                        + startingindex);
              }
          }
      }

    void MT2DForward::CalcBpol(const std::vector<double> &Periods)
      {
        const int nperiods = Periods.size();
        const long nx = XSizes.size();
        const long nz = ZSizes.size();
        const int nzearth = nz;
        const int modelsize = nx * nzearth;
        const int nelements = modelsize * nperiods;

        Hy_real.resize(boost::extents[nperiods][nx][nzearth]);
        Hy_imag.resize(boost::extents[nperiods][nx][nzearth]);
        Ex_real.resize(boost::extents[nperiods][nx][nzearth]);
        Ex_imag.resize(boost::extents[nperiods][nx][nzearth]);
        Ez_real.resize(boost::extents[nperiods][nx][nzearth]);
        Ez_imag.resize(boost::extents[nperiods][nx][nzearth]);
        std::fill_n(Hy_real.origin(), nelements, 0.0);
        std::fill_n(Hy_imag.origin(), nelements, 0.0);
        std::fill_n(Ex_real.origin(), nelements, 0.0);
        std::fill_n(Ex_imag.origin(), nelements, 0.0);
        std::fill_n(Ez_real.origin(), nelements, 0.0);
        std::fill_n(Ez_imag.origin(), nelements, 0.0);
#pragma omp parallel default(shared)
          {
            const t2DModelDim XS(XSizes);
            const t2DModelDim ZS(ZSizes);
            const t2DModelData Res(Resistivities);
#pragma omp for
            for (int i = 0; i < nperiods; ++i)
              {
                int startingindex = i * (nx * nzearth);
                hpol_(&Periods[i], &nx, &nz, XS.origin(), ZS.origin(),
                    Res.origin(), Hy_real.origin() + startingindex,
                    Hy_imag.origin() + startingindex, Ex_real.origin()
                        + startingindex, Ex_imag.origin() + startingindex,
                    Ez_real.origin() + startingindex, Ez_imag.origin()
                        + startingindex);
              }
          }
      }

  }
