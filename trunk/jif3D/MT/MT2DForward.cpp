//============================================================================
// Name        : MT2DForward.cpp
// Author      : Sep 24, 2008
// Version     :
// Copyright   : 2008, mmoorkamp
//============================================================================


#include "MT2DForward.h"
#include "Tarits2DMT.h"
namespace jif3D
  {

    MT2DForward::MT2DForward() :
      XSizes(), ZSizes(), Resistivities(), Hx_real(), Hx_imag(), Hy_real(),
          Hy_imag(), Ey_real(), Ey_imag(), Ex_real(), Ex_imag(), Ez_real(),
          Ez_imag()
      {
      }

    MT2DForward::~MT2DForward()
      {
      }

    void MT2DForward::CalcEpol(const std::vector<double> &Periods)
      {
        //we need a number of extra layers above the actual modeling domain
        const long nionos = 20;
        const long natmos = 20;

        //we define a number of constants that are used below
        const int nperiods = Periods.size();
        const long nx = XSizes.size();
        const int nzearth = ZSizes.size();
        const long nz = nzearth + nionos + natmos;

        const int modelsize = nx * nzearth;
        const int nelements = modelsize * nperiods;
        const double rionos = 1.0;

        //Initialize the field solutions to zero for
        //all periods and all model cells
        Hx_real.resize(boost::extents[nperiods][nx][nzearth]);
        Hx_imag.resize(boost::extents[nperiods][nx][nzearth]);
        Ey_real.resize(boost::extents[nperiods][nx][nzearth]);
        Ey_imag.resize(boost::extents[nperiods][nx][nzearth]);
        std::fill_n(Hx_real.origin(), nelements, 0.0);
        std::fill_n(Hx_imag.origin(), nelements, 0.0);
        std::fill_n(Ey_real.origin(), nelements, 0.0);
        std::fill_n(Ey_imag.origin(), nelements, 0.0);
        //we parallelize the calculation by frequency
        //this might not be the most effective way
        //but is easy to implement
#pragma omp parallel default(shared)
          {
            //extend the modeling domain in z-direction
            //by the extra layers
            const t2DModelDim XS(XSizes);
            t2DModelDim ZS(boost::extents[nz + nionos + natmos]);
            std::fill_n(ZS.begin(), nionos + natmos, 100.0);
            std::copy(ZSizes.begin(), ZSizes.end(), ZS.begin() + nionos
                + natmos);
#pragma omp for
            for (int i = 0; i < nperiods; ++i)
              {
                //calculate the index to store the solution for
                //the current period
                int startingindex = i * (nx * nzearth);
                //call the Fortran code
                epol_(&Periods[i], &nx, &nz, &nionos, &natmos, XS.origin(),
                    ZS.origin(), Resistivities.origin(), &rionos,
                    Hx_real.origin() + startingindex, Hx_imag.origin()
                        + startingindex, Ey_real.origin() + startingindex,
                    Ey_imag.origin() + startingindex);
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

        //Initialize the field solutions to zero for
          //all periods and all model cells
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
        //same as for CalcEPol, we parallelize by frequency
#pragma omp parallel default(shared)
          {
            const t2DModelDim XS(XSizes);
            const t2DModelDim ZS(ZSizes);
#pragma omp for
            for (int i = 0; i < nperiods; ++i)
              {
                int startingindex = i * (nx * nzearth);
                hpol_(&Periods[i], &nx, &nz, XS.origin(), ZS.origin(),
                    Resistivities.origin(), Hy_real.origin() + startingindex,
                    Hy_imag.origin() + startingindex, Ex_real.origin()
                        + startingindex, Ex_imag.origin() + startingindex,
                    Ez_real.origin() + startingindex, Ez_imag.origin()
                        + startingindex);
              }
          }
      }

  }
