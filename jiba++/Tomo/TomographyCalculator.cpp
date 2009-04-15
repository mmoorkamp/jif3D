//============================================================================
// Name        : TomographyCalculator.cpp
// Author      : Apr 14, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <boost/bind.hpp>
#include "TomographyCalculator.h"

namespace jiba
  {

    TomographyCalculator::TomographyCalculator()
      {
      }

    TomographyCalculator::~TomographyCalculator()
      {

      }

    void TomographyCalculator::Allocate(const size_t ngrid, const size_t ndata,
        const size_t npos)
      {
        const size_t oldngrid = (grid.nx + 1) * (grid.ny + 1) * (grid.nz + 1);
        if (oldngrid != ngrid)
          {
            if (grid.slow != NULL)
              {
                delete[] grid.slow;
              }
            grid.slow = new double[ngrid];
          }
        if (data.ndata_seis != ndata)
          {
            if (data.sno != NULL)
              {
                delete[] data.sno;
                delete[] data.rno;
                delete[] data.tcalc;
                delete[] raypath;
              }
            raypath = new jiba::RP_STRUCT[ndata];
            data.sno = new int[ndata];
            data.rno = new int[ndata];
            data.tcalc = new double[ndata];
          }
        if (geo.nrec + geo.nshot != npos)
          {
            if (geo.x != NULL)
              {
                delete[] geo.x;
                delete[] geo.y;
                delete[] geo.z;
              }
            geo.x = new float[npos];
            geo.y = new float[npos];
            geo.z = new float[npos];
          }

      }

    rvec TomographyCalculator::Calculate(const ThreeDSeismicModel &Model)
      {

        const size_t ngrid = (Model.GetXCellSizes().size() + 1)
            * (Model.GetYCellSizes().size() + 1)
            * (Model.GetZCellSizes().size() + 1);
        const size_t ndata = Model.GetSourceIndices().size();
        const size_t npos = Model.GetMeasPosX().size()
            + Model.GetSourcePosX().size();
        Allocate(ngrid, ndata, npos);
        grid.nx = Model.GetXCellSizes().size();
        grid.ny = Model.GetYCellSizes().size();
        grid.nz = Model.GetZCellSizes().size();
        grid.h = Model.GetXCellSizes()[0];
        grid.nborder = 0;

        //we take care of the origin in the model object
        std::fill_n(grid.org, 3, 0.0);

        std::fill_n(grid.slow, ngrid, 0.0);
        for (size_t i = 0; i < grid.nx; ++i)
          for (size_t j = 0; j < grid.ny; ++j)
            for (size_t k = 0; k < grid.nz; ++k)
              grid.slow[i * (grid.nz + 1) * (grid.ny + 1) + j * (grid.nz + 1)
                  + k] = Model.GetSlownesses()[i][j][k] * grid.h;

        data.ndata_seis = ndata;
        data.ndata_seis_act = ndata;
        std::transform(Model.GetReceiverIndices().begin(),
            Model.GetReceiverIndices().end(), data.sno, boost::bind(std::plus<
                int>(), _1, 1));
        std::transform(Model.GetSourceIndices().begin(),
            Model.GetSourceIndices().end(), data.rno, boost::bind(
                std::plus<int>(), _1, Model.GetMeasPosX().size() + 1));

        geo.nrec = Model.GetMeasPosX().size();
        geo.nshot = Model.GetSourcePosX().size();

        std::copy(Model.GetMeasPosX().begin(), Model.GetMeasPosX().end(), geo.x);
        std::copy(Model.GetMeasPosY().begin(), Model.GetMeasPosY().end(), geo.y);
        std::copy(Model.GetMeasPosZ().begin(), Model.GetMeasPosZ().end(), geo.z);
        std::copy(Model.GetSourcePosX().begin(), Model.GetSourcePosX().end(),
            geo.x + Model.GetMeasPosX().size());
        std::copy(Model.GetSourcePosY().begin(), Model.GetSourcePosY().end(),
            geo.y + Model.GetMeasPosY().size());
        std::copy(Model.GetSourcePosZ().begin(), Model.GetSourcePosZ().end(),
            geo.z + Model.GetMeasPosZ().size());

        ForwardModRay(geo, grid, &data, raypath, 0);
        jiba::rvec result(ndata);
        copy(data.tcalc, data.tcalc + ndata, result.begin());
        return result;
      }

    rvec TomographyCalculator::LQDerivative(const ThreeDSeismicModel &Model,
        const rvec &Misfit)
      {
        const size_t nmod = Model.GetSlownesses().num_elements();
        const size_t ndata = Misfit.size();
        assert(ndata == data.ndata_seis);
        jiba::rvec DerivMod(nmod);
        for (size_t i = 0; i < ndata; ++i)
          {
            const size_t nray = raypath[i].nray;
            for (size_t j = 0; j < nray; ++j)
              {
                DerivMod(raypath[i].ele[j]) = raypath[i].len[j] * Misfit(i);
            }
        }
        return DerivMod;
    }
}
