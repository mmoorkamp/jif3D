//============================================================================
// Name        : TomographyCalculator.cpp
// Author      : Apr 14, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <boost/bind.hpp>
#include "TomographyCalculator.h"
#include "ReadWriteTomographyData.h"
namespace jiba
  {

    TomographyCalculator::TomographyCalculator() :
      nairlayers(3)
      {
      }

    TomographyCalculator::~TomographyCalculator()
      {

      }

    void TomographyCalculator::Allocate(const size_t ngrid, const size_t ndata,
        const size_t npos)
      {
        //the actual size of the forward grid is always one cell extra in each direction
        const size_t oldngrid = (grid.nx + 1) * (grid.ny + 1) * (grid.nz + 1
            + nairlayers);
        //if the grid size changed we have to (re)allocate
        if (oldngrid != ngrid)
          {
            if (grid.slow != NULL)
              {
                delete[] grid.slow;
              }
            grid.slow = new float[ngrid];
          }
        //if the number of data changed we have to (re)allocate
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
        //if the geometry size changed we have to (re)allocate
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

        //first we calculate the size of the actual forward modeling grid
        //this needs to be one larger in each direction than our mdeol
        const size_t ngrid = (Model.GetXCellSizes().size() + 1)
            * (Model.GetYCellSizes().size() + 1)
            * (Model.GetZCellSizes().size() + 1 + nairlayers);
        const size_t ndata = Model.GetSourceIndices().size();
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nshot = Model.GetSourcePosX().size();
        const size_t npos = nmeas + Model.GetSourcePosX().size();
        //allocate the arrays
        Allocate(ngrid, ndata, npos);
        //now we can start to fill the different structures for the forward code
        //first we do the slowness grid
        grid.nx = Model.GetXCellSizes().size();
        grid.ny = Model.GetYCellSizes().size();
        grid.nz = Model.GetZCellSizes().size() + nairlayers;
        grid.h = Model.GetXCellSizes()[0];
        grid.nborder = 0;

        //we take care of the origin in the model object
        std::fill_n(grid.org, 3, 0.0);

        //the extra cells have to be zero, it is easier to initialize everything
        std::fill_n(grid.slow, ngrid, 0.0);
        //fill the real model with slowness values
        //we use physical slowness in s/m, the forward code requires
        //to multiply by the cell size
        //we have to skip the airlayers in the grid
        for (size_t i = 0; i < grid.nx; ++i)
          for (size_t j = 0; j < grid.ny; ++j)
            {
              const size_t layerindex = i * (grid.nz + 1) * (grid.ny + 1) + j * (grid.nz + 1);
              //fill the airlayers
              for (size_t k = 0; k < nairlayers; ++k)
                grid.slow[layerindex + k] = grid.h;
              //then copy the actual model
              for (size_t k = nairlayers; k < grid.nz; ++k)
                grid.slow[layerindex+ k] = Model.GetSlownesses()[i][j][k - nairlayers] * grid.h;
            }

        //fill the data structure
        data.ndata_seis = ndata;
        data.ndata_seis_act = ndata;
        //our indices are 0 based, the forward uses a base of 1
        std::transform(Model.GetSourceIndices().begin(),
            Model.GetSourceIndices().end(), data.sno, boost::bind(
                std::plus<int>(), _1, 1));

        //also we have different storage for source and receivers
        //but the forward stores the positions in a single array
        std::transform(Model.GetReceiverIndices().begin(),
            Model.GetReceiverIndices().end(), data.rno, boost::bind(std::plus<
                int>(), _1, nmeas + 1));

        geo.nrec = nmeas;
        geo.nshot = Model.GetSourcePosX().size();

        //we use the convention that first we store all shot positions
        std::copy(Model.GetSourcePosX().begin(), Model.GetSourcePosX().end(),
            geo.x);
        std::copy(Model.GetSourcePosY().begin(), Model.GetSourcePosY().end(),
            geo.y);
        //we also have to adjust for the offset by the airlayers
        std::transform(Model.GetSourcePosZ().begin(),
            Model.GetSourcePosZ().end(), geo.z, boost::bind(
                std::plus<double>(), _1, grid.h * nairlayers));
        //and then all measurement position
        std::copy(Model.GetMeasPosX().begin(), Model.GetMeasPosX().end(), geo.x
            + nshot);
        std::copy(Model.GetMeasPosY().begin(), Model.GetMeasPosY().end(), geo.y
            + nshot);
        std::transform(Model.GetMeasPosZ().begin(), Model.GetMeasPosZ().end(),
            geo.z + nshot, boost::bind(std::plus<double>(), _1, grid.h
                * nairlayers));

        //now we can do the forward modeling
        ForwardModRay(geo, grid, &data, raypath, 0);
        PlotRaypath("ray.vtk", raypath, ndata, grid.h, nairlayers);
        //and return the result as a vector
        jiba::rvec result(ndata);
        copy(data.tcalc, data.tcalc + ndata, result.begin());
        return result;
      }

    rvec TomographyCalculator::LQDerivative(const ThreeDSeismicModel &Model,
        const rvec &Misfit)
      {
        //in the forward modeling we also perform ray tracing
        //the length of the path in each cell is the sensitivity
        //so here we only have to sort them to the correct model parameter
        //and weight by the misfit
        const size_t nmod = Model.GetSlownesses().num_elements();
        const size_t ndata = Misfit.size();
        assert(ndata == data.ndata_seis);
        jiba::rvec DerivMod(nmod);
        std::fill(DerivMod.begin(), DerivMod.end(), 0.0);

        for (size_t i = 0; i < ndata; ++i)
          {
            const size_t nray = raypath[i].nray;
            for (size_t j = 0; j < nray; ++j)
              {
                const size_t offset = (grid.nz - nairlayers) * grid.ny * (int) floor(raypath[i].x[j])
                + (grid.nz - nairlayers) * (int) floor(raypath[i].y[j]) + (int) floor(raypath[i].z[j]) - nairlayers;

                DerivMod(offset) += raypath[i].len[j] * Misfit(i);
              }
          }
        return DerivMod;
      }
  }
