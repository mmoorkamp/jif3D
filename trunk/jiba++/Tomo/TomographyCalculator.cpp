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

    TomographyCalculator::TomographyCalculator(bool saverays) :
        writerays(saverays), nairlayers(3), geo(), grid(), data(), raypath()
      {
      }

    TomographyCalculator::~TomographyCalculator()
      {

      }

    void TomographyCalculator::Allocate(const size_t ngrid, const size_t ndata,
        const size_t npos)
      {
        //the actual size of the forward grid is always one cell extra in each direction
        const size_t oldngrid = (grid.nx + 1) * (grid.ny + 1)
            * (grid.nz + 1 + nairlayers);
        //if the grid size changed we have to (re)allocate
        if (oldngrid != ngrid)
          {
            grid.slow.resize(ngrid);
          }
        //if the number of data changed we have to (re)allocate
        if (data.ndata_seis != ndata)
          {
            raypath.resize(ndata);
            data.sno.resize(ndata);
            data.rno.resize(ndata);
            data.tcalc.resize(ndata);
          }
        //if the geometry size changed we have to (re)allocate
        if (geo.nrec + geo.nshot != npos)
          {
            geo.x.resize(npos);
            geo.y.resize(npos);
            geo.z.resize(npos);
          }

      }

    rvec TomographyCalculator::Calculate(const ThreeDSeismicModel &Model)
      {

        //first we calculate the size of the actual forward modeling grid
        //this needs to be one larger in each direction than our model
        const size_t ngrid = (Model.GetXCellSizes().size() + 1)
            * (Model.GetYCellSizes().size() + 1)
            * (Model.GetZCellSizes().size() + 1 + 2 * nairlayers);
        const size_t ndata = Model.GetSourceIndices().size();
        const size_t nmeas = Model.GetMeasPosX().size();
        const size_t nshot = Model.GetSourcePosX().size();
        const size_t npos = nmeas + nshot;
        //allocate the arrays
        Allocate(ngrid, ndata, npos);
        //now we can start to fill the different structures for the forward code
        //first we do the slowness grid
        grid.nx = Model.GetXCellSizes().size();
        grid.ny = Model.GetYCellSizes().size();
        grid.nz = Model.GetZCellSizes().size() + 2 * nairlayers;
        grid.h = Model.GetXCellSizes()[0];

        //the extra cells have to be zero, it is easier to initialize everything
        std::fill_n(grid.slow.begin(), ngrid, 0.0);
        //fill the real model with slowness values
        //we use physical slowness in s/m, the forward code requires
        //to multiply by the cell size
        //we have to skip the airlayers in the grid
        for (size_t i = 0; i < grid.nx; ++i)
          for (size_t j = 0; j < grid.ny; ++j)
            {
              const size_t layerindex = i * (grid.nz + 1) * (grid.ny + 1)
                  + j * (grid.nz + 1);
              //fill the airlayers
              std::fill_n(grid.slow.begin() + layerindex, nairlayers, grid.h);
              std::fill_n(grid.slow.begin() + layerindex + grid.nz - nairlayers,
                  nairlayers, grid.h);
              //then copy the actual model
              for (size_t k = nairlayers; k < grid.nz - nairlayers; ++k)
                grid.slow[layerindex + k] = Model.GetSlownesses()[i][j][k
                    - nairlayers] * grid.h;
            }

        //fill the data structure
        data.ndata_seis = ndata;
        data.ndata_seis_act = ndata;
        //our indices are 0 based, the forward uses a base of 1
        std::transform(Model.GetSourceIndices().begin(),
            Model.GetSourceIndices().end(), data.sno.begin(),
            boost::bind(std::plus<int>(), _1, 1));

        //also we have different storage for source and receivers
        //but the forward stores the positions in a single array
        //our indices are 0 based, the forward uses a base of 1
        //and we have to add the number of shots that we already stored
        std::transform(Model.GetReceiverIndices().begin(),
            Model.GetReceiverIndices().end(), data.rno.begin(),
            boost::bind(std::plus<int>(), _1, nshot + 1));

        geo.nrec = nmeas;
        geo.nshot = nshot;

        //we use the convention that first we store all shot positions
        std::copy(Model.GetSourcePosX().begin(), Model.GetSourcePosX().end(),
            geo.x.begin());
        std::copy(Model.GetSourcePosY().begin(), Model.GetSourcePosY().end(),
            geo.y.begin());
        //we also have to adjust for the offset by the airlayers
        std::transform(Model.GetSourcePosZ().begin(),
            Model.GetSourcePosZ().end(), geo.z.begin(),
            boost::bind(std::plus<double>(), _1, grid.h * nairlayers));
        //and then all measurement position
        std::copy(Model.GetMeasPosX().begin(), Model.GetMeasPosX().end(),
            geo.x.begin() + nshot);
        std::copy(Model.GetMeasPosY().begin(), Model.GetMeasPosY().end(),
            geo.y.begin() + nshot);
        std::transform(Model.GetMeasPosZ().begin(), Model.GetMeasPosZ().end(),
            geo.z.begin() + nshot,
            boost::bind(std::plus<double>(), _1, grid.h * nairlayers));

        //now we can do the forward modeling
        ForwardModRay(geo, grid, &data, &raypath[0]);
        //we only write out the file with the rays, if the corresponding
        //option is set to true
        if (writerays)
          {
            PlotRaypath("ray.vtk", &raypath[0], ndata, grid.h, nairlayers);
          }
        //and return the result as a vector
        jiba::rvec result(ndata, 0.0);
        /*for (size_t i = 0; i < ndata; ++i)
         {
         const size_t nray = raypath[i].nray;
         for (size_t j = 0; j < nray; ++j)
         {
         const size_t offset = (grid.nz - 2 * nairlayers) * grid.ny
         * floor(raypath[i].x[j]) + (grid.nz - 2 * nairlayers)
         * floor(raypath[i].y[j]) + floor(raypath[i].z[j])
         - nairlayers;
         result(i) += raypath[i].len[j] * grid.h * Model.GetSlownesses().data()[offset];
         }

         }*/
        std::copy(data.tcalc.begin(), data.tcalc.end(), result.begin());
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
        //the gradient is simply the length of the raypath
        //through each cell weighted by the data misfit
        for (size_t i = 0; i < ndata; ++i)
          {
            const size_t nray = raypath[i].nray;
            for (size_t j = 0; j < nray; ++j)
              {
                const size_t offset = (grid.nz - 2 * nairlayers) * grid.ny
                    * floor(raypath[i].x[j])
                    + (grid.nz - 2 * nairlayers) * floor(raypath[i].y[j])
                    + floor(raypath[i].z[j]) - nairlayers;

                DerivMod(offset) += 2.0 * raypath[i].len[j] * grid.h
                    * Misfit(i);
              }
          }
        return DerivMod;
      }
  }
