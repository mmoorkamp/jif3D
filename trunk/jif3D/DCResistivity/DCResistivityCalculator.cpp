//============================================================================
// Name        : DCResistivityCalculator.cpp
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#include "DCResistivityCalculator.h"

namespace jif3D
  {

    void DCResistivityCalculator::ModelToStruct(const ThreeDDCResistivityModel &Model,
        jif3D::GEOMETRY_RES &geo, jif3D::GRID_STRUCT_RES &grid)
      {
        const size_t ndata = Model.GetSourceIndices().size();
        const size_t nmeaspoint = Model.GetMeasPosX().size();
        const size_t nshot = Model.GetSourcePosPosX().size();

        grid.nx = Model.GetModelShape()[0];
        grid.ny = Model.GetModelShape()[1];
        grid.nz = Model.GetModelShape()[2];
        grid.dx = Model.GetXCellSizes()[0];
        grid.dy = Model.GetYCellSizes()[0];
        grid.dz.resize(Model.GetZCellSizes().size());
        std::copy(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
            grid.dz.begin());
        //grid.dz = Model.GetZCellSizes()[0];

        //we calculate the size of the actual forward modeling grid
        const size_t ngrid = grid.nx * grid.ny * grid.nz;
        //allocate the arrays
        Allocate(ngrid, nshot, nmeaspoint, ndata);
        //now we can start to fill the different structures for the forward code
        //first we do the resistivity grid
        //the extra cells have to be zero, it is easier to initialize everything
        //std::fill_n(grid.rho.begin(), ngrid, 0.0);
        //fill the real model with resistivity values
        //we use physical resistivity in ohm.m
        //The layerindex is given according to the sequence in DCResForwardBase code
        for (size_t i = 0; i < grid.nz; ++i)
          for (size_t j = 0; j < grid.ny; ++j)
            for (size_t k = 0; k < grid.nx; ++k)
              {
                const size_t layerindex = i * grid.nx * grid.ny + j * grid.nx + k;
                grid.rho[layerindex] = Model.GetResistivities()[k][j][i];
              }
        //grid.avg_cond = 1/grid.rho[0];
        //fill the geometry
        //In DCResForwardBase code, the coordinate system is base on zero-centred in x and y direction, but actual data's coordinate is not so,
        //we have to adjust source and receiver's coordinate to accord with coordinate system in DCResForwardBase.
        auto xmeasrange1 = std::minmax_element(Model.GetMeasPosX().begin(),
            Model.GetMeasPosX().end());
        auto xmeasrange2 = std::minmax_element(Model.GetMeasSecPosX().begin(),
            Model.GetMeasSecPosX().end());
        auto xpossourcerange = std::minmax_element(Model.GetSourcePosPosX().begin(),
            Model.GetSourcePosPosX().end());
        auto xnegsourcerange = std::minmax_element(Model.GetSourceNegPosX().begin(),
            Model.GetSourceNegPosX().end());
        auto ymeasrange1 = std::minmax_element(Model.GetMeasPosY().begin(),
            Model.GetMeasPosY().end());
        auto ymeasrange2 = std::minmax_element(Model.GetMeasSecPosY().begin(),
            Model.GetMeasSecPosY().end());
        auto ypossourcerange = std::minmax_element(Model.GetSourcePosPosY().begin(),
            Model.GetSourcePosPosY().end());
        auto ynegsourcerange = std::minmax_element(Model.GetSourceNegPosY().begin(),
            Model.GetSourceNegPosY().end());

        double minxmeas = std::min(*xmeasrange1.first, *xmeasrange2.first);
        double minxsource = std::min(*xpossourcerange.first, *xnegsourcerange.first);
        double minx = std::min(minxmeas, minxsource);
        double maxxmeas = std::max(*xmeasrange1.second, *xmeasrange2.second);
        double maxxsource = std::max(*xpossourcerange.second, *xnegsourcerange.second);
        double maxx = std::max(maxxmeas, maxxsource);
        double minymeas = std::min(*ymeasrange1.first, *ymeasrange2.first);
        double minysource = std::min(*ypossourcerange.first, *ynegsourcerange.first);
        double miny = std::min(minymeas, minysource);
        double maxymeas = std::max(*ymeasrange1.second, *ymeasrange2.second);
        double maxysource = std::max(*ypossourcerange.second, *ynegsourcerange.second);
        double maxy = std::max(maxymeas, maxysource);

        const double xshift = -minx - (maxx - minx) / 2.0;
        const double yshift = -miny - (maxy - miny) / 2.0;
        std::transform(Model.GetSourcePosPosX().begin(), Model.GetSourcePosPosX().end(),
            geo.PosSx.begin(), [xshift](double pos)
              { return pos+xshift;});
        std::transform(Model.GetSourcePosPosY().begin(), Model.GetSourcePosPosY().end(),
            geo.PosSy.begin(), [yshift](double pos)
              { return pos+yshift;});
        std::transform(Model.GetSourceNegPosX().begin(), Model.GetSourceNegPosX().end(),
            geo.NegSx.begin(), [xshift](double pos)
              { return pos+xshift;});
        std::transform(Model.GetSourceNegPosY().begin(), Model.GetSourceNegPosY().end(),
            geo.NegSy.begin(), [yshift](double pos)
              { return pos+yshift;});
        std::transform(Model.GetMeasPosX().begin(), Model.GetMeasPosX().end(),
            geo.rx1.begin(), [xshift](double pos)
              { return pos+xshift;});
        std::transform(Model.GetMeasPosY().begin(), Model.GetMeasPosY().end(),
            geo.ry1.begin(), [yshift](double pos)
              { return pos+yshift;});
        std::transform(Model.GetMeasSecPosX().begin(), Model.GetMeasSecPosX().end(),
            geo.rx2.begin(), [xshift](double pos)
              { return pos+xshift;});
        std::transform(Model.GetMeasSecPosY().begin(), Model.GetMeasSecPosY().end(),
            geo.ry2.begin(), [yshift](double pos)
              { return pos+yshift;});
        //we also need to change Z coordinate to the cell centre of the first layer
        const double zshift = grid.dz[0] / 2.0;
        if (Model.GetSourcePosPosZ()[0] == 0.0)
          {
            std::transform(Model.GetSourcePosPosZ().begin(),
                Model.GetSourcePosPosZ().end(), geo.PosSz.begin(), [zshift](double pos)
                  { return pos+zshift;});
            std::transform(Model.GetSourceNegPosZ().begin(),
                Model.GetSourceNegPosZ().end(), geo.NegSz.begin(), [zshift](double pos)
                  { return pos+zshift;});
            std::transform(Model.GetMeasPosZ().begin(), Model.GetMeasPosZ().end(),
                geo.rz1.begin(), [zshift](double pos)
                  { return pos+zshift;});
            std::transform(Model.GetMeasSecPosZ().begin(), Model.GetMeasSecPosZ().end(),
                geo.rz2.begin(), [zshift](double pos)
                  { return pos+zshift;});
          }
        else
          {
            geo.PosSz = Model.GetSourcePosPosZ();
            geo.NegSz = Model.GetSourceNegPosZ();
            geo.rz1 = Model.GetMeasPosZ();
            geo.rz2 = Model.GetMeasPosZ();

          }
        geo.nsource = nshot;
        std::copy(Model.GetSourceIndices().begin(), Model.GetSourceIndices().end(),
            geo.sno.begin());

      }

    DCResistivityCalculator::DCResistivityCalculator() :
        geo(), grid()
      {
      }

    DCResistivityCalculator::~DCResistivityCalculator()
      {

      }

    void DCResistivityCalculator::Allocate(const size_t ngrid, const size_t nshot,
        const size_t nmeaspoint, const size_t ndata)
      {
        //the actual size of the forward grid
        const size_t oldngrid = grid.rho.size();
        //if the grid size changed we have to (re)allocate
        if (oldngrid != ngrid)
          {
            grid.rho.resize(ngrid);
          }
        //if the geometry size changed we have to (re)allocate
        if (geo.nsource != nshot)
          {
            geo.PosSx.resize(nshot);
            geo.PosSy.resize(nshot);
            geo.PosSz.resize(nshot);
            geo.NegSx.resize(nshot);
            geo.NegSy.resize(nshot);
            geo.NegSz.resize(nshot);
          }
        if (geo.sno.size() != nmeaspoint)
          {
            geo.rx1.resize(nmeaspoint);
            geo.ry1.resize(nmeaspoint);
            geo.rz1.resize(nmeaspoint);
            geo.rx2.resize(nmeaspoint);
            geo.ry2.resize(nmeaspoint);
            geo.rz2.resize(nmeaspoint);
            geo.sno.resize(nmeaspoint);
          }

      }

    rvec DCResistivityCalculator::Calculate(const ModelType &Model)
      {
        const size_t ndata = Model.GetSourceIndices().size();
        ModelToStruct(Model, geo, grid);
        //now we can do the forward modeling
        jif3D::rvec result = ResForward(geo, grid, ndata);
        //return the result as a vector

        return result;

      }

    rvec DCResistivityCalculator::LQDerivative(const ModelType &Model, const rvec &Misfit)
      {

        ModelToStruct(Model, geo, grid);
        jif3D::rvec result = ResGradient(geo, grid, Misfit);

        return result;

      }
  }

