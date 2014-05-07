//============================================================================
// Name        : DCResistivityCalculator.cpp
// Author      : Zhanjie Shi and Richard.W Hobbs
// Version     : April 2014
// Copyright   : 2014, Zhanjie Shi and Richard.W Hobbs
//============================================================================

#include "DCResistivityCalculator.h"
#include <boost/bind.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace jif3D
  {

    using boost::numeric_cast;

    DCResistivityCalculator::DCResistivityCalculator() :
        geo(), grid(), data()
      {
      }

    DCResistivityCalculator::~DCResistivityCalculator()
      {

      }

    void DCResistivityCalculator::Allocate(const size_t ngrid, const size_t ndata,
        const size_t nshot, const size_t nmeaspoint)
      {
        //the actual size of the forward grid
        const size_t oldngrid = grid.nx * grid.ny * grid.nz;
        //if the grid size changed we have to (re)allocate
        if (oldngrid != ngrid)
          {
            grid.rho.resize(ngrid);
          }
        //if the number of data changed we have to (re)allocate
        if (data.ndata_res != ndata)
          {
            data.dcal.resize(ndata);
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
          }

      }

    rvec DCResistivityCalculator::Calculate(const ModelType &Model)
      {
        const size_t ndata = Model.GetSourceIndices().size();
        const size_t nmeaspoint = Model.GetMeasPosX().size();
        const size_t nshot = Model.GetSourcePosPosX().size();

        grid.nx = Model.GetModelShape()[0];
        grid.ny = Model.GetModelShape()[1];
        grid.nz = Model.GetModelShape()[2];
        grid.dx = Model.GetXCellSizes()[0];
        grid.dy = Model.GetYCellSizes()[0];
        std::copy(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
            grid.dz.begin());
        //grid.dz = Model.GetZCellSizes()[0];

        //we calculate the size of the actual forward modeling grid
        const size_t ngrid = grid.nx * grid.ny * grid.nz;
        //allocate the arrays
        Allocate(ngrid, ndata, nshot, nmeaspoint);
        //now we can start to fill the different structures for the forward code
        //first we do the resistivity grid
        //the extra cells have to be zero, it is easier to initialize everything
        std::fill_n(grid.rho.begin(), ngrid, 0.0);
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

        //fill the data structure
        data.ndata_res = ndata;

        //fill the geometry
        geo.nsource = nshot;
        std::transform(Model.GetSourceIndices().begin(), Model.GetSourceIndices().end(),
            geo.sno.begin(), [](int i)
              { return i;});
        //In DCResForwardBase code, the coordinate system is base on zero-centre in x and y direction, but actual coordinate zero is often on the upper left corner,
        //so we have to adjust source and receiver's coordinate to accord with that.
        //Moreover, for effective interpolation operation, there is at least one cell beyond the boundary of sources and receivers.
        const double xshift = grid.dx - (grid.nx * grid.dx) / 2.0;
        const double yshift = grid.dy - (grid.ny * grid.dy) / 2.0;
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
        std::transform(Model.GetSourcePosPosZ().begin(), Model.GetSourcePosPosZ().end(),
            geo.PosSz.begin(), [zshift](double pos)
              { return pos+zshift;});
        std::transform(Model.GetSourceNegPosZ().begin(), Model.GetSourceNegPosZ().end(),
            geo.NegSz.begin(), [zshift](double pos)
              { return pos+zshift;});
        std::transform(Model.GetMeasPosZ().begin(), Model.GetMeasPosZ().end(),
            geo.rz1.begin(), [zshift](double pos)
              { return pos+zshift;});
        std::transform(Model.GetMeasSecPosZ().begin(), Model.GetMeasSecPosZ().end(),
            geo.rz2.begin(), [zshift](double pos)
              { return pos+zshift;});

        //now we can do the forward modeling
        ResForward(geo, grid, &data);
        //return the result as a vector
        jif3D::rvec result(ndata, 0.0);
        std::copy(data.dcal.begin(), data.dcal.end(), result.begin());
        return result;
      }

    rvec DCResistivityCalculator::LQDerivative(const ModelType &Model, const rvec &Misfit)
      {

      }
  }
BOOST_CLASS_EXPORT_IMPLEMENT(jif3D::DCResistivityCalculator)

