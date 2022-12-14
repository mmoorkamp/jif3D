//============================================================================
// Name        : TomographyCalculator.cpp
// Author      : Apr 14, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Tomo/TomographyCalculator.h"
#include "../Tomo/modeling_seismic.h"
#include "../Tomo/ReadWriteTomographyData.h"

#include <boost/numeric/conversion/cast.hpp>
#include <algorithm>

namespace jif3D
  {

    using boost::numeric_cast;

    TomographyCalculator::TomographyCalculator() :
        minxindex(), minyindex(), nairlayers(3), geo(), grid(), data(), raypath()
      {
      }

    TomographyCalculator::~TomographyCalculator()
      {

      }

    void TomographyCalculator::WriteRays(const std::string &filename) const
      {

        PlotRaypath("ray.vtk", raypath, data.ndata_seis, grid.h, nairlayers, minxindex,
            minyindex);
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

    rvec TomographyCalculator::Calculate(const ThreeDSeismicModel &Model,
        const jif3D::TomographyData &Data)
      {
        const size_t ndata = Data.GetSourceIndices().size();
        const size_t nmeas = Data.GetMeasPosX().size();
        const size_t nshot = Data.GetSourcePosX().size();
        const size_t npos = nmeas + nshot;
        //we assume that all grid cells have the same size in all directions, so we just
        //read one of the values
        const double gridspacing = Model.GetXCellSizes()[0];

        const double orx = Model.GetXOrigin();
        const double ory = Model.GetYOrigin();
        const double orz = Model.GetZOrigin();

        //now check that this assumption holds
        if (!std::all_of(Model.GetXCellSizes().begin(), Model.GetXCellSizes().end(),
            [gridspacing](double d)
              { return ((gridspacing * 0.99 < d) && (gridspacing * 1.01 > d));}))
          {
            throw jif3D::FatalException(
                "Spacing of tomographic model not uniform in x-direction");
          }
        if (!std::all_of(Model.GetYCellSizes().begin(), Model.GetYCellSizes().end(),
            [gridspacing](double d)
              { return ((gridspacing * 0.99 < d) && (gridspacing * 1.01 > d));}))
          {
            throw jif3D::FatalException(
                "Spacing of tomographic model not uniform in y-direction");
          }
        if (!std::all_of(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
            [gridspacing](double d)
              { return ((gridspacing * 0.99 < d) && (gridspacing * 1.01 > d));}))
          {
            throw jif3D::FatalException(
                "Spacing of tomographic model not uniform in z-direction");
          }
        const int padding = 5;
        auto xmeasrange = std::minmax_element(Data.GetMeasPosX().begin(),
            Data.GetMeasPosX().end());
        auto ymeasrange = std::minmax_element(Data.GetMeasPosY().begin(),
            Data.GetMeasPosY().end());
        auto xsourcerange = std::minmax_element(Data.GetSourcePosX().begin(),
            Data.GetSourcePosX().end());
        auto ysourcerange = std::minmax_element(Data.GetSourcePosY().begin(),
            Data.GetSourcePosY().end());

        double minx = std::min(*xmeasrange.first, *xsourcerange.first) - orx;
        double maxx = std::max(*xmeasrange.second, *xsourcerange.second) - orx;
        double miny = std::min(*ymeasrange.first, *ysourcerange.first) - ory;
        double maxy = std::max(*ymeasrange.second, *ysourcerange.second) - ory;

        minxindex = std::max(numeric_cast<int>(std::floor(minx / gridspacing)) - padding,
            0);
        int maxxindex = std::min(
            numeric_cast<int>(std::ceil(maxx / gridspacing)) + padding,
            numeric_cast<int>(Model.GetModelShape()[0]));
        minyindex = std::max(numeric_cast<int>(std::floor(miny / gridspacing)) - padding,
            0);
        int maxyindex = std::min(
            numeric_cast<int>(std::ceil(maxy / gridspacing)) + padding,
            numeric_cast<int>(Model.GetModelShape()[1]));

        const double xshift = minxindex * gridspacing;
        const double yshift = minyindex * gridspacing;

        typedef boost::multi_array_types::index_range range;
        // OR typedef array_type::index_range range;
        ThreeDModelBase::t3DModelData::const_array_view<3>::type SubGrid =
            Model.GetSlownesses()[boost::indices[range(minxindex, maxxindex)][range(
                minyindex, maxyindex)][range(0, Model.GetModelShape()[2])]];

        grid.nx = SubGrid.shape()[0];
        grid.ny = SubGrid.shape()[1];
        grid.nz = SubGrid.shape()[2] + 2 * nairlayers;
        grid.h = gridspacing;

        //we calculate the size of the actual forward modeling grid
        //this needs to be one larger in each direction than our model
        const size_t ngrid = (grid.nx + 1) * (grid.ny + 1) * (grid.nz + 1);
        //allocate the arrays
        Allocate(ngrid, ndata, npos);
        //now we can start to fill the different structures for the forward code
        //first we do the slowness grid
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
              //fill the side padding areas
              std::fill_n(grid.slow.begin() + layerindex, nairlayers, grid.h);
              std::fill_n(grid.slow.begin() + layerindex + grid.nz - nairlayers,
                  nairlayers, grid.h);
              //then copy the actual model
              for (size_t k = nairlayers; k < grid.nz - nairlayers; ++k)
                grid.slow[layerindex + k] = SubGrid[i][j][k - nairlayers] * grid.h;
            }

        //fill the data structure
        data.ndata_seis = ndata;
        data.ndata_seis_act = ndata;
        //our indices are 0 based, the forward uses a base of 1
        std::transform(Data.GetSourceIndices().begin(), Data.GetSourceIndices().end(),
            data.sno.begin(), [](int i)
              { return i+1;});

        //also we have different storage for source and receivers
        //but the forward stores the positions in a single array
        //our indices are 0 based, the forward uses a base of 1
        //and we have to add the number of shots that we already stored
        std::transform(Data.GetReceiverIndices().begin(),
            Data.GetReceiverIndices().end(), data.rno.begin(), [nshot](int i)
              { return i+nshot+1;});

        geo.nrec = nmeas;
        geo.nshot = nshot;

        //we use the convention that first we store all shot positions
        //we have to adjust for the fact that we are only calculating on
        //a subgrid

        std::transform(Data.GetSourcePosX().begin(), Data.GetSourcePosX().end(),
            geo.x.begin(), [xshift,orx](double pos)
              { return pos-xshift-orx;});
        std::transform(Data.GetSourcePosY().begin(), Data.GetSourcePosY().end(),
            geo.y.begin(), [yshift,ory](double pos)
              { return pos-yshift -ory;});
        //we also have to adjust for the offset by the airlayers
        const double zoffset = grid.h * nairlayers;
        std::transform(Data.GetSourcePosZ().begin(), Data.GetSourcePosZ().end(),
            geo.z.begin(), [zoffset,orz] (double pos)
              { return pos + zoffset-orz;});
        //and then all measurement position
        std::transform(Data.GetMeasPosX().begin(), Data.GetMeasPosX().end(),
            geo.x.begin() + nshot, [xshift,orx](double pos)
              { return pos-xshift-orx;});
        std::transform(Data.GetMeasPosY().begin(), Data.GetMeasPosY().end(),
            geo.y.begin() + nshot, [yshift,ory](double pos)
              { return pos-yshift-ory;});
        std::transform(Data.GetMeasPosZ().begin(), Data.GetMeasPosZ().end(),
            geo.z.begin() + nshot, [zoffset,orz] (double pos)
              { return pos + zoffset-orz;});

        //now we can do the forward modeling
        ForwardModRay(geo, grid, data, raypath);

        //and return the result as a vector
        jif3D::rvec result(ndata, 0.0);
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
        const jif3D::TomographyData &Data, const rvec &Misfit)
      {
        //in the forward modeling we also perform ray tracing
        //the length of the path in each cell is the sensitivity
        //so here we only have to sort them to the correct model parameter
        //and weight by the misfit
        const size_t nmod = Model.GetSlownesses().num_elements();
        const size_t ndata = Misfit.size();
        assert(ndata == data.ndata_seis);
        jif3D::rvec DerivMod(nmod);
        std::fill(DerivMod.begin(), DerivMod.end(), 0.0);
        //the gradient is simply the length of the raypath
        //through each cell weighted by the data misfit
        const size_t yextent = Model.GetModelShape()[1];
        const size_t zextent = Model.GetModelShape()[2];
        for (size_t i = 0; i < ndata; ++i)
          {
            const size_t nray = raypath[i].nray;
            for (size_t j = 0; j < nray; ++j)
              {
                const size_t offset = (zextent) * yextent
                    * floor(raypath[i].x[j] + minxindex)
                    + zextent * floor(raypath[i].y[j] + minyindex)
                    + floor(raypath[i].z[j]) - nairlayers;

                DerivMod(offset) += 2.0 * raypath[i].len[j] * grid.h * Misfit(i);
              }
          }
        return DerivMod;
      }
  }

