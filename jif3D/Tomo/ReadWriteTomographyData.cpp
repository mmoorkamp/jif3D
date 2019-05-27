//============================================================================
// Name        : ReadWriteTomographyData.cpp
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Tomo/ReadWriteTomographyData.h"

#include "../Global/NumUtil.h"
#include "../Global/NetCDFTools.h"
#include "../Global/NetCDFPortHelper.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../Tomo/TomographyData.h"
#include<fstream>
#include <netcdf>


using netCDF::NcFile;
using netCDF::NcDim;
using netCDF::NcVar;

namespace jif3D
  {

    void PlotRaypath(const std::string &filename, const std::vector<RP_STRUCT> &raypath,
        const size_t nmeas, const double gridspacing, const size_t nairlayers,
        int minxindex, int minyindex)
      {
        std::ofstream outfile(filename.c_str());
        //write out the old style vtk header
        outfile << "# vtk DataFile Version 2.0\n";
        outfile << "Raypaths\n";
        outfile << "ASCII\n";
        outfile << "DATASET POLYDATA\n";
        size_t npoints = 0;
        size_t act_rays = 0;
        //count the number of points we need to plot the rays
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                ++act_rays;
                npoints += raypath[i].nray + 1;
              }

          }
        outfile << "POINTS " << npoints << " double\n";
        //write out the positions of each ray segment start and endpoint
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                for (size_t j = 0; j < raypath[i].nray + 1; ++j)
                  {
                    outfile << (raypath[i].x[j] + minxindex) * gridspacing << " "
                        << (raypath[i].y[j] + minyindex) * gridspacing << " "
                        << (raypath[i].z[j] - nairlayers) * gridspacing << "\n ";
                  }
              }
          }
        //now we connect the points by lines
        //by writing out the indices of the start and endpoint
        outfile << "\nLINES " << act_rays << " " << npoints + act_rays << std::endl;
        size_t index = 0;
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                outfile << raypath[i].nray + 1;
                for (size_t j = 0; j < raypath[i].nray + 1; ++j)
                  {
                    outfile << " " << index;
                    ++index;
                  }
                outfile << std::endl;
              }
          }
        //finally we write out the receiver positions
        //as an independent set of points for plotting
        outfile << "POINT_DATA " << npoints << std::endl;
        outfile << "SCALARS  Rays float\n";
        outfile << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < nmeas; ++i)
          {
            if (raypath[i].nray > 0)
              {
                std::fill_n(std::ostream_iterator<double>(outfile, "\n"),
                    raypath[i].nray + 1, double(i));
              }
          }
      }
  }
