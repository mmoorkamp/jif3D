//============================================================================
// Name        : makegravmesh.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file makemtmesh.cpp
 * Make a netcdf conductivity model file with a specified mesh. The conductivities in this file are all identical to the specified value.
 */

#include "../Global/Jif3DGlobal.h"
#include "../Global/FileUtil.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/mba.hpp"
#include "X3DModel.h"
#include "MTData.h"
#include "TipperData.h"
#include "ReadWriteImpedances.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <boost/multi_array.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <iostream>
#include <string>

using namespace std;

/*! \file makemtmesh.cpp
 * Make a forward modeling mesh for X3D. The program asks for the size in the three coordinate directions and
 * the cell size for each direction. All cells will have the same size for the two horizontal directions.
 * For the vertical direction we can specify the thickness of the top layer and a factor by which each cell
 * increases in thickness compared to the previous layer.
 * The program also asks for a conductivity value that the mesh will be filled with.
 * It will also create a layered background with the same number of layers as cells in x-direction and filled with
 * the same conductivity value as the mesh.
 */
namespace po = boost::program_options;

int main(int argc, char *argv[])
  {
    int incstart = 0;
    double rounding = 1.0;
    double dens = 0.0;
    double aircond = 1e-5;
    std::string DataName;
    std::string TopoName;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("rounding",
        po::value(&rounding)->default_value(1.0),
        "Round layer thicknesses to multiple of this number in meters.")("incstart",
        po::value(&incstart)->default_value(0),
        "Index of the layer where to start increasing the cell size in z-direction")(
        "center", po::value(&DataName),
        "Adjust the origin so the data in the specified file sits in the center of the mesh")(
        "writegrav", po::value(&dens)->default_value(0.0),
        "Write a gravity model with the same geometry. Specify the density value in the grid as a parameter")(
        "topofile", po::value(&TopoName), "File with topography information.")("aircond",
        po::value(&aircond)->default_value(1e-5), "Conductivity for the air");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }

    jif3D::X3DModel Model;
    int nx, ny, nz;
    double deltax, deltay, deltaz;
    //first find out the basic mesh parameters
    //the number of cells in each coordinate direction
    cout << "Number of cells in x-direction: ";
    cin >> nx;
    cout << "Number of cells in y-direction: ";
    cin >> ny;
    cout << "Number of cells in z-direction: ";
    cin >> nz;
    //and the size of the cells in each direction
    cout << "Cell size x [m]: ";
    cin >> deltax;
    cout << "Cell size y [m]: ";
    cin >> deltay;
    cout << "Cell size z [m]: ";
    cin >> deltaz;

    double factor = 1.0;
    cout << "Increase factor for z: ";
    cin >> factor;
    //set the cell sizes and allocate memory for the mesh
    Model.SetMeshSize(nx, ny, nz);
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(nz, deltaz);
    //the size of each cell in z-direction increases by the specified factor for each layer
    for (int i = 0; i < nz; ++i)
      {
        double thickness = deltaz;
        if (i >= incstart)
          thickness *= pow(factor, i - incstart);
        //x3d has some problems handling thicknesses over 10km with full meter precision
        //so if the thickness is > 10km we round to 100m
        thickness = floor(thickness / rounding) * rounding;
        ZCS[i] = thickness;
      }
    Model.SetZCellSizes(ZCS);
    double originx = 0.0;
    double originy = 0.0;
    double centerz = 0.0;
    jif3D::MTData Data;
    if (vm.count("center"))
      {
        Data.ReadNetCDF(DataName);

        auto mmx = boost::minmax_element(Data.GetMeasPosX().begin(),
            Data.GetMeasPosX().end());
        auto mmy = boost::minmax_element(Data.GetMeasPosY().begin(),
            Data.GetMeasPosY().end());
        auto mmz = boost::minmax_element(Data.GetMeasPosZ().begin(),
            Data.GetMeasPosZ().end());
        double centerx = (*mmx.second + *mmx.first) / 2.0;
        double centery = (*mmy.second + *mmy.first) / 2.0;
        centerz = (*mmz.second + *mmz.first) / 2.0;
        std::cout << "Center X: " << centerx << " Center Y: " << centery << " Center Z:"
            << centerz << std::endl;
        originx = centerx - (nx * deltax / 2.0) + deltax / 2.0;
        originy = centery - (ny * deltay / 2.0) + deltay / 2.0;
        std::cout << "Origin X: " << originx << " Origin  Y: " << originy << " Origin  Z:"
            << centerz << std::endl;
        Model.SetOrigin(originx, originy, centerz);
        if (originx > *mmx.first)
          {
            std::cout << " Warning, minimum model x: " << originx
                << " is larger than smallest X-Coordinate " << *mmx.first << std::endl;
          }
        if (originx + nx * deltax < *mmx.second)
          {
            std::cout << " Warning, maximum model x: " << originx + nx * deltax
                << " is smaller than largest X-Coordinate " << *mmx.second << std::endl;
          }
        if (originy > *mmy.first)
          {
            std::cout << " Warning, minimum model y: " << originy
                << " is larger than smallest y-Coordinate " << *mmy.first << std::endl;
          }
        if (originy + ny * deltay < *mmy.second)
          {
            std::cout << " Warning, maximum model x: " << originy + ny * deltay
                << " is smaller than largest y-Coordinate " << *mmy.second << std::endl;
          }
        Data.WriteMeasurementPoints(DataName + ".vtk");
      }
    //ask for a conductivity to fill the mesh with
    double defaultconductivity = 1.0;
    std::cout << "Conductivity: ";
    std::cin >> defaultconductivity;
    fill_n(Model.SetConductivities().origin(), Model.GetConductivities().num_elements(),
        defaultconductivity);
    //ask for a filename to write the mesh to
    std::string MeshFilename = jif3D::AskFilename("Meshfile name: ", false);

    //fill the background
    std::vector<double> bg_thicknesses(Model.GetZCellSizes().size());
    std::copy(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
        bg_thicknesses.begin());
    std::vector<double> bg_conductivities(Model.GetZCellSizes().size(),
        defaultconductivity);
    Model.SetBackgroundThicknesses(bg_thicknesses);
    Model.SetBackgroundConductivities(bg_conductivities);

    if (vm.count("topofile"))
      {
        std::ifstream topofile(TopoName.c_str());
        std::vector<double> topox, topoy, topoz;
        std::vector<mba::point<2>> points;
        while (topofile.good())
          {
            double currx, curry, currz;
            topofile >> currx >> curry >> currz;
            if (topofile.good())
              {
                topox.push_back(currx);
                topoy.push_back(curry);
                topoz.push_back(currz);
                points.push_back(
                  { currx, curry });
              }
          }
        boost::multi_array<double, 2> CellTopo(boost::extents[nx][ny]);
        // Bounding box containing the data points.
        auto mmx = boost::minmax_element(topox.begin(), topox.end());
        auto mmy = boost::minmax_element(topoy.begin(), topoy.end());
        auto mmz = boost::minmax_element(topoz.begin(), topoz.end());

        mba::point<2> lo =
          { *mmx.first - deltax, *mmy.first - deltay };
        mba::point<2> hi =
          { *mmx.second + deltax, *mmy.second + deltay };

        mba::index<2> grid =
          { 32, 32 };
        double minz = floor(*mmz.first / rounding) * rounding;
        Model.SetOrigin(originx, originy, minz);
        mba::MBA<2> interp(lo, hi, grid, points, topoz);
        jif3D::X3DModel Cov(Model);
        fill_n(Cov.SetConductivities().origin(), Cov.GetConductivities().num_elements(),
            1.0);

        for (size_t i = 0; i < nx; ++i)
          {
            double currposx = originx + (i + 0.5) * deltax;
            for (size_t j = 0; j < ny; ++j)
              {
                double currposy = originy + (j + 0.5) * deltay;
                double itopoz = interp(mba::point<2>
                  { currposx, currposy });
                //std::cout << currposx << " " << currposy << " " << itopoz << std::endl;
                int k = 0;
                double currposz = Model.GetZCoordinates().at(0);
                while (k < nz && currposz <= itopoz)
                  {
                    Model.SetConductivities()[i][j][k] = aircond;
                    Cov.SetConductivities()[i][j][k] = 1e-10;
                    ++k;
                    currposz = Model.GetZCoordinates().at(k);
                  }
                CellTopo[i][j] = Model.GetZCoordinates().at(k);
              }
          }

        std::vector<double> bg_cond(nz);
        for (size_t k = 0; k < nz; ++k)
          {

            using namespace boost::accumulators;
            accumulator_set<double, features<tag::mean, tag::median> > acc;
            for (size_t i = 0; i < nx; ++i)
              for (size_t j = 0; j < ny; ++j)
                {
                  acc(Model.GetConductivities()[i][j][k]);
                }
            bg_cond.at(k) = extract_result<tag::median>(acc);
          }
        Model.SetBackgroundConductivities(bg_cond);
        jif3D::X3DModel TearX(Cov), TearY(Cov), TearZ(Cov);
        for (size_t i = 0; i < nx - 1; ++i)
          for (size_t j = 0; j < ny; ++j)
            for (size_t k = 0; k < nz; ++k)
              TearX.SetConductivities()[i][j][k] = TearZ.GetConductivities()[i][j][k]
                  * TearZ.GetConductivities()[i + 1][j][k];
        for (size_t i = 0; i < nx; ++i)
          for (size_t j = 0; j < ny - 1; ++j)
            for (size_t k = 0; k < nz; ++k)
              TearY.SetConductivities()[i][j][k] = TearZ.GetConductivities()[i][j][k]
                  * TearZ.GetConductivities()[i][j + 1][k];
        Cov.WriteNetCDF(MeshFilename + ".cov.nc");
        Cov.WriteVTK(MeshFilename + ".cov.vtk");
        TearX.WriteNetCDF(MeshFilename + ".tearx.nc");
        TearX.WriteVTK(MeshFilename + ".tearx.vtk");
        TearY.WriteNetCDF(MeshFilename + ".teary.nc");
        TearY.WriteVTK(MeshFilename + ".teary.vtk");
        TearZ.WriteNetCDF(MeshFilename + ".tearz.nc");
        TearZ.WriteVTK(MeshFilename + ".tearz.vtk");
        if (vm.count("center"))
          {
            size_t ndata = Data.GetMeasPosX().size();
            std::vector<double> newz(ndata, 0.0);
            for (size_t i = 0; i < ndata; ++i)
              {
                auto indices = Model.FindAssociatedIndices(Data.GetMeasPosX().at(i),
                    Data.GetMeasPosY().at(i), Data.GetMeasPosZ().at(i));
                newz.at(i) = CellTopo[indices[0]][indices[1]];
              }
            Data.SetMeasurementPoints(Data.GetMeasPosX(), Data.GetMeasPosY(), newz);
            Data.WriteNetCDF(DataName + ".adjusted.nc");
            Data.WriteMeasurementPoints(DataName + ".adjusted.vtk");
            std::string TipperName = DataName + ".tip.nc";
            if (boost::filesystem::exists(TipperName))
              {
                jif3D::TipperData Tipper;
                Tipper.ReadNetCDF(TipperName);
                Tipper.SetMeasurementPoints(Data.GetMeasPosX(), Data.GetMeasPosY(), newz);
                Tipper.WriteNetCDF(TipperName+ ".adjusted.nc");
              }

          }

      }

    Model.WriteNetCDF(MeshFilename);
    Model.WriteVTK(MeshFilename + ".vtk");
    Model.WriteModEM(MeshFilename + ".dat");

    if (vm.count("writegrav"))
      {
        jif3D::ThreeDGravityModel GravModel;
        GravModel.SetMeshSize(nx, ny, nz);
        GravModel.SetXCoordinates(Model.GetXCoordinates());
        GravModel.SetYCoordinates(Model.GetYCoordinates());
        GravModel.SetZCoordinates(Model.GetZCoordinates());
        fill_n(GravModel.SetDensities().origin(), GravModel.SetDensities().num_elements(),
            dens);
        GravModel.WriteNetCDF(MeshFilename + ".grav.nc");
        GravModel.WriteVTK(MeshFilename + ".grav.vtk");
      }

  }
