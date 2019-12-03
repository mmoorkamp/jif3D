/*
 * mergemodels.cpp
 *
 *  Created on: 19 Nov 2019
 *      Author: moorkamp
 */

#include <vector>
#include <iostream>
#include <string>
#include <boost/algorithm/minmax.hpp>
#include "../MT/X3DModel.h"
#include "../MT/MTData.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"
#include "../Global/NumUtil.h"

double MinDist(double centerx, double centery, const std::vector<double> &XPos,
    const std::vector<double> &YPos)
  {
    const size_t npos = XPos.size();
    double MinDist = 1e44;
    for (size_t i = 0; i < npos; ++i)
      {
        //don't need sqrt because it affects all distances and we only need to know the minimum
        double dist = std::sqrt(
            jif3D::pow2(centerx - XPos.at(i)) + jif3D::pow2(centery - YPos.at(i)));
        if (dist < MinDist)
          MinDist = dist;
      }
    return MinDist;
  }

std::tuple<double, double> FindIndex(double centerx, double centery,
    const std::vector<double> &XCoord, const std::vector<double> &YCoord)
  {
    double xorigin = XCoord.front();
    double yorigin = YCoord.front();
    if (centerx < xorigin || centery < yorigin)
      {
        return
          { -1,-1};
      }
    if (centerx > XCoord.back() || centery > YCoord.back())
      {
        return
          { -2,-2};
      }
    int xindex = 0;
    int yindex = 0;
    while (XCoord.at(xindex + 1) < centerx && xindex < XCoord.size() - 1)
      ++xindex;
    while (YCoord.at(yindex + 1) < centery && yindex < YCoord.size() - 1)
      ++yindex;
    return
      { xindex,yindex};

  }

int main(int argc, char *argv[])
  {
    std::string modellistname = jif3D::AskFilename("Model list:");
    std::string datalistname = jif3D::AskFilename("data list:");

    std::ifstream modellist(modellistname.c_str());
    std::ifstream datalist(datalistname.c_str());

    std::vector<std::string> ModelNames;
    std::vector<std::string> DataNames;
    while (modellist.good())
      {
        std::string tmp;
        modellist >> tmp;
        if (modellist.good())
          {
            ModelNames.push_back(tmp);
          }
      }

    while (datalist.good())
      {
        std::string tmp;
        datalist >> tmp;
        if (datalist.good())
          {
            DataNames.push_back(tmp);
          }
      }

    if (ModelNames.size() != DataNames.size())
      {
        std::cerr << "Not the same number of models and data" << std::endl;
        return 100;
      }
    const size_t nregions = ModelNames.size();
    std::vector<jif3D::X3DModel> Models(nregions);
    std::vector<jif3D::MTData> Data(nregions);

    std::vector<double> MXCoords, MYCoords;

    for (size_t i = 0; i < nregions; ++i)
      {
        Models.at(i).ReadNetCDF(ModelNames.at(i));
        Data.at(i).ReadNetCDF(DataNames.at(i));
        Data.at(i).WriteMeasurementPoints(DataNames.at(i) + ".vtk");

        Models.at(i).WriteVTK(ModelNames.at(i) + ".vtk");
        MXCoords.push_back(Models.at(i).GetXCoordinates().front());
        MXCoords.push_back(Models.at(i).GetXCoordinates().back());
        MYCoords.push_back(Models.at(i).GetYCoordinates().front());
        MYCoords.push_back(Models.at(i).GetYCoordinates().back());
      }
    const double cellsize = Models.front().GetXCellSizes().front();
    auto minmaxx = std::minmax_element(MXCoords.begin(), MXCoords.end());
    auto minmaxy = std::minmax_element(MYCoords.begin(), MYCoords.end());
    size_t nx = (*(minmaxx.second) - *(minmaxx.first)) / cellsize + 1;
    size_t ny = (*(minmaxy.second) - *(minmaxy.first)) / cellsize + 1;
    size_t nz = Models.front().GetZCellSizes().size();
    jif3D::X3DModel MergeModel;
    MergeModel.SetMeshSize(nx, ny, nz);
    MergeModel.SetHorizontalCellSize(cellsize, cellsize, nx, ny);
    MergeModel.SetZCellSizes(Models.front().GetZCellSizes());
    MergeModel.SetOrigin(*(minmaxx.first), *(minmaxy.first), 0.0);
    std::fill_n(MergeModel.SetData().origin(), nx * ny * nz, 0.0);
    jif3D::X3DModel HitModel(MergeModel), MinDistModel(MergeModel);
    std::vector<jif3D::X3DModel> DistModels(Models);

    std::fill_n(MinDistModel.SetData().origin(), nx * ny * nz, 1e40);

    for (size_t i = 0; i < nx; ++i)
      {
        for (size_t j = 0; j < ny; ++j)
          {
            double centerx = MergeModel.GetXCoordinates().at(i)
                + MergeModel.GetXCellSizes().at(i) / 2.0;
            double centery = MergeModel.GetYCoordinates().at(j)
                + MergeModel.GetYCellSizes().at(j) / 2.0;
            std::vector<double> Distances(nregions);
            std::vector<int> XIndices(nregions), YIndices(nregions);
            for (size_t k = 0; k < nregions; ++k)
              {
                Distances.at(k) = MinDist(centerx, centery, Data.at(k).GetMeasPosX(),
                    Data.at(k).GetMeasPosY());
                auto Ind = FindIndex(centerx, centery, Models.at(k).GetXCoordinates(),
                    Models.at(k).GetYCoordinates());
                XIndices.at(k) = std::get < 0 > (Ind);
                YIndices.at(k) = std::get < 1 > (Ind);
                if (XIndices.at(k) >= 0 && YIndices.at(k) >= 0)
                  {
                    for (size_t l = 0; l < nz; ++l)
                      {
                        double dist = std::sqrt(jif3D::pow2(Distances.at(k))
                                                + jif3D::pow2(MergeModel.GetZCoordinates().at(l + 1)));
                        double weight = 1.0 / jif3D::pow2(dist);

                        DistModels.at(k).SetData()[XIndices.at(k)][YIndices.at(k)][l] =
                            weight;
                        MergeModel.SetData()[i][j][l] +=
                            weight
                                * std::log(
                                    Models.at(k).GetData()[XIndices.at(k)][YIndices.at(k)][l]);
                        HitModel.SetData()[i][j][l] += weight;
                        if (dist < MinDistModel.GetData()[i][j][l])
                          MinDistModel.SetData()[i][j][l] = dist;
                      }
                  }
              }
          }
      }

    for (size_t i = 0; i < nx; ++i)
      {
        for (size_t j = 0; j < ny; ++j)
          {
            for (size_t l = 0; l < nz; ++l)
              {
                MergeModel.SetData()[i][j][l] =
                    HitModel.GetData()[i][j][l] > 0 ?
                        std::exp(
                            MergeModel.SetData()[i][j][l] / HitModel.GetData()[i][j][l]) :
                        - 1.0;
              }
          }
      }
    MergeModel.WriteVTK("merged.vtk");
    MergeModel.WriteNetCDF("merge.nc");
    HitModel.WriteVTK("hit.vtk");
    HitModel.WriteNetCDF("hit.nc");
    MinDistModel.WriteVTK("mindist.vtk");
    MinDistModel.WriteNetCDF("mindist.nc");
    for (size_t k = 0; k < nregions; ++k)
      {
        DistModels.at(k).WriteVTK(DataNames.at(k) + ".dist.vtk");
      }
  }
