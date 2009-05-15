//============================================================================
// Name        : refinemodel.cpp
// Author      : May 12, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include "../ModelBase/ModelRefiner.h"
#include "ThreeDSeismicModel.h"

int main()
  {
    std::string modelfilename = jiba::AskFilename("Model Filename: ");

    jiba::ThreeDSeismicModel TomoModel;
    TomoModel.ReadNetCDF(modelfilename);

    double NewDelta = 0.0;
    std::cout << "New grid size: ";
    std::cin >> NewDelta;

    const double oldxmax =
        TomoModel.GetXCoordinates()[TomoModel.GetXCoordinates().size() - 1]
            + TomoModel.GetXCellSizes()[TomoModel.GetXCoordinates().size() - 1];

    const double oldymax =
        TomoModel.GetYCoordinates()[TomoModel.GetYCoordinates().size() - 1]
            + TomoModel.GetYCellSizes()[TomoModel.GetYCoordinates().size() - 1];

    const double oldzmax =
        TomoModel.GetZCoordinates()[TomoModel.GetZCoordinates().size() - 1]
            + TomoModel.GetZCellSizes()[TomoModel.GetZCoordinates().size() - 1];

    const size_t newnx = oldxmax / NewDelta;
    const size_t newny = oldymax / NewDelta;
    const size_t newnz = oldzmax / NewDelta;
    jiba::ThreeDModelBase::t3DModelDim XRefiner(boost::extents[newnx]);
    jiba::ThreeDModelBase::t3DModelDim YRefiner(boost::extents[newny]);
    jiba::ThreeDModelBase::t3DModelDim ZRefiner(boost::extents[newnz]);
    for (size_t i = 0; i < newnx; ++i)
      {
        XRefiner[i] = NewDelta * i;
      }
    for (size_t i = 0; i < newny; ++i)
      {
        YRefiner[i] = NewDelta * i;
      }
    for (size_t i = 0; i < newnz; ++i)
      {
        ZRefiner[i] = NewDelta * i;
      }
    jiba::ModelRefiner Refiner;
    Refiner.SetXCoordinates(XRefiner);
    Refiner.SetYCoordinates(YRefiner);
    Refiner.SetZCoordinates(ZRefiner);
    jiba::ThreeDSeismicModel FineModel;
    Refiner.RefineModel(TomoModel, FineModel);

    FineModel.WriteNetCDF(modelfilename+".fine.nc");
  }
