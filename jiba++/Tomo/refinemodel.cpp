//============================================================================
// Name        : refinemodel.cpp
// Author      : May 12, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include "../ModelBase/ModelRefiner.h"
#include "ThreeDSeismicModel.h"
#include <limits>

/*! \file refinemodel.cpp
 * Refine the grid of input model. You can specify a new grid size. The old
 * grid size should be an integer multiple of the new grid size in order
 * to ensure that the new model has a uniform grid size as required by
 * the tomography forward calculation object.
 */

int main()
  {
    //read in the name of the original model
    std::string modelfilename = jiba::AskFilename("Model Filename: ");
    //and then create an object containing the old model
    jiba::ThreeDSeismicModel TomoModel;
    TomoModel.ReadNetCDF(modelfilename);

    //seismic models have to have an equally spaced grid in all directions
    //so we only read in one grid size in m
    double NewDelta = 0.0;
    std::cout << "New grid size [m]: ";
    std::cin >> NewDelta;

    //check that grid size is positive
    if (NewDelta <= 0.0)
      {
        std::cerr << "Negative grid size " << NewDelta << " is not possible !" << std::endl;
        return 100;
      }
    //get the extend of the modeling domain in all three directions
    const double oldxmax =
        TomoModel.GetXCoordinates()[TomoModel.GetXCoordinates().size() - 1]
            + TomoModel.GetXCellSizes()[TomoModel.GetXCoordinates().size() - 1];

    const double oldymax =
        TomoModel.GetYCoordinates()[TomoModel.GetYCoordinates().size() - 1]
            + TomoModel.GetYCellSizes()[TomoModel.GetYCoordinates().size() - 1];

    const double oldzmax =
        TomoModel.GetZCoordinates()[TomoModel.GetZCoordinates().size() - 1]
            + TomoModel.GetZCellSizes()[TomoModel.GetZCoordinates().size() - 1];
    //calculate how many new cells we need to fill the modeling domain
    //this can be problematic if NewDelta does not divide the extend of the domain;
    const size_t newnx = oldxmax / NewDelta;
    const size_t newny = oldymax / NewDelta;
    const size_t newnz = oldzmax / NewDelta;
    //so we check that equally spaced cells fill the modeling domain to reasonable precision
    //we are guaranteed that a seismic grid has equal cell sizes in all directions, so
    //we only need to test one direction;
    if (std::fabs(newnx * NewDelta - oldxmax) > std::numeric_limits<float>::epsilon() )
      {
        std::cerr << "Refinement coordinates do not equally divide old grid !" << std::endl;
        return 100;
      }
    //setup the coordinates for the refined axes
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
    //now pass this information to the model refiner object
    jiba::ModelRefiner Refiner;
    Refiner.SetXCoordinates(XRefiner);
    Refiner.SetYCoordinates(YRefiner);
    Refiner.SetZCoordinates(ZRefiner);
    //and create the refined model
    jiba::ThreeDSeismicModel FineModel;
    Refiner.RefineModel(TomoModel, FineModel);
    //finally write out the refined model with an appropriate filename;
    FineModel.WriteNetCDF(modelfilename+".fine.nc");
  }
