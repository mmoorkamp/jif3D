//============================================================================
// Name        : refinemodel.cpp
// Author      : May 12, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include "../ModelBase/ModelRefiner.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../ModelBase/ReadAnyModel.h"

#include <limits>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

/*! \file refinemodel.cpp
 * Refine the grid of input model. You can specify a new grid size. The old
 * grid size should be an integer multiple of the new grid size in order
 * to ensure that the new model has a uniform grid size as required by
 * the tomography forward calculation object.
 */

int main(int argc, char *argv[])
  {
    bool refinez = true;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("refinez",
        po::value<bool>(&refinez)->default_value(false));

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }

    //read in the name of the original model
    std::string modelfilename = jif3D::AskFilename("Model Filename: ");
    //and then create an object containing the old model
    jif3D::ThreeDModelBase CoarseModel;
    CoarseModel = *jif3D::ReadAnyModel(modelfilename).get();

    //seismic models have to have an equally spaced grid in all directions
    //so we only read in one grid size in m
    double NewDelta = 0.0;
    std::cout << "New grid size [m]: ";
    std::cin >> NewDelta;

    //check that grid size is positive
    if (NewDelta <= 0.0)
      {
        std::cerr << "Negative grid size " << NewDelta << " is not possible !"
            << std::endl;
        return 100;
      }
    //get the extend of the modeling domain in all three directions
    const double oldxmax = CoarseModel.GetXCoordinates().back();
    const double oldymax = CoarseModel.GetYCoordinates().back();
    const double oldzmax = CoarseModel.GetZCoordinates().back();
    //calculate how many new cells we need to fill the modeling domain
    //this can be problematic if NewDelta does not divide the extend of the domain;
    const size_t newnx = oldxmax / NewDelta;
    const size_t newny = oldymax / NewDelta;
    const size_t newnz = oldzmax / NewDelta;
    //so we check that equally spaced cells fill the modeling domain to reasonable precision
    //we are guaranteed that a seismic grid has equal cell sizes in all directions, so
    //we only need to test one direction;
    if (std::fabs(newnx * NewDelta - oldxmax) > std::numeric_limits<float>::epsilon())
      {
        std::cerr << "Refinement coordinates do not equally divide old grid !"
            << std::endl;
        return 100;
      }
    //setup the coordinates for the refined axes
    jif3D::ThreeDModelBase::t3DModelDim XRefiner(newnx);
    for (size_t i = 0; i < newnx; ++i)
      {
        XRefiner[i] = NewDelta * i;
      }

    jif3D::ThreeDModelBase::t3DModelDim YRefiner(newny);
    for (size_t i = 0; i < newny; ++i)
      {
        YRefiner[i] = NewDelta * i;
      }
    jif3D::ThreeDModelBase::t3DModelDim ZRefiner;
    if (refinez)
      {
        ZRefiner.resize(newnz);
        for (size_t i = 0; i < newnz; ++i)
          {
            ZRefiner[i] = NewDelta * i;
          }
      }
    else
      {
        ZRefiner = CoarseModel.GetZCoordinates();
      }
    //now pass this information to the model refiner object
    jif3D::ModelRefiner Refiner;
    Refiner.SetXCoordinates(XRefiner);
    Refiner.SetYCoordinates(YRefiner);
    Refiner.SetZCoordinates(ZRefiner);
    //and create the refined model
    jif3D::ThreeDSeismicModel FineModel;
    Refiner.RefineModel(CoarseModel, FineModel);
    //finally write out the refined model with an appropriate filename;
    FineModel.WriteNetCDF(modelfilename + ".fine.nc");
  }
