//============================================================================
// Name        : test_ModelRefiner.cpp
// Author      : Feb 20, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE ModelRefiner test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <algorithm>
#include "ModelRefiner.h"
#include "ThreeDModelBase.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Global/NumUtil.h"

BOOST_AUTO_TEST_CASE(basic_axis_refinement)
  {
    jiba::ThreeDGravityModel CoarseModel;
    jiba::ThreeDGravityModel FineModel;
    size_t ncoarsecells = 2;
    CoarseModel.SetXCellSizes().resize(boost::extents[ncoarsecells]);
    CoarseModel.SetYCellSizes().resize(boost::extents[ncoarsecells]);
    CoarseModel.SetZCellSizes().resize(boost::extents[ncoarsecells]);
    for (size_t i = 0; i < ncoarsecells; ++i)
      {
        CoarseModel.SetXCellSizes()[i] = 50;
        CoarseModel.SetYCellSizes()[i] = 60;
        CoarseModel.SetZCellSizes()[i] = 70;
      }

    const size_t xrefinement = 10;
    jiba::ThreeDModelBase::t3DModelDim XRefiner(boost::extents[xrefinement]);
    for (size_t i = 0; i < xrefinement; ++i)
      {
        XRefiner[i] = 17 * i;
      }
    jiba::ModelRefiner Refiner;
    Refiner.SetXCoordinates(XRefiner);
    Refiner.RefineAxes(CoarseModel, FineModel);
    std::cout << "X-Axis: ";
    std::copy(FineModel.GetXCellSizes().origin(),
        FineModel.GetXCellSizes().origin() + FineModel.GetXCellSizes().size(),
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Y-Axis: ";
    std::copy(FineModel.GetYCellSizes().origin(),
        FineModel.GetYCellSizes().origin() + FineModel.GetYCellSizes().size(),
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "Z-Axis: ";
    std::copy(FineModel.GetZCellSizes().origin(),
        FineModel.GetZCellSizes().origin() + FineModel.GetZCellSizes().size(),
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    const size_t newnx = FineModel.GetXCellSizes().size();
    const size_t newny = FineModel.GetYCellSizes().size();
    const size_t newnz = FineModel.GetZCellSizes().size();
    BOOST_CHECK(FineModel.GetDensities().num_elements() == newnx * newny *newnz);
  }

BOOST_AUTO_TEST_CASE(model_projection)
  {
    jiba::ThreeDGravityModel CoarseModel;
    jiba::ThreeDGravityModel FineModel;
    size_t ncoarsecells = 2;
    CoarseModel.SetXCellSizes().resize(boost::extents[ncoarsecells]);
    CoarseModel.SetYCellSizes().resize(boost::extents[ncoarsecells]);
    CoarseModel.SetZCellSizes().resize(boost::extents[ncoarsecells]);
    for (size_t i = 0; i < ncoarsecells; ++i)
      {
        CoarseModel.SetXCellSizes()[i] = 50;
        CoarseModel.SetYCellSizes()[i] = 60;
        CoarseModel.SetZCellSizes()[i] = 70;
      }
    CoarseModel.SetDensities().resize(boost::extents[ncoarsecells][ncoarsecells][ncoarsecells]);
    std::generate_n(CoarseModel.SetDensities().origin(),
        CoarseModel.SetDensities().num_elements(), jiba::IntSequence(1.0));
    CoarseModel.WriteVTK("coarse.vtk");
    const size_t xrefinement = 10;
    jiba::ThreeDModelBase::t3DModelDim XRefiner(boost::extents[xrefinement]);
    jiba::ThreeDModelBase::t3DModelDim YRefiner(boost::extents[xrefinement]);
    jiba::ThreeDModelBase::t3DModelDim ZRefiner(boost::extents[xrefinement]);
    for (size_t i = 0; i < xrefinement; ++i)
      {
        XRefiner[i] = 17 * i;
        YRefiner[i] = 13 * i;
        ZRefiner[i] = 19 * i;
      }
    jiba::ModelRefiner Refiner;
    Refiner.SetXCoordinates(XRefiner);
    Refiner.SetYCoordinates(YRefiner);
    Refiner.SetZCoordinates(ZRefiner);
    Refiner.RefineAxes(CoarseModel, FineModel);
    const size_t newnx = FineModel.GetXCellSizes().size();
    const size_t newny = FineModel.GetYCellSizes().size();
    const size_t newnz = FineModel.GetZCellSizes().size();
    BOOST_CHECK(FineModel.GetDensities().num_elements() == newnx * newny *newnz);
    Refiner.ProjectValues(CoarseModel,FineModel);

    FineModel.WriteVTK("fine.vtk");
  }
