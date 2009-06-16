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

//check that one axis contains the coarse and the fine coordinates
void CheckAxis(const jiba::ThreeDModelBase::t3DModelDim &Refiner,
    const jiba::ThreeDModelBase::t3DModelDim &FineCoordinates,
    const jiba::ThreeDModelBase::t3DModelDim &CoarseCoordinates,
    const jiba::ThreeDModelBase::t3DModelDim &CoarseSizes)
  {
    const size_t ncoarsecells = CoarseCoordinates.size();
    std::vector<double> All;
    std::copy(Refiner.origin(), Refiner.origin() + Refiner.num_elements(),
        std::back_inserter(All));
    std::copy(CoarseCoordinates.begin(), CoarseCoordinates.end(),
        std::back_inserter(All));
    std::sort(All.begin(), All.end());
    All.erase(std::unique(All.begin(), All.end()), All.end());
    All.erase(std::remove_if(All.begin(), All.end(), boost::bind(std::greater<
        double>(), _1, CoarseCoordinates[ncoarsecells - 1]
        + CoarseSizes[ncoarsecells - 1])), All.end());

    BOOST_CHECK(std::equal(All.begin(),All.end(),FineCoordinates.origin()));
  }

void MakeModelandRefinement(jiba::ThreeDGravityModel &CoarseModel,
    jiba::ModelRefiner &Refiner, jiba::ThreeDModelBase::t3DModelDim &XRefiner,
    jiba::ThreeDModelBase::t3DModelDim &YRefiner,
    jiba::ThreeDModelBase::t3DModelDim &ZRefiner)
  {
    size_t ncoarsecellsx = 2;
    size_t ncoarsecellsy = 3;
    size_t ncoarsecellsz = 5;

    CoarseModel.SetXCellSizes().resize(boost::extents[ncoarsecellsx]);
    CoarseModel.SetYCellSizes().resize(boost::extents[ncoarsecellsy]);
    CoarseModel.SetZCellSizes().resize(boost::extents[ncoarsecellsz]);

    std::fill_n(CoarseModel.SetXCellSizes().origin(), ncoarsecellsx, 50.0);
    std::fill_n(CoarseModel.SetYCellSizes().origin(), ncoarsecellsy, 40.0);
    std::fill_n(CoarseModel.SetZCellSizes().origin(), ncoarsecellsz, 30.0);

    const size_t xrefinement = 10;
    const size_t yrefinement = 11;
    const size_t zrefinement = 13;
    XRefiner.resize(boost::extents[xrefinement]);
    YRefiner.resize(boost::extents[yrefinement]);
    ZRefiner.resize(boost::extents[zrefinement]);
    for (size_t i = 0; i < xrefinement; ++i)
      {
        XRefiner[i] = 17 * i;
      }
    for (size_t i = 0; i < yrefinement; ++i)
      {
        YRefiner[i] = 18 * i;
      }
    for (size_t i = 0; i < zrefinement; ++i)
      {
        ZRefiner[i] = 20 * i;
      }

    Refiner.SetXCoordinates(XRefiner);
    Refiner.SetYCoordinates(YRefiner);
    Refiner.SetZCoordinates(ZRefiner);

  }

BOOST_AUTO_TEST_CASE(basic_axis_refinement)
  {
    jiba::ThreeDGravityModel CoarseModel;
    jiba::ThreeDGravityModel FineModel;

    jiba::ModelRefiner Refiner;
    jiba::ThreeDModelBase::t3DModelDim XRefiner;
    jiba::ThreeDModelBase::t3DModelDim YRefiner;
    jiba::ThreeDModelBase::t3DModelDim ZRefiner;
    MakeModelandRefinement(CoarseModel, Refiner, XRefiner, YRefiner, ZRefiner);
    Refiner.RefineAxes(CoarseModel, FineModel);

    CheckAxis(XRefiner, FineModel.GetXCoordinates(),
        CoarseModel.GetXCoordinates(), CoarseModel.GetXCellSizes());
    CheckAxis(YRefiner, FineModel.GetYCoordinates(),
        CoarseModel.GetYCoordinates(), CoarseModel.GetYCellSizes());
    CheckAxis(ZRefiner, FineModel.GetZCoordinates(),
        CoarseModel.GetZCoordinates(), CoarseModel.GetZCellSizes());
    const size_t newnx = FineModel.GetXCellSizes().size();
    const size_t newny = FineModel.GetYCellSizes().size();
    const size_t newnz = FineModel.GetZCellSizes().size();
    const size_t ncoarsecellsx = CoarseModel.GetXCellSizes().size();
    const size_t ncoarsecellsy = CoarseModel.GetYCellSizes().size();
    const size_t ncoarsecellsz = CoarseModel.GetZCellSizes().size();

    BOOST_CHECK(FineModel.GetXCoordinates()[newnx-1] + FineModel.GetXCellSizes()[newnx-1] ==
        CoarseModel.GetXCoordinates()[ncoarsecellsx-1] + CoarseModel.GetXCellSizes()[ncoarsecellsx-1]);
    BOOST_CHECK(FineModel.GetYCoordinates()[newny-1] + FineModel.GetYCellSizes()[newny-1] ==
        CoarseModel.GetYCoordinates()[ncoarsecellsy-1] + CoarseModel.GetYCellSizes()[ncoarsecellsy-1]);
    BOOST_CHECK(FineModel.GetZCoordinates()[newnz-1] + FineModel.GetZCellSizes()[newnz-1] ==
        CoarseModel.GetZCoordinates()[ncoarsecellsz-1] + CoarseModel.GetZCellSizes()[ncoarsecellsz-1]);

    BOOST_CHECK(FineModel.GetDensities().num_elements() == newnx * newny *newnz);
  }

BOOST_AUTO_TEST_CASE(model_projection)
  {
    jiba::ThreeDGravityModel CoarseModel;
    jiba::ThreeDGravityModel FineModel;

    jiba::ModelRefiner Refiner;
    jiba::ThreeDModelBase::t3DModelDim XRefiner;
    jiba::ThreeDModelBase::t3DModelDim YRefiner;
    jiba::ThreeDModelBase::t3DModelDim ZRefiner;
    MakeModelandRefinement(CoarseModel, Refiner, XRefiner, YRefiner, ZRefiner);

    CoarseModel.SetDensities().resize(
        boost::extents[CoarseModel.GetXCellSizes().size()][CoarseModel.GetYCellSizes().size()][CoarseModel.GetZCellSizes().size()]);
    std::generate_n(CoarseModel.SetDensities().origin(),
        CoarseModel.SetDensities().num_elements(), jiba::IntSequence(1.0));
    CoarseModel.WriteVTK("coarse.vtk");

    Refiner.RefineAxes(CoarseModel, FineModel);
    const size_t newnx = FineModel.GetXCellSizes().size();
    const size_t newny = FineModel.GetYCellSizes().size();
    const size_t newnz = FineModel.GetZCellSizes().size();
    BOOST_CHECK(FineModel.GetDensities().num_elements() == newnx * newny *newnz);
    Refiner.ProjectValues(CoarseModel, FineModel);

    FineModel.WriteVTK("fine.vtk");
  }
