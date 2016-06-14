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
void CheckAxis(const jif3D::ThreeDModelBase::t3DModelDim &Refiner,
    const jif3D::ThreeDModelBase::t3DModelDim &FineCoordinates,
    const jif3D::ThreeDModelBase::t3DModelDim &CoarseCoordinates)
  {
    std::vector<double> All;
    std::copy(Refiner.begin(), Refiner.end(), std::back_inserter(All));
    std::copy(CoarseCoordinates.begin(), CoarseCoordinates.end(),
        std::back_inserter(All));
    std::sort(All.begin(), All.end());
    All.erase(std::unique(All.begin(), All.end()), All.end());

    //we have to do the comparison this way round
    //because all can have more elements at the end that are not
    //part of the original modeling domain
    BOOST_CHECK(std::equal(FineCoordinates.begin(), FineCoordinates.end(), All.begin()));
  }

void MakeModelandRefinement(jif3D::ThreeDGravityModel &CoarseModel,
    jif3D::ModelRefiner &Refiner, jif3D::ThreeDModelBase::t3DModelDim &XRefiner,
    jif3D::ThreeDModelBase::t3DModelDim &YRefiner,
    jif3D::ThreeDModelBase::t3DModelDim &ZRefiner)
  {
    size_t ncoarsecellsx = 2;
    size_t ncoarsecellsy = 3;
    size_t ncoarsecellsz = 5;

    CoarseModel.SetXCellSizes().resize(ncoarsecellsx);
    CoarseModel.SetYCellSizes().resize(ncoarsecellsy);
    CoarseModel.SetZCellSizes().resize(ncoarsecellsz);

    std::fill_n(CoarseModel.SetXCellSizes().begin(), ncoarsecellsx, 50.0);
    std::fill_n(CoarseModel.SetYCellSizes().begin(), ncoarsecellsy, 40.0);
    std::fill_n(CoarseModel.SetZCellSizes().begin(), ncoarsecellsz, 30.0);

    const size_t xrefinement = 20;
    const size_t yrefinement = 20;
    const size_t zrefinement = 50;
    XRefiner.resize(xrefinement);
    YRefiner.resize(yrefinement);
    ZRefiner.resize(zrefinement);
    for (size_t i = 0; i < xrefinement; ++i)
      {
        XRefiner[i] = 10 * i;
      }
    for (size_t i = 0; i < yrefinement; ++i)
      {
        YRefiner[i] = 20 * i;
      }
    for (size_t i = 0; i < zrefinement; ++i)
      {
        ZRefiner[i] = 5 * i;
      }

    Refiner.SetXCoordinates(XRefiner);
    Refiner.SetYCoordinates(YRefiner);
    Refiner.SetZCoordinates(ZRefiner);

  }

BOOST_AUTO_TEST_CASE(basic_axis_refinement)
  {
    jif3D::ThreeDGravityModel CoarseModel;
    jif3D::ThreeDGravityModel FineModel;

    jif3D::ModelRefiner Refiner;
    jif3D::ThreeDModelBase::t3DModelDim XRefiner;
    jif3D::ThreeDModelBase::t3DModelDim YRefiner;
    jif3D::ThreeDModelBase::t3DModelDim ZRefiner;
    MakeModelandRefinement(CoarseModel, Refiner, XRefiner, YRefiner, ZRefiner);
    Refiner.RefineAxes(CoarseModel, FineModel);

    CheckAxis(XRefiner, FineModel.GetXCoordinates(), CoarseModel.GetXCoordinates());
    CheckAxis(YRefiner, FineModel.GetYCoordinates(), CoarseModel.GetYCoordinates());
    CheckAxis(ZRefiner, FineModel.GetZCoordinates(), CoarseModel.GetZCoordinates());
    const size_t newnx = FineModel.GetXCellSizes().size();
    const size_t newny = FineModel.GetYCellSizes().size();
    const size_t newnz = FineModel.GetZCellSizes().size();
    const size_t ncoarsecellsx = CoarseModel.GetXCellSizes().size();
    const size_t ncoarsecellsy = CoarseModel.GetYCellSizes().size();
    const size_t ncoarsecellsz = CoarseModel.GetZCellSizes().size();

    BOOST_CHECK(
        FineModel.GetXCoordinates()[newnx - 1] + FineModel.GetXCellSizes()[newnx - 1]
            == CoarseModel.GetXCoordinates()[ncoarsecellsx - 1]
                + CoarseModel.GetXCellSizes()[ncoarsecellsx - 1]);
    BOOST_CHECK(
        FineModel.GetYCoordinates()[newny - 1] + FineModel.GetYCellSizes()[newny - 1]
            == CoarseModel.GetYCoordinates()[ncoarsecellsy - 1]
                + CoarseModel.GetYCellSizes()[ncoarsecellsy - 1]);
    BOOST_CHECK(
        FineModel.GetZCoordinates()[newnz - 1] + FineModel.GetZCellSizes()[newnz - 1]
            == CoarseModel.GetZCoordinates()[ncoarsecellsz - 1]
                + CoarseModel.GetZCellSizes()[ncoarsecellsz - 1]);

    BOOST_CHECK(FineModel.GetDensities().num_elements() == newnx * newny * newnz);
  }

BOOST_AUTO_TEST_CASE(model_projection)
  {
    jif3D::ThreeDGravityModel CoarseModel;
    jif3D::ThreeDGravityModel FineModel;

    jif3D::ModelRefiner Refiner;
    jif3D::ThreeDModelBase::t3DModelDim XRefiner;
    jif3D::ThreeDModelBase::t3DModelDim YRefiner;
    jif3D::ThreeDModelBase::t3DModelDim ZRefiner;
    MakeModelandRefinement(CoarseModel, Refiner, XRefiner, YRefiner, ZRefiner);

    CoarseModel.SetDensities().resize(
        boost::extents[CoarseModel.GetXCellSizes().size()][CoarseModel.GetYCellSizes().size()][CoarseModel.GetZCellSizes().size()]);
    std::iota(CoarseModel.SetDensities().origin(),
        CoarseModel.SetDensities().origin() + CoarseModel.SetDensities().num_elements(),
        0);
    CoarseModel.WriteVTK("coarse.vtk");

    Refiner.RefineAxes(CoarseModel, FineModel);
    const size_t newnx = FineModel.GetXCellSizes().size();
    const size_t newny = FineModel.GetYCellSizes().size();
    const size_t newnz = FineModel.GetZCellSizes().size();
    BOOST_CHECK(FineModel.GetDensities().num_elements() == newnx * newny * newnz);
    Refiner.ProjectValues(CoarseModel, FineModel);
    jif3D::rvec FineGradient(FineModel.GetDensities().num_elements());
    std::fill(FineGradient.begin(), FineGradient.end(), 1.0);
    jif3D::rvec CoarseGradient = Refiner.CombineGradient(FineGradient, CoarseModel,
        FineModel);
    BOOST_CHECK(
        std::count(CoarseGradient.begin(), CoarseGradient.end(), 60)
            == int(CoarseGradient.size()));
    FineModel.WriteVTK("fine.vtk");
  }
