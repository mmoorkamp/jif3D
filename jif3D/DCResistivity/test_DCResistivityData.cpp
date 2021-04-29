//============================================================================
// Name        : test_DCResistivityData.cpp
// Author      : Mar, 2021
// Version     :
// Copyright   : 2021, zhanjie
//============================================================================

#define BOOST_TEST_MODULE DCResistivityData test
#define BOOST_TEST_MAIN ...
#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#endif
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include "DCResistivityData.h"
#include "ThreeDDCResistivityModel.h"
#include "../Global/VecMat.h"

BOOST_AUTO_TEST_SUITE( DCResistivity_Suite )

    BOOST_AUTO_TEST_CASE (read_write_DCResistivityModel_test)
      {
        //create a random number of cells and background layers
        const size_t xsize = rand() % 10 + 2;
        const size_t ysize = rand() % 10 + 2;
        const size_t zsize = rand() % 10 + 2;

        const double xcellsize = drand48() + 1e-3;
        const double ycellsize = drand48() + 1e-3;
        jif3D::ThreeDDCResistivityModel Model;
        Model.SetMeshSize(xsize, ysize, zsize);
        Model.SetHorizontalCellSize(xcellsize, ycellsize, xsize, ysize);
        jif3D::ThreeDModelBase::t3DModelDim ZCS(zsize, 1.0);

        std::generate(ZCS.begin(), ZCS.end(), drand48);
        Model.SetZCellSizes(ZCS);
        std::generate_n(Model.SetData().origin(), Model.SetData().num_elements(),
            drand48);

        Model.WriteNetCDF("dcres_out.nc");

        jif3D::ThreeDDCResistivityModel InModel;
        InModel.ReadNetCDF("dcres_out.nc");

        BOOST_CHECK(Model.GetXCellSizes().size() == InModel.GetXCellSizes().size());
        BOOST_CHECK(Model.GetYCellSizes().size() == InModel.GetYCellSizes().size());
        BOOST_CHECK(Model.GetZCellSizes().size() == InModel.GetZCellSizes().size());
        BOOST_CHECK(Model.GetData().num_elements() == InModel.GetData().num_elements());

        double precision = 0.001;
        for (size_t i = 0; i < Model.GetXCellSizes().size(); ++i)
          {
            BOOST_CHECK_CLOSE(Model.GetXCellSizes()[i], InModel.GetXCellSizes()[i],
                precision);
          }
        for (size_t i = 0; i < Model.GetYCellSizes().size(); ++i)
          {
            BOOST_CHECK_CLOSE(Model.GetYCellSizes()[i], InModel.GetYCellSizes()[i],
                precision);
          }
        for (size_t i = 0; i < Model.GetZCellSizes().size(); ++i)
          {
            BOOST_CHECK_CLOSE(Model.GetZCellSizes()[i], InModel.GetZCellSizes()[i],
                precision);
          }

        for (size_t i = 0; i < Model.GetData().num_elements(); ++i)
          {
            BOOST_CHECK_CLOSE(*(Model.GetData().origin() + i),
                *(InModel.GetData().origin() + i), precision);
          }

      }

    BOOST_AUTO_TEST_CASE (read_write_DCResistivityData_test)
      {
        std::string filename = "dcresdata.nc";

        const size_t nsources = rand() % 10 + 2;
        const size_t ndata = (rand() % 10 + 2) * nsources;

        jif3D::rvec Data(ndata), Error(ndata);
        std::generate(Data.begin(), Data.end(), drand48);
        std::generate(Error.begin(), Error.end(), drand48);
        jif3D::DCResistivityData DcData;

        for (size_t i = 0; i < nsources; ++i)
          {
        	DcData.AddSource(drand48(), drand48(), drand48(), drand48(), drand48(),
                drand48());
          }

        for (size_t i = 0; i < ndata; ++i)
          {
        	DcData.AddAllMeasurementPoint(drand48(), drand48(), drand48(), drand48(),
                drand48(), drand48(), rand() % nsources);
          }

        DcData.SetDataAndErrors(std::vector<double>(Data.begin(),Data.end()),std::vector<double>(Error.begin(),Error.end()));

        DcData.WriteNetCDF(filename);

        jif3D::DCResistivityData InDcData;
        InDcData.ReadNetCDF(filename);

        BOOST_CHECK(Data.size() == DcData.GetData().size());
        BOOST_CHECK(Error.size() == DcData.GetErrors().size());
        double precision = 0.001;
        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Data(i), DcData.GetData()[i], precision);
            BOOST_CHECK_CLOSE(Error(i), DcData.GetErrors()[i], precision);
            BOOST_CHECK_CLOSE(DcData.GetMeasPosX()[i], InDcData.GetMeasPosX()[i],
                precision);
            BOOST_CHECK_CLOSE(DcData.GetMeasPosY()[i], InDcData.GetMeasPosY()[i],
                precision);
            BOOST_CHECK_CLOSE(DcData.GetMeasPosZ()[i], InDcData.GetMeasPosZ()[i],
                precision);

            BOOST_CHECK_CLOSE(DcData.GetMeasSecPosX()[i], InDcData.GetMeasSecPosX()[i],
                precision);
            BOOST_CHECK_CLOSE(DcData.GetMeasSecPosY()[i], InDcData.GetMeasSecPosY()[i],
                precision);
            BOOST_CHECK_CLOSE(DcData.GetMeasSecPosZ()[i], InDcData.GetMeasSecPosZ()[i],
                precision);
            BOOST_CHECK_EQUAL(DcData.GetSourceIndices()[i], InDcData.GetSourceIndices()[i]);
          }

        for (size_t i = 0; i < nsources; ++i)
          {
            BOOST_CHECK_CLOSE(DcData.GetSourcePosPosX()[i], InDcData.GetSourcePosPosX()[i],
                precision);
            BOOST_CHECK_CLOSE(DcData.GetSourcePosPosY()[i], InDcData.GetSourcePosPosY()[i],
                precision);
            BOOST_CHECK_CLOSE(DcData.GetSourcePosPosZ()[i], InDcData.GetSourcePosPosZ()[i],
                precision);

            BOOST_CHECK_CLOSE(DcData.GetSourceNegPosX()[i], InDcData.GetSourceNegPosX()[i],
                precision);
            BOOST_CHECK_CLOSE(DcData.GetSourceNegPosY()[i], InDcData.GetSourceNegPosY()[i],
                precision);
            BOOST_CHECK_CLOSE(DcData.GetSourceNegPosZ()[i], InDcData.GetSourceNegPosZ()[i],
                precision);
          }

      }
    BOOST_AUTO_TEST_SUITE_END()
