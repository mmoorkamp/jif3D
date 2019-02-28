//============================================================================
// Name        : test_ReadWriteDCResistivity.cpp
// Author      : May, 2014
// Version     :
// Copyright   : 2014, zhanjie shi
//============================================================================

#define BOOST_TEST_MODULE ReadWriteDCResistivity test
#define BOOST_TEST_MAIN ...
#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#endif
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "ReadWriteDCResistivityData.h"
#include "ThreeDDCResistivityModel.h"

BOOST_AUTO_TEST_SUITE( ReadWriteDCResistivity_Suite )

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
        jif3D::ThreeDDCResistivityModel Model;

        for (size_t i = 0; i < nsources; ++i)
          {
            Model.AddSource(drand48(), drand48(), drand48(), drand48(), drand48(),
                drand48());
          }

        for (size_t i = 0; i < ndata; ++i)
          {
            Model.AddMeasurementPoint(drand48(), drand48(), drand48(), drand48(),
                drand48(), drand48(), rand() % nsources);
          }

        SaveApparentResistivity(filename, Data, Error, Model);

        jif3D::rvec InData, InError;
        jif3D::ThreeDDCResistivityModel InModel;
        ReadApparentResistivity(filename, InData, InError, InModel);

        BOOST_CHECK(Data.size() == InData.size());
        BOOST_CHECK(Error.size() == InError.size());
        double precision = 0.001;
        for (size_t i = 0; i < ndata; ++i)
          {
            BOOST_CHECK_CLOSE(Data(i), InData(i), precision);
            BOOST_CHECK_CLOSE(Error(i), InError(i), precision);
            BOOST_CHECK_CLOSE(Model.GetMeasPosX()[i], InModel.GetMeasPosX()[i],
                precision);
            BOOST_CHECK_CLOSE(Model.GetMeasPosY()[i], InModel.GetMeasPosY()[i],
                precision);
            BOOST_CHECK_CLOSE(Model.GetMeasPosZ()[i], InModel.GetMeasPosZ()[i],
                precision);

            BOOST_CHECK_CLOSE(Model.GetMeasSecPosX()[i], InModel.GetMeasSecPosX()[i],
                precision);
            BOOST_CHECK_CLOSE(Model.GetMeasSecPosY()[i], InModel.GetMeasSecPosY()[i],
                precision);
            BOOST_CHECK_CLOSE(Model.GetMeasSecPosZ()[i], InModel.GetMeasSecPosZ()[i],
                precision);
            BOOST_CHECK_EQUAL(Model.GetSourceIndices()[i], InModel.GetSourceIndices()[i]);
          }

        for (size_t i = 0; i < nsources; ++i)
          {
            BOOST_CHECK_CLOSE(Model.GetSourcePosPosX()[i], InModel.GetSourcePosPosX()[i],
                precision);
            BOOST_CHECK_CLOSE(Model.GetSourcePosPosY()[i], InModel.GetSourcePosPosY()[i],
                precision);
            BOOST_CHECK_CLOSE(Model.GetSourcePosPosZ()[i], InModel.GetSourcePosPosZ()[i],
                precision);

            BOOST_CHECK_CLOSE(Model.GetSourceNegPosX()[i], InModel.GetSourceNegPosX()[i],
                precision);
            BOOST_CHECK_CLOSE(Model.GetSourceNegPosY()[i], InModel.GetSourceNegPosY()[i],
                precision);
            BOOST_CHECK_CLOSE(Model.GetSourceNegPosZ()[i], InModel.GetSourceNegPosZ()[i],
                precision);
          }

      }
    BOOST_AUTO_TEST_SUITE_END()
