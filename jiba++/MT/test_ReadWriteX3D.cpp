//============================================================================
// Name        : test_ReadWriteX3D.cpp
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#define BOOST_TEST_MODULE ReadWriteX3D test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/lambda/lambda.hpp>
#include "ReadWriteX3D.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"

BOOST_AUTO_TEST_SUITE( ReadWriteX3D_Suite )

BOOST_AUTO_TEST_CASE (readX3D_test)
    {

      jiba::ThreeDModelBase::t3DModelDim XCellSizes, YCellSizes, ZCellSizes;
      jiba::ThreeDModelBase::t3DModelData Data;
      const std::string filename("x3d.in");
      jiba::Read3DModelFromX3D(filename, XCellSizes, YCellSizes, ZCellSizes,
          Data);
      jiba::Write3DModelToVTK("x3d.vtk","Conductivity",XCellSizes,YCellSizes,ZCellSizes,Data);
      jiba::Write3DModelToNetCDF("x3d.nc","Conductivity","Siemens",XCellSizes,YCellSizes,ZCellSizes,Data);

    }

BOOST_AUTO_TEST_SUITE_END()
