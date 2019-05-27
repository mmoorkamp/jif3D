//============================================================================
// Name        : test_InterpolateFields.cpp
// Author      : Jul 4, 2013
// Version     :
// Copyright   : 2013, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE OneDMTCalculator test
#define BOOST_TEST_MAIN ...
#include "../Global/Jif3DTesting.h"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <cstdlib>
#include "InterpolateField.h"
#include "../ModelBase/CellBoundaries.h"
#include "../Global/NumUtil.h"
#include "../Global/Jif3DPlatformHelper.h"
#include "X3DModel.h"
#include "MTData.h"

BOOST_AUTO_TEST_SUITE (InterpolateFields_Suite)

bool Between(const double limit1, const double limit2, const double value)
  {
    const double upper = std::max(limit1, limit2);
    const double lower = std::min(limit1, limit2);
    return (lower <= value) && (upper >= value);
  }

//std::complex<double> InterpolateField(std::vector<std::complex<double> > &Field,
//   const X3DModel &Model, size_t MeasIndex,
//   const std::vector<size_t> &MeasDepthIndices)

BOOST_AUTO_TEST_CASE (simpleinter_test)
  {
    const size_t nx = 3, ny = 4, nz = 5;
    jif3D::X3DModel Model;
    jif3D::MTData Data;
    Model.SetMeshSize(nx, ny, nz);
    const double deltax = 100, deltay = 70, deltaz = 80;
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(nz, deltaz);
    Model.SetZCellSizes(ZCS);

    const size_t nmeas = 2;
    for (size_t i = 0; i < nmeas; ++i)
      {
        Data.AddMeasurementPoint(3.0 * deltax / 2, 3.0 * deltay / 2, i * deltaz);
      }
    Data.AddMeasurementPoint(3.0 * deltax / 2, deltay, 0.0);
    Data.AddMeasurementPoint(3.0 * deltax / 2, deltay, 80.0);
    Data.AddMeasurementPoint(deltax, 3.0 * deltay / 2.0, 0.0);
    Data.AddMeasurementPoint(deltax, 3.0 * deltay / 2.0, 80.0);
    std::vector<size_t> MeasDepthIndices;
    std::vector<double> ShiftDepth;
    size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model,
        Data.GetMeasPosZ());

    std::vector<std::complex<double> > Field(nx * ny * nlevels);
    std::iota(Field.begin(), Field.end(), 0);

    std::complex<double> Real = InterpolateField(Field, Model, Data, 0, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Real.real(), 5);
    BOOST_CHECK_EQUAL(Real.imag(), 0);
    Real = InterpolateField(Field, Model, Data,1, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Real.real(), 17);
    BOOST_CHECK_EQUAL(Real.imag(), 0);

    Real = InterpolateField(Field, Model, Data,2, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Real.real(), 4.5);
    BOOST_CHECK_EQUAL(Real.imag(), 0);
    Real = InterpolateField(Field, Model, Data,3, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Real.real(), 16.5);
    BOOST_CHECK_EQUAL(Real.imag(), 0);

    Real = InterpolateField(Field, Model, Data,4, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Real.real(), 3.0);
    BOOST_CHECK_EQUAL(Real.imag(), 0);
    Real = InterpolateField(Field, Model, Data,5, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Real.real(), 15);
    BOOST_CHECK_EQUAL(Real.imag(), 0);

    std::transform(Field.begin(), Field.end(), Field.begin(),
        [] (std::complex<double> val)
          { return val * std::complex<double>(0.0, 1.0);});
    std::complex<double> Imag = InterpolateField(Field, Model, Data, 0, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Imag.imag(), 5);
    BOOST_CHECK_EQUAL(Imag.real(), 0);
    Imag = InterpolateField(Field, Model, Data,1, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Imag.imag(), 17);
    BOOST_CHECK_EQUAL(Imag.real(), 0);
    Imag = InterpolateField(Field, Model, Data,2, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Imag.imag(), 4.5);
    BOOST_CHECK_EQUAL(Imag.real(), 0);
    Imag = InterpolateField(Field, Model, Data,3, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Imag.imag(), 16.5);
    BOOST_CHECK_EQUAL(Imag.real(), 0);

    Imag = InterpolateField(Field, Model, Data,4, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Imag.imag(), 3.0);
    BOOST_CHECK_EQUAL(Imag.real(), 0);
    Imag = InterpolateField(Field, Model, Data,5, MeasDepthIndices);
    BOOST_CHECK_EQUAL(Imag.imag(), 15.0);
    BOOST_CHECK_EQUAL(Imag.real(), 0);
  }

BOOST_AUTO_TEST_CASE (functioninter_test)
  {
    jif3D::platform::srand48((int) time(nullptr));
    const size_t nx = 11, ny = 12, nz = 5;
    jif3D::X3DModel Model;
    jif3D::MTData Data;
    Model.SetMeshSize(nx, ny, nz);
    const double deltax = 100, deltay = 70, deltaz = 80;
    Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
    jif3D::ThreeDModelBase::t3DModelDim ZCS(nz, deltaz);
    Model.SetZCellSizes(ZCS);
    const double xcoeff = jif3D::platform::drand48();
    const double ycoeff = jif3D::platform::drand48();

    std::vector<std::complex<double> > Field(nx * ny);
    for (size_t i = 0; i < nx; ++i)
      {
        for (size_t j = 0; j < ny; ++j)
          {

            double xpos = deltax / 2.0 + i * deltax;
            double ypos = deltay / 2.0 + j * deltay;
            Field.at(ny * i + j) = xpos * xcoeff + ypos * ycoeff;
          }
      }
    const size_t ntest = 50;
    for (size_t i = 0; i < ntest; ++i)
      {
        Data.ClearMeasurementPoints();
        double xpos = deltax / 2 + jif3D::platform::drand48() * (deltax * (nx - 1));
        double ypos = deltay / 2 + jif3D::platform::drand48() * (deltay * (ny - 1));
        double zpos = 0;
        Data.AddMeasurementPoint(xpos, ypos, zpos);
        std::vector<size_t> MeasDepthIndices;
        std::vector<double> ShiftDepth;
        size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model,Data.GetMeasPosZ());
        std::complex<double> value = InterpolateField(Field, Model, Data,0,
            MeasDepthIndices);
        double trueval = xpos * xcoeff + ypos * ycoeff;
        BOOST_CHECK_CLOSE(value.real(), trueval, 0.01);
      }
  }
BOOST_AUTO_TEST_SUITE_END()
