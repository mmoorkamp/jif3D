//============================================================================
// Name        : test_InterpolateFields.cpp
// Author      : Jul 4, 2013
// Version     :
// Copyright   : 2013, mmoorkamp
//============================================================================

#define BOOST_TEST_MODULE OneDMTCalculator test
#define BOOST_TEST_MAIN ...
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>
#include "InterpolateField.h"
#include "../ModelBase/CellBoundaries.h"
#include "../Global/NumUtil.h"
#include "X3DModel.h"

BOOST_AUTO_TEST_SUITE( InterpolateFields_Suite )

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
        Model.SetMeshSize(nx, ny, nz);
        const double deltax = 100, deltay = 70, deltaz = 80;
        Model.SetHorizontalCellSize(deltax, deltay, nx, ny);
        std::fill_n(Model.SetZCellSizes().origin(), nz, deltaz);
        const size_t nmeas = 2;
        for (size_t i = 0; i < nmeas; ++i)
          {
            Model.AddMeasurementPoint(3.0 * deltax / 2, 3.0 * deltay / 2, i * deltaz);
          }
        Model.AddMeasurementPoint(3.0 * deltax / 2, deltay, 0.0);
        Model.AddMeasurementPoint(3.0 * deltax / 2, deltay, 80.0);
        Model.AddMeasurementPoint(deltax, 3.0 * deltay / 2.0, 0.0);
        Model.AddMeasurementPoint(deltax, 3.0 * deltay / 2.0, 80.0);
        std::vector<size_t> MeasDepthIndices;
        std::vector<double> ShiftDepth;
        size_t nlevels = ConstructDepthIndices(MeasDepthIndices, ShiftDepth, Model);

        std::vector<std::complex<double> > Field(nx * ny * nlevels);
        std::generate(Field.begin(), Field.end(), jif3D::IntSequence(0));

        std::complex<double> Real = InterpolateField(Field, Model, 0, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Real.real(), 5);
        BOOST_CHECK_EQUAL(Real.imag(), 0);
        Real = InterpolateField(Field, Model, 1, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Real.real(), 17);
        BOOST_CHECK_EQUAL(Real.imag(), 0);

        Real = InterpolateField(Field, Model, 2, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Real.real(), 4.5);
        BOOST_CHECK_EQUAL(Real.imag(), 0);
        Real = InterpolateField(Field, Model, 3, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Real.real(), 16.5);
        BOOST_CHECK_EQUAL(Real.imag(), 0);

        Real = InterpolateField(Field, Model, 4, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Real.real(), 3.0);
        BOOST_CHECK_EQUAL(Real.imag(), 0);
        Real = InterpolateField(Field, Model, 5, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Real.real(), 15);
        BOOST_CHECK_EQUAL(Real.imag(), 0);

        std::transform(Field.begin(), Field.end(), Field.begin(),
            boost::bind(std::multiplies<std::complex<double> >(), _1,
                std::complex<double>(0.0, 1.0)));
        std::complex<double> Imag = InterpolateField(Field, Model, 0, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Imag.imag(), 5);
        BOOST_CHECK_EQUAL(Imag.real(), 0);
        Imag = InterpolateField(Field, Model, 1, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Imag.imag(), 17);
        BOOST_CHECK_EQUAL(Imag.real(), 0);
        Imag = InterpolateField(Field, Model, 2, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Imag.imag(), 4.5);
        BOOST_CHECK_EQUAL(Imag.real(), 0);
        Imag = InterpolateField(Field, Model, 3, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Imag.imag(), 16.5);
        BOOST_CHECK_EQUAL(Imag.real(), 0);

        Imag = InterpolateField(Field, Model, 4, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Imag.imag(), 3.0);
        BOOST_CHECK_EQUAL(Imag.real(), 0);
        Imag = InterpolateField(Field, Model, 5, MeasDepthIndices);
        BOOST_CHECK_EQUAL(Imag.imag(), 15.0);
        BOOST_CHECK_EQUAL(Imag.real(), 0);
      }
    BOOST_AUTO_TEST_SUITE_END()
