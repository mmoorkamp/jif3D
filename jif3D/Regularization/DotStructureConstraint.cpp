/*
 * DotStructureConstraint.cpp
 *
 *  Created on: Feb 5, 2013
 *      Author: mmoorkamp
 */

#include "DotStructureConstraint.h"
#include "../Global/NumUtil.h"
#include <boost/math/special_functions/pow.hpp>
namespace jif3D
  {
    namespace bm = boost::math;
    double GradientElement(const GradientRegularization &GradA,
        const GradientRegularization &GradB, size_t origin, size_t index, size_t nmod)
      {
        const size_t halfmod = nmod / 2;
        double result = GradB.GetDataDifference()(origin)
            * (GradA.GetDataDifference()(index) * GradB.GetDataDifference()(index)
                + GradA.GetDataDifference()(index + halfmod)
                    * GradB.GetDataDifference()(index + halfmod)
                + GradA.GetDataDifference()(index + nmod)
                    * GradB.GetDataDifference()(index + nmod))
            - (bm::pow<2>(GradB.GetDataDifference()(index))
                + bm::pow<2>(GradB.GetDataDifference()(index + halfmod))
                + bm::pow<2>(GradB.GetDataDifference()(index + nmod)))
                * GradA.GetDataDifference()(origin);
        return result;
      }

    void DotStructureConstraint::ImplDataDifference(const jif3D::rvec &Model,
        jif3D::rvec &Diff)
      {
        const size_t nmod = Model.size();
        //We have two complete models and three components of the cross gradient
        //we calculate the cross-gradient for each cell of the grid
        const size_t halfmod = nmod / 2;
        Diff.resize(halfmod);
        //we need the spatial gradients of both models, so we use the
        //appropriate objective function objects for this
        FirstGradient.CalcMisfit(
            ublas::vector_range<const jif3D::rvec>(Model, ublas::range(0, halfmod)));
        SecondGradient.CalcMisfit(
            ublas::vector_range<const jif3D::rvec>(Model, ublas::range(halfmod, nmod)));
        for (size_t i = 0; i < halfmod; ++i)
          {
            Diff(i) = (bm::pow<2>(FirstGradient.GetDataDifference()(i))
                + bm::pow<2>(FirstGradient.GetDataDifference()(i + halfmod))
                + bm::pow<2>(FirstGradient.GetDataDifference()(i + nmod)))
                * (bm::pow<2>(SecondGradient.GetDataDifference()(i))
                    + bm::pow<2>(SecondGradient.GetDataDifference()(i + halfmod))
                    + bm::pow<2>(SecondGradient.GetDataDifference()(i + nmod)))
                - bm::pow<2>(
                    FirstGradient.GetDataDifference()(i)
                        * SecondGradient.GetDataDifference()(i)
                        + FirstGradient.GetDataDifference()(i + halfmod)
                            * SecondGradient.GetDataDifference()(i + halfmod)
                        + FirstGradient.GetDataDifference()(i + nmod)
                            * SecondGradient.GetDataDifference()(i + nmod));
            Diff(i) = sign(Diff(i)) * sqrt(std::abs(Diff(i)));

          }
      }

    jif3D::rvec DotStructureConstraint::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {
        const size_t nmod = Model.size();
        const size_t halfmod = nmod / 2;
        FirstGradient.CalcGradient(
            ublas::vector_range<const jif3D::rvec>(Model, ublas::range(0, halfmod)));
        SecondGradient.CalcGradient(
            ublas::vector_range<const jif3D::rvec>(Model, ublas::range(halfmod, nmod)));
        jif3D::rvec Gradient(nmod, 0.0);
        for (size_t i = 0; i < halfmod; ++i)
          {
            Gradient(i) = (SecondGradient.GetDataDifference()(i)
                + SecondGradient.GetDataDifference()(i + halfmod)
                + SecondGradient.GetDataDifference()(i + nmod))
                * (FirstGradient.GetDataDifference()(i)
                    * SecondGradient.GetDataDifference()(i)
                    + FirstGradient.GetDataDifference()(i + halfmod)
                        * SecondGradient.GetDataDifference()(i + halfmod)
                    + FirstGradient.GetDataDifference()(i + nmod)
                        * SecondGradient.GetDataDifference()(i + nmod))
                - (bm::pow<2>(SecondGradient.GetDataDifference()(i))
                    + bm::pow<2>(SecondGradient.GetDataDifference()(i + halfmod))
                    + bm::pow<2>(SecondGradient.GetDataDifference()(i + nmod)))
                    * (FirstGradient.GetDataDifference()(i)
                        + FirstGradient.GetDataDifference()(i + halfmod)
                        + FirstGradient.GetDataDifference()(i + nmod));

            Gradient(i + halfmod) = (FirstGradient.GetDataDifference()(i)
                + FirstGradient.GetDataDifference()(i + halfmod)
                + FirstGradient.GetDataDifference()(i + nmod))
                * (FirstGradient.GetDataDifference()(i)
                    * SecondGradient.GetDataDifference()(i)
                    + FirstGradient.GetDataDifference()(i + halfmod)
                        * SecondGradient.GetDataDifference()(i + halfmod)
                    + FirstGradient.GetDataDifference()(i + nmod)
                        * SecondGradient.GetDataDifference()(i + nmod))
                - (bm::pow<2>(FirstGradient.GetDataDifference()(i))
                    + bm::pow<2>(FirstGradient.GetDataDifference()(i + halfmod))
                    + bm::pow<2>(FirstGradient.GetDataDifference()(i + nmod)))
                    * (SecondGradient.GetDataDifference()(i)
                        + SecondGradient.GetDataDifference()(i + halfmod)
                        + SecondGradient.GetDataDifference()(i + nmod));

          }
        const size_t xsize = ModelGeometry.GetModelShape()[0];
        const size_t ysize = ModelGeometry.GetModelShape()[1];
        const size_t zsize = ModelGeometry.GetModelShape()[2];
        for (size_t i = 1; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {

                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    const size_t frontindex = ModelGeometry.IndexToOffset(i - 1, j, k);
                    Gradient(index) -= GradientElement(FirstGradient, SecondGradient,
                        frontindex, frontindex, nmod);
                    Gradient(index + halfmod) -= GradientElement(SecondGradient,
                        FirstGradient, frontindex, frontindex, nmod);
                  }
              }
          }
        //now we exclude the first cell in y-direction
        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 1; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    const size_t leftindex = ModelGeometry.IndexToOffset(i, j - 1, k);
                    Gradient(index) -= GradientElement(FirstGradient, SecondGradient,
                        leftindex + halfmod, leftindex, nmod);
                    Gradient(index + halfmod) -= GradientElement(SecondGradient,
                        FirstGradient, leftindex + halfmod, leftindex, nmod);
                  }
              }
          }
        //and now in z-direction
        for (size_t i = 0; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 1; k < zsize; ++k)
                  {
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    const size_t topindex = ModelGeometry.IndexToOffset(i, j, k - 1);
                    Gradient(index) -= GradientElement(FirstGradient, SecondGradient,
                        topindex + nmod, topindex, nmod);
                    Gradient(index + halfmod) -= GradientElement(SecondGradient,
                        FirstGradient, topindex + nmod, topindex, nmod);
                  }
              }
          }

        return 2.0 * Gradient;
      }

    DotStructureConstraint::~DotStructureConstraint()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
