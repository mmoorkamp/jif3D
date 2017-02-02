//============================================================================
// Name        : CrossGradient.cpp
// Author      : Oct 22, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "CrossGradient.h"

namespace jif3D
  {

    void CrossGradient::ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff)
      {
        const size_t nmod = Model.size();
        //We have two complete models and three components of the cross gradient
        //we calculate the cross-gradient for each cell of the grid
        Diff.resize(nmod / 2 * 3);
        //we need the spatial gradients of both models, so we use the
        //appropriate objective function objects for this
        FirstGradient.CalcMisfit(
            ublas::vector_range<const jif3D::rvec>(Model, ublas::range(0, nmod / 2)));
        SecondGradient.CalcMisfit(
            ublas::vector_range<const jif3D::rvec>(Model, ublas::range(nmod / 2, nmod)));
        //the gradient objective function stores the three directional components
        //as a vector of size 3*m, we create a range for each component
        //that makes the cross gradient calculation more simple below
        ublas::vector_range<const jif3D::rvec> FirstXGrad(
            FirstGradient.GetDataDifference(), ublas::range(0, nmod / 2));
        ublas::vector_range<const jif3D::rvec> FirstYGrad(
            FirstGradient.GetDataDifference(), ublas::range(nmod / 2, nmod));
        ublas::vector_range<const jif3D::rvec> FirstZGrad(
            FirstGradient.GetDataDifference(), ublas::range(nmod, Diff.size()));
        ublas::vector_range<const jif3D::rvec> SecondXGrad(
            SecondGradient.GetDataDifference(), ublas::range(0, nmod / 2));
        ublas::vector_range<const jif3D::rvec> SecondYGrad(
            SecondGradient.GetDataDifference(), ublas::range(nmod / 2, nmod));
        ublas::vector_range<const jif3D::rvec> SecondZGrad(
            SecondGradient.GetDataDifference(), ublas::range(nmod, Diff.size()));
        //now come the actual calculations of the cross-gradient objective
        //function, in the usual form of a difference vector that will be
        //summed and squared in the base class, we have three components
        //of the cross-gradient vector that we store consecutively
        ublas::vector_range<jif3D::rvec>(Diff, ublas::range(0, nmod / 2)) =
            ublas::element_prod(FirstYGrad, SecondZGrad)
                - ublas::element_prod(FirstZGrad, SecondYGrad);
        ublas::vector_range<jif3D::rvec>(Diff, ublas::range(nmod / 2, nmod)) =
            ublas::element_prod(FirstZGrad, SecondXGrad)
                - ublas::element_prod(FirstXGrad, SecondZGrad);
        ublas::vector_range<jif3D::rvec>(Diff, ublas::range(nmod, Diff.size())) =
            ublas::element_prod(FirstXGrad, SecondYGrad)
                - ublas::element_prod(FirstYGrad, SecondXGrad);
      }

    jif3D::rvec CrossGradient::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {
        const size_t nmod = Model.size();
        const size_t halfmod = nmod / 2;
        FirstGradient.CalcGradient(
            ublas::vector_range<const jif3D::rvec>(Model, ublas::range(0, halfmod)));
        SecondGradient.CalcGradient(
            ublas::vector_range<const jif3D::rvec>(Model, ublas::range(halfmod, nmod)));
        jif3D::rvec Gradient(nmod);
        Gradient.clear();
        const size_t xsize = ModelGeometry.GetModelShape()[0];
        const size_t ysize = ModelGeometry.GetModelShape()[1];
        const size_t zsize = ModelGeometry.GetModelShape()[2];
        //first we perform some calculations for all model cells
        //that do not involve explicit differences to neighboring cells
        for (size_t index = 0; index < halfmod; ++index)
          {
            //the gradient of the x-component of the cross-gradient vector
            Gradient(index) += Diff(index)
                * (SecondGradient.GetDataDifference()(index + halfmod)
                    - SecondGradient.GetDataDifference()(index + nmod));
            Gradient(index + halfmod) += Diff(index)
                * (FirstGradient.GetDataDifference()(index + nmod)
                    - FirstGradient.GetDataDifference()(index + halfmod));
            //the gradient of the y-component of the cross-gradient vector
            Gradient(index) += Diff(index + halfmod)
                * (SecondGradient.GetDataDifference()(index + nmod)
                    - SecondGradient.GetDataDifference()(index));
            Gradient(index + halfmod) += Diff(index + halfmod)
                * (FirstGradient.GetDataDifference()(index)
                    - FirstGradient.GetDataDifference()(index + nmod));
            //the gradient of the z-component of the cross-gradient vector
            Gradient(index) += Diff(index + nmod)
                * (SecondGradient.GetDataDifference()(index)
                    - SecondGradient.GetDataDifference()(index + halfmod));
            Gradient(index + halfmod) += Diff(index + nmod)
                * (FirstGradient.GetDataDifference()(index + halfmod)
                    - FirstGradient.GetDataDifference()(index));
          }
        //now there is a always a pair of components that involves
        //the values of neighboring cells, typically for a rotation
        //the y and z components involve the x neighbor etc
        //as we do forward differences this means excluding the first
        //cell in x-direction and is a consequence of the derivatives
        //of eq. 4-6 in Linde et al. 2006
        for (size_t i = 1; i < xsize; ++i)
          {
            for (size_t j = 0; j < ysize; ++j)
              {
                for (size_t k = 0; k < zsize; ++k)
                  {
                    const size_t index = ModelGeometry.IndexToOffset(i, j, k);
                    const size_t frontindex = ModelGeometry.IndexToOffset(i - 1, j, k);
                    Gradient(index) -= Diff(frontindex + halfmod)
                        * SecondGradient.GetDataDifference()(frontindex + nmod);
                    Gradient(index + halfmod) += Diff(frontindex + halfmod)
                        * FirstGradient.GetDataDifference()(frontindex + nmod);

                    Gradient(index) += Diff(frontindex + nmod)
                        * SecondGradient.GetDataDifference()(frontindex + halfmod);
                    Gradient(index + halfmod) -= Diff(frontindex + nmod)
                        * FirstGradient.GetDataDifference()(frontindex + halfmod);
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
                    //for the x-component
                    Gradient(index) += Diff(leftindex)
                        * SecondGradient.GetDataDifference()(leftindex + nmod);
                    Gradient(index + halfmod) -= Diff(leftindex)
                        * FirstGradient.GetDataDifference()(leftindex + nmod);
                    //for the z-component
                    Gradient(index) -= Diff(leftindex + nmod)
                        * SecondGradient.GetDataDifference()(leftindex);
                    Gradient(index + halfmod) += Diff(leftindex + nmod)
                        * FirstGradient.GetDataDifference()(leftindex);
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
                    //this is the x-component of the cross-gradient
                    Gradient(index) -= Diff(topindex)
                        * SecondGradient.GetDataDifference()(topindex + halfmod);
                    Gradient(index + halfmod) += Diff(topindex)
                        * FirstGradient.GetDataDifference()(topindex + halfmod);
                    //this is the y-component of the cross-gradient
                    Gradient(index) += Diff(topindex + halfmod)
                        * SecondGradient.GetDataDifference()(topindex);
                    Gradient(index + halfmod) -= Diff(topindex + halfmod)
                        * FirstGradient.GetDataDifference()(topindex);
                  }
              }
          }

        return 2.0 * Gradient;
      }

  }
