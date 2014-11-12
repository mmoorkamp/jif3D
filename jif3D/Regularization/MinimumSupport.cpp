//============================================================================
// Name        : MatOpRegularization.cpp
// Author      : Jan 11, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#include <boost/math/special_functions/pow.hpp>
#include "MinimumSupport.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/convert.h"

namespace jif3D
  {
    namespace bm = boost::math;
    void MinimumSupport::ImplDataDifference(const jif3D::rvec &Model, jif3D::rvec &Diff)
      {
        const size_t nmod = Model.size();
        if (GetReferenceModel().size() != nmod)
          {
            jif3D::rvec RefMod(nmod);
            RefMod.clear();
            SetReferenceModel(RefMod);
          }
        jif3D::rvec X(Model - GetReferenceModel());
        RegFunc->CalcMisfit(X);


        RegDiff = RegFunc->GetDataDifference();
/*        Write3DModelToVTK("regdiffx" + stringify(ModelNumber) + ".vtk", "WaveCoeff",
            RegFunc->ModelGeo.GetXCellSizes(), RegFunc->ModelGeo.GetYCellSizes(),
            RegFunc->ModelGeo.GetZCellSizes(),
            boost::const_multi_array_ref<double, 3>(&RegDiff(0),
                boost::extents[nx][ny][nz]));
        Write3DModelToVTK("regdiffy" + stringify(ModelNumber) + ".vtk", "WaveCoeff",
            RegFunc->ModelGeo.GetXCellSizes(), RegFunc->ModelGeo.GetYCellSizes(),
            RegFunc->ModelGeo.GetZCellSizes(),
            boost::const_multi_array_ref<double, 3>(&RegDiff(nmod),
                boost::extents[nx][ny][nz]));
        Write3DModelToVTK("regdiffz" + stringify(ModelNumber) + ".vtk", "WaveCoeff",
            RegFunc->ModelGeo.GetXCellSizes(), RegFunc->ModelGeo.GetYCellSizes(),
            RegFunc->ModelGeo.GetZCellSizes(),
            boost::const_multi_array_ref<double, 3>(&RegDiff(2*nmod),
                boost::extents[nx][ny][nz]));
        Write3DModelToVTK("x" + stringify(ModelNumber) + ".vtk", "WaveCoeff",
            RegFunc->ModelGeo.GetXCellSizes(), RegFunc->ModelGeo.GetYCellSizes(),
            RegFunc->ModelGeo.GetZCellSizes(),
            boost::const_multi_array_ref<double, 3>(&X(0), boost::extents[nx][ny][nz]));*/
        Diff.resize(3 * nmod, false);
        const double b = beta * beta;

        for (size_t i = 0; i < 3 * nmod; ++i)
          {
            //double mag = bm::pow<2>(RegDiff(i)) + bm::pow<2>(RegDiff(nmod + i))
            //    + bm::pow<2>(RegDiff(2 * nmod + i));
            //Diff(i) = sqrt(Geometry.GetData().data()[i] * mag / (b + mag));
            Diff(i) = RegDiff(i) / sqrt(bm::pow<2>(RegDiff(i)) + b);
          }

/*        Write3DModelToVTK("supp" + stringify(ModelNumber) + ".vtk", "Support",
            RegFunc->ModelGeo.GetXCellSizes(), RegFunc->ModelGeo.GetYCellSizes(),
            RegFunc->ModelGeo.GetZCellSizes(),
            boost::const_multi_array_ref<double, 3>(&Diff(0),
                boost::extents[nx][ny][nz]));*/
        ++ModelNumber;
      }

    jif3D::rvec MinimumSupport::ImplGradient(const jif3D::rvec &Model,
        const jif3D::rvec &Diff)
      {
        const size_t nmod = Model.size();
        jif3D::rvec XGrad(nmod);
        jif3D::rvec YGrad(nmod);
        jif3D::rvec ZGrad(nmod);
        jif3D::rvec Grad(nmod);
        Grad.clear();
        const double b = beta * beta;

        for (size_t i = 0; i < nmod; ++i)
          {
            XGrad(i) = 2 * b * RegDiff(i) / bm::pow<2>(bm::pow<2>(RegDiff(i)) + b);
            YGrad(i) = 2 * b * RegDiff(i + nmod)
                / bm::pow<2>(bm::pow<2>(RegDiff(i + nmod)) + b);
            ZGrad(i) = 2 * b * RegDiff(i + 2 * nmod)
                / bm::pow<2>(bm::pow<2>(RegDiff(i + 2 * nmod)) + b);
          }
        ublas::axpy_prod(ublas::trans(RegFunc->GetXOperator()), XGrad, Grad);
        ublas::axpy_prod(ublas::trans(RegFunc->GetYOperator()), YGrad, Grad, false);
        ublas::axpy_prod(ublas::trans(RegFunc->GetZOperator()), ZGrad, Grad, false);
        jif3D::comp_mat Mat(
            RegFunc->GetXOperator() + RegFunc->GetYOperator() + RegFunc->GetZOperator());
        return Grad;
      }

    MinimumSupport::MinimumSupport(boost::shared_ptr<jif3D::MatOpRegularization> RF,
        double b) :
        ModelNumber(0), beta(b), RegFunc(RF), Geometry(RF->ModelGeo)
      {

        /*
         std::fill_n(Geometry.SetData().origin(), Geometry.SetData().num_elements(), 1.0);
         const size_t nx = RegFunc->ModelGeo.GetXCellSizes().size();
         const size_t ny = RegFunc->ModelGeo.GetYCellSizes().size();
         const size_t nz = RegFunc->ModelGeo.GetZCellSizes().size();
         for (size_t i = 3*nx/4; i < nx; ++i)
         {
         for (size_t j = 3*ny/4; j < ny; ++j)
         {
         for (size_t k = 3*nz/4; k < nz; ++k)
         {
         Geometry.SetData()[i][j][k] *= 100;
         }
         }
         }
         */

      }

    MinimumSupport::~MinimumSupport()
      {
      }

  } /* namespace jif3D */
