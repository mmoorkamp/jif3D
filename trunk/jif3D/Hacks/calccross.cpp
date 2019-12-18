//============================================================================
// Name        : calccross.cpp
// Author      : Sep 29, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/ReadAnyModel.h"
#include "../Regularization/CrossGradient.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Inversion/StochasticCovariance.h"

/*! \file calccross.cpp
 * Calculate the cross-gradient constraint for each model cell between two models.
 * For simplicity the two input models have to be seismic models and we do not
 * consider any model covariances.
 */

int main()
  {
    std::string ModelName1 = jif3D::AskFilename("Model file 1: ");
    std::string ModelName2 = jif3D::AskFilename("Model file 2: ");

    boost::shared_ptr<jif3D::ThreeDModelBase> Model1(jif3D::ReadAnyModel(ModelName1));
    boost::shared_ptr<jif3D::ThreeDModelBase> Model2(jif3D::ReadAnyModel(ModelName2));

    const size_t nparam = Model1->GetData().num_elements();

    if (nparam != Model2->GetData().num_elements())
      {
        std::cerr << "Model sizes do not match !";
        return 100;
      }
    jif3D::rvec ModelVec(2 * nparam);
    std::copy(Model1->GetData().origin(), Model1->GetData().origin() + nparam,
        ModelVec.begin());
    std::copy(Model2->GetData().origin(), Model2->GetData().origin() + nparam,
        ModelVec.begin() + nparam);
    bool considersize = true;
    jif3D::CrossGradient CGObjective(*Model1,considersize);

    double cg = CGObjective.CalcMisfit(ModelVec);

    jif3D::rvec DiffVec(nparam);
    for (size_t i = 0; i < nparam; ++i)
      {
        DiffVec(i) = std::pow(CGObjective.GetDataDifference()(i), 2)
            + std::pow(CGObjective.GetDataDifference()(i + nparam), 2)
            + std::pow(CGObjective.GetDataDifference()(i + 2 * nparam), 2);
      }

    std::copy(DiffVec.begin(), DiffVec.end(), Model1->SetData().origin());

    std::cout << "CG: " << cg << std::endl;
    //std::string outfilename = jif3D::AskFilename("Outfile name: ", false);
    //Model1->WriteVTK(outfilename+"cg_mag.vtk","CG Mag");

    jif3D::rvec CG(CGObjective.GetIndividualMisfit());
    const size_t nx = Model1->GetData().shape()[0];
    const size_t ny = Model1->GetData().shape()[1];
    const size_t nz = Model1->GetData().shape()[2];
    const size_t nmod = nx * ny * nz;
    jif3D::ThreeDModelBase::t3DModelData XGrad(boost::extents[nx][ny][nz]);
    jif3D::ThreeDModelBase::t3DModelData YGrad(boost::extents[nx][ny][nz]);
    jif3D::ThreeDModelBase::t3DModelData ZGrad(boost::extents[nx][ny][nz]);
    std::copy(CG.begin(), CG.begin() + nmod, XGrad.origin());
    std::copy(CG.begin() + nmod, CG.begin() + 2 * nmod, YGrad.origin());
    std::copy(CG.begin() + 2 * nmod, CG.begin() + 3 * nmod, ZGrad.origin());
    jif3D::Write3DVectorModelToVTK("crossgrad.vtk", "CrossGrad",
        Model1->GetXCoordinates(), Model1->GetYCoordinates(), Model1->GetZCoordinates(),
        XGrad, YGrad, ZGrad);

    jif3D::rvec CGGrad(CGObjective.CalcGradient(ModelVec));
    jif3D::ThreeDModelBase::t3DModelData GradMod1(boost::extents[nx][ny][nz]);
    jif3D::ThreeDModelBase::t3DModelData GradMod2(boost::extents[nx][ny][nz]);

    std::copy(CGGrad.begin(), CGGrad.begin() + nmod, GradMod1.origin());
    std::copy(CGGrad.begin() + nmod, CGGrad.begin() + 2 * nmod, GradMod2.origin());
    jif3D::Write3DModelToVTK("crossgrad_grad1.vtk", "CrossGrad",
        Model1->GetXCoordinates(), Model1->GetYCoordinates(), Model1->GetZCoordinates(),
        GradMod1);
    jif3D::Write3DModelToVTK("crossgrad_grad2.vtk", "CrossGrad",
        Model1->GetXCoordinates(), Model1->GetYCoordinates(), Model1->GetZCoordinates(),
        GradMod2);

    double a = 2.0;
    double nu = 1.0;
    double sigma = 1.0;

    jif3D::StochasticCovariance Cov(nx, ny, nz, a, nu, sigma);
    jif3D::rvec mCm = Cov.ApplyCovar(ublas::subrange(CGGrad, 0, nmod));
    std::copy(mCm.begin(), mCm.end(), GradMod1.origin());

    jif3D::Write3DModelToVTK("crossgrad_grad1_cov.vtk", "CrossGrad",
        Model1->GetXCoordinates(), Model1->GetYCoordinates(), Model1->GetZCoordinates(),
        GradMod1);
    mCm = Cov.ApplyCovar(ublas::subrange(CGGrad, nmod, 2 * nmod));
    std::copy(mCm.begin(), mCm.end(), GradMod1.origin());
    jif3D::Write3DModelToVTK("crossgrad_grad2_cov.vtk", "CrossGrad",
        Model1->GetXCoordinates(), Model1->GetYCoordinates(), Model1->GetZCoordinates(),
        GradMod1);

    jif3D::rvec m(nx * ny * nz, 0.0);
    //std::generate(m.begin(), m.end(), drand48);
    //std::iota( m.begin(), m.end(),1);
    m((nx * ny * nz) / 2  + (ny * nz) / 2 + nz / 2) = 1.0;
    mCm = Cov.ApplyCovar(m);
    std::copy(mCm.begin(), mCm.end(), GradMod1.origin());
    jif3D::Write3DModelToVTK("covtest.vtk", "Cov",
        Model1->GetXCoordinates(), Model1->GetYCoordinates(), Model1->GetZCoordinates(),
        GradMod1);

  }

