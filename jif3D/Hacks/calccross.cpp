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


/*! \file calccross.cpp
 * Calculate the cross-gradient constraint for each model cell between two models.
 * For simplicity the two input models have to be seismic models and we do not
 * consider any model covariances.
 */

int main()
  {
    std::string ModelName1 = jif3D::AskFilename("Model file 1: ");
    std::string ModelName2 = jif3D::AskFilename("Model file 2: ");



    boost::shared_ptr<jif3D::ThreeDModelBase> Model1(
        jif3D::ReadAnyModel(ModelName1));
    boost::shared_ptr<jif3D::ThreeDModelBase> Model2(
            jif3D::ReadAnyModel(ModelName2));


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

    jif3D::CrossGradient CGObjective(*Model1);

    double cg  = CGObjective.CalcMisfit(ModelVec);

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
    const size_t nx = Model1->GetXCoordinates().size();
    const size_t ny = Model1->GetYCoordinates().size();
    const size_t nz = Model1->GetZCoordinates().size();
    const size_t nmod = nx * ny * nz;
    jif3D::ThreeDModelBase::t3DModelData XGrad(boost::extents[nx][ny][nz]);
    jif3D::ThreeDModelBase::t3DModelData YGrad(boost::extents[nx][ny][nz]);
    jif3D::ThreeDModelBase::t3DModelData ZGrad(boost::extents[nx][ny][nz]);
    std::copy(CG.begin(), CG.begin() + nmod, XGrad.origin());
    std::copy(CG.begin() + nmod, CG.begin() + 2 * nmod, YGrad.origin());
    std::copy(CG.begin() + 2 * nmod, CG.begin() + 3 * nmod, ZGrad.origin());
    jif3D::Write3DVectorModelToVTK("crossgrad.vtk", "CroosGrad", Model1->GetXCellSizes(),
        Model1->GetYCellSizes(), Model1->GetZCellSizes(), XGrad, YGrad, ZGrad);
  }

