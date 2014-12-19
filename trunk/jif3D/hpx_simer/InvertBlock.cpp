//============================================================================
// Name        : InvertBlock.cpp
// Author      : 19 Dec 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================



#include "InvertBlock.h"
#include "../Global/Noise.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/MTEquations.h"
#include "../ModelTransforms/LogTransform.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/JointObjective.h"

double InvertBlock(jif3D::X3DModel Model, jif3D::rvec Data, std::string tempdir, std::string x3dname)
  {
    jif3D::X3DMTCalculator Calculator(tempdir, x3dname, false);
    boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> > X3DObjective(
        new jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator>(Calculator));


    X3DObjective->SetObservedData(Data);
    jif3D::X3DModel CoarseModel;

    const size_t nz = 2;
    CoarseModel.SetMeshSize(1, 1, nz);
    double maxx = Model.GetXCellSizes()[0] * Model.GetXCoordinates().num_elements();
    double maxy = Model.GetYCellSizes()[0] * Model.GetYCoordinates().num_elements();
    const size_t nzold = Model.GetZCoordinates().num_elements();
    double maxz = Model.GetZCoordinates()[nzold - 1] + Model.GetZCellSizes()[nzold - 1];
    CoarseModel.SetHorizontalCellSize(maxx, maxy, 1, 1);
    CoarseModel.SetZCellSizes()[0] = Model.GetZCellSizes()[0];
    CoarseModel.SetZCellSizes()[1] = maxz - Model.GetZCellSizes()[0];
    CoarseModel.SetBackgroundConductivities(Model.GetBackgroundConductivities());
    CoarseModel.SetBackgroundThicknesses(Model.GetBackgroundThicknesses());
    CoarseModel.CopyMeasurementConfigurations(Model);

    jif3D::rvec InvVec(2, CoarseModel.GetBackgroundConductivities()[0] * 1.01);

    typedef std::complex<double> cd;
    cd Zdet = std::sqrt(
        cd(Data(0), Data(1)) * cd(Data(6), Data(7))
            - cd(Data(2), Data(3)) * cd(Data(4), Data(5)));
    double AppRho = jif3D::AppRes(Zdet, Model.GetFrequencies()[0]);

    InvVec(1) = 1.0 / AppRho;
    std::cout << "Initial guess: " << InvVec(1) << std::endl;
    X3DObjective->SetCoarseModelGeometry(CoarseModel);
    X3DObjective->SetFineModelGeometry(Model);
    jif3D::rvec Error(jif3D::ConstructMTError(Data, 0.001));
    X3DObjective->SetDataError(Error);

    boost::shared_ptr<jif3D::JointObjective> Joint(new jif3D::JointObjective(false));
    jif3D::rvec RefModel(Data.size(), 1.0);
    boost::shared_ptr<jif3D::GeneralModelTransform> ConductivityTransform(
        new jif3D::LogTransform(RefModel));
    jif3D::rvec LogVec = ConductivityTransform->PhysicalToGeneralized(InvVec);
    Joint->AddObjective(X3DObjective, ConductivityTransform);

    double chi = 1e10;
    double cond = 1.0;
    size_t iteration = 0;
    while (chi > 100 && iteration < 5)
      {
        chi = Joint->CalcMisfit(LogVec);
        std::cout << "Chi: " << chi << " " << InvVec(1) << std::endl;
        if (chi > 100)
          {
            jif3D::rvec Synth1(X3DObjective->GetDataDifference());

            double delta = 0.01;
            jif3D::rvec PertVec(LogVec);
            PertVec(1) += delta;
            Joint->CalcMisfit(PertVec);
            jif3D::rvec Synth2(X3DObjective->GetDataDifference());
            jif3D::rvec J = (Synth2 - Synth1) / delta;
            double y = ublas::inner_prod(J, Synth1);

            double JtJ = ublas::inner_prod(J, J);
            cond = LogVec(1) - y / JtJ;
            InvVec(1) = exp(cond);
            LogVec(1) = cond;
            iteration++;
          }
      }
    return InvVec(1);
  }
