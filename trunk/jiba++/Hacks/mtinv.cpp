//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <boost/bind.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/GradientRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DObjective.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/ReadWriteImpedances.h"

namespace ublas = boost::numeric::ublas;

int main(int argc, char *argv[])
  {

    //first we read in the starting model and the measured data
    std::string modelfilename, datafilename;
    std::cout << "Starting model Filename: ";
    std::cin >> modelfilename;
    //we read in the starting modelfile
    jiba::X3DModel Model;
    Model.ReadNetCDF(modelfilename);
    Model.SetFrequencies().assign(1, 1.0);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;

    //read in data
    jiba::rvec Data, Frequencies;
    std::vector<double> XCoord, YCoord, ZCoord;
    jiba::ReadImpedancesFromNetCDF(datafilename, Frequencies, XCoord, YCoord,
        ZCoord, Data);
    //if we don't have data inversion doesn't make sense;
    if (Data.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    //we define a few constants that are used throughout the inversion

    const size_t ndata = Data.size();
    const double errorlevel = 0.02;
    const double maxdata = *std::max_element(Data.begin(), Data.end(),
        jiba::absLess<double, double>());
    //create objects for the misfit and a very basic error estimate
    jiba::rvec DataError(ndata);
    for (size_t i = 0; i < ndata; ++i)
      {
        DataError(i) = maxdata
            * errorlevel;
      }
    for (size_t i = 0; i < Model.GetConductivities().shape()[2]; ++i)
      {
        Model.SetConductivities()[0][0][i] *= (1 + 0.005 * (i + 1));
      }
    jiba::rvec InvModel(Model.GetConductivities().num_elements());
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin()
            + Model.GetConductivities().num_elements(), InvModel.begin());
    Model.WriteVTK("start.vtk");
    const double mincond = 5e-2;
    const double maxcond = 5;
    boost::shared_ptr<jiba::GeneralModelTransform> ConductivityTransform(
        new jiba::TanhTransform(mincond, maxcond));

    jiba::rvec RefModel(InvModel);
    InvModel = ConductivityTransform->PhysicalToGeneralized(InvModel);

    boost::shared_ptr<jiba::X3DObjective>
        X3DObjective(new jiba::X3DObjective());
    X3DObjective->SetObservedData(Data);
    X3DObjective->SetModelGeometry(Model);
    X3DObjective->SetDataCovar(DataError);

    boost::shared_ptr<jiba::JointObjective> Objective(
        new jiba::JointObjective());
    boost::shared_ptr<jiba::GradientRegularization> Regularization(
        new jiba::GradientRegularization(Model));

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    Objective->AddObjective(X3DObjective, ConductivityTransform);
    Objective->AddObjective(Regularization, boost::shared_ptr<
        jiba::GeneralModelTransform>(new jiba::ModelCopyTransform()), lambda);

    std::cout << "Performing inversion." << std::endl;

    jiba::LimitedMemoryQuasiNewton LBFGS(Objective, 5);
    std::ofstream misfitfile("misfit.out");
    misfitfile << "0 " << Objective->CalcMisfit(InvModel) << " ";
    std::copy(Objective->GetIndividualFits().begin(),
        Objective->GetIndividualFits().end(), std::ostream_iterator<double>(
            misfitfile, " "));
    misfitfile << std::endl;
    for (size_t i = 0; i < 5; ++i)
      {
        std::cout << "Iteration: " << i + 1 << std::endl;
        LBFGS.MakeStep(InvModel);
        jiba::rvec CondInvModel = ConductivityTransform->GeneralizedToPhysical(
            InvModel);
        std::copy(CondInvModel.begin(), CondInvModel.end(),
            Model.SetConductivities().origin());
        Model.WriteVTK(modelfilename + jiba::stringify(i)
                        + ".mt.inv.vtk");
        std::cout << "Currrent Misfit: " << LBFGS.GetMisfit() << std::endl;
        std::cout << "Currrent Gradient: " << LBFGS.GetGradNorm() << std::endl;
        std::cout << std::endl;

        misfitfile << i + 1 << " " << LBFGS.GetMisfit() << " ";
        std::copy(Objective->GetIndividualFits().begin(),
            Objective->GetIndividualFits().end(),
            std::ostream_iterator<double>(misfitfile, " "));
        misfitfile << " " << Objective->GetNEval();
        misfitfile << std::endl;
      }
    InvModel = ConductivityTransform->GeneralizedToPhysical(InvModel);

    std::copy(InvModel.begin(), InvModel.begin()
        + Model.GetConductivities().num_elements(),
        Model.SetConductivities().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec InvData(jiba::X3DMTCalculator().Calculate(Model));
    jiba::WriteImpedancesToNetCDF(modelfilename + ".inv_imp.nc", Frequencies,
        XCoord, YCoord, ZCoord, InvData);

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
