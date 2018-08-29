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
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Regularization/GradientRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyCalculator.h"

namespace ublas = boost::numeric::ublas;

int main()
  {
    //these objects hold information about the measurements and their geometry
    jif3D::rvec Data, Errors;

    //first we read in the starting model and the measured data
    std::string modelfilename, datafilename;
    std::cout << "Starting model Filename: ";
    std::cin >> modelfilename;
    //we read in the starting modelfile
    jif3D::ThreeDSeismicModel Model;
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;

    //read in data
    jif3D::ReadTraveltimes(datafilename, Data, Errors, Model);
    //if we don't have data inversion doesn't make sense;
    if (Data.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    //we define a few constants that are used throughout the inversion

    const size_t ndata = Data.size();


    jif3D::rvec InvModel(Model.GetSlownesses().num_elements());
    std::copy(Model.GetSlownesses().origin(), Model.GetSlownesses().origin()
        + Model.GetSlownesses().num_elements(), InvModel.begin());

    jif3D::rvec RefModel(InvModel);
    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel(i) = log(InvModel(i) / RefModel(i));
    jif3D::TomographyCalculator Calculator;
    boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::TomographyCalculator> >
        TomoObjective(
            new jif3D::ThreeDModelObjective<jif3D::TomographyCalculator>(
                Calculator));
    TomoObjective->SetObservedData(Data);
    TomoObjective->SetFineModelGeometry(Model);
    TomoObjective->SetCoarseModelGeometry(Model);
    TomoObjective->SetDataError(Errors);

    boost::shared_ptr<jif3D::JointObjective> Objective(
        new jif3D::JointObjective());
    boost::shared_ptr<jif3D::ObjectiveFunction> Regularization(
        new jif3D::GradientRegularization(Model));

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    boost::shared_ptr<jif3D::GeneralModelTransform> ModelTransform(
        new jif3D::LogTransform(RefModel));
    Objective->AddObjective(TomoObjective, ModelTransform);
    Objective->AddObjective(Regularization, ModelTransform, lambda);

    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;
    std::cout << "Performing inversion." << std::endl;

    jif3D::LimitedMemoryQuasiNewton LBFGS(Objective);
    std::ofstream misfitfile("misfit.out");
    misfitfile << "0 " << Objective->CalcMisfit(InvModel) << " ";
    size_t iteration = 0;

    do
      {

        try
          {
            misfitfile << iteration + 1 << " " << LBFGS.GetMisfit() << " ";
            std::copy(Objective->GetIndividualFits().begin(),
                Objective->GetIndividualFits().end(), std::ostream_iterator<
                    double>(misfitfile, " "));
            misfitfile << " " << Objective->GetNEval();
            misfitfile << std::endl;
            LBFGS.MakeStep(InvModel);

            jif3D::rvec TomoModel = ModelTransform->GeneralizedToPhysical(
                InvModel);
            std::copy(TomoModel.begin(), TomoModel.end(),
                Model.SetSlownesses().origin());
            Model.WriteVTK(modelfilename + jif3D::stringify(iteration)
                + ".tomo.inv.vtk");

            std::cout << std::endl;
          } catch (jif3D::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            iteration = maxiter;
          }
        ++iteration;
      } while (iteration < maxiter && LBFGS.GetMisfit() > ndata
        && LBFGS.GetGradNorm() > 1e-6);
    InvModel = ModelTransform->GeneralizedToPhysical(InvModel);

    std::copy(InvModel.begin(), InvModel.end(), Model.SetSlownesses().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jif3D::rvec InvData(jif3D::TomographyCalculator().Calculate(Model));
    jif3D::SaveTraveltimes(modelfilename + ".inv_tt.nc", InvData, Errors, Model);
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
