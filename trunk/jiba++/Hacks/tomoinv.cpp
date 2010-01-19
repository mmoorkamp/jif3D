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
#include "../Regularization/MinDiffRegularization.h"
#include "../Regularization/GradientRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyObjective.h"
#include "../Tomo/TomographyCalculator.h"

namespace ublas = boost::numeric::ublas;

int main(int argc, char *argv[])
  {
    //these objects hold information about the measurements and their geometry
    jiba::rvec Data;

    //first we read in the starting model and the measured data
    std::string modelfilename, datafilename;
    std::cout << "Starting model Filename: ";
    std::cin >> modelfilename;
    //we read in the starting modelfile
    jiba::ThreeDSeismicModel Model;
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;

    //read in data
    jiba::ReadTraveltimes(datafilename, Data, Model);
    //if we don't have data inversion doesn't make sense;
    if (Data.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    //we define a few constants that are used throughout the inversion

    const size_t ndata = Data.size();
    const double maxdata = *std::max_element(Data.begin(), Data.end(),
        jiba::absLess<double, double>());
    //create objects for the misfit and a very basic error estimate
    jiba::rvec DataError(ndata);
    std::fill_n(DataError.begin(), ndata, 5.0);

    jiba::rvec InvModel(Model.GetSlownesses().num_elements());
    std::copy(Model.GetSlownesses().origin(), Model.GetSlownesses().origin()
        + Model.GetSlownesses().num_elements(), InvModel.begin());

    jiba::rvec RefModel(InvModel);
    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel(i) = log(InvModel(i) / RefModel(i));
    boost::shared_ptr<jiba::TomographyObjective> TomoObjective(
        new jiba::TomographyObjective());
    TomoObjective->SetObservedData(Data);
    TomoObjective->SetFineModelGeometry(Model);
    TomoObjective->SetCoarseModelGeometry(Model);
    TomoObjective->SetDataCovar(DataError);

    boost::shared_ptr<jiba::JointObjective> Objective(
        new jiba::JointObjective());
    boost::shared_ptr<jiba::ObjectiveFunction> Regularization(
        new jiba::GradientRegularization(Model));

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    boost::shared_ptr<jiba::GeneralModelTransform> ModelTransform(
        new jiba::LogTransform(RefModel));
    Objective->AddObjective(TomoObjective, ModelTransform);
    Objective->AddObjective(Regularization, ModelTransform, lambda);

    size_t maxiter = 1;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;
    std::cout << "Performing inversion." << std::endl;

    jiba::LimitedMemoryQuasiNewton LBFGS(Objective, 5);
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

            jiba::rvec TomoModel = ModelTransform->GeneralizedToPhysical(
                InvModel);
            std::copy(TomoModel.begin(), TomoModel.end(),
                Model.SetSlownesses().origin());
            Model.WriteVTK(modelfilename + jiba::stringify(iteration)
                + ".tomo.inv.vtk");

            std::cout << std::endl;
          } catch (jiba::FatalException &e)
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
    jiba::rvec InvData(jiba::TomographyCalculator().Calculate(Model));
    jiba::SaveTraveltimes(modelfilename + ".inv_tt.nc", InvData, Model);
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
