//============================================================================
// Name        : onedmtinv.cpp
// Author      : 13 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../Global/Noise.h"
#include "../MT/X3DModel.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/OneDMTObjective.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ModelTransforms.h"
#include "../Regularization/OneDRegularization.h"

namespace po = boost::program_options;

int main(int argc, char *argv[])
  {
    double mincond = 1e-6;
    double maxcond = 10;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("covmod", po::value<std::string>(),
        "A file containing the model covariance")("mincond",
        po::value(&mincond)->default_value(1e-6),
        "The minimum value for conductivity in S/m")("maxcond",
        po::value(&maxcond)->default_value(10),
        "The maximum value for conductivity in S/m");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    std::string modelfilename = jiba::AskFilename("Model file: ");
    std::string datafilename = jiba::AskFilename("Data file: ");
    jiba::X3DModel Model;
    Model.ReadNetCDF(modelfilename);
    jiba::rvec Impedances, ImpError;
    std::vector<double> Frequencies, StatXCoord, StatYCoord, StatZCoord;
    jiba::ReadImpedancesFromMTT(datafilename, Frequencies, Impedances, ImpError);

    boost::shared_ptr<jiba::OneDMTObjective> MTObjective(new jiba::OneDMTObjective);
    const size_t nfreq = Frequencies.size();
    jiba::rvec OneDImp(nfreq * 2);
    for (size_t i = 0; i < nfreq; ++i)
      {
        OneDImp(i * 2) = Impedances(i * 8 + 2);
        OneDImp(i * 2 + 1) = Impedances(i * 8 + 3);
      }
    std::cout << "Imp: " << OneDImp << std::endl;
    const size_t nlayers = Model.GetBackgroundConductivities().size();
    Model.SetFrequencies() = Frequencies;
    MTObjective->SetModelGeometry(Model);
    MTObjective->SetObservedData(OneDImp);
    jiba::rvec AvgImp(OneDImp.size());
    for (size_t i =0; i < nfreq; ++i)
      {
        AvgImp(2*i) =(OneDImp(2*i) + OneDImp(2*i+1))/2.0;
        AvgImp(2*i+1) =AvgImp(2*i);
      }
    jiba::rvec DataError(jiba::ConstructError(AvgImp, 0.05, 0.0));
    MTObjective->SetDataError(DataError);

    jiba::rvec InvModel(nlayers);
    std::copy(Model.GetBackgroundConductivities().begin(),
        Model.GetBackgroundConductivities().end(), InvModel.begin());
    boost::shared_ptr<jiba::OneDRegularization> Regularization(
        new jiba::OneDRegularization(nlayers));
    //Regularization->SetReferenceModel(InvModel);
    double reglambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> reglambda;

    boost::shared_ptr<jiba::ChainedTransform> ConductivityTransform(
        new jiba::ChainedTransform);

    jiba::rvec RefModel(InvModel);
    std::fill(RefModel.begin(), RefModel.end(), 1.0);
    //because the tanh transform is used inside a logarithmic transform
    //we need to take the natural logarithm of the actual minimum and maximum
    ConductivityTransform->AppendTransform(
        boost::shared_ptr<jiba::GeneralModelTransform>(
            new jiba::TanhTransform(std::log(mincond), std::log(maxcond))));
    ConductivityTransform->AppendTransform(
        boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::LogTransform(RefModel)));
    // ConductivityTransform->AddTransform(
    //    boost::shared_ptr<jiba::GeneralModelTransform>(new jiba::ModelCopyTransform));
    InvModel = ConductivityTransform->PhysicalToGeneralized(InvModel);

    boost::shared_ptr<jiba::JointObjective> Objective(new jiba::JointObjective(true));
    Objective->AddObjective(MTObjective, ConductivityTransform, 1.0, "MT");
    Objective->AddObjective(Regularization, ConductivityTransform, reglambda, "Reg");

    jiba::LimitedMemoryQuasiNewton Optimizer(Objective);
    Optimizer.SetModelCovDiag(InvModel);
    std::ofstream misfitfile("misfit.out");
    misfitfile << "0 " << Objective->CalcMisfit(InvModel) << " ";

    std::copy(Objective->GetIndividualFits().begin(),
        Objective->GetIndividualFits().end(),
        std::ostream_iterator<double>(misfitfile, " "));
    misfitfile << std::endl;

    size_t iteration = 0;
    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    double maxiter = 100.0;
    std::cout << "Maximum iterations: ";
    std::cin >> maxiter;
    const size_t ndata = OneDImp.size();
    bool terminate = false;
    while (iteration < maxiter && !terminate && Objective->GetIndividualFits()[0] > ndata)
      {
        try
          {
            std::cout << "Iteration: " << iteration + 1 << std::endl;
            std::cout << "InvModel:  " << InvModel << std::endl;
            Optimizer.MakeStep(InvModel);
            std::cout << "Current Misfit: " << Optimizer.GetMisfit() << std::endl;
            std::cout << "Current Gradient: " << Optimizer.GetGradNorm() << std::endl;
            std::cout << std::endl;
            jiba::rvec CondInvModel = ConductivityTransform->GeneralizedToPhysical(
                InvModel);
            std::vector<double> Tmp(CondInvModel.begin(), CondInvModel.end());
            Model.SetBackgroundConductivities(Tmp);

            Model.WriteNetCDF(modelfilename + jiba::stringify(iteration) + ".mt.inv.nc");
            misfitfile << iteration + 1 << " " << Optimizer.GetMisfit() << " ";
            std::copy(Objective->GetIndividualFits().begin(),
                Objective->GetIndividualFits().end(),
                std::ostream_iterator<double>(misfitfile, " "));
            misfitfile << " " << Objective->GetNEval();
            misfitfile << std::endl;
            iteration++;
            terminate = jiba::WantAbort();
          } catch (jiba::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            terminate = true;
          }
      }

    boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;
    InvModel = ConductivityTransform->GeneralizedToPhysical(InvModel);

    std::copy(InvModel.begin(),
        InvModel.begin() + Model.GetConductivities().num_elements(),
        Model.SetConductivities().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec InvData(jiba::OneDMTCalculator().Calculate(Model));
    jiba::rvec FullImp(nfreq * 8, 0.0), FullErr(nfreq * 8, 0.0);
    for (size_t i = 0; i < nfreq; ++i)
      {
        FullImp(i * 8 + 2) = InvData(i * 2);
        FullImp(i * 8 + 3) = InvData(i * 2 + 1);
        FullImp(i * 8 + 4) = -InvData(i * 2);
        FullImp(i * 8 + 5) = -InvData(i * 2 + 1);
        FullErr(i * 8 + 2) = DataError(i * 2);
        FullErr(i * 8 + 3) = DataError(i * 2);
        FullErr(i * 8 + 4) = DataError(i * 2);
        FullErr(i * 8 + 5) = DataError(i * 2);
      }
    jiba::WriteImpedancesToMtt(modelfilename + ".inv_imp", Frequencies, FullImp, FullErr);

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
