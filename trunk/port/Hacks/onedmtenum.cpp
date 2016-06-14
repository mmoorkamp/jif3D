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

    std::string modelfilename = jif3D::AskFilename("Model file: ");
    std::string datafilename = jif3D::AskFilename("Data file: ");
    jif3D::X3DModel Model;
    Model.ReadNetCDF(modelfilename);
    jif3D::rvec Impedances, ImpError;
    std::vector<double> Frequencies, StatXCoord, StatYCoord, StatZCoord;
    jif3D::ReadImpedancesFromMTT(datafilename, Frequencies, Impedances, ImpError);

    boost::shared_ptr<jif3D::OneDMTObjective> MTObjective(new jif3D::OneDMTObjective);
    const size_t nfreq = Frequencies.size();
    jif3D::rvec OneDImp(nfreq * 2);
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
    jif3D::rvec AvgImp(OneDImp.size());
    for (size_t i =0; i < nfreq; ++i)
      {
        AvgImp(2*i) =(OneDImp(2*i) + OneDImp(2*i+1))/2.0;
        AvgImp(2*i+1) =AvgImp(2*i);
      }
    jif3D::rvec DataError(jif3D::ConstructError(AvgImp, 0.05, 0.0));
    MTObjective->SetDataError(DataError);

    jif3D::rvec InvModel(nlayers);
    std::copy(Model.GetBackgroundConductivities().begin(),
        Model.GetBackgroundConductivities().end(), InvModel.begin());
    boost::shared_ptr<jif3D::OneDRegularization> Regularization(
        new jif3D::OneDRegularization(nlayers));
    //Regularization->SetReferenceModel(InvModel);
    double reglambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> reglambda;

    boost::shared_ptr<jif3D::ChainedTransform> ConductivityTransform(
        new jif3D::ChainedTransform);

    jif3D::rvec RefModel(InvModel);
    std::fill(RefModel.begin(), RefModel.end(), 1.0);
    //because the tanh transform is used inside a logarithmic transform
    //we need to take the natural logarithm of the actual minimum and maximum
    ConductivityTransform->AddTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
    ConductivityTransform->AddTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(new jif3D::LogTransform(RefModel)));
    // ConductivityTransform->AddTransform(
    //    boost::shared_ptr<jif3D::GeneralModelTransform>(new jif3D::ModelCopyTransform));
    InvModel = ConductivityTransform->PhysicalToGeneralized(InvModel);

    boost::shared_ptr<jif3D::JointObjective> Objective(new jif3D::JointObjective(true));
    Objective->AddObjective(MTObjective, ConductivityTransform, 1.0, "MT");
    Objective->AddObjective(Regularization, ConductivityTransform, reglambda, "Reg");

    jif3D::LimitedMemoryQuasiNewton Optimizer(Objective);
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
            jif3D::rvec CondInvModel = ConductivityTransform->GeneralizedToPhysical(
                InvModel);
            std::vector<double> Tmp(CondInvModel.begin(), CondInvModel.end());
            Model.SetBackgroundConductivities(Tmp);

            Model.WriteNetCDF(modelfilename + jif3D::stringify(iteration) + ".mt.inv.nc");
            misfitfile << iteration + 1 << " " << Optimizer.GetMisfit() << " ";
            std::copy(Objective->GetIndividualFits().begin(),
                Objective->GetIndividualFits().end(),
                std::ostream_iterator<double>(misfitfile, " "));
            misfitfile << " " << Objective->GetNEval();
            misfitfile << std::endl;
            iteration++;
            terminate = jif3D::WantAbort();
          } catch (jif3D::FatalException &e)
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
    jif3D::rvec InvData(jif3D::OneDMTCalculator().Calculate(Model));
    jif3D::rvec FullImp(nfreq * 8, 0.0), FullErr(nfreq * 8, 0.0);
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
    jif3D::WriteImpedancesToMtt(modelfilename + ".inv_imp", Frequencies, FullImp, FullErr);

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
