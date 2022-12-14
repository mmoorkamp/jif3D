//============================================================================
// Name        : onedmtinv.cpp
// Author      : 13 Jun 2012
// Version     : 
// Copyright   : 2012, mm489
//============================================================================

#include <iostream>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../Global/Noise.h"
#include "../Global/VecMat.h"
#include "../MT/X3DModel.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/OneDMTObjective.h"
#include "../MT/MTEquations.h"
#include "../MT/MTData.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ModelTransforms.h"
#include "../Regularization/OneDRegularization.h"

namespace po = boost::program_options;

void MakeInvCovar(jif3D::OneDMTObjective &MTObjective, const std::vector<double> Imp,
    const std::vector<double> &Frequencies)
  {
    double reserr, phaseerr;
    std::cout << "Relative resistivity error: ";
    std::cin >> reserr;
    std::cout << "Relative phase error: ";
    std::cin >> phaseerr;
    const size_t nImp = Imp.size();
    jif3D::comp_mat InvCov(nImp, nImp);
    for (size_t i = 0; i < nImp - 1; i += 2)
      {
        std::string impfilename = "imp_conf_" + jif3D::stringify(i) + ".out";
        std::string appresfilename = "ap_conf_" + jif3D::stringify(i) + ".out";
        std::ofstream impfile(impfilename);
        std::ofstream appresfile(appresfilename);
        std::complex<double> Z(Imp.at(i), Imp.at(i + 1));
        double freq = Frequencies[i / 2];
        std::cout << "Frequency " << freq << std::endl;
        double appres = jif3D::AppRes(Z, freq);
        double phase = jif3D::ImpedancePhase(Z) / 180.0
            * boost::math::constants::pi<double>();
        double sigma_p2 = std::pow(phase * phaseerr, 2);
        double sigma_r2 = std::pow(appres * reserr, 2);
        std::cout << "Z: " << Z << " " << appres << " " << phase << "\n";
        const double factor = 1.0
            / (freq * jif3D::twopimu * appres * sigma_p2 * sigma_r2);
        InvCov(i, i) =
            factor
                * (pow(sin(phase), 2) * sigma_r2
                    + 4 * pow(cos(phase) * appres, 2) * sigma_p2);
        InvCov(i + 1, i + 1) =
            factor
                * (pow(cos(phase), 2) * sigma_r2
                    + 4 * pow(sin(phase) * appres, 2) * sigma_p2);
        InvCov(i + 1, i) = factor * cos(phase) * sin(phase)
            * (4 * pow(appres, 2) * sigma_p2 - sigma_r2);
        InvCov(i, i + 1) = InvCov(i + 1, i);
        std::cout << "InvCov: " << InvCov(i, i) << " " << InvCov(i + 1, i) << " "
            << InvCov(i + 1, i + 1) << "\n";
        const size_t nsamples = 100;
        double length1 = sqrt(freq * jif3D::twopimu * sigma_r2 / (4.) * appres);
        double length2 = sqrt(jif3D::twopimu * freq * appres * sigma_p2);
        double range = std::max(length1, length2);
        double step = 2 * range / (nsamples);
        std::cout << "Length1: " << length1 << "Length2: " << length2 << "Step: " << step
            << " " << range << "\n";
        for (double zreal = Z.real() - range; zreal < Z.real() + range; zreal += step)
          {
            for (double zimag = Z.imag() - range; zimag < Z.imag() + range; zimag += step)
              {
                double distance2 = ((zreal - Z.real()) * InvCov(i, i)
                    + (zimag - Z.imag()) * InvCov(i, i + 1)) * (zreal - Z.real())
                    + ((zreal - Z.real()) * InvCov(i + 1, i)
                        + (zimag - Z.imag()) * InvCov(i + 1, i + 1)) * (zimag - Z.imag());
                if (distance2 < 1)
                  {
                    impfile << zreal << " " << zimag << "\n";
                    appresfile
                        << jif3D::ImpedancePhase(std::complex<double>(zreal, zimag))
                        << " " << jif3D::AppRes(std::complex<double>(zreal, zimag), freq)
                        << "\n";
                  }
              }
          }
      }
    MTObjective.SetInvCovMat(InvCov);
  }

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
        "The maximum value for conductivity in S/m")("appres",
        "Emulate inversion of apparent resistivity and phase");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    std::string modelfilename = jif3D::AskFilename("Model file: ");
    std::string datafilename = jif3D::AskFilename("Data file: ");
    std::string dataextension = jif3D::GetFileExtension(modelfilename);
    jif3D::X3DModel Model;
    Model.ReadNetCDF(modelfilename);
    jif3D::MTData Data;

    Data.ReadNetCDF(datafilename);
    std::vector<double> Impedances(Data.GetData()), ImpError(Data.GetErrors());
    boost::shared_ptr<jif3D::OneDMTObjective> MTObjective(new jif3D::OneDMTObjective);
    const size_t nfreq = Data.GetFrequencies().size();
    std::vector<double> OneDImp(nfreq * 2, 0.0), OneDErr(nfreq * 2, 1.0);
    std::vector<size_t> nest(nfreq, 0);
    const size_t nsites = Data.GetMeasPosX().size();

    for (size_t i = 0; i < nfreq; ++i)
      {
        for (size_t j = 0; j < nsites; ++j)
          {
            size_t siteindex = (i * nsites + j) * 8;
            double absImp1 = std::sqrt(
                std::pow(Impedances.at(siteindex + 2), 2)
                    + std::pow(Impedances.at(siteindex + 3), 2));
            double absImp2 = std::sqrt(
                std::pow(Impedances.at(siteindex + 4), 2)
                    + std::pow(Impedances.at(siteindex + 5), 2));
            if ((ImpError.at(siteindex + 2) / absImp1 < 0.5)
                && (ImpError.at(siteindex + 4) / absImp2 < 0.5))
              {
                OneDImp.at(i * 2) += Impedances.at(siteindex + 2)
                    - Impedances.at(siteindex + 4);
                OneDImp.at(i * 2 + 1) += Impedances.at(siteindex + 3)
                    - Impedances.at(siteindex + 5);
                ++nest.at(i);
              }
          }
        if (nest.at(i) > 0)
          {
            OneDImp.at(i * 2) /= 2 * nest.at(i);
            OneDImp.at(i * 2 + 1) /= 2 * nest.at(i);
          }
        else
          {
            OneDImp.at(i * 2) = 1.0;
            OneDImp.at(i * 2 + 1) = 1.0;
          }
      }
    jif3D::MTData OneDData;
    OneDData.SetFrequencies(Data.GetFrequencies());
    OneDData.SetMeasurementPoints(
      { 0 },
      { 0 },
      { 0 });
    OneDData.SetDataAndErrors(OneDImp, OneDErr);
    OneDData.CompleteObject();

    std::cout << "Considered: ";
    std::copy(nest.begin(), nest.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout << " out of " << nsites << std::endl;
    const size_t nlayers = Model.GetBackgroundConductivities().size();
    MTObjective->SetModelGeometry(Model);
    MTObjective->SetObservedData(OneDData);
    std::vector<double> AvgImp(OneDImp.size(), 0.0);

    for (size_t i = 0; i < nfreq; ++i)
      {
        AvgImp.at(2 * i) = (OneDImp.at(2 * i) + OneDImp.at(2 * i + 1)) / 2.0;
        AvgImp.at(2 * i + 1) = AvgImp.at(2 * i);
      }
    std::vector<double> DataError(AvgImp.size(), 0.0);
    for (size_t i = 0; i < nfreq; ++i)
      {
        if (nest.at(i) > 0)
          {
            DataError.at(2 * i) = AvgImp.at(2 * i) * 0.05;
            DataError.at(2 * i + 1) = AvgImp.at(2 * i) * 0.05;
          }
        else
          {
            DataError.at(2 * i) = 10.0;
            DataError.at(2 * i + 1) = 10.0;
          }
      }
    if (vm.count("appres"))
      {
        MakeInvCovar(*MTObjective.get(), OneDImp, Data.GetFrequencies());
      }
    else
      {

        MTObjective->SetDataError(DataError);
      }

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
    ConductivityTransform->AppendTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
    ConductivityTransform->AppendTransform(
        boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::LogTransform(RefModel)));
    // ConductivityTransform->AddTransform(
    //    boost::shared_ptr<jif3D::GeneralModelTransform>(new jif3D::ModelCopyTransform));
    InvModel = ConductivityTransform->PhysicalToGeneralized(InvModel);

    boost::shared_ptr<jif3D::ModelCopyTransform> RegTrans(new jif3D::ModelCopyTransform);
    boost::shared_ptr<jif3D::JointObjective> Objective(new jif3D::JointObjective(true));
    Objective->AddObjective(MTObjective, ConductivityTransform, 1.0, "MT");
    Objective->AddObjective(Regularization, RegTrans, reglambda, "Reg");

    jif3D::LimitedMemoryQuasiNewton Optimizer(Objective);
    //Optimizer.SetModelCovDiag(InvModel);
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
    Model.SetBackgroundConductivities(
        std::vector<double>(InvModel.begin(), InvModel.end()));

    for (size_t i = 0; i < Model.GetXCellSizes().size(); ++i)
      for (size_t j = 0; j < Model.GetYCellSizes().size(); ++j)
        for (size_t k = 0; k < Model.GetZCellSizes().size(); ++k)
          Model.SetData()[i][j][k] = InvModel(k);
    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    std::vector<double> InvData(jif3D::OneDMTCalculator().Calculate(Model, Data));
    std::vector<double> FullImp(nfreq * 8, 0.0), FullErr(nfreq * 8, 0.0), Averaged(
        nfreq * 8, 1.0), AvgErr(nfreq * 8, 1.0);
    for (size_t i = 0; i < nfreq; ++i)
      {
        FullImp.at(i * 8 + 2) = InvData.at(i * 2);
        FullImp.at(i * 8 + 3) = InvData.at(i * 2 + 1);
        FullImp.at(i * 8 + 4) = -InvData.at(i * 2);
        FullImp.at(i * 8 + 5) = -InvData.at(i * 2 + 1);
        FullErr.at(i * 8 + 2) = DataError.at(i * 2);
        FullErr.at(i * 8 + 3) = DataError.at(i * 2);
        FullErr.at(i * 8 + 4) = DataError.at(i * 2);
        FullErr.at(i * 8 + 5) = DataError.at(i * 2);
        Averaged.at(i * 8 + 2) = OneDImp.at(i * 2);
        Averaged.at(i * 8 + 3) = OneDImp.at(i * 2 + 1);
        Averaged.at(i * 8 + 4) = -OneDImp.at(i * 2);
        Averaged.at(i * 8 + 5) = -OneDImp.at(i * 2 + 1);
        AvgErr.at(i * 8 + 2) = DataError.at(i * 2);
        AvgErr.at(i * 8 + 3) = DataError.at(i * 2 + 1);
        AvgErr.at(i * 8 + 4) = DataError.at(i * 2);
        AvgErr.at(i * 8 + 5) = DataError.at(i * 2 + 1);
      }
    std::vector<double> Coord =
      { 0.0 };
    std::vector<double> C =
      { 1, 0, 0, 1 };
    std::vector<std::string> Names;
    std::vector<double> Angles(Coord.size(),0.0);
    jif3D::WriteImpedancesToNetCDF(modelfilename + ".inv_imp.nc", Data.GetFrequencies(),
        Coord, Coord, Coord, FullImp, FullErr, C, Names, Angles);

    jif3D::WriteImpedancesToNetCDF(modelfilename + ".avg_imp.nc", Data.GetFrequencies(),
        Coord, Coord, Coord, Averaged, AvgErr, C, Names, Angles);

    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
