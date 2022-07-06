//============================================================================
// Name        : jointinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/hpx_init_params.hpp>
#endif

#ifdef HAVEOPENMP
#include <omp.h>
#endif

#include "../Global/Serialization.h"
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/NonLinearConjugateGradient.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/DiagonalCovariance.h"
#include "../Inversion/StochasticCovariance.h"
#include "../Inversion/MultiSectionCovariance.h"
#include "../ModelBase/ReadAnyModel.h"
#include "../ModelBase/VTKTools.h"
#include "../Regularization/CrossGradient.h"
#include "../MI/MutualInformationConstraint.h"
#include "../Joint/SetupRegularization.h"
#include "../Joint/SetupInversion.h"
#include "../Joint/SetupGravity.h"
#include "../Joint/SetupMagnetics.h"
#include "../Joint/SetupMagnetization.h"
#include "../Joint/SetupMT.h"
#include "../Joint/SetupTomo.h"
#include "../Joint/SetupDCResistivity.h"
#include "../Joint/SetupSW.h"
#include "../Joint/InversionOutput.h"
#include "../Joint/GeneralDataSetup.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <boost/program_options.hpp>

namespace ublas = boost::numeric::ublas;
#ifdef HAVEHPX
#include <hpx/modules/program_options.hpp>
#endif
#ifdef HAVEOPENMP
#include <boost/program_options.hpp>
#endif

#ifdef HAVEHPX
namespace po = hpx::program_options;
#endif
#ifdef HAVEOPENMP
namespace po = boost::program_options;
#endif

double CovWidth = 3.0;
double Cov_nu = 1.0;
double Cov_sigma = 1.0;
jif3D::SetupRegularization RegSetup;
jif3D::SetupInversion InversionSetup;
std::vector<boost::shared_ptr<jif3D::GeneralDataSetup>> DataSetups;
int mibins = 10;

int hpx_main(po::variables_map &vm)
  {

    const size_t nDataSetups = DataSetups.size();

#ifdef HAVEOPENMP
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int>());
      }
#endif

    //some objects accept a directory name as a path to store temporary files
    //we check that this directory actually exists to avoid double checking
    //and early detection of problems
    boost::filesystem::path TempDir = boost::filesystem::current_path();
    if (vm.count("tempdir"))
      {
        TempDir = vm["tempdir"].as<std::string>();
        if (!boost::filesystem::is_directory(TempDir))
          {
            std::cerr << TempDir.string() << " is not a directory or does not exist ! \n";
            return 500;
          }
      }

    std::string meshfilename = jif3D::AskFilename("Mesh filename: ");
    boost::shared_ptr<jif3D::ThreeDModelBase> Mesh = jif3D::ReadAnyModel(meshfilename);

    std::vector<size_t> startindices(
      { 0 });

    std::vector<std::string> SegmentNames;
    std::vector<jif3D::GeneralDataSetup::parametertype> SegmentTypes;

    auto Objective = boost::make_shared<jif3D::JointObjective>(true);
    std::vector<jif3D::rvec> Covs;

    for (size_t i = 0; i < nDataSetups; ++i)
      {
        jif3D::rvec CurrCov;
        DataSetups.at(i)->SetupObjective(vm, *Objective.get(), *Mesh.get(), CurrCov,
            startindices, SegmentNames, SegmentTypes, TempDir);
        Covs.push_back(CurrCov);
      }

    boost::shared_ptr<jif3D::GeneralModelTransform> Copy = boost::make_shared<
        jif3D::ModelCopyTransform>();

    const size_t ninv = startindices.back();
    jif3D::rvec InvModel(ninv, 0.0);
    jif3D::rvec CovModVec(ninv, 1.0);

    boost::shared_ptr<jif3D::ObjectiveFunction> Regularization = RegSetup.SetupObjective(
        vm, *Mesh);

    size_t currsegment = 0;
    for (size_t i = 0; i < nDataSetups; ++i)
      {
        if ((DataSetups.at(i)->GetTransform() != nullptr)
            && (DataSetups.at(i)->GetStartingParameters().size() > 0))
          {
            DataSetups.at(i)->GetTransform()->SetLength(ninv);
            std::copy(Covs.at(i).begin(), Covs.at(i).end(),
                CovModVec.begin() + startindices.at(currsegment));
            jif3D::rvec CurrParm =
                DataSetups.at(i)->GetTransform()->PhysicalToGeneralized(
                    DataSetups.at(i)->GetStartingParameters());
            currsegment++;
            InvModel += CurrParm;
          }
      }

    const size_t ngrid = Mesh->GetNModelElements();
    jif3D::rvec MiValVec(ngrid, 1.0);
    boost::shared_ptr<jif3D::ThreeDModelBase> ValMod;
    if (vm.count("coupling_validity"))
      {
        std::string ValModName = vm["coupling_validity"].as<std::string>();
        ValMod = jif3D::ReadAnyModel(ValModName);
        if (ValMod->GetNModelElements() != ngrid)
          {
            std::cerr << "Size of model to specify coupling validity "
                << ValMod->GetNModelElements() << " does not match inversion grid size: "
                << ngrid << std::endl;
            return 100;
          }
        std::copy(ValMod->GetData().origin(), ValMod->GetData().origin() + ngrid,
            MiValVec.begin());
      }

    boost::shared_ptr<jif3D::ObjectiveFunction> Coupling;
    std::string CouplingName;
    if (vm.count("mutual_information"))
      {
        //The mutual information constraint needs a range over which the histogram is constructed
        //all physical model parameters are transferred to general model parameters to limit variability
        //while the result of the mapping can range from -infinity to +infinity, the vast majority of values
        //will be in the +/-2 range and only values very close to the boundary will be outside
        const double minmi = -2.0;
        const double maxmi = 2.0;
        Coupling = boost::make_shared<jif3D::MutualInformationConstraint>(minmi, maxmi,
            minmi, maxmi, mibins, MiValVec);
        CouplingName = "MI";
      }
    else
      {
        Coupling = boost::make_shared<jif3D::CrossGradient>(*Mesh);
        CouplingName = "CG";
      }

    for (size_t i = 0; i < SegmentNames.size(); ++i)
      {
        for (size_t j = i + 1; j < SegmentNames.size(); ++j)
          {
            if ((SegmentTypes.at(i) == jif3D::GeneralDataSetup::gridparameter)
                && (SegmentTypes.at(j) == jif3D::GeneralDataSetup::gridparameter))
              {
                double couplinglambda = 10.0;
                std::cout << CouplingName << " " << SegmentNames.at(i) << "-"
                    << SegmentNames.at(j) << " weight:";
                std::cin >> couplinglambda;
                if (couplinglambda > 0)
                  {
                    boost::shared_ptr<jif3D::MultiSectionTransform> CouplingTrans =
                        boost::make_shared<jif3D::MultiSectionTransform>(ninv);
                    CouplingTrans->AddSection(startindices.at(i), startindices.at(i + 1),
                        Copy);
                    CouplingTrans->AddSection(startindices.at(j), startindices.at(j + 1),
                        Copy);
                    auto CurrCoupling = boost::shared_ptr<jif3D::ObjectiveFunction>(
                        Coupling->clone());
                    std::string SegmentDisplay = CouplingName + "-" + SegmentNames.at(i)
                        + "-" + SegmentNames.at(j);
                    Objective->AddObjective(CurrCoupling, CouplingTrans, couplinglambda,
                        SegmentDisplay, jif3D::JointObjective::coupling);
                  }
              }
          }
      }

    auto CovObj = boost::make_shared<jif3D::MultiSectionCovariance>(InvModel.size());

    for (size_t i = 0; i < SegmentNames.size(); ++i)
      {
        if (SegmentTypes.at(i) == jif3D::GeneralDataSetup::gridparameter)
          {
            double lambda = 0.0;
            std::cout << SegmentNames.at(i) << "  Regularization Lambda: ";
            std::cin >> lambda;
            if (lambda > 0.0)
              {
                std::string DispName = SegmentNames.at(i) + "Reg";
                boost::shared_ptr<jif3D::GeneralModelTransform> RegTransform =
                    boost::make_shared<jif3D::MultiSectionTransform>(ninv,
                        startindices.at(i), startindices.at(i + 1), Copy);
                auto CurrReg = boost::shared_ptr<jif3D::ObjectiveFunction>(
                    Regularization->clone());
                Objective->AddObjective(CurrReg, RegTransform, lambda, DispName,
                    jif3D::JointObjective::regularization);
              }
            if (CovWidth != 0.0)
              {
                jif3D::rvec CurrCov = ublas::subrange(CovModVec, startindices.at(i),
                    startindices.at(i + 1));
                boost::shared_ptr<jif3D::GeneralCovariance> StochCov = boost::make_shared<
                    jif3D::StochasticCovariance>(CurrCov, Mesh->GetModelShape()[0],
                    Mesh->GetModelShape()[1], Mesh->GetModelShape()[2], CovWidth, Cov_nu,
                    Cov_sigma);
                CovObj->AddSection(startindices.at(i), startindices.at(i + 1), StochCov);
              }
          }
      }
    if (CovWidth == 0)
      {
        double min = *std::min_element(CovModVec.begin(), CovModVec.end());
        if (min <= 0.0)
          {
            std::cerr << "Covariance vector contains non-positive element: " << min
                << " Aborting!" << std::endl;
            return 100;
          }
        CovObj->AddSection(0, InvModel.size(),
            boost::make_shared<jif3D::DiagonalCovariance>(CovModVec));
      }

    boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovObj);

    size_t iteration = 0;
    size_t maxiter = 30;
    std::cout << "Maximum number of iterations: ";
    std::cin >> maxiter;
    std::string modelfilename("result");
    std::ofstream misfitfile("misfit.out");
    std::ofstream rmsfile("rms.out");

    std::cout << "Performing inversion." << std::endl;
    double InitialMisfit = Objective->CalcMisfit(InvModel);
    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);
    StoreRMS(rmsfile, iteration, *Objective);

    bool terminate = false;
    while (iteration < maxiter && !terminate)
      {
        try
          {
            std::cout << "Iteration: " << iteration << std::endl;
            Optimizer->MakeStep(InvModel);
            const std::string currname = modelfilename + jif3D::stringify(iteration);
            for (auto DS : DataSetups)
              DS->IterationOutput(currname, InvModel);

            ++iteration;

            StoreMisfit(misfitfile, iteration, Optimizer->GetMisfit(), *Objective);
            StoreRMS(rmsfile, iteration, *Objective);
            std::cout << "\n\n";
          } catch (jif3D::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            iteration = maxiter;
          }
        //we stop when either we do not make any improvement any more
        terminate = CheckConvergence(*Objective);
        //or the file abort exists in the current directory
        terminate = terminate || jif3D::WantAbort();
      }

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;

    for (auto DS : DataSetups)
      DS->FinalOutput(modelfilename, InvModel);

    if (vm.count("mutual_information"))
      {

      }
    else
      {

        auto ObjectiveTypes = Objective->GetObjectiveTypes();
        std::vector<std::string> Names = Objective->GetObjectiveNames();
        for (size_t i = 0; i < ObjectiveTypes.size(); ++i)
          {
            if (ObjectiveTypes.at(i) == jif3D::JointObjective::coupling)
              {
                jif3D::rvec CG(Objective->GetObjective(i).GetIndividualMisfit());
                const size_t nx = Mesh->GetData().shape()[0];
                const size_t ny = Mesh->GetData().shape()[1];
                const size_t nz = Mesh->GetData().shape()[2];
                const size_t nmod = nx * ny * nz;
                jif3D::ThreeDModelBase::t3DModelData XGrad(boost::extents[nx][ny][nz]);
                jif3D::ThreeDModelBase::t3DModelData YGrad(boost::extents[nx][ny][nz]);
                jif3D::ThreeDModelBase::t3DModelData ZGrad(boost::extents[nx][ny][nz]);
                std::copy(CG.begin(), CG.begin() + nmod, XGrad.origin());
                std::copy(CG.begin() + nmod, CG.begin() + 2 * nmod, YGrad.origin());
                std::copy(CG.begin() + 2 * nmod, CG.begin() + 3 * nmod, ZGrad.origin());
                std::string Name = Names.at(i);
                jif3D::Write3DVectorModelToVTK(Name + ".vtk", Name,
                    Mesh->GetXCoordinates(), Mesh->GetYCoordinates(),
                    Mesh->GetZCoordinates(), XGrad, YGrad, ZGrad);
                jif3D::ThreeDGravityModel CGModel;
                CGModel.SetMeshSize(nx, ny, nz);
                CGModel.SetXCoordinates(Mesh->GetXCoordinates());
                CGModel.SetYCoordinates(Mesh->GetYCoordinates());
                CGModel.SetZCoordinates(Mesh->GetZCoordinates());
                jif3D::ThreeDModelBase::t3DModelData AbsCG(boost::extents[nx][ny][nz]);

                std::transform(CG.begin(), CG.begin() + nmod, AbsCG.origin(),
                    [](double val)
                      { return val*val;});
                std::transform(CG.begin() + nmod, CG.begin() + 2 * nmod, AbsCG.origin(),
                    AbsCG.origin(), [](double val1, double val2)
                      { return val2 + val1*val1;});
                std::transform(CG.begin() + 2 * nmod, CG.begin() + 3 * nmod,
                    AbsCG.origin(), AbsCG.origin(), [](double val1, double val2)
                      { return val2 + val1*val1;});
                std::transform(AbsCG.origin(), AbsCG.origin() + nmod, AbsCG.origin(),
                    [](double val)
                      { return std::sqrt(val);});
                CGModel.SetDensities() = AbsCG;
                CGModel.WriteNetCDF(Name + ".cg_abs.nc");
                CGModel.WriteVTK(Name + ".cg_abs.vtk");
                std::vector<double> CG_Cov(Objective->GetObjective(i).GetDataError());
                std::copy(CG_Cov.begin(), CG_Cov.begin() + nmod, XGrad.origin());
                std::copy(CG_Cov.begin() + nmod, CG_Cov.begin() + 2 * nmod,
                    YGrad.origin());
                std::copy(CG_Cov.begin() + 2 * nmod, CG_Cov.begin() + 3 * nmod,
                    ZGrad.origin());
                jif3D::Write3DVectorModelToVTK(Name + "_cov.vtk", Name,
                    Mesh->GetXCoordinates(), Mesh->GetYCoordinates(),
                    Mesh->GetZCoordinates(), XGrad, YGrad, ZGrad);
              }
          }

      }

#ifdef HAVEHPX
    return hpx::finalize();
#endif
    return 0;

  }

int main(int argc, char *argv[])
  {

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("iterations", po::value<int>(),
        "The maximum number of iterations")("geometrygrid", po::value<std::string>(),
        "A file containing the model geometry of the inversion grid")("threads",
        po::value<int>(), "The number of openmp threads")("mutual_information",
        po::value(&mibins),
        "Use mutual information coupling, specify number of bins in histogram")(
        "stochcov", po::value(&CovWidth)->default_value(0),
        "Width of stochastic regularization, enabled if > 0")("tempdir",
        po::value<std::string>(), "The name of the directory to store temporary files in")(
        "stoch_nu", po::value(&Cov_nu)->default_value(1.00),
        "The value of nu for the stochastic covariance")("stoch_a",
        po::value(&Cov_sigma)->default_value(1.0),
        "The value of a for the stochastic covariance")("coupling_validity",
        po::value<std::string>(),
        "A model file that specifies where the coupling should be applied (1) and where not (0)");

    DataSetups.push_back(boost::make_shared<jif3D::SetupGravity>());
    DataSetups.push_back(boost::make_shared<jif3D::SetupMagnetics>());
    DataSetups.push_back(boost::make_shared<jif3D::SetupMagnetization>());
    DataSetups.push_back(boost::make_shared<jif3D::SetupMT>());
    DataSetups.push_back(boost::make_shared<jif3D::SetupDCResistivity>());
    DataSetups.push_back(boost::make_shared<jif3D::SetupTomo>());
    DataSetups.push_back(boost::make_shared<jif3D::SetupSW>());

    desc.add(RegSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());
    const size_t nDataSetups = DataSetups.size();

    for (size_t i = 0; i < nDataSetups; ++i)
      desc.add(DataSetups.at(i)->SetupOptions());



#ifdef HAVEHPX
    hpx::init_params initparms;
    initparms.desc_cmdline = desc;
    return hpx::init(argc, argv, initparms);
#endif

    //set up the command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }

    return hpx_main(vm);

  }

