//============================================================================
// Name        : mtinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================
#ifdef HAVEHPX
#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#endif

#ifdef HAVEOPENMP
#include <omp.h>
#endif

#include <boost/program_options.hpp>
#include <boost/program_options/config.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/make_shared.hpp>

#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../Global/ReadWriteSparseMatrix.h"
#include "../ModelBase/ReadAnyModel.h"
#include "../ModelBase/VTKTools.h"
#include "../Regularization/GradientRegularization.h"
#include "../Regularization/CrossGradient.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/DiagonalCovariance.h"
#include "../Inversion/StochasticCovariance.h"
#include "../Inversion/MultiSectionCovariance.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/X3DTipperCalculator.h"
#include "../MT/X3DFieldCalculator.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTTransforms.h"
#include "MTData.h"
#include "TipperData.h"
#include "../Joint/SetupRegularization.h"
#include "../Joint/SetupInversion.h"
#include "../Joint/InversionOutput.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;

jif3D::SetupRegularization RegSetup;
jif3D::SetupInversion InversionSetup;
double mincond = 1e-6;
double maxcond = 10;
double relerr = 0.02;
double coolingfactor = 1.0;
double xorigin = 0.0;
double yorigin = 0.0;
double conddelta = 0.001;
std::string X3DName = "x3d";
std::string MTInvCovarName;
std::string RefModelName;
std::string CrossModelName;
std::string TipperName;
double TipperWeight;
double MTWeight;
double tipperr, tiprelerr;
double DistCorr = 0;
bool CleanFiles = true;
double CovWidth = 3.0;
double zshift = 0.0;

int hpx_main(boost::program_options::variables_map& vm)
  {

    boost::shared_ptr<jif3D::JointObjective> Objective(new jif3D::JointObjective(true));

#ifdef HAVEOPENMP
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int>());
      }
#endif

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

    if (!boost::filesystem::exists(X3DName) && !vm.count("regcheck"))
      {
        std::cerr << X3DName << " is not accessible or  does not exist ! \n";
        return 500;
      }

//first we read in the starting model and the measured data
    std::string modelfilename = jif3D::AskFilename("Starting model Filename: ");
    std::string extension = jif3D::GetFileExtension(modelfilename);
    std::cout << "Extension: " << extension << std::endl;
//we read in the starting modelfile
    jif3D::X3DModel Model;
    Model.SetOrigin(xorigin, yorigin, 0.0);
    if (extension.compare(".mod") == 0)
      {
        Model.ReadModEM(modelfilename);
        Model.WriteNetCDF(modelfilename + ".nc");
      }
    else
      {
        Model.ReadNetCDF(modelfilename);
        Model.WriteModEM(modelfilename + ".dat");
      }

    Model.WriteUBC(modelfilename + ".ubc");
    for (size_t i = 0; i < Model.GetConductivities().shape()[2]; ++i)
      {
        Model.SetConductivities()[0][0][i] *= (1 + conddelta * (i + 1));
      }

    const size_t ngrid = Model.GetConductivities().num_elements();
    jif3D::rvec InvModel(ngrid);
    std::copy(Model.GetConductivities().origin(),
        Model.GetConductivities().origin() + ngrid, InvModel.begin());
    Model.WriteVTK("start.vtk");
    auto Copier = boost::make_shared<jif3D::ModelCopyTransform>();

    boost::shared_ptr<jif3D::ChainedTransform> ConductivityTransform(
        new jif3D::ChainedTransform);

    if (vm.count("loglim"))
      {
        ConductivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::LogLimTransform(mincond, maxcond)));
      }
    else
      {
        jif3D::rvec RefModel(ngrid, 1.0);
//because the tanh transform is used inside a logarithmic transform
//we need to take the natural logarithm of the actual minimum and maximum
        ConductivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
        ConductivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::LogTransform(RefModel)));
      }

    boost::shared_ptr<jif3D::MultiSectionTransform> MTTransform(
        new jif3D::MultiSectionTransform(InvModel.size(), 0, ngrid,
            ConductivityTransform));
    subrange(InvModel, 0, ngrid) = subrange(MTTransform->PhysicalToGeneralized(InvModel),
        0, ngrid);
    jif3D::rvec GridCov(ngrid, 1.0);
    //std::copy(CovModVec.begin(), CovModVec.begin() + ngrid, GridCov.begin());
    //now we setup the regularization
    jif3D::ThreeDMTModel TearModX, TearModY, TearModZ;
    boost::shared_ptr<jif3D::RegularizationFunction> Regularization =
        RegSetup.SetupObjective(vm, Model, GridCov, TearModX, TearModY, TearModZ);
    boost::shared_ptr<jif3D::MultiSectionTransform> ModRegTrans(
        new jif3D::MultiSectionTransform(InvModel.size(), 0, ngrid, Copier));
    if (vm.count("refmodel"))
      {
        jif3D::X3DModel RefModel;
        RefModel.ReadNetCDF(RefModelName);
        jif3D::rvec RefVec(RefModel.GetNModelElements());
        std::copy(RefModel.GetConductivities().origin(),
            RefModel.GetConductivities().origin() + RefModel.GetNModelElements(),
            RefVec.begin());
        Regularization->SetReferenceModel(ModRegTrans->PhysicalToGeneralized(RefVec));
      }

    if (vm.count("regcheck"))
      {
        std::cout << " Regularization: " << Regularization->CalcMisfit(InvModel)
            << std::endl;
        jif3D::rvec reggrad = Regularization->CalcGradient(InvModel);
        jif3D::rvec RegVals(Regularization->GetDataDifference());
        const size_t nmod = InvModel.size();
        if (RegVals.size() != nmod * 3)
          {
            std::cerr << "Number of model parameters " << nmod << " is not 3* "
                << RegVals.size() << std::endl;
            return 100;
          }
        std::copy(RegVals.begin(), RegVals.begin() + nmod, Model.SetData().origin());
        Model.WriteVTK(modelfilename + ".regx.vtk");
        std::copy(RegVals.begin() + nmod, RegVals.begin() + 2 * nmod,
            Model.SetData().origin());
        Model.WriteVTK(modelfilename + ".regy.vtk");
        std::copy(RegVals.begin() + 2 * nmod, RegVals.begin() + 3 * nmod,
            Model.SetData().origin());
        Model.WriteVTK(modelfilename + ".regz.vtk");
        std::copy(reggrad.begin(), reggrad.end(), Model.SetData().origin());
        Model.WriteVTK(modelfilename + ".reggrad.vtk");
        return 0;
      }

    //if we have not specified a background in the model file
    if (Model.GetBackgroundConductivities().empty())
      {
        jif3D::rvec Cond(Model.GetConductivities().num_elements());
        std::copy(Model.GetConductivities().origin(),
            Model.GetConductivities().origin() + Model.GetConductivities().num_elements(),
            Cond.begin());
        std::nth_element(Cond.begin(), Cond.begin() + Cond.size() / 2, Cond.end());
        double MedCond = Cond(Cond.size() / 2);
        std::vector<double> bg_thick(2), bg_cond(2, MedCond);
        bg_thick.at(0) = Model.GetZCoordinates().at(Model.GetZCoordinates().size() - 1);
        bg_thick.at(1) = bg_thick.at(0);
        Model.SetBackgroundConductivities(bg_cond);
        Model.SetBackgroundThicknesses(bg_thick);
      }

//get the name of the file containing the data and read it in
    std::string datafilename = jif3D::AskFilename("Data Filename: ");
    extension = jif3D::GetFileExtension(datafilename);
//read in data
    jif3D::MTData DataMT;
    if (extension.compare(".dat") == 0)
      {
        DataMT.ReadModEM(datafilename);
        DataMT.WriteNetCDF("start_data.nc");
      }
    else
      {
        DataMT.ReadNetCDF(datafilename);
        DataMT.WriteModEM("start_data.dat");
      }
    DataMT.PlotMeasurementConfiguration(datafilename + ".vtk");
    std::vector<double> ZError = DataMT.GetErrors();
    std::vector<double> DataError = ZError;
    std::vector<double> MinErr(ZError.size());
    if (vm.count("inderrors"))
      {
        MinErr = jif3D::ConstructError(DataMT.GetData(), ZError, relerr);
      }
    else
      {
        if (vm.count("rowerrors"))
          {
            MinErr = jif3D::ConstructMTRowError(DataMT.GetData(), relerr);
          }
        else
          {
            MinErr = jif3D::ConstructMTError(DataMT.GetData(), relerr);
          }
      }

    if (vm.count("forcecommon"))
      {
        std::copy(MinErr.begin(), MinErr.end(), ZError.begin());
      }
    else
      {
        std::transform(ZError.begin(), ZError.end(), MinErr.begin(), ZError.begin(),
            [](double a, double b)
              { return std::max(a,b);});
      }

    auto negerr = std::find_if(ZError.begin(), ZError.end(), [](double val)
      { return val <= 0.0;});
    if (negerr != ZError.end())
      {
        std::cerr << "Negative data error in element "
            << std::distance(ZError.begin(), negerr) << " value " << *negerr
            << " will cause problem in inversion, exiting " << std::endl;
        return 100;
      }

//if we don't have data inversion doesn't make sense;
    if (DataMT.GetData().empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

//we define a few constants that are used throughout the inversion
//    const size_t ndata = Data.size(); // unused

    jif3D::rvec CovModVec(ngrid, 1.0);
    if (vm.count("covmod"))
      {
        jif3D::ThreeDModelBase CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        std::string Filename(vm["covmod"].as<std::string>());
        CovModel = *jif3D::ReadAnyModel(Filename).get();
        const size_t ncovmod = CovModel.GetData().num_elements();
        if (ncovmod == 0)
          {
            std::cerr << "Trying to read in covariance model, but no values read !"
                << std::endl;
            return 200;
          }
        if (ncovmod != ngrid)
          {
            std::cerr
                << "Covariance model, does not have the same number of cells as inversion model !"
                << std::endl;
            return 200;
          }
        CovModVec.resize(ncovmod);
        std::copy(CovModel.GetData().origin(), CovModel.GetData().origin() + ncovmod,
            CovModVec.begin());
        jif3D::X3DModel WriteModel;
        WriteModel = CovModel;
        WriteModel.WriteVTK("cov.vtk");
      }

    std::vector<double> C = DataMT.GetDistortion();

    bool WantDistCorr = (DistCorr > 0);
    //if we want to correct for distortion within the inversion
    if (WantDistCorr)
      {
        //we need to expand the model vector to hold the
        //elements of the distortion matrix
        jif3D::rvec Grid(InvModel);
        InvModel.resize(ngrid + C.size());
        std::copy(Grid.begin(), Grid.end(), InvModel.begin());
        std::copy(C.begin(), C.end(), InvModel.begin() + ngrid);
        //also the diagonal of the model covariance needs to
        //accommodate the new parameters
        jif3D::rvec OldCov(CovModVec);
        CovModVec.resize(InvModel.size());
        std::fill(CovModVec.begin(), CovModVec.end(), 1.0);
        std::copy(OldCov.begin(), OldCov.end(), CovModVec.begin());
      }

    if (vm.count("crossmodel"))
      {
        boost::shared_ptr<jif3D::ThreeDModelBase> CrossModel(
            jif3D::ReadAnyModel(CrossModelName));
        if (CrossModel->GetNModelElements() != ngrid)
          {
            std::cerr
                << "Cross gradient model does not have the same size as inversion grid"
                << std::endl;
            return 100;
          }
        double cgweight = 1.0;
        std::cout << "Cross-gradient weight: ";
        std::cin >> cgweight;

        jif3D::rvec OldInv(InvModel);
        InvModel.resize(OldInv.size() + ngrid, true);
        CovModVec.resize(InvModel.size(), true);
        std::copy(CrossModel->GetData().origin(), CrossModel->GetData().origin() + ngrid,
            InvModel.begin() + OldInv.size());
        std::fill_n(CovModVec.begin() + OldInv.size(), ngrid, 1e-32);

        std::cout << "Inversion vector length: " << InvModel.size() << " " << ngrid << " "
            << OldInv.size() << std::endl;
        boost::shared_ptr<jif3D::MultiSectionTransform> CrossTrans(
            new jif3D::MultiSectionTransform(InvModel.size(), 0, ngrid, Copier));
        CrossTrans->AddSection(OldInv.size(), InvModel.size(), Copier);

        boost::shared_ptr<jif3D::CrossGradient> CrossConstr(
            new jif3D::CrossGradient(Model, TearModX, TearModY, TearModZ));
        Objective->AddObjective(CrossConstr, CrossTrans, cgweight, "Cross",
            jif3D::JointObjective::coupling);
      }

    auto TipperTransform = boost::shared_ptr<jif3D::MultiSectionTransform>(
        MTTransform->clone());
    if (WantDistCorr)
      {
        MTTransform->AddSection(ngrid, ngrid + DataMT.GetDistortion().size(), Copier);
        MTTransform->SetLength(ngrid + DataMT.GetDistortion().size());
      }

    boost::shared_ptr<jif3D::X3DFieldCalculator> FC = boost::make_shared<
        jif3D::X3DFieldCalculator>(TempDir, X3DName, zshift);

    jif3D::X3DMTCalculator Calculator(TempDir, X3DName, WantDistCorr, CleanFiles, FC);
    if (vm.count("opt"))
      {
        Calculator.SetGreenType1(jif3D::GreenCalcType::opt);
        Calculator.SetGreenType4(jif3D::GreenCalcType::opt);
      }
    boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> > X3DObjective(
        new jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator>(Calculator));
    if (vm.count("rhophi"))
      {
        X3DObjective->SetDataTransform(boost::make_shared<jif3D::ComplexLogTransform>());
        std::transform(ZError.begin(), ZError.end(), DataMT.GetData().begin(),
            ZError.begin(),
            [](double err, double dat)
              { return err/std::max(std::numeric_limits<double>::epsilon(),std::abs(dat));});
      }
    X3DObjective->SetObservedData(DataMT);
    X3DObjective->SetCoarseModelGeometry(Model);

//create objects for the misfit and error estimates
    if (vm.count("mtinvcovar"))
      {
        jif3D::comp_mat InvCov(DataMT.GetData().size(), DataMT.GetData().size());
        jif3D::ReadSparseMatrixFromNetcdf(MTInvCovarName, InvCov, "InvCovariance");
        X3DObjective->SetInvCovMat(InvCov);
      }
    else
      {
        X3DObjective->SetDataError(ZError);
      }

    if (MTWeight > 0.0)
      {
        Objective->AddObjective(X3DObjective, MTTransform, MTWeight, "MT",
            jif3D::JointObjective::datafit);
      }
    jif3D::TipperData DataTip;

    jif3D::X3DTipperCalculator TipCalc(TempDir, X3DName, true, FC);
    boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::X3DTipperCalculator> > TipperObjective(
        new jif3D::ThreeDModelObjective<jif3D::X3DTipperCalculator>(TipCalc));
    if (vm.count("tipperdata"))
      {

        std::string tipext = jif3D::GetFileExtension(TipperName);
        if (tipext == ".nc")
          {
            DataTip.ReadNetCDF(TipperName);
          }
        else
          {
            DataTip.ReadModEM(TipperName);
            DataTip.WriteNetCDF(TipperName + ".conv.nc");
          }
        size_t ntip = DataTip.GetData().size();
        jif3D::rvec TE(ntip, 0.0);
        std::vector<double> TError = DataTip.GetErrors();
        for (size_t i = 0; i < ntip / 4; ++i)
          {
            double abs_re = std::sqrt(
                std::pow(DataTip.GetData().at(4 * i), 2)
                    + std::pow(DataTip.GetData().at(4 * i + 2), 2));
            double abs_im = std::sqrt(
                std::pow(DataTip.GetData().at(4 * i + 1), 2)
                    + std::pow(DataTip.GetData().at(4 * i + 3), 2));
            TE(4 * i) = std::max(std::max(abs_re * tiprelerr, TError.at(4 * i)),
                TError.at(4 * i + 2));
            TE(4 * i + 2) = TE(4 * i);
            TE(4 * i + 1) = std::max(std::max(abs_im * tiprelerr, TError.at(4 * i)),
                TError.at(4 * i + 2));
            TE(4 * i + 3) = TE(4 * i + 1);
          }
        std::copy(TE.begin(), TE.end(), TError.begin());

        std::transform(TError.begin(), TError.end(), TError.begin(), [](double a)
          { return std::max(a,tipperr);});
        TipperObjective->SetObservedData(DataTip);
        TipperObjective->SetCoarseModelGeometry(Model);
        TipperObjective->SetDataError(TError);
        Objective->AddObjective(TipperObjective, TipperTransform, TipperWeight, "Tipper",
            jif3D::JointObjective::datafit);
      }

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    Objective->AddObjective(Regularization, ModRegTrans, lambda, "Regularization",
        jif3D::JointObjective::regularization);

    if (WantDistCorr)
      {
        size_t nstats = DataMT.GetExIndices().size() / DataMT.GetFrequencies().size();
        jif3D::rvec CRef(nstats * 4);

        for (size_t i = 0; i < nstats; ++i)
          {
            CRef(i * 4) = 1.0;
            CRef(i * 4 + 1) = 0.0;
            CRef(i * 4 + 2) = 0.0;
            CRef(i * 4 + 3) = 1.0;
          }
        jif3D::X3DModel DistModel;

        DistModel.SetMeshSize(nstats * 4, 1, 1);
        boost::shared_ptr<jif3D::RegularizationFunction> DistReg(
            new jif3D::MinDiffRegularization(DistModel));
        DistReg->SetReferenceModel(CRef);
        boost::shared_ptr<jif3D::MultiSectionTransform> DistRegTrans(
            new jif3D::MultiSectionTransform(InvModel.size(), ngrid, ngrid + CRef.size(),
                Copier));

        Objective->AddObjective(DistReg, DistRegTrans, DistCorr, "DistReg",
            jif3D::JointObjective::regularization);
      }

    auto CovObj = boost::make_shared<jif3D::MultiSectionCovariance>(InvModel.size());
    if (CovWidth != 0.0)
      {
        boost::shared_ptr<jif3D::GeneralCovariance> StochCov = boost::make_shared<
            jif3D::StochasticCovariance>(CovModVec, Model.GetModelShape()[0],
            Model.GetModelShape()[1], Model.GetModelShape()[2], CovWidth, 1.0, 1.0);
        CovObj->AddSection(0, ngrid, StochCov);
        if (WantDistCorr)
          {
            boost::shared_ptr<jif3D::GeneralCovariance> DistCov = boost::make_shared<
                jif3D::DiagonalCovariance>();
            CovObj->AddSection(ngrid, InvModel.size(), DistCov);
          }
      }
    else
      {
        CovObj->AddSection(0, InvModel.size(),
            boost::make_shared<jif3D::DiagonalCovariance>(CovModVec));
      }

    boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovObj);

    if (vm.count("writegrad"))
      {
        jif3D::rvec GM = MTTransform->GeneralizedToPhysical(InvModel);
        std::cout << " MT Data Misfit: " << X3DObjective->CalcMisfit(GM) << std::endl;
        jif3D::rvec grad = X3DObjective->CalcGradient(GM);

        const size_t nval = Model.GetNModelElements();
        std::copy(grad.begin(), grad.begin() + nval, Model.SetData().origin());
        Model.WriteVTK(modelfilename + ".mtdatagrad.vtk");
        jif3D::rvec covgrad = CovObj->ApplyCovar(grad);
        std::copy(covgrad.begin(), covgrad.begin() + nval, Model.SetData().origin());
        Model.WriteVTK(modelfilename + ".mtdatacovgrad.vtk");

        if (vm.count("tipperdata"))
          {
            GM = TipperTransform->GeneralizedToPhysical(InvModel);
            std::cout << " Tipper Data Misfit: " << TipperObjective->CalcMisfit(GM)
                << std::endl;
            jif3D::rvec tipgrad = TipperObjective->CalcGradient(GM);

            ublas::subrange(grad, 0, nval) = tipgrad;
            std::copy(tipgrad.begin(), tipgrad.begin() + nval, Model.SetData().origin());
            Model.WriteVTK(modelfilename + ".tipdatagrad.vtk");

            jif3D::rvec tipcovgrad = CovObj->ApplyCovar(grad);

            std::copy(tipcovgrad.begin(), tipcovgrad.begin() + nval,
                Model.SetData().origin());
            Model.WriteVTK(modelfilename + ".tipdatacovgrad.vtk");
          }

        GM = ModRegTrans->GeneralizedToPhysical(InvModel);
        std::cout << " Regularization: " << Regularization->CalcMisfit(GM) << std::endl;
        //the full gradient can also contain values for distortion correction
        //this becomes important when applying the covariance below
        //for this output we are not interested in the distortion part
        //so we pad the vector and only look at the conductivity values
        jif3D::rvec reggrad(grad.size(), 0.0);
        ublas::subrange(grad, 0, nval) = Regularization->CalcGradient(GM);

        std::copy(grad.begin(), grad.begin() + nval, Model.SetData().origin());
        Model.WriteVTK(modelfilename + ".reggrad.vtk");
        covgrad = CovObj->ApplyCovar(grad);
        std::copy(covgrad.begin(), covgrad.begin() + nval, Model.SetData().origin());
        Model.WriteVTK(modelfilename + ".regcovgrad.vtk");

        return 0;
      }

    size_t maxiter = 30;
    std::cout << "Maximum number of iterations: ";
    std::cin >> maxiter;
    std::cout << "Performing inversion." << std::endl;

    std::ofstream misfitfile("misfit.out");
    std::ofstream rmsfile("rms.out");
    std::ofstream weightfile("weights.out");

    size_t iteration = 0;
    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();
    double InitialMisfit = Objective->CalcMisfit(InvModel);
    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);
    StoreRMS(rmsfile, 0, *Objective);

    bool terminate = false;
    while (iteration < maxiter && !terminate)
      {
        try
          {

            std::cout << "\n\n Iteration: " << iteration << std::endl;

            //update the inversion model
            Optimizer->MakeStep(InvModel);
            Objective->MultiplyWeights(jif3D::JointObjective::regularization,
                coolingfactor);
            ++iteration;

            //write out some information about misfit to the screen
            jif3D::rvec DistModel = MTTransform->GeneralizedToPhysical(InvModel);
            std::copy(DistModel.begin(), DistModel.begin() + Model.GetNModelElements(),
                Model.SetData().origin());
            Model.WriteVTK(modelfilename + jif3D::stringify(iteration) + ".mt.inv.vtk");
            Model.WriteNetCDF(modelfilename + jif3D::stringify(iteration) + ".mt.inv.nc");
            if (DistCorr > 0)
              {
                std::copy(DistModel.begin() + Model.GetNModelElements(),
                    DistModel.begin() + Model.GetNModelElements() + C.size(), C.begin());
                DataMT.SetDistortion(C);
                DataMT.WriteNetCDF(
                    modelfilename + jif3D::stringify(iteration) + ".dist_imp.nc");
              }
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm() << std::endl;
            //and write the current misfit for all objectives to a misfit file
            StoreMisfit(misfitfile, iteration, Optimizer->GetMisfit(), *Objective);
            StoreRMS(rmsfile, iteration, *Objective);
            StoreWeights(weightfile, iteration, *Objective);
            std::cout << "\n\n";

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

    if (vm.count("crossmodel"))
      {
        int index = std::distance(Objective->GetObjectiveTypes().begin(),
            std::find(Objective->GetObjectiveTypes().begin(),
                Objective->GetObjectiveTypes().end(), jif3D::JointObjective::coupling));
        std::cout << "Index for CG: " << index << std::endl;
        //write out the reference model, this is mainly for debugging to see
        //if it has been modified at all (which it should not)
        std::copy(InvModel.end() - Model.GetNModelElements(), InvModel.end(),
            Model.SetConductivities().origin());
        Model.WriteVTK("crossref.final.vtk");
        //write out the final cross gradient model
        jif3D::rvec CG(Objective->GetObjective(index).GetIndividualMisfit());
        const size_t nx = Model.GetData().shape()[0];
        const size_t ny = Model.GetData().shape()[1];
        const size_t nz = Model.GetData().shape()[2];
        const size_t nmod = nx * ny * nz;
        jif3D::ThreeDModelBase::t3DModelData XGrad(boost::extents[nx][ny][nz]);
        jif3D::ThreeDModelBase::t3DModelData YGrad(boost::extents[nx][ny][nz]);
        jif3D::ThreeDModelBase::t3DModelData ZGrad(boost::extents[nx][ny][nz]);
        std::copy(CG.begin(), CG.begin() + nmod, XGrad.origin());
        std::copy(CG.begin() + nmod, CG.begin() + 2 * nmod, YGrad.origin());
        std::copy(CG.begin() + 2 * nmod, CG.begin() + 3 * nmod, ZGrad.origin());
        jif3D::Write3DVectorModelToVTK("crossgrad.vtk", "CrossGrad",
            Model.GetXCoordinates(), Model.GetYCoordinates(), Model.GetZCoordinates(),
            XGrad, YGrad, ZGrad);
      }

    InvModel = MTTransform->GeneralizedToPhysical(InvModel);

    std::copy(InvModel.begin(),
        InvModel.begin() + Model.GetConductivities().num_elements(),
        Model.SetConductivities().origin());
    std::copy(InvModel.begin() + Model.GetNModelElements(), InvModel.end(), C.begin());
//calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;

    if (MTWeight > 0.0)
      {
        DataMT.SetDistortion(C);
        DataMT.WriteNetCDF(modelfilename + ".dist_imp.nc");
        std::vector<double> SData(X3DObjective->GetSyntheticData().begin(),
            X3DObjective->GetSyntheticData().end());
        DataMT.SetDataAndErrors(SData, X3DObjective->GetDataError());
        DataMT.WriteNetCDF(modelfilename + ".inv_imp.nc");
        std::vector<double> Misfit(X3DObjective->GetIndividualMisfit().begin(),
            X3DObjective->GetIndividualMisfit().end());
        DataMT.SetDataAndErrors(Misfit, X3DObjective->GetDataError());
        DataMT.WriteNetCDF(modelfilename + ".diff_imp.nc");
      }

    if (vm.count("tipperdata"))
      {
        std::vector<double> STipper(TipperObjective->GetSyntheticData().begin(),
            TipperObjective->GetSyntheticData().end());
        std::vector<double> TippE = TipperObjective->GetDataError();
        DataTip.SetDataAndErrors(STipper, TippE);
        DataTip.WriteNetCDF(modelfilename + ".inv_tip.nc");
        DataTip.PlotMeasurementConfiguration(modelfilename + ".tip.");
        std::vector<double> TippFit(TipperObjective->GetIndividualMisfit().begin(),
            TipperObjective->GetIndividualMisfit().end());
        DataTip.SetDataAndErrors(TippFit, TippE);
        DataTip.WriteNetCDF(modelfilename + ".diff_tip.nc");

      }
//and write out the data and model
//here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    Model.WriteXYZ(modelfilename + ".inv.xyz");
    Model.WriteUBC(modelfilename + ".inv.ubc");
    std::cout << std::endl;
#ifdef HAVEHPX
    return hpx::finalize();
#endif
    return 0;
  }

int main(int argc, char* argv[])
  {

    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("xorigin",
        po::value(&xorigin)->default_value(0.0),
        "The origin for the inversion grid in x-direction")("yorigin",
        po::value(&yorigin)->default_value(0.0),
        "The origin for the inversion grid in y-direction")("covmod",
        po::value<std::string>(), "A file containing the model covariance")("mincond",
        po::value(&mincond)->default_value(1e-6),
        "The minimum value for conductivity in S/m")("maxcond",
        po::value(&maxcond)->default_value(5),
        "The maximum value for conductivity in S/m")("tempdir", po::value<std::string>(),
        "The name of the directory to store temporary files in")("distcorr",
        po::value(&DistCorr)->default_value(0),
        "Correct for distortion within inversion, value is regularization factor")(
        "x3dname", po::value(&X3DName)->default_value("x3d"),
        "The name of the executable for x3d")("mtinvcovar",
        po::value<std::string>(&MTInvCovarName),
        "Inverse covariance matrix to use in MT misfit calculation.")("inderrors",
        "Use the individual errors for each element instead of the same error floor for all elements")(
        "rowerrors",
        "Use rowwise errors for each element instead of the same error floor for all elements")(
        "forcecommon", "Use exactly the same error mtrelerr * Berd(Z) for all elements.")(
        "mtrelerr", po::value(&relerr)->default_value(0.02),
        "Error floor for impedance estimates.")("coolingfactor",
        po::value(&coolingfactor)->default_value(1.0),
        "The factor to multiply the weight for the regularization at each iteration EXPERIMENTAL")(
        "debug", "Show debugging output.")("opt",
        "Use opt for Green's function calculation in x3d.")("zshift",
        po::value(&zshift)->default_value(0.0),
        "Shift the depth of all sites by this distance in m, EXPERIMENTAL")("refmodel",
        po::value(&RefModelName),
        "The name of the reference model to substract before calculating smoothness")(
        "crossmodel", po::value(&CrossModelName),
        "The name of a model to use as a cross-gradient constraint")("regcheck",
        "Only perform a regularization calculation")("writegrad",
        "Write out the gradient for the data and quit")("conddelta",
        po::value(&conddelta)->default_value(0.001),
        "The relative amount by which the conductivities in the first row of cells is disturbed to ensure proper gradient calculation")(
        "rhophi", "Use apparent resistivity and phase instead of impedance")("cleanfiles",
        po::value(&CleanFiles)->default_value(true),
        "Clean up all temporary files at the end of the program.")("titan",
        "Read in Titan24 data.")("loglim",
        "Use logarithmic transform to bound model parameters")("tipperdata",
        po::value(&TipperName), "The name for the tipper data")("tipperlambda",
        po::value(&TipperWeight)->default_value(1.0), "The weight for the tipper data")(
        "tippererr", po::value(&tipperr)->default_value(0.005),
        "The absolute error for the tipper data")("tiprelerr",
        po::value(&tiprelerr)->default_value(0.02),
        "The relative error for the tipper data")("mtlambda",
        po::value(&MTWeight)->default_value(1.0), "The weight for the MT data")(
        "stochcov", po::value(&CovWidth)->default_value(0),
        "Width of stochastic regularization, enabled if > 0, EXPERIMENTAL");

    desc.add(RegSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());

#ifdef HAVEHPX
    return hpx::init(desc, argc, argv);
#else
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
#endif

  }
