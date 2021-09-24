//============================================================================
// Name        : jointinv.cpp
// Author      : May 12, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
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
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../ModelBase/EqualGeometry.h"
#include "../ModelBase/ReadAnyModel.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/MultiSectionCovariance.h"
#include "../Inversion/StochasticCovariance.h"

#include "../Regularization/MinDiffRegularization.h"
#include "../Inversion/ModelTransforms.h"
#include "../SurfaceWaves/SurfaceWaveModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../SurfaceWaves/SurfaceWaveCalculator.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../MT/ReadWriteImpedances.h"
#include "SetupSW.h"
#include "SetupTomo.h"
#include "SetupGravity.h"
#include "SetupMT.h"
#include "SetupInversion.h"
#include "SetupRegularization.h"
#include "SetupCoupling.h"
#include "InversionOutput.h"

namespace ublas = boost::numeric::ublas;
namespace po = boost::program_options;

/** \addtogroup joint Joint inversion routines */
/* @{ */

/*! \file jointinv.cpp
 * The main joint inversion program. The main task of the program is to read in the appropriate files
 * and options and from these settings assemble the objects that perform the actual work. Also,
 * the main program manages the output of inversion results and data among other statistics.
 */

//first we create objects that manage the setup of individual parts
//that way we can reuse these objects so that other programs use
//exactly the same options
#ifdef HAVETRAVELTIME
jif3D::SetupTomo TomoSetup;
#else
jif3D::SetupSW TomoSetup;
#endif

jif3D::SetupGravity GravitySetup;
jif3D::SetupMT MTSetup;
jif3D::SetupInversion InversionSetup;
jif3D::SetupRegularization RegSetup;
jif3D::SetupCoupling CouplingSetup;
bool WaveletParm = false;
double xorigin, yorigin;
double coolingfactor = 1.0;
int saveinterval = 1;
double CovWidth = 3.0;
double Cov_nu = 1.0;
double Cov_sigma = 1.0;

int hpx_main(boost::program_options::variables_map &vm)
  {

    if (vm.count("wavelet"))
      {
        WaveletParm = true;
      }
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

    //we want some output so we set Verbose in the constructor to true
    boost::shared_ptr<jif3D::JointObjective> Objective;
    Objective = boost::shared_ptr<jif3D::JointObjective>(new jif3D::JointObjective(true));

    //we need a number of transformation objects to translate the generalized model parameters
    //used by the inversion algorithm to physically meaningful parameters for the forward
    //calculation. The exact type of transformation depends on the chosen coupling
    //and is assigned below
    boost::shared_ptr<jif3D::GeneralModelTransform> TomoTransform, MTTransform,
        GravityTransform;
    //coupling setup is responsible to set the appropriate transformation
    //as we have to use different ones depending on the chose coupling mechanism

    //we need the geometry of the starting model to setup
    //the transformations
    std::string geometryfilename;
    if (vm.count("geometrygrid"))
      {
        geometryfilename = vm["geometrygrid"].as<std::string>();
      }
    else
      {
        geometryfilename = jif3D::AskFilename("Inversion Model Geometry: ");
      }
    boost::shared_ptr<jif3D::ThreeDModelBase> StartModel = jif3D::ReadAnyModel(
        geometryfilename);
    const size_t ngrid = StartModel->GetData().num_elements();

    jif3D::rvec CovModVec(3 * ngrid,1.0);
    if (vm.count("covmod"))
      {
        jif3D::ThreeDModelBase CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        std::string Filename(vm["covmod"].as<std::string>());
        if (!boost::filesystem::exists(Filename))
          {
            std::cerr << "Covariance model has been specified but file: " << Filename
                << " does not exist " << std::endl;
            return 100;
          }
        CovModel = *jif3D::ReadAnyModel(Filename).get();
        const size_t ncovmod = CovModel.GetData().num_elements();
        if (ncovmod != ngrid)
          {
            std::cerr << "Covariance grid not the same size as inversion grid "
                << std::endl;
            return 100;
          }
        std::copy(CovModel.GetData().origin(), CovModel.GetData().origin() + ncovmod,
            CovModVec.begin());
        std::copy(CovModel.GetData().origin(), CovModel.GetData().origin() + ncovmod,
            CovModVec.begin() + ngrid);
        std::copy(CovModel.GetData().origin(), CovModel.GetData().origin() + ncovmod,
            CovModVec.begin() + 2*ngrid);
      }

    if (vm.count("tomocov") || vm.count("gravcov") || vm.count("mtcov"))
      {
        std::cout
            << "Setting individual covariances currently only works with cross-gradient or mutual information coupling !"
            << std::endl;
        std::fill(CovModVec.begin(), CovModVec.end(), 1.0);
      }

    if (vm.count("tomocov"))
      {
        jif3D::ThreeDModelBase CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        std::string Filename(vm["tomocov"].as<std::string>());
        CovModel = *jif3D::ReadAnyModel(Filename).get();
        const size_t ncovmod = CovModel.GetData().num_elements();
        if (ncovmod != ngrid)
          {
            std::cerr << "Tomography covariance grid not the same size as inversion grid "
                << std::endl;
            return 100;
          }
        std::copy(CovModel.GetData().origin(), CovModel.GetData().origin() + ncovmod,
            CovModVec.begin());
      }

    if (vm.count("gravcov"))
      {
        jif3D::ThreeDModelBase CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        std::string Filename(vm["gravcov"].as<std::string>());
        CovModel = *jif3D::ReadAnyModel(Filename).get();
        const size_t ncovmod = CovModel.GetData().num_elements();
        if (ncovmod != ngrid)
          {
            std::cerr << "Gravity covariance grid not the same size as inversion grid "
                << std::endl;
            return 100;
          }
        std::copy(CovModel.GetData().origin(), CovModel.GetData().origin() + ncovmod,
            CovModVec.begin() + ngrid);
      }

    if (vm.count("mtcov"))
      {
        jif3D::ThreeDModelBase CovModel;
        //we store the covariances in a seismic model file
        //but we do not have to have an equidistant grid
        std::string Filename(vm["mtcov"].as<std::string>());
        CovModel = *jif3D::ReadAnyModel(Filename).get();
        const size_t ncovmod = CovModel.GetData().num_elements();
        if (ncovmod != ngrid)
          {
            std::cerr << "MT covariance grid not the same size as inversion grid "
                << std::endl;
            return 100;
          }
        std::copy(CovModel.GetData().origin(), CovModel.GetData().origin() + ncovmod,
            CovModVec.begin() + 2 * ngrid);
      }

    CouplingSetup.SetupTransforms(vm, *StartModel, TomoTransform, GravityTransform,
        MTTransform, WaveletParm);

    //read in the tomography model and setup the options that are applicable
    //for the seismic tomography part of the inversion
    bool havetomo = TomoSetup.SetupObjective(vm, *Objective.get(), TomoTransform, xorigin,
        yorigin);
    if (havetomo && !EqualGridGeometry(TomoSetup.GetModel(), *StartModel))
      {
        throw jif3D::FatalException(
            "Tomography model does not have the same geometry as starting model");
      }
    //setup the gravity part of the joint inversion
    bool havegrav = GravitySetup.SetupObjective(vm, *Objective.get(), GravityTransform,
        xorigin, yorigin, TempDir);
    //if we have a seismic and a gravity objective function, we have
    //to make sure that the starting models have the same geometry (not considering refinement)
    if (havegrav && !EqualGridGeometry(*StartModel, GravitySetup.GetScalModel()))
      {
        throw jif3D::FatalException(
            "Gravity model does not have the same geometry as starting model");
      }
    //setup the MT part of the joint inversion
    bool havemt = MTSetup.SetupObjective(vm, *Objective.get(), MTTransform, xorigin,
        yorigin, TempDir);
    //if we have a seismic and a MT objective function, we have
    //to make sure that the starting models have the same geometry (not considering refinement)
    if (havemt && !EqualGridGeometry(MTSetup.GetModel(), *StartModel))
      {
        throw jif3D::FatalException(
            "MT model does not have the same geometry as starting model");
      }
    //now we setup the regularization
    jif3D::ThreeDSeismicModel TearModX, TearModY, TearModZ;
    boost::shared_ptr<jif3D::RegularizationFunction> Regularization =
        RegSetup.SetupObjective(vm, *StartModel, jif3D::rvec(ngrid, 1.0), TearModX,
            TearModY, TearModZ);

    //the vector InvModel will hold the current inversion model
    //depending on the chosen coupling mechanism it will have different size
    //so we fill its content in the object Coupling setup
    jif3D::rvec InvModel;
    CouplingSetup.SetupModelVector(vm, InvModel, *StartModel, TomoSetup.GetModel(),
        GravitySetup.GetScalModel(), MTSetup.GetModel(), *Objective.get(), Regularization,
        RegSetup.GetSubStart(), TearModX, TearModY, TearModZ, CovModVec);

    //finally ask for the maximum number of iterations
    size_t maxiter = 1;
    if (!vm.count("iterations"))
      {
        std::cout << "Maximum iterations: ";
        std::cin >> maxiter;
      }
    else
      {
        maxiter = vm["iterations"].as<int>();
      }
    //note the start time of the core calculations for statistics
    //and output some status information
    boost::posix_time::ptime starttime = boost::posix_time::microsec_clock::local_time();

    std::string modelfilename = "result";
    //write out the seismic source and receiver positions for plotting
    //and general quality control
    if (havetomo)
      {
#ifdef HAVETRAVELTIME
        auto TomoData(TomoSetup.GetTomoObjective().GetObservedData());
#else
        auto TomoData(TomoSetup.GetSurfaceWaveObjective().GetObservedData());
#endif
        TomoData.WriteMeasurementPoints(modelfilename + ".rec.vtk");
        //TomoData.WriteSourcePoints(modelfilename + ".sor.vtk");
      }

    const size_t nparm = InvModel.size();
    const size_t ncovmod = CovModVec.size();
    std::cout << nparm << " Inversion parameters " << ncovmod << " Covariance values "
        << std::endl;


    boost::shared_ptr<jif3D::MultiSectionTransform> DistRegTrans;
    if (havemt)
      {
        const size_t nmtsites =
            MTSetup.GetMTObjective().GetObservedData().GetMeasPosX().size();
        jif3D::rvec CRef(nmtsites * 4);
        for (size_t i = 0; i < nmtsites; ++i)
          {
            CRef(i * 4) = 1.0;
            CRef(i * 4 + 1) = 0.0;
            CRef(i * 4 + 2) = 0.0;
            CRef(i * 4 + 3) = 1.0;
          }
        boost::shared_ptr<jif3D::GeneralModelTransform> Copier(
            new jif3D::ModelCopyTransform);
        //this transformation only becomes active if we use distortion correction with MT data
        //in this case we add extra inversion parameters beyond the current ones
        //so we set it up here that we can access it later, but with a parameter setting that only
        //works if we actually have distortion correction, so we have to be careful later
        DistRegTrans = boost::make_shared<jif3D::MultiSectionTransform>(
            InvModel.size() + CRef.size(), InvModel.size(), InvModel.size() + CRef.size(),
            Copier);
        auto MTData = MTSetup.GetMTObjective().GetObservedData();
        MTData.WriteMeasurementPoints(modelfilename + ".mt_sites.vtk");

        //if we want to correct for distortion within the inversion
        if (MTSetup.GetDistCorr() > 0)
          {
            std::vector<double> C =
                MTSetup.GetMTObjective().GetObservedData().GetDistortion();

            //we need to expand the model vector to hold the
            //elements of the distortion matrix
            jif3D::rvec Grid(InvModel);
            InvModel.resize(InvModel.size() + C.size());
            std::copy(Grid.begin(), Grid.end(), InvModel.begin());
            std::copy(C.begin(), C.end(), InvModel.begin() + Grid.size());
            //also the diagonal of the model covariance needs to
            //accommodate the new parameters
            jif3D::rvec OldCov(CovModVec);
            CovModVec.resize(InvModel.size());
            std::fill(CovModVec.begin(), CovModVec.end(), 1.0);
            std::copy(OldCov.begin(), OldCov.end(), CovModVec.begin());

            jif3D::X3DModel DistModel;
            DistModel.SetMeshSize(nmtsites * 4, 1, 1);
            boost::shared_ptr<jif3D::RegularizationFunction> DistReg(
                new jif3D::MinDiffRegularization(DistModel));
            DistReg->SetReferenceModel(CRef);

            dynamic_cast<jif3D::MultiSectionTransform*>(MTTransform.get())->SetLength(
                InvModel.size());
            //DistRegTrans->SetLength(InvModel.size());
            dynamic_cast<jif3D::MultiSectionTransform*>(MTTransform.get())->AddSection(
                Grid.size(), InvModel.size(), Copier);
            Objective->AddObjective(DistReg, DistRegTrans, MTSetup.GetDistCorr(),
                "DistReg", jif3D::JointObjective::regularization);
          }
      }

    std::cout << "Calculating initial misfit." << std::endl;
    size_t iteration = 0;
    std::ofstream misfitfile("misfit.out");
    std::ofstream rmsfile("rms.out");
    std::ofstream weightfile("weights.out");
    //calculate initial misfit
    double InitialMisfit = Objective->CalcMisfit(InvModel);
    StoreMisfit(misfitfile, 0, InitialMisfit, *Objective);
    StoreRMS(rmsfile, 0, *Objective);
    StoreWeights(weightfile, 0, *Objective);
    jif3D::ThreeDGravityModel GravModel(GravitySetup.GetScalModel());
    jif3D::X3DModel MTModel(MTSetup.GetModel());
    auto SWModel(TomoSetup.GetModel());

    if (!havetomo)
      SWModel = *StartModel;
    if (!havemt)
      MTModel = *StartModel;
    if (!havegrav)
      GravModel = *StartModel;

    std::cout << "Performing inversion." << std::endl;

    boost::shared_ptr<jif3D::MultiSectionCovariance> CovObj = boost::make_shared<
        jif3D::MultiSectionCovariance>(InvModel.size());
    if (CovWidth != 0.0)
      {
        ublas::vector_range<jif3D::rvec> TomoCov(CovModVec, ublas::range (0, ngrid));
        ublas::vector_range<jif3D::rvec> GravCov(CovModVec, ublas::range (ngrid, 2*ngrid));
        ublas::vector_range<jif3D::rvec> MTCov(CovModVec, ublas::range (2*ngrid, 3*ngrid));
        boost::shared_ptr<jif3D::GeneralCovariance> StochCovTomo = boost::make_shared<
            jif3D::StochasticCovariance>(TomoCov, StartModel->GetModelShape()[0],
            StartModel->GetModelShape()[1], StartModel->GetModelShape()[2], CovWidth,
            Cov_nu, Cov_sigma);
        boost::shared_ptr<jif3D::GeneralCovariance> StochCovGrav = boost::make_shared<
                    jif3D::StochasticCovariance>(GravCov, StartModel->GetModelShape()[0],
                    StartModel->GetModelShape()[1], StartModel->GetModelShape()[2], CovWidth,
                    Cov_nu, Cov_sigma);
        boost::shared_ptr<jif3D::GeneralCovariance> StochCovMT = boost::make_shared<
                    jif3D::StochasticCovariance>(MTCov, StartModel->GetModelShape()[0],
                    StartModel->GetModelShape()[1], StartModel->GetModelShape()[2], CovWidth,
                    Cov_nu, Cov_sigma);

        CovObj->AddSection(0, ngrid, StochCovTomo);
        CovObj->AddSection(ngrid, 2 * ngrid, StochCovGrav);
        CovObj->AddSection(2 * ngrid, 3 * ngrid, StochCovMT);
        if (MTSetup.GetDistCorr() > 0)
          {
            boost::shared_ptr<jif3D::GeneralCovariance> DistCov = boost::make_shared<
                jif3D::DiagonalCovariance>();
            CovObj->AddSection(3 * ngrid, InvModel.size(), DistCov);
          }
      }
    else
      {
        CovObj->AddSection(0, InvModel.size(),
            boost::make_shared<jif3D::DiagonalCovariance>(CovModVec));
      }

    //auto CovObj = boost::make_shared<jif3D::DiagonalCovariance>(CovVec);

    boost::shared_ptr<jif3D::GradientBasedOptimization> Optimizer =
        InversionSetup.ConfigureInversion(vm, Objective, InvModel, CovObj);

    bool terminate = false;
    jif3D::rvec OldModel(InvModel);
    //this is the core inversion loop, we make optimization steps
    //until either we reach the maximum number of iterations
    //or fulfill a termination criterion
    while (iteration < maxiter && !terminate)
      {
        terminate = true;
        //we catch all jif3D internal exceptions so that we can graciously
        //exit and write out some final information before stopping the program
        try
          {
            std::cout << "\n\n Iteration: " << iteration << std::endl;
            //we save the current model so we can go back to it
            //in case the optimization step fails
            OldModel = InvModel;
            //update the inversion model
            Optimizer->MakeStep(InvModel);

            Objective->MultiplyWeights(jif3D::JointObjective::regularization,
                coolingfactor);
            ++iteration;
            //we save all models at the iterations the user has selected (default is every iteration)
            //this way  we can look at the development
            // and use intermediate models in case something goes wrong
            if (iteration % saveinterval == 0)
              {
                if (havetomo)
                  {
#ifdef HAVETRAVELTIME
                    SaveModel(InvModel, *TomoTransform.get(), SWModel,
                        modelfilename + jif3D::stringify(iteration) + ".tomo.inv");
#else
                    jif3D::rvec TransModel = TomoTransform->GeneralizedToPhysical(
                        InvModel);
                    const size_t ncells = SWModel.GetNModelElements();
                    std::copy(TransModel.begin(), TransModel.begin() + ncells,
                        SWModel.SetData().origin());
                    jif3D::SurfaceWaveModel::t3DModelData values(SWModel.GetVp());
                    std::copy(TransModel.begin() + ncells,
                        TransModel.begin() + 2 * ncells, values.origin());
                    SWModel.SetVp(values);
                    std::copy(TransModel.begin() + 2 * ncells, TransModel.end(),
                        values.origin());
                    SWModel.SetDensAnomaly(values);
                    SWModel.WriteVTK(
                        modelfilename + jif3D::stringify(iteration) + ".tomo.inv.vtk");
                    SWModel.WriteNetCDF(
                        modelfilename + jif3D::stringify(iteration) + ".tomo.inv.nc");
#endif
                  }
                if (havemt)
                  {
                    SaveModel(InvModel, *MTTransform.get(), MTModel,
                        modelfilename + jif3D::stringify(iteration) + ".mt.inv");
                  }
                if (havegrav)
                  {
                    SaveModel(InvModel, *GravityTransform.get(), GravModel,
                        modelfilename + jif3D::stringify(iteration) + ".grav.inv");
                  }
              }
            //write out some information about misfit to the screen
            std::cout << "Currrent Misfit: " << Optimizer->GetMisfit() << std::endl;
            std::cout << "Currrent Gradient: " << Optimizer->GetGradNorm() << std::endl;
            //and write the current misfit for all objectives to a misfit file
            StoreMisfit(misfitfile, iteration, Optimizer->GetMisfit(), *Objective);
            StoreRMS(rmsfile, iteration, *Objective);
            StoreWeights(weightfile, iteration, *Objective);
            std::cout << "\n\n";
          } catch (jif3D::FatalException &e)
          {
            std::cerr << e.what() << std::endl;
            InvModel = OldModel;
            iteration = maxiter;
          }
        //we stop when either we do not make any improvement any more
        terminate = CheckConvergence(*Objective);
        //or the file abort exists in the current directory
        terminate = terminate || jif3D::WantAbort();
      }

    //calculate the predicted refraction data
    std::cout << "Calculating response of inversion model." << std::endl;
    //during the last iteration we might have performed steps in the line search
    //so we update the forward calculation for the last proper inversion model
    //we use the results from this calculation to save the final inversion synthetic data
    if (iteration > 0)
      {
        Objective->CalcMisfit(InvModel);
      }
    if (havetomo)
      {
#ifdef HAVETRAVELTIME
                    SaveModel(InvModel, *TomoTransform.get(), SWModel,
                        modelfilename + ".tomo.inv");
                    if (vm.count("writerays"))
                     {
                     TomoSetup.GetTomoObjective().GetCalculator().WriteRays("rays.vtk");
                     }
#else
        jif3D::rvec TransModel = TomoTransform->GeneralizedToPhysical(InvModel);
        const size_t ncells = SWModel.GetNModelElements();
        std::copy(TransModel.begin(), TransModel.begin() + ncells,
            SWModel.SetData().origin());
        jif3D::SurfaceWaveModel::t3DModelData values(SWModel.GetVp());
        std::copy(TransModel.begin() + ncells, TransModel.begin() + 2 * ncells,
            values.origin());
        SWModel.SetVp(values);
        std::copy(TransModel.begin() + 2 * ncells, TransModel.end(), values.origin());
        SWModel.SetDensAnomaly(values);
        SWModel.WriteVTK(modelfilename + ".tomo.inv.vtk");
        SWModel.WriteNetCDF(modelfilename + ".tomo.inv.nc");

        auto TomoData = TomoSetup.GetSurfaceWaveObjective().GetObservedData();
        jif3D::rvec TomoInvData(TomoSetup.GetSurfaceWaveObjective().GetSyntheticData());
        std::vector<double> TomoError(TomoSetup.GetSurfaceWaveObjective().GetDataError());
        TomoData.SetDataAndErrors(
            std::vector<double>(TomoInvData.begin(), TomoInvData.end()), TomoError);
        TomoData.WriteNetCDF(modelfilename + ".inv_tt.nc");
        jif3D::rvec TomoDiff(TomoSetup.GetSurfaceWaveObjective().GetIndividualMisfit());
        TomoData.SetDataAndErrors(std::vector<double>(TomoDiff.begin(), TomoDiff.end()),
            TomoError);
        TomoData.WriteNetCDF(modelfilename + ".diff_tt.nc");
#endif

      }
    //if we are inverting gravity data and have specified site locations
    if (havegrav)
      {

        SaveModel(InvModel, *GravityTransform.get(), GravModel,
            modelfilename + ".grav.inv");
        //and write out the data and model
        //here we have to distinguish again between scalar and ftg data
        if (GravitySetup.GetHaveScal())
          {
            GravModel = GravitySetup.GetScalModel();
            auto ScalGravData = GravitySetup.GetScalGravObjective().GetObservedData();
            jif3D::rvec ScalGravInvData(
                GravitySetup.GetScalGravObjective().GetSyntheticData());
            jif3D::SaveScalarGravityMeasurements(modelfilename + ".inv_sgd.nc",
                std::vector<double>(ScalGravInvData.begin(), ScalGravInvData.end()),
                ScalGravData.GetMeasPosX(), ScalGravData.GetMeasPosY(),
                ScalGravData.GetMeasPosZ(),
                GravitySetup.GetScalGravObjective().GetDataError());
            jif3D::rvec ScalDiff(
                GravitySetup.GetScalGravObjective().GetIndividualMisfit());
            jif3D::SaveScalarGravityMeasurements(modelfilename + ".diff_sgd.nc",
                std::vector<double>(ScalDiff.begin(), ScalDiff.end()),
                ScalGravData.GetMeasPosX(), ScalGravData.GetMeasPosY(),
                ScalGravData.GetMeasPosZ(),
                GravitySetup.GetScalGravObjective().GetDataError());
            jif3D::Write3DDataToVTK(modelfilename + ".inv_sgd.vtk", "grav_accel",
                std::vector<double>(ScalGravInvData.begin(), ScalGravInvData.end()),
                ScalGravData.GetMeasPosX(), ScalGravData.GetMeasPosY(),
                ScalGravData.GetMeasPosZ());

          }
        if (GravitySetup.GetHaveFTG())
          {
            auto FTGData = GravitySetup.GetFTGObjective().GetObservedData();
            GravModel = GravitySetup.GetFTGModel();
            jif3D::rvec FTGInvData(GravitySetup.GetFTGObjective().GetSyntheticData());
            jif3D::SaveTensorGravityMeasurements(modelfilename + ".inv_ftg.nc",
                std::vector<double>(FTGInvData.begin(), FTGInvData.end()),
                FTGData.GetMeasPosX(), FTGData.GetMeasPosY(), FTGData.GetMeasPosZ(),
                GravitySetup.GetFTGObjective().GetDataError());
            jif3D::rvec FTGDiff(GravitySetup.GetFTGObjective().GetIndividualMisfit());
            jif3D::SaveTensorGravityMeasurements(modelfilename + ".diff_ftg.nc",
                std::vector<double>(FTGDiff.begin(), FTGDiff.end()),
                FTGData.GetMeasPosX(), FTGData.GetMeasPosY(), FTGData.GetMeasPosZ(),
                GravitySetup.GetFTGObjective().GetDataError());
            jif3D::Write3DTensorDataToVTK(modelfilename + ".inv_ftg.vtk", "U",
                std::vector<double>(FTGInvData.begin(), FTGInvData.end()),
                FTGData.GetMeasPosX(), FTGData.GetMeasPosY(), FTGData.GetMeasPosZ());
          }
      }

    //if we are inverting MT data and have specified site locations
    if (havemt)
      {
        SaveModel(InvModel, *MTTransform.get(), MTModel, modelfilename + ".mt.inv");
        auto MTData = MTSetup.GetMTObjective().GetObservedData();
        if (MTSetup.GetDistCorr() > 0.0)
          {
            std::vector<double> C;
            jif3D::rvec tmp = DistRegTrans->GeneralizedToPhysical(InvModel);
            std::copy(tmp.begin(), tmp.end(), std::back_inserter(C));
            MTData.SetDistortion(C);
          }
        //calculate MT inversion result
        MTData.WriteNetCDF(modelfilename + ".dist_imp.nc");
        jif3D::rvec MTInvData(MTSetup.GetMTObjective().GetSyntheticData());
        MTData.SetDataAndErrors(std::vector<double>(MTInvData.begin(), MTInvData.end()),
            MTSetup.GetMTObjective().GetDataError());
        MTData.WriteNetCDF(modelfilename + ".inv_imp.nc");
        jif3D::rvec MTDiff(MTSetup.GetMTObjective().GetIndividualMisfit());
        MTData.SetDataAndErrors(std::vector<double>(MTDiff.begin(), MTDiff.end()),
            MTSetup.GetMTObjective().GetDataError());
        MTData.WriteNetCDF(modelfilename + ".diff_imp.nc");

        if (MTSetup.GetTipLambda() > 0)
          {
            auto TipData = MTSetup.GetTipObjective().GetObservedData();
            jif3D::rvec TipInvData(MTSetup.GetTipObjective().GetSyntheticData());
            TipData.SetDataAndErrors(
                std::vector<double>(TipInvData.begin(), TipInvData.end()),
                MTSetup.GetTipObjective().GetDataError());
            TipData.WriteNetCDF(modelfilename + ".inv_tip.nc");
            jif3D::rvec TipDiff(MTSetup.GetTipObjective().GetIndividualMisfit());
            TipData.SetDataAndErrors(std::vector<double>(TipDiff.begin(), TipDiff.end()),
                MTSetup.GetTipObjective().GetDataError());
            TipData.WriteNetCDF(modelfilename + ".diff_tip.nc");
          }

      }
    //if we are using cross gradient coupling, we want to output the final
    //values of the cross-gradient
    if (vm.count("crossgrad"))
      {
        auto ObjectiveTypes = Objective->GetObjectiveTypes();
        std::vector<std::string> Names = Objective->GetObjectiveNames();
        for (size_t i = 0; i < ObjectiveTypes.size(); ++i)
          {
            if (ObjectiveTypes.at(i) == jif3D::JointObjective::coupling)
              {
                jif3D::rvec CG(Objective->GetObjective(i).GetIndividualMisfit());
                const size_t nx = StartModel->GetData().shape()[0];
                const size_t ny = StartModel->GetData().shape()[1];
                const size_t nz = StartModel->GetData().shape()[2];
                const size_t nmod = nx * ny * nz;
                jif3D::ThreeDModelBase::t3DModelData XGrad(boost::extents[nx][ny][nz]);
                jif3D::ThreeDModelBase::t3DModelData YGrad(boost::extents[nx][ny][nz]);
                jif3D::ThreeDModelBase::t3DModelData ZGrad(boost::extents[nx][ny][nz]);
                std::copy(CG.begin(), CG.begin() + nmod, XGrad.origin());
                std::copy(CG.begin() + nmod, CG.begin() + 2 * nmod, YGrad.origin());
                std::copy(CG.begin() + 2 * nmod, CG.begin() + 3 * nmod, ZGrad.origin());
                std::string Name = Names.at(i);
                jif3D::Write3DVectorModelToVTK(Name + ".vtk", Name,
                    StartModel->GetXCoordinates(), StartModel->GetYCoordinates(),
                    StartModel->GetZCoordinates(), XGrad, YGrad, ZGrad);
                std::vector<double> CG_Cov(Objective->GetObjective(i).GetDataError());
                std::copy(CG_Cov.begin(), CG_Cov.begin() + nmod, XGrad.origin());
                std::copy(CG_Cov.begin() + nmod, CG_Cov.begin() + 2 * nmod,
                    YGrad.origin());
                std::copy(CG_Cov.begin() + 2 * nmod, CG_Cov.begin() + 3 * nmod,
                    ZGrad.origin());
                jif3D::Write3DVectorModelToVTK(Name + "_cov.vtk", Name,
                    StartModel->GetXCoordinates(), StartModel->GetYCoordinates(),
                    StartModel->GetZCoordinates(), XGrad, YGrad, ZGrad);
              }
          }

      }
    std::ofstream datadiffile("data.diff");
    std::copy(Objective->GetDataDifference().begin(),
        Objective->GetDataDifference().end(),
        std::ostream_iterator<double>(datadiffile, "\n"));
    boost::posix_time::ptime endtime = boost::posix_time::microsec_clock::local_time();
    double cachedruntime = (endtime - starttime).total_seconds();
    std::cout << "Runtime: " << cachedruntime << " s" << std::endl;
    std::cout << std::endl;

#ifdef HAVEHPX
    return hpx::finalize();
#endif
    return 0;
  }

int main(int argc, char *argv[])
  {
    //we also create a number of options that are specific to our joint inversion
    //or act globally so that they cannot be associated with one subsystem
    po::options_description desc("General options");

    desc.add_options()("help", "produce help message")("debug",
        "Write debugging information")("iterations", po::value<int>(),
        "The maximum number of iterations")("threads", po::value<int>(),
        "The number of openmp threads")("covmod", po::value<std::string>(),
        "A file containing the model covariance for all methods")("tomocov",
        po::value<std::string>(),
        "A file containing the model covariance only for tomography")("gravcov",
        po::value<std::string>(),
        "A file containing the model covariance only for gravity")("mtcov",
        po::value<std::string>(), "A file containing the model covariance only for MT")(
        "geometrygrid", po::value<std::string>(),
        "A file containing the model geometry of the inversion grid")("tempdir",
        po::value<std::string>(), "The name of the directory to store temporary files in")(
        "wavelet", "Parametrize inversion by wavelet coefficients")("xorigin",
        po::value(&xorigin)->default_value(0.0),
        "The origin for the inversion grid in x-direction")("yorigin",
        po::value(&yorigin)->default_value(0.0),
        "The origin for the inversion grid in y-direction")("coolingfactor",
        po::value(&coolingfactor)->default_value(1.0),
        "The factor to multiply the weight for the regularization at each iteration EXPERIMENTAL")(
        "sequential",
        "Do not create a single objective function, but split into on OF per method EXPERIMENTAL")(
        "saveinterval", po::value(&saveinterval)->default_value(1),
        "The interval in iterations at which intermediate models are saved.")("stochcov",
        po::value(&CovWidth)->default_value(0),
        "Width of stochastic regularization, enabled if > 0")("stoch_nu",
        po::value(&Cov_nu)->default_value(1.00),
        "The value of nu for the stochastic covariance")("stoch_a",
        po::value(&Cov_sigma)->default_value(1.0),
        "The value of a for the stochastic covariance");
//we need to add the description for each part to the boost program options object
//that way the user can get a help output and the parser object recongnizes these options
    desc.add(TomoSetup.SetupOptions());
    desc.add(GravitySetup.SetupOptions());
    desc.add(MTSetup.SetupOptions());
    desc.add(InversionSetup.SetupOptions());
    desc.add(RegSetup.SetupOptions());
    desc.add(CouplingSetup.SetupOptions());

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
/* @} */
