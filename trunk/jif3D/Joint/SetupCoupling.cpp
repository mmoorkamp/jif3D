//============================================================================
// Name        : SetupCoupling.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupCoupling.h"
#include "../Global/FileUtil.h"
#include "../Regularization/CrossGradient.h"
#include "../Inversion/ModelTransforms.h"
#include "../MI/MutualInformationConstraint.h"
#include "SaltRelConstraint.h"
#include <boost/make_shared.hpp>

namespace ublas = boost::numeric::ublas;

namespace jif3D
  {
    void SetupModelCovar(jif3D::rvec &Covar, const jif3D::rvec &InvModel,
        const std::vector<double> &OldCov, size_t ngrid)
      {
        assert(Covar.size() == 3 * ngrid);
        assert(InvModel.size() == ngrid);
        if (OldCov.empty())
          {
            ublas::subrange(Covar, 0, ngrid) = InvModel;
            ublas::subrange(Covar, ngrid, 2 * ngrid) = InvModel;
            ublas::subrange(Covar, 2 * ngrid, 3 * ngrid) = InvModel;
          }
        else
          {
            assert(OldCov.size() == 3 * ngrid);
            std::transform(InvModel.begin(), InvModel.end(), OldCov.begin(),
                Covar.begin(), std::multiplies<double>());
            std::transform(InvModel.begin(), InvModel.end(), OldCov.begin() + ngrid,
                Covar.begin() + ngrid, std::multiplies<double>());
            std::transform(InvModel.begin(), InvModel.end(), OldCov.begin() + 2 * ngrid,
                Covar.begin() + 2 * ngrid, std::multiplies<double>());

          }
      }

    SetupCoupling::SetupCoupling() :
        mincond(0.0), maxcond(0.0), minslow(0.0), maxslow(0.0), mindens(0.0), maxdens(
            0.0), density_a(0.0), density_b(0.0), cond_a(0.0), cond_b(0.0), cond_c(0.0), vpvs(
            std::sqrt(3.0)), DensReplace(0.0), CondReplace(3.0)
      {
      }

    SetupCoupling::~SetupCoupling()
      {
      }

    po::options_description SetupCoupling::SetupOptions()
      {
        po::options_description desc("Coupling options");
        desc.add_options()("crossgrad", "Use cross-gradient coupling")(
            "mutual_information", po::value(&mibins)->default_value(50),
            "Use a MI based coupling constraint")("minslow",
            po::value(&minslow)->default_value(1e-4))("maxslow",
            po::value(&maxslow)->default_value(0.005))("minvel",
            po::value(&minvel)->default_value(1000))("maxvel",
            po::value(&maxvel)->default_value(8000))("mincond",
            po::value(&mincond)->default_value(1e-6))("maxcond",
            po::value(&maxcond)->default_value(5))("mindens",
            po::value(&mindens)->default_value(-500.0))("maxdens",
            po::value(&maxdens)->default_value(500.0))("vpvs",
            po::value(&vpvs)->default_value(std::sqrt(3.0)))("density_a",
            po::value(&density_a)->default_value(5000),
            "The slope of the velocity-density relationship")("density_b",
            po::value(&density_b)->default_value(8500),
            "The offset of the velocity-density relationship")("cond_a",
            po::value(&cond_a)->default_value(2.31e-7),
            "The quadratic term of the velocity-conductivity relationship")("cond_b",
            po::value(&cond_b)->default_value(-5.79e-4),
            "The linear term of the velocity-conductivity relationship")("cond_c",
            po::value(&cond_c)->default_value(0.124),
            "The constant term of the velocity-conductivity relationship")("relmodel",
            po::value(&RelModelName),
            "Name of a model file that specifies where to apply the parameter relationship")(
            "DensReplace", po::value(&DensReplace)->default_value(0.0),
            "Value to use for Density where relationship is not valid")("CondReplace",
            po::value(&CondReplace)->default_value(3.3),
            "Value to use for Conductivity where relationship is not valid")("seisano",
            "Use only the seismic anomaly for coupling");

        return desc;
      }

    void SetupCoupling::SetupTransforms(const po::variables_map &vm,
        ThreeDModelBase &GeometryModel,
        boost::shared_ptr<jif3D::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jif3D::GeneralModelTransform> &GravityTransform,
        boost::shared_ptr<jif3D::GeneralModelTransform> &MTTransform, bool Wavelet)
      {
        seisano = vm.count("seisano");

        const std::size_t ngrid = GeometryModel.GetNModelElements();

        //we declare WaveletTrans here so that we can use identical
        //transform for the different physical properties below
        boost::shared_ptr<jif3D::GeneralModelTransform> WaveletTrans;
        //but we only assign something to it if we actually want to
        //work in the wavelet domain as the constructor checks that
        //the grid size is a power of two, but this is not required otherwise
        if (Wavelet)
          {
            WaveletTrans = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::WaveletModelTransform(GeometryModel.GetData().shape()[0],
                    GeometryModel.GetData().shape()[1],
                    GeometryModel.GetData().shape()[2]));
          }
        //if we want to do a cross-gradient type joint inversion
        //we need to set transformations for each data type
        //that extract the right part of the model vector
        //and then transform from generalized to physical parameters

        boost::shared_ptr<jif3D::GeneralModelTransform> Copier(
            new jif3D::ModelCopyTransform);

        //each set of transformations is chained together in a similar way
        //the section transform takes the right part of the model vector
        //the TanHTransform sorts out the constraints
        //and the ExpansionTransform puts the derivatives in the right part of the whole model vector

        //we declare a local variable, so that we can use all properties
        //of ChainedTransform before we assign it to the function parameter
        //of type GeneralTransform
        //we start with density, because we also need it for surface waves
        boost::shared_ptr<jif3D::ChainedTransform> DensityTransform(
            new jif3D::ChainedTransform);
        DensityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(mindens, maxdens)));
        DensCrossTrans = Copier;
        GravityTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::MultiSectionTransform(3 * ngrid, ngrid, 2 * ngrid,
                DensityTransform));
        DensRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
            new jif3D::ChainedTransform);
        DensRegTrans->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::MultiSectionTransform(3 * ngrid, ngrid, 2 * ngrid, Copier)));

        //then we do velocity, for surface waves we need both velocity and density
        boost::shared_ptr<jif3D::GeneralModelTransform> SVelTrans = boost::make_shared<
            jif3D::TanhTransform>(minvel, maxvel);
        boost::shared_ptr<jif3D::ChainedTransform> SurfTrans(new jif3D::ChainedTransform);
        boost::shared_ptr<jif3D::MultiSectionTransform> SDensTrans = boost::make_shared<
            jif3D::MultiSectionTransform>(3 * ngrid);
        SDensTrans->AddSection(0, ngrid, SVelTrans);
        SDensTrans->AddSection(ngrid, 2 * ngrid, DensityTransform);
        SurfTrans->AppendTransform(SDensTrans);
        boost::shared_ptr<jif3D::DensPVelTransform> DensPVel = boost::make_shared<
            jif3D::DensPVelTransform>(std::vector<double>(), vpvs);
        SurfTrans->AppendTransform(DensPVel);

        boost::shared_ptr<jif3D::ChainedTransform> SlownessTransform(
            new jif3D::ChainedTransform);
        SlownessTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(minslow, maxslow)));
        SlowCrossTrans = Copier;
        //TomoTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
        //    new jif3D::MultiSectionTransform(3 * ngrid, 0, ngrid, SlownessTransform));
        //we regularize on the raw model parameters as these are more evenly spaced than slowness values
        SlowRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
            new jif3D::ChainedTransform);
        SlowRegTrans->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::MultiSectionTransform(3 * ngrid, 0, ngrid, Copier)));
#ifdef HAVETRAVELTIME
            TomoTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                            new jif3D::MultiSectionTransform(3 * ngrid, 0, ngrid, SlownessTransform));
#else
        TomoTransform = SurfTrans;

#endif
        //and then conductivity, conductivity has a LogTransform in addition
        //to reduce the range of inversion parameters
        boost::shared_ptr<jif3D::ChainedTransform> ConductivityTransform(
            new jif3D::ChainedTransform);
        jif3D::rvec RefModel(ngrid);
        std::fill(RefModel.begin(), RefModel.end(), 1.0);
        ConductivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(std::log(mincond), std::log(maxcond))));
        ConductivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::LogTransform(RefModel)));
        CondCrossTrans = Copier;
        MTTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::MultiSectionTransform(3 * ngrid, 2 * ngrid, 3 * ngrid,
                ConductivityTransform));
        //we regularize on the raw model parameters as these are more evenly spaced than conductivities
        CondRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
            new jif3D::ChainedTransform);
        CondRegTrans->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::MultiSectionTransform(3 * ngrid, 2 * ngrid, 3 * ngrid,
                    Copier)));
        //if we want to regularize in the wavelet  domain
        //we need to add a wavelet transform to the regularization
        if (Wavelet)
          {
            SlowRegTrans->AppendTransform(WaveletTrans);
            DensRegTrans->AppendTransform(WaveletTrans);
            CondRegTrans->AppendTransform(WaveletTrans);
          }

        //finished setting up cross-gradient coupling

        SlowTrans = TomoTransform;
        CondTrans = MTTransform;
        DensTrans = GravityTransform;
      }

    void SetupCoupling::SetupCrossGradModel(jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry, const SeisModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart,
        const jif3D::ThreeDModelBase &TearModelX,
        const jif3D::ThreeDModelBase &TearModelY,
        const jif3D::ThreeDModelBase &TearModelZ, const jif3D::rvec &CovVec)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);
#ifdef HAVETRAVELTIME
        jif3D::rvec SeisModel(3 * ngrid, 0.0);
#else
        jif3D::rvec SeisModel(3 * ngrid, 0.0);
#endif

        if (SeisMod.GetNModelElements() == ngrid)
          {
            std::copy(SeisMod.GetData().origin(), SeisMod.GetData().origin() + ngrid,
                SeisModel.begin());
            std::cout << "Transforming slowness model. " << std::endl;
            ublas::subrange(InvModel, 0, ngrid) = ublas::subrange(
                SlowTrans->PhysicalToGeneralized(SeisModel), 0, ngrid);
          }

        jif3D::rvec GravModel(ngrid, 1.0);
        if (GravMod.GetNModelElements() == ngrid)
          {
            std::copy(GravMod.GetDensities().origin(),
                GravMod.GetDensities().origin() + ngrid, GravModel.begin());
            std::cout << "Transforming Density model. " << std::endl;
            ublas::subrange(InvModel, ngrid, 2 * ngrid) = ublas::subrange(
                DensTrans->PhysicalToGeneralized(GravModel), ngrid, 2 * ngrid);
          }

        jif3D::rvec MTModel(ngrid, 1.0);
        if (MTMod.GetNModelElements() == ngrid)
          {

            std::copy(MTMod.GetConductivities().origin(),
                MTMod.GetConductivities().origin() + ngrid, MTModel.begin());
            //make sure that each layer has a slightly different conductivity
            //at least in one cell, otherwise x3d optimizes by joining layers
            //and messes up the gradient calculation
            for (size_t i = 0; i < MTMod.GetConductivities().shape()[2]; ++i)
              {
                MTModel(MTMod.IndexToOffset(0, 0, i)) *= (1 + 0.0001 * (i + 1));
              }
            std::cout << "Transforming conductivity model. " << std::endl;
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(
                CondTrans->PhysicalToGeneralized(MTModel), 2 * ngrid, 3 * ngrid);
          }
        //then we construct the three cross gradient terms
        //the double section transform takes two sections of the model
        //and feeds them to the objective function
        boost::shared_ptr<jif3D::CrossGradient> SeisGravCross(
            new jif3D::CrossGradient(ModelGeometry, TearModelX, TearModelY, TearModelZ));
        boost::shared_ptr<jif3D::MultiSectionTransform> SeisGravTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        SeisGravTrans->AddSection(0, ngrid, SlowCrossTrans);
        SeisGravTrans->AddSection(ngrid, 2 * ngrid, DensCrossTrans);

        //for each cross-gradient term we ask for a weight
        double seisgravlambda = 1.0;
        std::cout << "Weight for seismic-gravity cross-gradient term: ";
        std::cin >> seisgravlambda;
        if (seisgravlambda > 0.0)
          {
            Objective.AddObjective(SeisGravCross, SeisGravTrans, seisgravlambda,
                "SeisGrav", JointObjective::coupling);
          }
        boost::shared_ptr<jif3D::CrossGradient> SeisMTCross(
            new jif3D::CrossGradient(ModelGeometry, TearModelX, TearModelY, TearModelZ));
        boost::shared_ptr<jif3D::MultiSectionTransform> SeisMTTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        SeisMTTrans->AddSection(0, ngrid, SlowCrossTrans);
        SeisMTTrans->AddSection(2 * ngrid, 3 * ngrid, CondCrossTrans);

        double seismtlambda = 1.0;
        std::cout << "Weight for seismic-MT cross-gradient term: ";
        std::cin >> seismtlambda;
        if (seismtlambda > 0.0)
          {
            Objective.AddObjective(SeisMTCross, SeisMTTrans, seismtlambda, "SeisMT",
                JointObjective::coupling);
          }
        boost::shared_ptr<jif3D::CrossGradient> GravMTCross(
            new jif3D::CrossGradient(ModelGeometry, TearModelX, TearModelY, TearModelZ));
        boost::shared_ptr<jif3D::MultiSectionTransform> GravMTTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        GravMTTrans->AddSection(ngrid, 2 * ngrid, DensCrossTrans);
        GravMTTrans->AddSection(2 * ngrid, 3 * ngrid, CondCrossTrans);

        double gravmtlambda = 1.0;
        std::cout << "Weight for gravity-MT cross-gradient term: ";
        std::cin >> gravmtlambda;
        if (gravmtlambda > 0.0)
          {
            Objective.AddObjective(GravMTCross, GravMTTrans, gravmtlambda, "GravMT",
                JointObjective::coupling);
          }
        //finally we construct the regularization terms
        //we ask for a weight and construct a regularization object
        //for each type of physical parameter separately
        //first we set up seismic tomography

        jif3D::rvec TomoCovVec(ngrid, 1.0);
        if (!CovVec.empty())
          {
            if (CovVec.size() == ngrid)
              {
                TomoCovVec = CovVec;
              }
            else
              {
                TomoCovVec = ublas::subrange(CovVec, 0, ngrid);
              }
          }
        double seisreglambda = 1.0;
        std::cout << " Weight for seismic regularization: ";
        std::cin >> seisreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> SeisReg(Regularization->clone());
        jif3D::rvec TomoCovar(3 * ngrid);
        SetupModelCovar(TomoCovar, TomoCovVec, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(std::vector<double>(TomoCovar.begin(), TomoCovar.end()));

        //then the regularization of densities

        jif3D::rvec GravCovVec(ngrid, 1.0);
        if (!CovVec.empty())
          {
            if (CovVec.size() == ngrid)
              {
                GravCovVec = CovVec;
              }
            else
              {
                GravCovVec = ublas::subrange(CovVec, ngrid, 2 * ngrid);
              }
          }
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> GravReg(Regularization->clone());
        jif3D::rvec GravCovar(3 * ngrid);

        SetupModelCovar(GravCovar, GravCovVec, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(std::vector<double>(GravCovar.begin(), GravCovar.end()));

        //and finally conductivities
        jif3D::rvec CondCovVec(ngrid, 1.0);
        if (!CovVec.empty())
          {
            if (CovVec.size() == ngrid)
              {
                CondCovVec = CovVec;
              }
            else
              {
                CondCovVec = ublas::subrange(CovVec, 2 * ngrid, 3 * ngrid);
              }
          }
        double mtreglambda = 1.0;
        std::cout << " Weight for MT regularization: ";
        std::cin >> mtreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> MTReg(Regularization->clone());
        jif3D::rvec MTCovar(3 * ngrid);
        SetupModelCovar(MTCovar, CondCovVec, MTReg->GetDataError(), ngrid);
        MTReg->SetDataError(std::vector<double>(MTCovar.begin(), MTCovar.end()));
        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
        if (substart)
          {
            SeisReg->SetReferenceModel(SlowRegTrans->GeneralizedToPhysical(InvModel));
            GravReg->SetReferenceModel(DensRegTrans->GeneralizedToPhysical(InvModel));
            MTReg->SetReferenceModel(CondRegTrans->GeneralizedToPhysical(InvModel));
          }
        if (seisreglambda > 0.0)
          {
            Objective.AddObjective(SeisReg, SlowRegTrans, seisreglambda, "SeisReg",
                JointObjective::regularization);
          }
        if (gravreglambda > 0.0)
          {
            Objective.AddObjective(GravReg, DensRegTrans, gravreglambda, "GravReg",
                JointObjective::regularization);
          }
        if (mtreglambda > 0.0)
          {
            Objective.AddObjective(MTReg, CondRegTrans, mtreglambda, "MTReg",
                JointObjective::regularization);
          }
      }

    void SetupCoupling::SetupMIModel(jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry, const SeisModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);

        jif3D::rvec SeisModel(3 * ngrid, 0.0);
        if (SeisMod.GetNModelElements() == ngrid)
          {
            std::copy(SeisMod.GetData().origin(), SeisMod.GetData().origin() + ngrid,
                SeisModel.begin());
            std::cout << "Transforming slowness model. " << std::endl;
            ublas::subrange(InvModel, 0, ngrid) = ublas::subrange(
                SlowTrans->PhysicalToGeneralized(SeisModel), 0, ngrid);
          }

        jif3D::rvec GravModel(ngrid, 1.0);
        if (GravMod.GetNModelElements() == ngrid)
          {
            std::copy(GravMod.GetDensities().origin(),
                GravMod.GetDensities().origin() + ngrid, GravModel.begin());
            std::cout << "Transforming Density model. " << std::endl;
            ublas::subrange(InvModel, ngrid, 2 * ngrid) = ublas::subrange(
                DensTrans->PhysicalToGeneralized(GravModel), ngrid, 2 * ngrid);
          }

        jif3D::rvec MTModel(ngrid, 1.0);
        if (MTMod.GetNModelElements() == ngrid)
          {

            std::copy(MTMod.GetConductivities().origin(),
                MTMod.GetConductivities().origin() + ngrid, MTModel.begin());
            //make sure that each layer has a slightly different conductivity
            //at least in one cell, otherwise x3d optimizes by joining layers
            //and messes up the gradient calculation
            for (size_t i = 0; i < MTMod.GetConductivities().shape()[2]; ++i)
              {
                MTModel(MTMod.IndexToOffset(0, 0, i)) *= (1 + 0.0001 * (i + 1));
              }
            std::cout << "Transforming conductivity model. " << std::endl;
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(
                CondTrans->PhysicalToGeneralized(MTModel), 2 * ngrid, 3 * ngrid);
          }
        boost::shared_ptr<jif3D::GeneralModelTransform> Copier(
            new jif3D::ModelCopyTransform);

        boost::shared_ptr<jif3D::GeneralModelTransform> Anomaly(
            new jif3D::Anomalator(InvModel));

        //then we construct the three cross gradient terms
        //the double section transform takes two sections of the model
        //and feeds them to the objective function
        boost::shared_ptr<jif3D::MutualInformationConstraint> SeisGravMI =
            boost::make_shared<jif3D::MutualInformationConstraint>(-2.0, 2.0, -2.0, 2.0,
                mibins);
        boost::shared_ptr<jif3D::MultiSectionTransform> SeisGravTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        if (seisano)
          {
            SeisGravTrans->AddSection(0, ngrid, Anomaly);
          }
        else
          {
            SeisGravTrans->AddSection(0, ngrid, Copier);
          }
        SeisGravTrans->AddSection(ngrid, 2 * ngrid, Copier);

        //for each cross-gradient term we ask for a weight
        double seisgravlambda = 1.0;
        std::cout << "Weight for seismic-gravity MI term: ";
        std::cin >> seisgravlambda;
        if (seisgravlambda > 0.0)
          {
            Objective.AddObjective(SeisGravMI, SeisGravTrans, seisgravlambda, "SeisGrav",
                JointObjective::coupling);
          }
        boost::shared_ptr<jif3D::MutualInformationConstraint> SeisMTMI =
            boost::make_shared<jif3D::MutualInformationConstraint>(-2.0, 2.0, -2.0, 2.0,
                mibins);

        boost::shared_ptr<jif3D::MultiSectionTransform> SeisMTTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        if (seisano)
          {
            SeisGravTrans->AddSection(0, ngrid, Anomaly);
          }
        else
          {
            SeisGravTrans->AddSection(0, ngrid, Copier);
          }
        SeisMTTrans->AddSection(2 * ngrid, 3 * ngrid, Copier);

        double seismtlambda = 1.0;
        std::cout << "Weight for seismic-MT MI term: ";
        std::cin >> seismtlambda;
        if (seismtlambda > 0.0)
          {
            Objective.AddObjective(SeisMTMI, SeisMTTrans, seismtlambda, "SeisMT",
                JointObjective::coupling);
          }
        boost::shared_ptr<jif3D::MutualInformationConstraint> GravMTMI =
            boost::make_shared<jif3D::MutualInformationConstraint>(-2.0, 2.0, -2.0, 2.0,
                mibins);
        boost::shared_ptr<jif3D::MultiSectionTransform> GravMTTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        GravMTTrans->AddSection(ngrid, 2 * ngrid, Copier);
        GravMTTrans->AddSection(2 * ngrid, 3 * ngrid, Copier);

        double gravmtlambda = 1.0;
        std::cout << "Weight for gravity-MT MI term: ";
        std::cin >> gravmtlambda;
        if (gravmtlambda > 0.0)
          {
            Objective.AddObjective(GravMTMI, GravMTTrans, gravmtlambda, "GravMT",
                JointObjective::coupling);
          }
        //finally we construct the regularization terms
        //we ask for a weight and construct a regularization object
        //for each type of physical parameter separately
        //first we set up seismic tomography

        jif3D::rvec TomoCovVec(ngrid, 1.0);
        double seisreglambda = 1.0;
        std::cout << " Weight for seismic regularization: ";
        std::cin >> seisreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> SeisReg(Regularization->clone());

        //then the regularization of densities

        jif3D::rvec GravCovVec(ngrid, 1.0);
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> GravReg(Regularization->clone());

        //and finally conductivities
        jif3D::rvec CondCovVec(ngrid, 1.0);
        double mtreglambda = 1.0;
        std::cout << " Weight for MT regularization: ";
        std::cin >> mtreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> MTReg(Regularization->clone());

        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
        if (substart)
          {
            SeisReg->SetReferenceModel(SlowRegTrans->GeneralizedToPhysical(InvModel));
            GravReg->SetReferenceModel(DensRegTrans->GeneralizedToPhysical(InvModel));
            MTReg->SetReferenceModel(CondRegTrans->GeneralizedToPhysical(InvModel));
          }
        if (seisreglambda > 0.0)
          {
            Objective.AddObjective(SeisReg, SlowRegTrans, seisreglambda, "SeisReg",
                JointObjective::regularization);
          }
        if (gravreglambda > 0.0)
          {
            Objective.AddObjective(GravReg, DensRegTrans, gravreglambda, "GravReg",
                JointObjective::regularization);
          }
        if (mtreglambda > 0.0)
          {
            Objective.AddObjective(MTReg, CondRegTrans, mtreglambda, "MTReg",
                JointObjective::regularization);
          }
      }

    void SetupCoupling::SetupModelVector(const po::variables_map &vm,
        jif3D::rvec &InvModel, const jif3D::ThreeDModelBase &ModelGeometry,
        const SeisModel &SeisMod, const jif3D::ThreeDGravityModel &GravMod,
        const jif3D::ThreeDMTModel &MTMod, jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart,
        const jif3D::ThreeDModelBase &TearModelX,
        const jif3D::ThreeDModelBase &TearModelY,
        const jif3D::ThreeDModelBase &TearModelZ, const jif3D::rvec &CovVec)
      {
        //depending on the type of coupling the model vector looks quite different
        //and we have to setup a number of different things
        //so we separate it into different functions even though it means
        //passing around many parameters

        if (vm.count("crossgrad"))
          {
            SetupCrossGradModel(InvModel, ModelGeometry, SeisMod, GravMod, MTMod,
                Objective, Regularization, substart, TearModelX, TearModelY, TearModelZ,
                CovVec);
          }
        else
          {
            SetupMIModel(InvModel, ModelGeometry, SeisMod, GravMod, MTMod, Objective,
                Regularization, substart);

          }

      }

  }
