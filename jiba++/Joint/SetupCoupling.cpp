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
#include "SaltRelConstraint.h"

namespace ublas = boost::numeric::ublas;

namespace jiba
  {
    void SetupModelCovar(jiba::rvec &Covar, const jiba::rvec &InvModel,
        const jiba::rvec &OldCov, size_t ngrid)
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
            ublas::subrange(Covar, 0, ngrid) = ublas::element_prod(InvModel,
                ublas::subrange(OldCov, 0, ngrid));
            ublas::subrange(Covar, ngrid, 2 * ngrid) = ublas::element_prod(InvModel,
                ublas::subrange(OldCov, ngrid, 2 * ngrid));
            ublas::subrange(Covar, 2 * ngrid, 3 * ngrid) = ublas::element_prod(InvModel,
                ublas::subrange(OldCov, 2 * ngrid, 3 * ngrid));
          }
      }

    SetupCoupling::SetupCoupling() :
        mincond(0.0), maxcond(0.0), minslow(0.0), maxslow(0.0), mindens(0.0), maxdens(
            0.0), density_a(0.0), density_b(0.0), cond_a(0.0), cond_b(0.0), cond_c(0.0), DensReplace(
            0.0), CondReplace(3.0)
      {
      }

    SetupCoupling::~SetupCoupling()
      {
      }

    po::options_description SetupCoupling::SetupOptions()
      {
        po::options_description desc("Coupling options");
        desc.add_options()("crossgrad", "Use cross-gradient coupling")("saltrel",
            "Use a parameter constraint designed for salt")("minslow",
            po::value(&minslow)->default_value(1e-4))("maxslow",
            po::value(&maxslow)->default_value(0.005))("mincond",
            po::value(&mincond)->default_value(1e-4))("maxcond",
            po::value(&maxcond)->default_value(5))("mindens",
            po::value(&mindens)->default_value(0.5))("maxdens",
            po::value(&maxdens)->default_value(4.0))("density_a",
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
            "Value to use for Conductivity where relationship is not valid");

        return desc;
      }

    void SetupCoupling::SetupTransforms(const po::variables_map &vm,
        ThreeDSeismicModel &StartModel,
        boost::shared_ptr<jiba::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &GravityTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &MTTransform, bool Wavelet)
      {

        //we need the geometry of the starting model to setup
        //the transformations
        std::string modelfilename = jiba::AskFilename("Inversion Model Geometry: ");
        //although we specify the starting model as a seismic model
        //it does not need to have the same grid specifications
        StartModel.ReadNetCDF(modelfilename, false);

        const std::size_t ngrid = StartModel.GetNModelElements();

        //we declare WaveletTrans here so that we can use identical
        //transform for the different physical properties below
        boost::shared_ptr<jiba::GeneralModelTransform> WaveletTrans;
        //but we only assign something to it if we actually want to
        //work in the wavelet domain as the constructor checks that
        //the grid size is a power of two, but this is not required otherwise
        if (Wavelet)
          {
            WaveletTrans = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::WaveletModelTransform(StartModel.GetData().shape()[0],
                    StartModel.GetData().shape()[1], StartModel.GetData().shape()[2]));
          }
        //if we want to do a cross-gradient type joint inversion
        //we need to set transformations for each data type
        //that extract the right part of the model vector
        //and then transform from generalized to physical parameters

        if (vm.count("crossgrad") || vm.count("saltrel"))
          {
            //each set of transformations is chained together in a similar way
            //the section transform takes the right part of the model vector
            //the TanHTransform sorts out the constraints
            //and the ExpansionTransform puts the derivatives in the right part of the whole model vector
            //first we do slowness

            //we declare a local variable, so that we can use all properties
            //of ChainedTransform before we assign it to the function parameter
            //of type GeneralTransform
            //we start with tomography
            boost::shared_ptr<jiba::ChainedTransform> SlownessTransform(
                new jiba::ChainedTransform);

            SlownessTransform->AppendTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(minslow, maxslow)));
            SlowCrossTrans = boost::shared_ptr<jiba::GeneralModelTransform>(
                SlownessTransform->clone());
            SlowRegTrans = boost::shared_ptr<jiba::ChainedTransform>(
                SlownessTransform->clone());
            TomoTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::MultiSectionTransform(3 * ngrid, 0, ngrid, SlownessTransform));

            //then we do density
            boost::shared_ptr<jiba::ChainedTransform> DensityTransform(
                new jiba::ChainedTransform);
            DensityTransform->AppendTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(mindens, maxdens)));
            DensCrossTrans = boost::shared_ptr<jiba::GeneralModelTransform>(
                DensityTransform->clone());
            DensRegTrans = boost::shared_ptr<jiba::ChainedTransform>(
                DensityTransform->clone());
            GravityTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::MultiSectionTransform(3 * ngrid, ngrid, 2 * ngrid,
                    DensityTransform));

            //and then conductivity, conductivity has a LogTransform in addition
            //to reduce the range of inversion parameters
            boost::shared_ptr<jiba::ChainedTransform> ConductivityTransform(
                new jiba::ChainedTransform);
            jiba::rvec RefModel(ngrid);
            std::fill(RefModel.begin(), RefModel.end(), 1.0);
            ConductivityTransform->AppendTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(std::log(mincond), std::log(maxcond))));
            ConductivityTransform->AppendTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::LogTransform(RefModel)));
            CondCrossTrans = boost::shared_ptr<jiba::GeneralModelTransform>(
                ConductivityTransform->clone());
            CondRegTrans = boost::shared_ptr<jiba::ChainedTransform>(
                ConductivityTransform->clone());
            MTTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::MultiSectionTransform(3 * ngrid, 2 * ngrid, 3 * ngrid,
                    ConductivityTransform));
            //if we want to regularize in the wavelet  domain
            //we need to add a wavelet transform to the regularization
            if (Wavelet)
              {
                SlowRegTrans->AppendTransform(WaveletTrans);
                DensRegTrans->AppendTransform(WaveletTrans);
                CondRegTrans->AppendTransform(WaveletTrans);
              }
            //finished setting up cross-gradient coupling

          }
        else
          {
            // We can specify a model file where we want
            // the parameter relationship to be valid
            jiba::ThreeDSeismicModel RelModel;
            if (vm.count("relmodel"))
              {
                RelModel.ReadNetCDF(RelModelName, false);
              }
            else
              {
                const size_t nx = StartModel.GetSlownesses().shape()[0];
                const size_t ny = StartModel.GetSlownesses().shape()[1];
                const size_t nz = StartModel.GetSlownesses().shape()[2];
                RelModel.SetSlownesses().resize(boost::extents[nx][ny][nz]);
                std::fill_n(RelModel.SetSlownesses().origin(), nx * ny * nz, 1.0);

              }
            //if we want direct parameter coupling we do not need the section transforms
            //but we directly calculate conductivity and density from slowness

            TomoTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::TanhTransform(minslow, maxslow));

            GravityTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::DensityTransform(TomoTransform, RelModel, DensReplace,
                    density_a, density_b));
            MTTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::ConductivityTransform(TomoTransform, RelModel, CondReplace,
                    cond_a, cond_b, cond_c));

            if (Wavelet)
              {

                boost::shared_ptr<jiba::ChainedTransform> SlownessTransform(
                    new jiba::ChainedTransform);

                SlownessTransform->AppendTransform(
                    boost::shared_ptr<jiba::GeneralModelTransform>(
                        new jiba::TanhTransform(minslow, maxslow)));
                SlownessTransform->AppendTransform(WaveletTrans);

                SlowRegTrans = SlownessTransform;

              }

          }
        SlowTrans = TomoTransform;
        CondTrans = MTTransform;
        DensTrans = GravityTransform;
      }

    void SetupCoupling::SetupCrossGradModel(jiba::rvec &InvModel,
        const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDSeismicModel &SeisMod, const jiba::ThreeDGravityModel &GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::ObjectiveFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);

        jiba::rvec SeisModel(ngrid, 0.0);
        if (SeisMod.GetNModelElements() == ngrid)
          {
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, SeisModel.begin());
            std::cout << "Transforming slowness model. " << std::endl;
            ublas::subrange(InvModel, 0, ngrid) = ublas::subrange(
                SlowTrans->PhysicalToGeneralized(SeisModel), 0, ngrid);
          }

        jiba::rvec GravModel(ngrid, 1.0);
        if (GravMod.GetNModelElements() == ngrid)
          {
            std::copy(GravMod.GetDensities().origin(),
                GravMod.GetDensities().origin() + ngrid, GravModel.begin());
            std::cout << "Transforming Density model. " << std::endl;
            ublas::subrange(InvModel, ngrid, 2 * ngrid) = ublas::subrange(
                DensTrans->PhysicalToGeneralized(GravModel), ngrid, 2 * ngrid);
          }

        jiba::rvec MTModel(ngrid, 1.0);
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
        boost::shared_ptr<jiba::CrossGradient> SeisGravCross(
            new jiba::CrossGradient(ModelGeometry));
        boost::shared_ptr<jiba::MultiSectionTransform> SeisGravTrans(
            new jiba::MultiSectionTransform(3 * ngrid));
        SeisGravTrans->AddSection(0, ngrid, SlowCrossTrans);
        SeisGravTrans->AddSection(ngrid, 2 * ngrid, DensCrossTrans);

        //for each cross-gradient term we ask for a weight
        double seisgravlambda = 1.0;
        std::cout << "Weight for seismic-gravity cross-gradient term: ";
        std::cin >> seisgravlambda;
        if (seisgravlambda > 0.0)
          {
            Objective.AddObjective(SeisGravCross, SeisGravTrans, seisgravlambda,
                "SeisGrav");
          }
        boost::shared_ptr<jiba::CrossGradient> SeisMTCross(
            new jiba::CrossGradient(ModelGeometry));
        boost::shared_ptr<jiba::MultiSectionTransform> SeisMTTrans(
            new jiba::MultiSectionTransform(3 * ngrid));
        SeisMTTrans->AddSection(0, ngrid, SlowCrossTrans);
        SeisMTTrans->AddSection(2 * ngrid, 3 * ngrid, CondCrossTrans);

        double seismtlambda = 1.0;
        std::cout << "Weight for seismic-MT cross-gradient term: ";
        std::cin >> seismtlambda;
        if (seismtlambda > 0.0)
          {
            Objective.AddObjective(SeisMTCross, SeisMTTrans, seismtlambda, "SeisMT");
          }
        boost::shared_ptr<jiba::CrossGradient> GravMTCross(
            new jiba::CrossGradient(ModelGeometry));
        boost::shared_ptr<jiba::MultiSectionTransform> GravMTTrans(
            new jiba::MultiSectionTransform(3 * ngrid));
        GravMTTrans->AddSection(ngrid, 2 * ngrid, DensCrossTrans);
        GravMTTrans->AddSection(2 * ngrid, 3 * ngrid, CondCrossTrans);

        double gravmtlambda = 1.0;
        std::cout << "Weight for gravity-MT cross-gradient term: ";
        std::cin >> gravmtlambda;
        if (gravmtlambda > 0.0)
          {
            Objective.AddObjective(GravMTCross, GravMTTrans, seismtlambda, "GravMT");
          }
        //finally we construct the regularization terms
        //we ask for a weight and construct a regularization object
        //for each type of physical parameter separately
        //first we set up seismic tomography
        double seisreglambda = 1.0;
        std::cout << " Weight for seismic regularization: ";
        std::cin >> seisreglambda;
        boost::shared_ptr<jiba::ObjectiveFunction> SeisReg(Regularization->clone());
        jiba::rvec TomoCovar(3 * ngrid);
        SetupModelCovar(TomoCovar, SeisModel, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(TomoCovar);

        //then the regularization of densities
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jiba::ObjectiveFunction> GravReg(Regularization->clone());
        jiba::rvec GravCovar(3 * ngrid);
        jiba::rvec Ones(GravModel.size(), 1.0);
        SetupModelCovar(GravCovar, Ones, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(GravCovar);

        //and finally conductivities
        double mtreglambda = 1.0;
        std::cout << " Weight for MT regularization: ";
        std::cin >> mtreglambda;
        boost::shared_ptr<jiba::ObjectiveFunction> MTReg(Regularization->clone());
        jiba::rvec MTCovar(3 * ngrid);
        SetupModelCovar(MTCovar, MTModel, MTReg->GetDataError(), ngrid);
        MTReg->SetDataError(MTCovar);
        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
//        if (substart)
//          {
//            SeisReg->SetReferenceModel(SeisModel);
//            GravReg->SetReferenceModel(GravModel);
//            MTReg->SetReferenceModel(MTModel);
//          }
        if (seisreglambda > 0.0)
          {
            Objective.AddObjective(SeisReg, SlowRegTrans, seisreglambda, "SeisReg");
          }
        if (gravreglambda > 0.0)
          {
            Objective.AddObjective(GravReg, DensRegTrans, gravreglambda, "GravReg");
          }
        if (mtreglambda > 0.0)
          {
            Objective.AddObjective(MTReg, CondRegTrans, mtreglambda, "MTReg");
          }
      }

    void SetupCoupling::SetupSaltModel(const po::variables_map &vm, jiba::rvec &InvModel,
        const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDSeismicModel &SeisMod, const jiba::ThreeDGravityModel &GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::ObjectiveFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);

        jiba::rvec SeisModel(ngrid, 0.0);
        if (SeisMod.GetNModelElements() == ngrid)
          {
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, SeisModel.begin());
            std::cout << "Transforming slowness model. " << std::endl;
            ublas::subrange(InvModel, 0, ngrid) = ublas::subrange(
                SlowTrans->PhysicalToGeneralized(SeisModel), 0, ngrid);
          }

        jiba::rvec GravModel(ngrid, 1.0);
        if (GravMod.GetNModelElements() == ngrid)
          {
            std::copy(GravMod.GetDensities().origin(),
                GravMod.GetDensities().origin() + ngrid, GravModel.begin());
            std::cout << "Transforming Density model. " << std::endl;
            ublas::subrange(InvModel, ngrid, 2 * ngrid) = ublas::subrange(
                DensTrans->PhysicalToGeneralized(GravModel), ngrid, 2 * ngrid);
          }
        else
          {
            std::fill_n(ublas::subrange(InvModel, ngrid, 2 * ngrid).begin(), ngrid, 0.0);
          }

        jiba::rvec MTModel(ngrid, 1.0);
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
        else
          {
            std::fill_n(ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid).begin(), ngrid,
                0.0);
          }

        jiba::ThreeDSeismicModel RelModel;
        if (vm.count("relmodel"))
          {
            RelModel.ReadNetCDF(RelModelName, false);
          }
        else
          {
            const size_t nx = ModelGeometry.GetData().shape()[0];
            const size_t ny = ModelGeometry.GetData().shape()[1];
            const size_t nz = ModelGeometry.GetData().shape()[2];
            RelModel.SetSlownesses().resize(boost::extents[nx][ny][nz]);
            std::fill_n(RelModel.SetSlownesses().origin(), nx * ny * nz, 1.0);

          }

        //we currently hacked out the density transform
        boost::shared_ptr<jiba::SaltRelConstraint> SaltRel(
            new jiba::SaltRelConstraint(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::ModelCopyTransform()),
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::ConductivityTransform(
                        boost::shared_ptr<jiba::GeneralModelTransform>(
                            new jiba::ModelCopyTransform()), RelModel))));

        SaltRel->SetExcludeCells(RelModel);
        boost::shared_ptr<jiba::MultiSectionTransform> SaltRelTrans(
            new jiba::MultiSectionTransform(3 * ngrid));
        SaltRelTrans->AddSection(0, ngrid, SlowCrossTrans);
        //evil HACK !
        SaltRelTrans->AddSection(0, ngrid, SlowCrossTrans);
        //SaltRelTrans->AddSection(ngrid, 2 * ngrid, DensCrossTrans);
        SaltRelTrans->AddSection(2 * ngrid, 3 * ngrid, CondCrossTrans);

        //for each cross-gradient term we ask for a weight
        double saltrellambda = 1.0;
        std::cout << "Weight for salt relationship term: ";
        std::cin >> saltrellambda;
        if (saltrellambda > 0.0)
          {
            Objective.AddObjective(SaltRel, SaltRelTrans, saltrellambda, "SaltRel");
          }
        //finally we construct the regularization terms
        //we ask for a weight and construct a regularization object
        //for each type of physical parameter separately
        //first we set up seismic tomography
        double seisreglambda = 1.0;
        std::cout << " Weight for seismic regularization: ";
        std::cin >> seisreglambda;
        boost::shared_ptr<jiba::ObjectiveFunction> SeisReg(Regularization->clone());
        jiba::rvec TomoCovar(3 * ngrid);
        SetupModelCovar(TomoCovar, SeisModel, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(TomoCovar);

        //then the regularization of densities
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jiba::ObjectiveFunction> GravReg(Regularization->clone());
        jiba::rvec GravCovar(3 * ngrid);
        jiba::rvec Ones(GravModel.size(), 1.0);
        SetupModelCovar(GravCovar, Ones, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(GravCovar);

        //and finally conductivities
        double mtreglambda = 1.0;
        std::cout << " Weight for MT regularization: ";
        std::cin >> mtreglambda;
        boost::shared_ptr<jiba::ObjectiveFunction> MTReg(Regularization->clone());
        jiba::rvec MTCovar(3 * ngrid);
        SetupModelCovar(MTCovar, MTModel, MTReg->GetDataError(), ngrid);
        MTReg->SetDataError(MTCovar);

        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
/*        if (substart)
          {
            SeisReg->SetReferenceModel(SeisModel);
            GravReg->SetReferenceModel(GravModel);
            MTReg->SetReferenceModel(MTModel);
          }*/
        if (seisreglambda > 0.0)
          {
            Objective.AddObjective(SeisReg, SlowRegTrans, seisreglambda, "SeisReg");

          }
        if (gravreglambda > 0.0)
          {
            Objective.AddObjective(GravReg, DensRegTrans, gravreglambda, "GravReg");

          }
        if (mtreglambda > 0.0)
          {
            Objective.AddObjective(MTReg, CondRegTrans, mtreglambda, "MTReg");
          }
      }

    void SetupCoupling::SetupFixedCouplingModel(jiba::rvec &InvModel,
        const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDSeismicModel &SeisMod, const jiba::ThreeDGravityModel &GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::ObjectiveFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();

        //if we want direct parameter coupling
        //we just need a single regularization term
        //as we invert for slowness, the regularization also works on slowness
        if (SeisMod.GetNModelElements() != ngrid)
          {
            throw jiba::FatalException(
                "Size of seismic model does not match model geometry !");
          }
        InvModel.resize(ngrid);
        std::copy(SeisMod.GetSlownesses().origin(),
            SeisMod.GetSlownesses().origin() + ngrid, InvModel.begin());
        double reglambda = 1.0;
        std::cout << " Weight for regularization: ";
        std::cin >> reglambda;
        if (reglambda > 0.0)
          {
            jiba::rvec TomoCovar(3 * ngrid);
            SetupModelCovar(TomoCovar, InvModel, Regularization->GetDataError(), ngrid);
            Regularization->SetDataError(TomoCovar);
          }
/*        if (substart)
          {
            Regularization->SetReferenceModel(InvModel);
          }*/
        Objective.AddObjective(Regularization, SlowTrans, reglambda, "Reg");
        InvModel = SlowTrans->PhysicalToGeneralized(InvModel);
      }

    void SetupCoupling::SetupModelVector(const po::variables_map &vm,
        jiba::rvec &InvModel, const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDSeismicModel &SeisMod, const jiba::ThreeDGravityModel &GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::ObjectiveFunction> Regularization, bool substart)
      {
        //if we want to do corss-gradient type inversion
        //we need three cross-gradient objective functions
        //and three regularization terms
        if (vm.count("crossgrad"))
          {
            SetupCrossGradModel(InvModel, ModelGeometry, SeisMod, GravMod, MTMod,
                Objective, Regularization, substart);
          }
        else
          {
            if (vm.count("saltrel"))
              {
                SetupSaltModel(vm, InvModel, ModelGeometry, SeisMod, GravMod, MTMod,
                    Objective, Regularization, substart);
              }
            else
              {
                SetupFixedCouplingModel(InvModel, ModelGeometry, SeisMod, GravMod, MTMod,
                    Objective, Regularization, substart);
              }
          }
      }
  }
