//============================================================================
// Name        : SetupCoupling.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupCoupling.h"
#include "../Global/FileUtil.h"
#include "../Regularization/CrossGradient.h"
#include "../Inversion/WaveletModelTransform.h"
#include "SaltRelConstraint.h"
#include "SaltRelTrans.h"

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

    SetupCoupling::SetupCoupling()
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
            po::value(&maxdens)->default_value(4.0));

        return desc;
      }

    void SetupCoupling::SetupTransforms(const po::variables_map &vm,
        ThreeDSeismicModel &StartModel,
        boost::shared_ptr<jiba::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &GravityTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &MTTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &RegTransform, bool Wavelet)
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
            boost::shared_ptr<jiba::ChainedTransform> SlownessTransform(
                new jiba::ChainedTransform);
            boost::shared_ptr<jiba::ChainedTransform> SlownessCrossTransform(
                new jiba::ChainedTransform);
            SlownessTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::SectionTransform(0, ngrid)));
            if (Wavelet)
              {
                SlownessTransform->AddTransform(WaveletTrans);
                SlownessCrossTransform->AddTransform(WaveletTrans);
              }
            SlownessTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(minslow, maxslow)));
            SlownessCrossTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(minslow, maxslow)));
            SlownessTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::ExpansionTransform(3 * ngrid, 0, ngrid)));
            //then we do density
            boost::shared_ptr<jiba::ChainedTransform> DensityTransform(
                new jiba::ChainedTransform);
            boost::shared_ptr<jiba::ChainedTransform> DensityCrossTransform(
                new jiba::ChainedTransform);
            DensityTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::SectionTransform(ngrid, 2 * ngrid)));
            if (Wavelet)
              {
                DensityTransform->AddTransform(WaveletTrans);
                DensityCrossTransform->AddTransform(WaveletTrans);
              }
            DensityTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(mindens, maxdens)));
            DensityCrossTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(mindens, maxdens)));

            DensityTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::ExpansionTransform(3 * ngrid, ngrid, 2 * ngrid)));
            //and then conductivity, conductivity has a LogTransform in addition
            //to reduce the range of inversion parameters
            boost::shared_ptr<jiba::ChainedTransform> ConductivityTransform(
                new jiba::ChainedTransform);
            boost::shared_ptr<jiba::ChainedTransform> ConductivityCrossTransform(
                new jiba::ChainedTransform);
            jiba::rvec RefModel(ngrid);
            std::fill(RefModel.begin(), RefModel.end(), 1.0);
            ConductivityTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::SectionTransform(2 * ngrid, 3 * ngrid)));
            if (Wavelet)
              {
                ConductivityTransform->AddTransform(WaveletTrans);
                ConductivityCrossTransform->AddTransform(WaveletTrans);
              }
            ConductivityTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(std::log(mincond), std::log(maxcond))));
            ConductivityCrossTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(std::log(mincond), std::log(maxcond))));
            ConductivityTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::LogTransform(RefModel)));
            ConductivityCrossTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::LogTransform(RefModel)));
            ConductivityTransform->AddTransform(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::ExpansionTransform(3 * ngrid, 2 * ngrid, 3 * ngrid)));

            TomoTransform = SlownessTransform;
            GravityTransform = DensityTransform;
            MTTransform = ConductivityTransform;
            RegTransform = TomoTransform;
            SlowCrossTrans = SlownessCrossTransform;
            CondCrossTrans = ConductivityCrossTransform;
            DensCrossTrans = DensityCrossTransform;
          }
        else
          {
            //if we want direct parameter coupling we do not need the section transforms
            //but we directly calculate conductivity and density from slowness
            if (!Wavelet)
              {
                TomoTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::TanhTransform(minslow, maxslow));
                GravityTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::DensityTransform(TomoTransform));
                MTTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::ConductivityTransform(TomoTransform));

                RegTransform = TomoTransform;
              }
            else
              {
                boost::shared_ptr<jiba::ChainedTransform> SlownessTransform(
                    new jiba::ChainedTransform);
                SlownessTransform->AddTransform(WaveletTrans);
                SlownessTransform->AddTransform(
                    boost::shared_ptr<jiba::GeneralModelTransform>(
                        new jiba::TanhTransform(minslow, maxslow)));

                boost::shared_ptr<jiba::ChainedTransform> DensityTransform(
                    new jiba::ChainedTransform);
                DensityTransform->AddTransform(WaveletTrans);
                DensityTransform->AddTransform(
                    boost::shared_ptr<jiba::GeneralModelTransform>(
                        new jiba::TanhTransform(mindens, maxdens)));

                boost::shared_ptr<jiba::ChainedTransform> ConductivityTransform(
                    new jiba::ChainedTransform);
                ConductivityTransform->AddTransform(WaveletTrans);
                ConductivityTransform->AddTransform(
                    boost::shared_ptr<jiba::GeneralModelTransform>(
                        new jiba::TanhTransform(std::log(mincond), std::log(maxcond))));

                RegTransform = boost::shared_ptr<jiba::ModelCopyTransform>(
                    new jiba::ModelCopyTransform());

                TomoTransform = SlownessTransform;
                GravityTransform = DensityTransform;
                MTTransform = ConductivityTransform;
                RegTransform = TomoTransform;
              }

          }
        SlowTrans = TomoTransform;
        CondTrans = MTTransform;
        DensTrans = GravityTransform;
        RegTrans = RegTransform;

      }

    void SetupCoupling::SetupCrossGradModel(jiba::rvec &InvModel,
        const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDSeismicModel &SeisMod, const jiba::ThreeDGravityModel &GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::MatOpRegularization> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);

        jiba::rvec SeisModel(ngrid, 0.0);
        if (SeisMod.GetNModelElements() == ngrid)
          {
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, SeisModel.begin());
            std::cout << "Transforming slowness model. " << std::endl;
            ublas::subrange(InvModel, 0, ngrid) = SlowTrans->PhysicalToGeneralized(
                SeisModel);
          }

        jiba::rvec GravModel(ngrid, 1.0);
        if (GravMod.GetNModelElements() == ngrid)
          {
            std::copy(GravMod.GetDensities().origin(),
                GravMod.GetDensities().origin() + ngrid, GravModel.begin());
            std::cout << "Transforming Density model. " << std::endl;
            ublas::subrange(InvModel, ngrid, 2 * ngrid) =
                DensTrans->PhysicalToGeneralized(GravModel);
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
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) =
                CondTrans->PhysicalToGeneralized(MTModel);
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
        boost::shared_ptr<jiba::MatOpRegularization> SeisReg(Regularization->clone());
        jiba::rvec TomoCovar(3 * ngrid);
        SetupModelCovar(TomoCovar, SeisModel, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(TomoCovar);

        //then the regularization of densities
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jiba::MatOpRegularization> GravReg(Regularization->clone());
        jiba::rvec GravCovar(3 * ngrid);
        jiba::rvec Ones(GravModel.size(), 1.0);
        SetupModelCovar(GravCovar, Ones, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(GravCovar);

        //and finally conductivities
        double mtreglambda = 1.0;
        std::cout << " Weight for MT regularization: ";
        std::cin >> mtreglambda;
        boost::shared_ptr<jiba::MatOpRegularization> MTReg(Regularization->clone());
        jiba::rvec MTCovar(3 * ngrid);
        SetupModelCovar(MTCovar, MTModel, MTReg->GetDataError(), ngrid);
        MTReg->SetDataError(MTCovar);
        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
        if (substart)
          {
            SeisReg->SetReferenceModel(SeisModel);
            GravReg->SetReferenceModel(GravModel);
            MTReg->SetReferenceModel(MTModel);
          }
        if (seisreglambda > 0.0)
          {
            Objective.AddObjective(SeisReg, SlowTrans, seisreglambda, "SeisReg");
          }
        if (gravreglambda > 0.0)
          {
            Objective.AddObjective(GravReg, DensTrans, gravreglambda, "GravReg");
          }
        if (mtreglambda > 0.0)
          {
            Objective.AddObjective(MTReg, CondTrans, mtreglambda, "MTReg");
          }
      }

    void SetupCoupling::SetupSaltModel(jiba::rvec &InvModel,
        const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDSeismicModel &SeisMod, const jiba::ThreeDGravityModel &GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::MatOpRegularization> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);

        jiba::rvec SeisModel(ngrid, 0.0);
        if (SeisMod.GetNModelElements() == ngrid)
          {
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, SeisModel.begin());
            std::cout << "Transforming slowness model. " << std::endl;
            ublas::subrange(InvModel, 0, ngrid) = SlowTrans->PhysicalToGeneralized(
                SeisModel);
          }

        jiba::rvec GravModel(ngrid, 1.0);
        if (GravMod.GetNModelElements() == ngrid)
          {
            std::copy(GravMod.GetDensities().origin(),
                GravMod.GetDensities().origin() + ngrid, GravModel.begin());
            std::cout << "Transforming Density model. " << std::endl;
            ublas::subrange(InvModel, ngrid, 2 * ngrid) =
                DensTrans->PhysicalToGeneralized(GravModel);
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
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) =
                CondTrans->PhysicalToGeneralized(MTModel);
          }
        else
          {
            std::fill_n(ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid).begin(), ngrid,
                0.0);
          }
        //we currently hacked out the density transform
        boost::shared_ptr<jiba::SaltRelConstraint> SaltRel(
            new jiba::SaltRelConstraint(
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::ModelCopyTransform()),
                boost::shared_ptr<jiba::GeneralModelTransform>(
                    new jiba::SlowCondTrans(
                        boost::shared_ptr<jiba::GeneralModelTransform>(
                            new jiba::ModelCopyTransform())))));

        jiba::ThreeDSeismicModel SaltExclude(SeisMod);
        //this is an EVIL HACK
        for (size_t i = 0; i < SaltExclude.GetSlownesses().shape()[0]; ++i)
          {
            for (size_t j = 0; j < SaltExclude.GetSlownesses().shape()[1]; ++j)
              {
                SaltExclude.SetSlownesses()[i][j][0] = -1.0;
              }
          }
        SaltRel->SetExcludeCells(SaltExclude);
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
        boost::shared_ptr<jiba::MatOpRegularization> SeisReg(Regularization->clone());
        jiba::rvec TomoCovar(3 * ngrid);
        SetupModelCovar(TomoCovar, SeisModel, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(TomoCovar);

        //then the regularization of densities
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jiba::MatOpRegularization> GravReg(Regularization->clone());
        jiba::rvec GravCovar(3 * ngrid);
        jiba::rvec Ones(GravModel.size(), 1.0);
        SetupModelCovar(GravCovar, Ones, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(GravCovar);

        //and finally conductivities
        double mtreglambda = 1.0;
        std::cout << " Weight for MT regularization: ";
        std::cin >> mtreglambda;
        boost::shared_ptr<jiba::MatOpRegularization> MTReg(Regularization->clone());
        jiba::rvec MTCovar(3 * ngrid);
        SetupModelCovar(MTCovar, MTModel, MTReg->GetDataError(), ngrid);
        MTReg->SetDataError(MTCovar);
        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
        if (substart)
          {
            SeisReg->SetReferenceModel(SeisModel);
            GravReg->SetReferenceModel(GravModel);
            MTReg->SetReferenceModel(MTModel);
          }
        if (seisreglambda > 0.0)
          {
            Objective.AddObjective(SeisReg, SlowTrans, seisreglambda, "SeisReg");
          }
        if (gravreglambda > 0.0)
          {
            Objective.AddObjective(GravReg, DensTrans, gravreglambda, "GravReg");
          }
        if (mtreglambda > 0.0)
          {
            Objective.AddObjective(MTReg, CondTrans, mtreglambda, "MTReg");
          }
      }

    void SetupCoupling::SetupFixedCouplingModel(jiba::rvec &InvModel,
        const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDSeismicModel &SeisMod, const jiba::ThreeDGravityModel &GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::MatOpRegularization> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();

        //if we want direct parameter coupling
        //we just need a single regularization term
        //as we invert for slowness, the regularization also works on slowness
        if (SeisMod.GetNModelElements() != ngrid)
          {
            throw jiba::FatalException("Size of seismic model does not match model geometry !");
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
        if (substart)
          {
            Regularization->SetReferenceModel(InvModel);
          }
        Objective.AddObjective(Regularization, RegTrans, reglambda, "Reg");
        InvModel = SlowTrans->PhysicalToGeneralized(InvModel);
      }

    void SetupCoupling::SetupModelVector(const po::variables_map &vm,
        jiba::rvec &InvModel, const jiba::ThreeDModelBase &ModelGeometry,
        const jiba::ThreeDSeismicModel &SeisMod, const jiba::ThreeDGravityModel &GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::MatOpRegularization> Regularization, bool substart)
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
                SetupSaltModel(InvModel, ModelGeometry, SeisMod, GravMod, MTMod,
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
