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

namespace jif3D
  {
    void SetupModelCovar(jif3D::rvec &Covar, const jif3D::rvec &InvModel,
        const jif3D::rvec &OldCov, size_t ngrid)
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
            "Value to use for Conductivity where relationship is not valid")(
            "usetransformed", "Use transformed values for calculating cross-gradient");

        return desc;
      }

    void SetupCoupling::SetupTransforms(const po::variables_map &vm,
        ThreeDSeismicModel &GeometryModel,
        boost::shared_ptr<jif3D::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jif3D::GeneralModelTransform> &GravityTransform,
        boost::shared_ptr<jif3D::GeneralModelTransform> &MTTransform, bool Wavelet)
      {

        //we need the geometry of the starting model to setup
        //the transformations
        std::string modelfilename = jif3D::AskFilename("Inversion Model Geometry: ");
        //although we specify the starting model as a seismic model
        //it does not need to have the same grid specifications
        GeometryModel.ReadNetCDF(modelfilename, false);

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
            boost::shared_ptr<jif3D::ChainedTransform> SlownessTransform(
                new jif3D::ChainedTransform);
            boost::shared_ptr<jif3D::GeneralModelTransform> Copier(
                new jif3D::ModelCopyTransform);
            SlownessTransform->AppendTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::TanhTransform(minslow, maxslow)));
            SlowCrossTrans = boost::shared_ptr<jif3D::GeneralModelTransform>(
                SlownessTransform->clone());
            TomoTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::MultiSectionTransform(3 * ngrid, 0, ngrid, SlownessTransform));
            //we regularize on the raw model parameters as these are more evenly spaced than slowness values
            SlowRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
                new jif3D::ChainedTransform);
            SlowRegTrans->AppendTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::MultiSectionTransform(3 * ngrid, 0, ngrid, Copier)));

            //then we do density
            boost::shared_ptr<jif3D::ChainedTransform> DensityTransform(
                new jif3D::ChainedTransform);
            DensityTransform->AppendTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::TanhTransform(mindens, maxdens)));
            DensCrossTrans = boost::shared_ptr<jif3D::GeneralModelTransform>(
                DensityTransform->clone());
            GravityTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::MultiSectionTransform(3 * ngrid, ngrid, 2 * ngrid,
                    DensityTransform));
            DensRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
                new jif3D::ChainedTransform);
            DensRegTrans->AppendTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::MultiSectionTransform(3 * ngrid, ngrid, 2 * ngrid,
                        Copier)));

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
            CondCrossTrans = boost::shared_ptr<jif3D::GeneralModelTransform>(
                ConductivityTransform->clone());
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
            if (vm.count("usetransformed"))
              {
                CondCrossTrans = Copier;
                DensCrossTrans = Copier;
                SlowCrossTrans = Copier;
              }
            //finished setting up cross-gradient coupling

          }
        else
          {
            // We can specify a model file where we want
            // the parameter relationship to be valid
            jif3D::ThreeDSeismicModel RelModel;
            if (vm.count("relmodel"))
              {
                RelModel.ReadNetCDF(RelModelName, false);
              }
            else
              {
                const size_t nx = GeometryModel.GetSlownesses().shape()[0];
                const size_t ny = GeometryModel.GetSlownesses().shape()[1];
                const size_t nz = GeometryModel.GetSlownesses().shape()[2];
                RelModel.SetSlownesses().resize(boost::extents[nx][ny][nz]);
                std::fill_n(RelModel.SetSlownesses().origin(), nx * ny * nz, 1.0);

              }
            //if we want direct parameter coupling we do not need the section transforms
            //but we directly calculate conductivity and density from slowness

            TomoTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(minslow, maxslow));

            GravityTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::DensityTransform(TomoTransform, RelModel, DensReplace,
                    density_a, density_b));
            MTTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::ConductivityTransform(TomoTransform, RelModel, CondReplace,
                    cond_a, cond_b, cond_c));
            SlowRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
                new jif3D::ChainedTransform);
            SlowRegTrans->AppendTransform(TomoTransform);
            if (Wavelet)
              {
                std::cout << " Using wavelet parametrization " << std::endl;
                SlowRegTrans->AppendTransform(WaveletTrans);
              }

          }
        SlowTrans = TomoTransform;
        CondTrans = MTTransform;
        DensTrans = GravityTransform;
      }

    void SetupCoupling::SetupCrossGradModel(jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDSeismicModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);

        jif3D::rvec SeisModel(ngrid, 0.0);
        if (SeisMod.GetNModelElements() == ngrid)
          {
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, SeisModel.begin());
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
            new jif3D::CrossGradient(ModelGeometry));
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
            new jif3D::CrossGradient(ModelGeometry));
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
            new jif3D::CrossGradient(ModelGeometry));
        boost::shared_ptr<jif3D::MultiSectionTransform> GravMTTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        GravMTTrans->AddSection(ngrid, 2 * ngrid, DensCrossTrans);
        GravMTTrans->AddSection(2 * ngrid, 3 * ngrid, CondCrossTrans);

        double gravmtlambda = 1.0;
        std::cout << "Weight for gravity-MT cross-gradient term: ";
        std::cin >> gravmtlambda;
        if (gravmtlambda > 0.0)
          {
            Objective.AddObjective(GravMTCross, GravMTTrans, seismtlambda, "GravMT",
                JointObjective::coupling);
          }
        //finally we construct the regularization terms
        //we ask for a weight and construct a regularization object
        //for each type of physical parameter separately
        //first we set up seismic tomography
        jif3D::rvec Ones(GravModel.size(), 1.0);
        double seisreglambda = 1.0;
        std::cout << " Weight for seismic regularization: ";
        std::cin >> seisreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> SeisReg(Regularization->clone());
        jif3D::rvec TomoCovar(3 * ngrid);
        SetupModelCovar(TomoCovar, Ones, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(TomoCovar);

        //then the regularization of densities
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> GravReg(Regularization->clone());
        jif3D::rvec GravCovar(3 * ngrid);

        SetupModelCovar(GravCovar, Ones, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(GravCovar);

        //and finally conductivities
        double mtreglambda = 1.0;
        std::cout << " Weight for MT regularization: ";
        std::cin >> mtreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> MTReg(Regularization->clone());
        jif3D::rvec MTCovar(3 * ngrid);
        SetupModelCovar(MTCovar, Ones, MTReg->GetDataError(), ngrid);
        MTReg->SetDataError(MTCovar);
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

    void SetupCoupling::SetupSaltModel(const po::variables_map &vm, jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDSeismicModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);

        jif3D::rvec SeisModel(ngrid, 0.0);
        if (SeisMod.GetNModelElements() == ngrid)
          {
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, SeisModel.begin());
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
        else
          {
            std::fill_n(ublas::subrange(InvModel, ngrid, 2 * ngrid).begin(), ngrid, 0.0);
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
        else
          {
            std::fill_n(ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid).begin(), ngrid,
                0.0);
          }

        jif3D::ThreeDSeismicModel RelModel;
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
        boost::shared_ptr<jif3D::SaltRelConstraint> SaltRel(
            new jif3D::SaltRelConstraint(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::ModelCopyTransform()),
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::ConductivityTransform(
                        boost::shared_ptr<jif3D::GeneralModelTransform>(
                            new jif3D::ModelCopyTransform()), RelModel))));

        SaltRel->SetExcludeCells(RelModel);
        boost::shared_ptr<jif3D::MultiSectionTransform> SaltRelTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
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
            Objective.AddObjective(SaltRel, SaltRelTrans, saltrellambda, "SaltRel",
                JointObjective::coupling);
          }
        //finally we construct the regularization terms
        //we ask for a weight and construct a regularization object
        //for each type of physical parameter separately
        //first we set up seismic tomography
        jif3D::rvec Ones(ngrid, 1.0);
        double seisreglambda = 1.0;
        std::cout << " Weight for seismic regularization: ";
        std::cin >> seisreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> SeisReg(Regularization->clone());
        jif3D::rvec TomoCovar(3 * ngrid);
        SetupModelCovar(TomoCovar, Ones, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(TomoCovar);

        //then the regularization of densities
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> GravReg(Regularization->clone());
        jif3D::rvec GravCovar(3 * ngrid);
        SetupModelCovar(GravCovar, Ones, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(GravCovar);

        //and finally conductivities
        double mtreglambda = 1.0;
        std::cout << " Weight for MT regularization: ";
        std::cin >> mtreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> MTReg(Regularization->clone());
        jif3D::rvec MTCovar(3 * ngrid);
        SetupModelCovar(MTCovar, Ones, MTReg->GetDataError(), ngrid);
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

    void SetupCoupling::SetupFixedCouplingModel(jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDSeismicModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();

        //if we want direct parameter coupling
        //we just need a single regularization term
        //as we invert for slowness, the regularization also works on slowness
        if (SeisMod.GetNModelElements() != ngrid)
          {
            throw jif3D::FatalException(
                "Size of seismic model does not match model geometry !", __FILE__, __LINE__);
          }
        InvModel.resize(ngrid);
        std::copy(SeisMod.GetSlownesses().origin(),
            SeisMod.GetSlownesses().origin() + ngrid, InvModel.begin());
        double reglambda = 1.0;
        std::cout << " Weight for regularization: ";
        std::cin >> reglambda;
        if (reglambda > 0.0)
          {
            jif3D::rvec TomoCovar(3 * ngrid);
            SetupModelCovar(TomoCovar, InvModel, Regularization->GetDataError(), ngrid);
            Regularization->SetDataError(TomoCovar);
          }
        if (substart)
          {

            jif3D::rvec RefModel = SlowTrans->PhysicalToGeneralized(InvModel);
            RefModel = SlowRegTrans->GeneralizedToPhysical(RefModel);
            Regularization->SetReferenceModel(RefModel);
          }
        Objective.AddObjective(Regularization, SlowRegTrans, reglambda, "Reg",
            JointObjective::regularization);
        InvModel = SlowTrans->PhysicalToGeneralized(InvModel);
      }

    void SetupCoupling::SetupModelVector(const po::variables_map &vm,
        jif3D::rvec &InvModel, const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDSeismicModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {
        //depending on the type of coupling the model vector looks quite different
        //and we have to setup a number of different things
        //so we separate it into different functions even though it means
        //passing around many parameters
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
