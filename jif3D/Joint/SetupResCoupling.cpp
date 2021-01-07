/*
 * SetupResCoupling.cpp
 *
 *  Created on: Nov 18, 2020
 *      Author: zhanjie and mmoorkamp
 */


#include "../Global/FileUtil.h"
#include "../Regularization/CrossGradient.h"
#include "../Inversion/ModelTransforms.h"
#include "../MI/MutualInformationConstraint.h"
#include "../Joint/SaltRelConstraint.h"
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include "SetupResCoupling.h"

namespace ublas = boost::numeric::ublas;

namespace jif3D
  {
    void SetupResCModelCovar(jif3D::rvec &Covar, const jif3D::rvec &InvModel,
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

    SetupResCoupling::SetupResCoupling() :
        minres(0.0), maxres(0.0), minslow(0.0), maxslow(0.0), mindens(0.0), maxdens(
            0.0), density_a(0.0), density_b(0.0), res_a(0.0), res_b(0.0), res_c(0.0), vpvs(
            std::sqrt(3.0)), DensReplace(0.0), ResReplace(3.0)
      {
      }

    SetupResCoupling::~SetupResCoupling()
      {
      }

    po::options_description SetupResCoupling::SetupOptions()
      {
        po::options_description desc("Coupling options");
        desc.add_options()("crossgrad", "Use cross-gradient coupling")("saltrel",
            "Use a parameter constraint designed for salt")("mutual_information",
            po::value(&mibins)->default_value(50), "Use a MI based coupling constraint")(
            "minslow", po::value(&minslow)->default_value(1e-4))("maxslow",
            po::value(&maxslow)->default_value(0.005))("minvel",
            po::value(&minvel)->default_value(1000))("maxvel",
            po::value(&maxvel)->default_value(8000))("minres",
            po::value(&minres)->default_value(1e-6))("maxres",
            po::value(&maxres)->default_value(1000))("mindens",
            po::value(&mindens)->default_value(-1000.0))("maxdens",
            po::value(&maxdens)->default_value(1000.0))("vpvs",
            po::value(&vpvs)->default_value(std::sqrt(3.0)))("density_a",
            po::value(&density_a)->default_value(5000),
            "The slope of the velocity-density relationship")("density_b",
            po::value(&density_b)->default_value(8500),
            "The offset of the velocity-density relationship")("res_a",
            po::value(&res_a)->default_value(2.31e-7),
            "The quadratic term of the velocity-conductivity relationship")("res_b",
            po::value(&res_b)->default_value(-5.79e-4),
            "The linear term of the velocity-conductivity relationship")("res_c",
            po::value(&res_c)->default_value(0.124),
            "The constant term of the velocity-conductivity relationship")("relmodel",
            po::value(&RelModelName),
            "Name of a model file that specifies where to apply the parameter relationship")(
            "DensReplace", po::value(&DensReplace)->default_value(0.0),
            "Value to use for Density where relationship is not valid")("ResReplace",
            po::value(&ResReplace)->default_value(3.3),
            "Value to use for Conductivity where relationship is not valid");

        return desc;
      }

    void SetupResCoupling::SetupTransforms(const po::variables_map &vm,
        ThreeDModelBase &GeometryModel,
        boost::shared_ptr<jif3D::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jif3D::GeneralModelTransform> &GravityTransform,
        boost::shared_ptr<jif3D::GeneralModelTransform> &DCResTransform, bool Wavelet)
      {

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

        if (vm.count("crossgrad") || vm.count("saltrel")
            || vm.count("mutual_information"))
          {
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
                    new jif3D::MultiSectionTransform(3 * ngrid, ngrid, 2 * ngrid,
                        Copier)));

            //then we do velocity, for surface waves we need both velocity and density
            boost::shared_ptr<jif3D::GeneralModelTransform> SVelTrans =
                boost::make_shared<jif3D::TanhTransform>(minvel, maxvel);
            boost::shared_ptr<jif3D::ChainedTransform> SurfTrans(
                new jif3D::ChainedTransform);
            boost::shared_ptr<jif3D::MultiSectionTransform> SDensTrans =
                boost::make_shared<jif3D::MultiSectionTransform>(3 * ngrid);
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
            //and then resistivity, resistivity has a LogTransform in addition
            //to reduce the range of inversion parameters
            boost::shared_ptr<jif3D::ChainedTransform> ResistivityTransform(
                new jif3D::ChainedTransform);
            jif3D::rvec RefModel(ngrid);
            std::fill(RefModel.begin(), RefModel.end(), 1.0);
            ResistivityTransform->AppendTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::TanhTransform(std::log(minres), std::log(maxres))));
            ResistivityTransform->AppendTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::LogTransform(RefModel)));
            ResCrossTrans = Copier;
            DCResTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::MultiSectionTransform(3 * ngrid, 2 * ngrid, 3 * ngrid,
                		ResistivityTransform));
            //we regularize on the raw model parameters as these are more evenly spaced than resistivities
            ResRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
                new jif3D::ChainedTransform);
            ResRegTrans->AppendTransform(
                boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::MultiSectionTransform(3 * ngrid, 2 * ngrid, 3 * ngrid,
                        Copier)));
            //if we want to regularize in the wavelet  domain
            //we need to add a wavelet transform to the regularization
            if (Wavelet)
              {
                SlowRegTrans->AppendTransform(WaveletTrans);
                DensRegTrans->AppendTransform(WaveletTrans);
                ResRegTrans->AppendTransform(WaveletTrans);
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
                const size_t nx = GeometryModel.GetData().shape()[0];
                const size_t ny = GeometryModel.GetData().shape()[1];
                const size_t nz = GeometryModel.GetData().shape()[2];
                RelModel.SetSlownesses().resize(boost::extents[nx][ny][nz]);
                std::fill_n(RelModel.SetSlownesses().origin(), nx * ny * nz, 1.0);

              }
            //if we want direct parameter coupling we do not need the section transforms
            //but we directly calculate resistivity and density from slowness

            TomoTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(minslow, maxslow));

            GravityTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::DensityTransform(TomoTransform, RelModel, DensReplace,
                    density_a, density_b));
            DCResTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::ResistivityTransform(TomoTransform, RelModel, ResReplace,
                    res_a, res_b, res_c));
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
        ResTrans = DCResTransform;
        DensTrans = GravityTransform;
      }

    void SetupResCoupling::SetupCrossGradModel(jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry, const SeisModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDDCResistivityModel &DCResMod,
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

        jif3D::rvec DCResModel(ngrid, 1.0);
        if (DCResMod.GetNModelElements() == ngrid)
          {

            std::copy(DCResMod.GetResistivities().origin(),
                DCResMod.GetResistivities().origin() + ngrid, DCResModel.begin());
            std::cout << "Transforming resistivity model. " << std::endl;
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(
                ResTrans->PhysicalToGeneralized(DCResModel), 2 * ngrid, 3 * ngrid);
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
        boost::shared_ptr<jif3D::CrossGradient> SeisDCResCross(
            new jif3D::CrossGradient(ModelGeometry, TearModelX, TearModelY, TearModelZ));
        boost::shared_ptr<jif3D::MultiSectionTransform> SeisDCResTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        SeisDCResTrans->AddSection(0, ngrid, SlowCrossTrans);
        SeisDCResTrans->AddSection(2 * ngrid, 3 * ngrid, ResCrossTrans);

        double seisdclambda = 1.0;
        std::cout << "Weight for seismic-DC cross-gradient term: ";
        std::cin >> seisdclambda;
        if (seisdclambda > 0.0)
          {
            Objective.AddObjective(SeisDCResCross, SeisDCResTrans, seisdclambda, "SeisDC",
                JointObjective::coupling);
          }
        boost::shared_ptr<jif3D::CrossGradient> GravDCResCross(
            new jif3D::CrossGradient(ModelGeometry, TearModelX, TearModelY, TearModelZ));
        boost::shared_ptr<jif3D::MultiSectionTransform> GravDCResTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        GravDCResTrans->AddSection(ngrid, 2 * ngrid, DensCrossTrans);
        GravDCResTrans->AddSection(2 * ngrid, 3 * ngrid, ResCrossTrans);

        double gravdclambda = 1.0;
        std::cout << "Weight for gravity-DC cross-gradient term: ";
        std::cin >> gravdclambda;
        if (gravdclambda > 0.0)
          {
            Objective.AddObjective(GravDCResCross, GravDCResTrans, gravdclambda, "GravDC",
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
        SetupResCModelCovar(TomoCovar, TomoCovVec, SeisReg->GetDataError(), ngrid);
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

        SetupResCModelCovar(GravCovar, GravCovVec, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(std::vector<double>(GravCovar.begin(), GravCovar.end()));

        //and finally resistivities
        jif3D::rvec ResCovVec(ngrid, 1.0);
        if (!CovVec.empty())
          {
            if (CovVec.size() == ngrid)
              {
                ResCovVec = CovVec;
              }
            else
              {
                ResCovVec = ublas::subrange(CovVec, 2 * ngrid, 3 * ngrid);
              }
          }
        double dcreglambda = 1.0;
        std::cout << " Weight for DC regularization: ";
        std::cin >> dcreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> DCReg(Regularization->clone());
        jif3D::rvec ResCovar(3 * ngrid);
        SetupResCModelCovar(ResCovar, ResCovVec, DCReg->GetDataError(), ngrid);
        DCReg->SetDataError(std::vector<double>(ResCovar.begin(), ResCovar.end()));
        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
        if (substart)
          {
            SeisReg->SetReferenceModel(SlowRegTrans->GeneralizedToPhysical(InvModel));
            GravReg->SetReferenceModel(DensRegTrans->GeneralizedToPhysical(InvModel));
            DCReg->SetReferenceModel(ResRegTrans->GeneralizedToPhysical(InvModel));
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
        if (dcreglambda > 0.0)
          {
            Objective.AddObjective(DCReg, ResRegTrans, dcreglambda, "DCReg",
                JointObjective::regularization);
          }
      }

    void SetupResCoupling::SetupMIModel(jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry, const SeisModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDDCResistivityModel &DCResMod,
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

        jif3D::rvec DCResModel(ngrid, 1.0);
        if (DCResMod.GetNModelElements() == ngrid)
          {

            std::copy(DCResMod.GetResistivities().origin(),
                DCResMod.GetResistivities().origin() + ngrid, DCResModel.begin());
            //make sure that each layer has a slightly different conductivity
            //at least in one cell, otherwise x3d optimizes by joining layers
            //and messes up the gradient calculation
            std::cout << "Transforming resistivity model. " << std::endl;
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(
                ResTrans->PhysicalToGeneralized(DCResModel), 2 * ngrid, 3 * ngrid);
          }
        boost::shared_ptr<jif3D::GeneralModelTransform> Copier(
            new jif3D::ModelCopyTransform);

        //then we construct the three cross gradient terms
        //the double section transform takes two sections of the model
        //and feeds them to the objective function
        boost::shared_ptr<jif3D::MutualInformationConstraint> SeisGravMI =
            boost::make_shared<jif3D::MutualInformationConstraint>(-2.0, 2.0, -2.0, 2.0,
                mibins);
        boost::shared_ptr<jif3D::MultiSectionTransform> SeisGravTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        SeisGravTrans->AddSection(0, ngrid, Copier);
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
        boost::shared_ptr<jif3D::MutualInformationConstraint> SeisDCResMI =
            boost::make_shared<jif3D::MutualInformationConstraint>(-2.0, 2.0, -2.0, 2.0,
                mibins);

        boost::shared_ptr<jif3D::MultiSectionTransform> SeisDCResTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        SeisDCResTrans->AddSection(0, ngrid, Copier);
        SeisDCResTrans->AddSection(2 * ngrid, 3 * ngrid, Copier);

        double seisdclambda = 1.0;
        std::cout << "Weight for seismic-DC MI term: ";
        std::cin >> seisdclambda;
        if (seisdclambda > 0.0)
          {
            Objective.AddObjective(SeisDCResMI, SeisDCResTrans, seisdclambda, "SeisDC",
                JointObjective::coupling);
          }
        boost::shared_ptr<jif3D::MutualInformationConstraint> GravDCResMI =
            boost::make_shared<jif3D::MutualInformationConstraint>(-2.0, 2.0, -2.0, 2.0,
                mibins);
        boost::shared_ptr<jif3D::MultiSectionTransform> GravDCResTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        GravDCResTrans->AddSection(ngrid, 2 * ngrid, Copier);
        GravDCResTrans->AddSection(2 * ngrid, 3 * ngrid, Copier);

        double gravdclambda = 1.0;
        std::cout << "Weight for gravity-DC MI term: ";
        std::cin >> gravdclambda;
        if (gravdclambda > 0.0)
          {
            Objective.AddObjective(GravDCResMI, GravDCResTrans, gravdclambda, "GravDC",
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

        //and finally resistivities
        jif3D::rvec ResCovVec(ngrid, 1.0);
        double dcreglambda = 1.0;
        std::cout << " Weight for DC regularization: ";
        std::cin >> dcreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> DCReg(Regularization->clone());

        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
        if (substart)
          {
            SeisReg->SetReferenceModel(SlowRegTrans->GeneralizedToPhysical(InvModel));
            GravReg->SetReferenceModel(DensRegTrans->GeneralizedToPhysical(InvModel));
            DCReg->SetReferenceModel(ResRegTrans->GeneralizedToPhysical(InvModel));
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
        if (dcreglambda > 0.0)
          {
            Objective.AddObjective(DCReg, ResRegTrans, dcreglambda, "DCReg",
                JointObjective::regularization);
          }
      }

    void SetupResCoupling::SetupSaltModel(const po::variables_map &vm, jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry, const SeisModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDDCResistivityModel &DCResMod,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();
        InvModel.resize(3 * ngrid, 0.0);

        jif3D::rvec SeisModel(ngrid, 0.0);
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
        else
          {
            std::fill_n(ublas::subrange(InvModel, ngrid, 2 * ngrid).begin(), ngrid, 0.0);
          }

        jif3D::rvec DCResModel(ngrid, 1.0);
        if (DCResMod.GetNModelElements() == ngrid)
          {

            std::copy(DCResMod.GetResistivities().origin(),
                DCResMod.GetResistivities().origin() + ngrid, DCResModel.begin());
            //make sure that each layer has a slightly different conductivity
            //at least in one cell, otherwise x3d optimizes by joining layers
            //and messes up the gradient calculation
            std::cout << "Transforming resistivity model. " << std::endl;
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(
                ResTrans->PhysicalToGeneralized(DCResModel), 2 * ngrid, 3 * ngrid);
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
                    new jif3D::ResistivityTransform(
                        boost::shared_ptr<jif3D::GeneralModelTransform>(
                            new jif3D::ModelCopyTransform()), RelModel))));

        SaltRel->SetExcludeCells(RelModel);
        boost::shared_ptr<jif3D::MultiSectionTransform> SaltRelTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        SaltRelTrans->AddSection(0, ngrid, SlowCrossTrans);
        //evil HACK !
        SaltRelTrans->AddSection(0, ngrid, SlowCrossTrans);
        //SaltRelTrans->AddSection(ngrid, 2 * ngrid, DensCrossTrans);
        SaltRelTrans->AddSection(2 * ngrid, 3 * ngrid, ResCrossTrans);

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
        SetupResCModelCovar(TomoCovar, Ones, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(std::vector<double>(TomoCovar.begin(), TomoCovar.end()));

        //then the regularization of densities
        double gravreglambda = 1.0;
        std::cout << " Weight for gravity regularization: ";
        std::cin >> gravreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> GravReg(Regularization->clone());
        jif3D::rvec GravCovar(3 * ngrid);
        SetupResCModelCovar(GravCovar, Ones, GravReg->GetDataError(), ngrid);
        GravReg->SetDataError(std::vector<double>(GravCovar.begin(), GravCovar.end()));

        //and finally conductivities
        double dcreglambda = 1.0;
        std::cout << " Weight for DC regularization: ";
        std::cin >> dcreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> DCReg(Regularization->clone());
        jif3D::rvec DCResCovar(3 * ngrid);
        SetupResCModelCovar(DCResCovar, Ones, DCReg->GetDataError(), ngrid);
        DCReg->SetDataError(std::vector<double>(DCResCovar.begin(), DCResCovar.end()));

        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
        if (substart)
          {
            SeisReg->SetReferenceModel(SeisModel);
            GravReg->SetReferenceModel(GravModel);
            DCReg->SetReferenceModel(DCResModel);
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
        if (dcreglambda > 0.0)
          {
            Objective.AddObjective(DCReg, ResRegTrans, dcreglambda, "DCReg",
                JointObjective::regularization);
          }
      }

    void SetupResCoupling::SetupFixedCouplingModel(jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry, const SeisModel &SeisMod,
        const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDDCResistivityModel &DCResMod,
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
                "Size of seismic model does not match model geometry !", __FILE__,
                __LINE__);
          }
        InvModel.resize(ngrid);
        std::copy(SeisMod.GetData().origin(), SeisMod.GetData().origin() + ngrid,
            InvModel.begin());
        double reglambda = 1.0;
        std::cout << " Weight for regularization: ";
        std::cin >> reglambda;
        if (reglambda > 0.0)
          {
            jif3D::rvec TomoCovar(3 * ngrid);
            SetupResCModelCovar(TomoCovar, InvModel, Regularization->GetDataError(), ngrid);
            Regularization->SetDataError(
                std::vector<double>(TomoCovar.begin(), TomoCovar.end()));
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


    void SetupResCoupling::SetupModelVector(const po::variables_map &vm,
        jif3D::rvec &InvModel, const jif3D::ThreeDModelBase &ModelGeometry,
        const SeisModel &SeisMod, const jif3D::ThreeDGravityModel &GravMod,
        const jif3D::ThreeDDCResistivityModel &DCResMod, jif3D::JointObjective &Objective,
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
            SetupCrossGradModel(InvModel, ModelGeometry, SeisMod, GravMod, DCResMod,
                Objective, Regularization, substart, TearModelX, TearModelY, TearModelZ,
                CovVec);
          }
        else
          {
            if (vm.count("saltrel"))
              {
                SetupSaltModel(vm, InvModel, ModelGeometry, SeisMod, GravMod, DCResMod,
                    Objective, Regularization, substart);
              }
            else
              {
                if (vm.count("mutual_information"))
                  {
                    SetupMIModel(InvModel, ModelGeometry, SeisMod, GravMod, DCResMod,
                        Objective, Regularization, substart);
                  }
                else
                  {
                    SetupFixedCouplingModel(InvModel, ModelGeometry, SeisMod, GravMod,
                        DCResMod, Objective, Regularization, substart);
                  }
              }
          }
      }
  }
