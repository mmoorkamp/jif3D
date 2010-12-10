//============================================================================
// Name        : SetupCoupling.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include "SetupCoupling.h"
#include "../Global/FileUtil.h"
#include "../Regularization/CrossGradient.h"

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
            ublas::subrange(Covar, ngrid, 2 * ngrid) = ublas::element_prod(
                InvModel, ublas::subrange(OldCov, ngrid, 2 * ngrid));
            ublas::subrange(Covar, 2 * ngrid, 3 * ngrid) = ublas::element_prod(
                InvModel, ublas::subrange(OldCov, 2 * ngrid, 3 * ngrid));
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
        desc.add_options()("crossgrad", "Use cross-gradient coupling")(
            "minslow", po::value(&minslow)->default_value(1e-4))("maxslow",
            po::value(&maxslow)->default_value(0.005))("mincond", po::value(
            &mincond)->default_value(1e-4))("maxcond",
            po::value(&maxcond)->default_value(5))("mindens", po::value(
            &mindens)->default_value(0.5))("maxdens",
            po::value(&maxdens)->default_value(4.0));

        return desc;
      }

    void SetupCoupling::SetupTransforms(const po::variables_map &vm,
        boost::shared_ptr<jiba::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &GravityTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &MTTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &RegTransform)
      {

        //we need the geometry of the starting model to setup
        //the transformations
        jiba::ThreeDSeismicModel StartModel;
        std::string modelfilename = jiba::AskFilename(
            "Inversion Model Geometry: ");
        //although we specify the starting model as a seismic model
        //it does not need to have the same grid specifications
        StartModel.ReadNetCDF(modelfilename, false);

        const size_t ngrid = StartModel.GetNModelElements();
        //if we want to do a cross-gradient type joint inversion
        //we need to set transformations for each data type
        //that extract the right part of the model vector
        //and then transform from generalized to physical parameters
        if (vm.count("crossgrad"))
          {
            //each set of transformations is chained together in a similar way
            //the section transform takes the right part of the model vector
            //the TanHTransform sorts out the constraints
            //and the ExpansionTransform puts the derivatives in the right part of the whole model vector
            //first we do slowness
            boost::shared_ptr<jiba::ChainedTransform> SlownessTransform(
                new jiba::ChainedTransform);
            SlownessTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::SectionTransform(0,
                ngrid)));
            SlownessTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::TanhTransform(minslow,
                maxslow)));
            SlownessTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::ExpansionTransform(3
                * ngrid, 0, ngrid)));
            //then we do density
            boost::shared_ptr<jiba::ChainedTransform> DensityTransform(
                new jiba::ChainedTransform);
            DensityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::SectionTransform(ngrid,
                2 * ngrid)));
            DensityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::TanhTransform(mindens,
                maxdens)));
            DensityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::ExpansionTransform(3
                * ngrid, ngrid, 2 * ngrid)));
            //and then conductivity, conductivity has a LogTransform in addition
            //to reduce the range of inversion parameters
            boost::shared_ptr<jiba::ChainedTransform> ConductivityTransform(
                new jiba::ChainedTransform);
            jiba::rvec RefModel(ngrid);
            std::fill(RefModel.begin(), RefModel.end(), 1.0);
            ConductivityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::SectionTransform(2
                * ngrid, 3 * ngrid)));
            ConductivityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::TanhTransform(std::log(
                mincond), std::log(maxcond))));
            ConductivityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::LogTransform(RefModel)));
            ConductivityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::ExpansionTransform(3
                * ngrid, 2 * ngrid, 3 * ngrid)));

            TomoTransform = SlownessTransform;
            GravityTransform = DensityTransform;
            MTTransform = ConductivityTransform;
            RegTransform = TomoTransform;
          }
        else
          {
            //if we want direct parameter coupling we do not need the section transforms
            //but we directly calculate conductivity and density from slowness
            TomoTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::TanhTransform(minslow, maxslow));
            GravityTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::DensityTransform(TomoTransform));
            MTTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::ConductivityTransform(TomoTransform));

            RegTransform = TomoTransform;
          }
        SlowTrans = TomoTransform;
        CondTrans = MTTransform;
        DensTrans = GravityTransform;

      }

    void SetupCoupling::SetupModelVector(const po::variables_map &vm,
        jiba::rvec &InvModel, const jiba::ThreeDSeismicModel &SeisMod,
        const jiba::ThreeDGravityModel GravMod,
        const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::MatOpRegularization> Regularization,
        bool substart)
      {

        const size_t ngrid = SeisMod.GetSlownesses().num_elements();

        //if we want to do corss-gradient type inversion
        //we need three cross-gradient objective functions
        //and three regularization terms
        if (vm.count("crossgrad"))
          {
            //first we check whether the size of the three inversion grid matches
            //if (MTMod.GetConductivities().num_elements() != ngrid
            //    || GravMod.GetDensities().num_elements() != ngrid)
            //  {
            //    throw jiba::FatalException(" Grids have different sizes !");
            //  }
            // we construct the model vector for the starting model
            //from the different starting models and the transformations
            InvModel.resize(3 * ngrid, 0.0);

            jiba::rvec SeisModel(ngrid, 0.0);
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, SeisModel.begin());
            ublas::subrange(InvModel, 0, ngrid)
                = SlowTrans->PhysicalToGeneralized(SeisModel);

            jiba::rvec GravModel(ngrid, 1.0);
            if (GravMod.GetDensities().num_elements() == ngrid)
              {
                std::copy(GravMod.GetDensities().origin(),
                    GravMod.GetDensities().origin() + ngrid, GravModel.begin());
                ublas::subrange(InvModel, ngrid, 2 * ngrid)
                    = DensTrans->PhysicalToGeneralized(GravModel);
              }

            jiba::rvec MTModel(ngrid, 1.0);
            if (MTMod.GetConductivities().num_elements() == ngrid)
              {
                std::copy(MTMod.GetConductivities().origin(),
                    MTMod.GetConductivities().origin() + ngrid, MTModel.begin());
                ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid)
                    = CondTrans->PhysicalToGeneralized(MTModel);
              }
            //then we constrcut the three cross gradient terms
            //the double section transform takes two sections of the model
            //and feeds them to the objective function
            boost::shared_ptr<jiba::CrossGradient> SeisGravCross(
                new jiba::CrossGradient(SeisMod));
            boost::shared_ptr<jiba::GeneralModelTransform> SeisGravTrans(
                new jiba::DoubleSectionTransform(3 * ngrid, 0, ngrid, ngrid, 2
                    * ngrid));
            //for each cross-gradient term we ask for a weight
            double seisgravlambda = 1.0;
            std::cout << "Weight for seismic-gravity cross-gradient term: ";
            std::cin >> seisgravlambda;
            if (seisgravlambda > 0.0)
              {
                Objective.AddObjective(SeisGravCross, SeisGravTrans,
                    seisgravlambda, "SeisGrav");
              }
            boost::shared_ptr<jiba::CrossGradient> SeisMTCross(
                new jiba::CrossGradient(SeisMod));
            boost::shared_ptr<jiba::GeneralModelTransform> SeisMTTrans(
                new jiba::DoubleSectionTransform(3 * ngrid, 0, ngrid,
                    2 * ngrid, 3 * ngrid));

            double seismtlambda = 1.0;
            std::cout << "Weight for seismic-MT cross-gradient term: ";
            std::cin >> seismtlambda;
            if (seismtlambda > 0.0)
              {
                Objective.AddObjective(SeisMTCross, SeisMTTrans, seismtlambda,
                    "SeisMT");
              }
            boost::shared_ptr<jiba::CrossGradient> GravMTCross(
                new jiba::CrossGradient(SeisMod));
            boost::shared_ptr<jiba::GeneralModelTransform> GravMTTrans(
                new jiba::DoubleSectionTransform(3 * ngrid, ngrid, 2 * ngrid, 2
                    * ngrid, 3 * ngrid));

            double gravmtlambda = 1.0;
            std::cout << "Weight for gravity-MT cross-gradient term: ";
            std::cin >> gravmtlambda;
            if (gravmtlambda > 0.0)
              {
                Objective.AddObjective(GravMTCross, GravMTTrans, seismtlambda,
                    "GravMT");
              }
            //finally we construct the regularization terms
            //we ask for a weight and construct a regularization object
            //for each type of physical parameter separately
            //first we set up seismic tomography
            double seisreglambda = 1.0;
            std::cout << " Weight for seismic regularization: ";
            std::cin >> seisreglambda;
            boost::shared_ptr<jiba::MatOpRegularization> SeisReg(
                Regularization->clone());
            jiba::rvec TomoCovar(3 * ngrid);
            SetupModelCovar(TomoCovar, SeisModel, SeisReg->GetDataCovar(),
                ngrid);
            SeisReg->SetDataCovar(TomoCovar);

            //then the regularization of densities
            double gravreglambda = 1.0;
            std::cout << " Weight for gravity regularization: ";
            std::cin >> gravreglambda;
            boost::shared_ptr<jiba::MatOpRegularization> GravReg(
                Regularization->clone());
            jiba::rvec GravCovar(3 * ngrid);
            SetupModelCovar(GravCovar, GravModel, GravReg->GetDataCovar(),
                ngrid);
            GravReg->SetDataCovar(GravCovar);

            //and finally conductivities
            double mtreglambda = 1.0;
            std::cout << " Weight for MT regularization: ";
            std::cin >> mtreglambda;
            boost::shared_ptr<jiba::MatOpRegularization> MTReg(
                Regularization->clone());
            jiba::rvec MTCovar(3 * ngrid);
            SetupModelCovar(MTCovar, MTModel, MTReg->GetDataCovar(), ngrid);
            MTReg->SetDataCovar(MTCovar);
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
                Objective.AddObjective(SeisReg, SlowTrans, seisreglambda,
                    "SeisReg");
              }
            if (gravreglambda > 0.0)
              {
                Objective.AddObjective(GravReg, DensTrans, gravreglambda,
                    "GravReg");
              }
            if (mtreglambda > 0.0)
              {
                Objective.AddObjective(MTReg, CondTrans, mtreglambda, "MTReg");
              }
          }
        else
          {
            //if we want direct parameter coupling
            //we just need a single regularization term
            //as we invert for slowness, the regularization also works on slowness
            InvModel.resize(ngrid);
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, InvModel.begin());
            double reglambda = 1.0;
            std::cout << " Weight for regularization: ";
            std::cin >> reglambda;
            jiba::rvec TomoCovar(3 * ngrid);
            SetupModelCovar(TomoCovar, InvModel,
                Regularization->GetDataCovar(), ngrid);
            Regularization->SetDataCovar(TomoCovar);
            if (substart)
              {
                Regularization->SetReferenceModel(InvModel);
              }
            Objective.AddObjective(Regularization, SlowTrans, reglambda, "Reg");
            InvModel = SlowTrans->PhysicalToGeneralized(InvModel);

          }
      }
  }
