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

    SetupCoupling::SetupCoupling()
      {
      }

    SetupCoupling::~SetupCoupling()
      {
      }

    po::options_description SetupCoupling::SetupOptions()
      {
        po::options_description desc("Coupling options");
        desc.add_options()("crossgrad", "Use cross-gradient coupling");

        return desc;
      }

    void SetupCoupling::SetupTransforms(const po::variables_map &vm,
        boost::shared_ptr<jiba::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &GravityTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &MTTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &RegTransform)
      {
        const double minslow = 1e-4;
        const double maxslow = 0.005;
        const double mincond = -10;
        const double maxcond = 3;
        const double mindens = 0.5;
        const double maxdens = 10.0;

        jiba::ThreeDSeismicModel StartModel;
        std::string modelfilename = jiba::AskFilename(
            "Inversion Model Geometry: ");
        StartModel.ReadNetCDF(modelfilename);

        const size_t ngrid = StartModel.GetNModelElements();
        if (vm.count("crossgrad"))
          {
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

            boost::shared_ptr<jiba::ChainedTransform> ConductivityTransform(
                new jiba::ChainedTransform);
            jiba::rvec RefModel(ngrid);
            std::fill(RefModel.begin(), RefModel.end(), 1.0);
            ConductivityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::SectionTransform(2
                * ngrid, 3 * ngrid)));
            ConductivityTransform->AddTransform(boost::shared_ptr<
                jiba::GeneralModelTransform>(new jiba::TanhTransform(mincond,
                maxcond)));
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
        boost::shared_ptr<jiba::MatOpRegularization> Regularization)
      {

        const size_t ngrid = SeisMod.GetSlownesses().num_elements();
        if (MTMod.GetConductivities().num_elements() != ngrid
            || GravMod.GetDensities().num_elements() != ngrid)
          {
            throw jiba::FatalException(" Grids have different sizes !");
          }
        if (vm.count("crossgrad"))
          {
            InvModel.resize(3 * ngrid);
            jiba::rvec TempModel(ngrid);
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, TempModel.begin());
            ublas::subrange(InvModel, 0, ngrid)
                = SlowTrans->PhysicalToGeneralized(TempModel);

            std::copy(GravMod.GetDensities().origin(),
                GravMod.GetDensities().origin() + ngrid, TempModel.begin());
            ublas::subrange(InvModel, ngrid, 2 * ngrid)
                = DensTrans->PhysicalToGeneralized(TempModel);

            std::copy(MTMod.GetConductivities().origin(),
                MTMod.GetConductivities().origin() + ngrid, TempModel.begin());
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid)
                = CondTrans->PhysicalToGeneralized(TempModel);

            boost::shared_ptr<jiba::CrossGradient> SeisGravCross(
                new jiba::CrossGradient(SeisMod));
            boost::shared_ptr<jiba::GeneralModelTransform> SeisGravTrans(
                new jiba::DoubleSectionTransform(3 * ngrid, 0, ngrid, ngrid, 2
                    * ngrid));

            double seisgravlambda = 1.0;
            std::cout << "Weight for seismic-gravity cross-gradient term: ";
            std::cin >> seisgravlambda;
            Objective.AddObjective(SeisGravCross, SeisGravTrans,
                seisgravlambda, "SeisGrav");

            boost::shared_ptr<jiba::CrossGradient> SeisMTCross(
                new jiba::CrossGradient(SeisMod));
            boost::shared_ptr<jiba::GeneralModelTransform> SeisMTTrans(
                new jiba::DoubleSectionTransform(3 * ngrid, 0, ngrid,
                    2 * ngrid, 3 * ngrid));

            double seismtlambda = 1.0;
            std::cout << "Weight for seismic-MT cross-gradient term: ";
            std::cin >> seismtlambda;
            Objective.AddObjective(SeisMTCross, SeisMTTrans, seismtlambda,
                "SeisMT");

            boost::shared_ptr<jiba::CrossGradient> GravMTCross(
                new jiba::CrossGradient(SeisMod));
            boost::shared_ptr<jiba::GeneralModelTransform> GravMTTrans(
                new jiba::DoubleSectionTransform(3 * ngrid, ngrid, 2 * ngrid, 2
                    * ngrid, 3 * ngrid));

            double gravmtlambda = 1.0;
            std::cout << "Weight for gravity-MT cross-gradient term: ";
            std::cin >> gravmtlambda;
            Objective.AddObjective(GravMTCross, GravMTTrans, seismtlambda,
                "GravMT");

            double seisreglambda = 1.0;
            std::cout << " Weight for seismic regularization: ";
            std::cin >> seisreglambda;
            Objective.AddObjective(
                boost::shared_ptr<jiba::MatOpRegularization>(
                    Regularization->clone()), SlowTrans, seisreglambda,
                "SeisReg");

            double gravreglambda = 1.0;
            std::cout << " Weight for gravity regularization: ";
            std::cin >> gravreglambda;
            Objective.AddObjective(
                boost::shared_ptr<jiba::MatOpRegularization>(
                    Regularization->clone()), DensTrans, gravreglambda,
                "GravReg");

            double mtreglambda = 1.0;
            std::cout << " Weight for MT regularization: ";
            std::cin >> mtreglambda;
            Objective.AddObjective(
                boost::shared_ptr<jiba::MatOpRegularization>(
                    Regularization->clone()), CondTrans, mtreglambda, "MTReg");
          }
        //if we want direct parameter coupling
        else
          {
            InvModel.resize(ngrid);
            std::copy(SeisMod.GetSlownesses().origin(),
                SeisMod.GetSlownesses().origin() + ngrid, InvModel.begin());
            InvModel = SlowTrans->PhysicalToGeneralized(InvModel);

            double reglambda = 1.0;
            std::cout << " Weight for regularization: ";
            std::cin >> reglambda;
            Objective.AddObjective(Regularization, SlowTrans, reglambda, "Reg");
          }
      }
  }
