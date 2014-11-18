//============================================================================
// Name        : SetupTomoDCMagCoupling.cpp
// Author      : August 3, 2014
// Version     : 
// Copyright   : 2014, zhanjie
//============================================================================

#include "SetupTomoDCMagCoupling.h"
#include "../Global/FileUtil.h"
#include "../Regularization/CrossGradient.h"
#include "../Inversion/ModelTransforms.h"
#include "SaltRelConstraint.h"

namespace ublas = boost::numeric::ublas;

namespace jif3D
  {
    void SetupTomoDCMagModelCovar(jif3D::rvec &Covar, const jif3D::rvec &InvModel,
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

    SetupTomoDCMagCoupling::SetupTomoDCMagCoupling() :
        minres(0.0), maxres(0.0), minslow(0.0), maxslow(0.0), minsuscep(0.0), maxsuscep(
            0.0)
      {
      }

    SetupTomoDCMagCoupling::~SetupTomoDCMagCoupling()
      {
      }

    po::options_description SetupTomoDCMagCoupling::SetupOptions()
      {
        po::options_description desc("Coupling options");
        desc.add_options()("minslow",
                po::value(&minslow)->default_value(1e-4))("maxslow",
            po::value(&maxslow)->default_value(0.005))("minres",
            po::value(&minres)->default_value(0.01))("maxres",
            po::value(&maxres)->default_value(1e6))("minsuscep",
            po::value(&minsuscep)->default_value(0.01))("maxsuscep",
            po::value(&maxsuscep)->default_value(1e6));

        return desc;
      }

    void SetupTomoDCMagCoupling::SetupTransforms(const po::variables_map &vm, ThreeDSeismicModel &GeometryModel,
            boost::shared_ptr<jif3D::GeneralModelTransform> &TomoTransform,
            boost::shared_ptr<jif3D::GeneralModelTransform> &DCTransform,
            boost::shared_ptr<jif3D::GeneralModelTransform> &MagTransform, bool Wavelet)
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
                    GeometryModel.GetData().shape()[1], GeometryModel.GetData().shape()[2]));
          }
        //if we want to do a cross-gradient type joint inversion
        //we need to set transformations for each data type
        //that extract the right part of the model vector
        //and then transform from generalized to physical parameters


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
        boost::shared_ptr<jif3D::GeneralModelTransform> Copier(new jif3D::ModelCopyTransform);
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
        SlowRegTrans->AppendTransform(boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::MultiSectionTransform(3 * ngrid, 0, ngrid, Copier)));

        //then we do resistivity, resistivity has a LogTransform in addition
        //to reduce the range of inversion parameters
        boost::shared_ptr<jif3D::ChainedTransform> ResistivityTransform(
            new jif3D::ChainedTransform);
        jif3D::rvec ResRefModel(ngrid);
        std::fill(ResRefModel.begin(), ResRefModel.end(), 1.0);
        ResistivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(std::log(minres), std::log(maxres))));
        ResistivityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::LogTransform(ResRefModel)));
        ResCrossTrans = boost::shared_ptr<jif3D::GeneralModelTransform>(
        		ResistivityTransform->clone());
        DCTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::MultiSectionTransform(3 * ngrid, ngrid, 2 * ngrid,
            		ResistivityTransform));
        //we regularize on the raw model parameters as these are more evenly spaced than resistivity
        ResRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
            new jif3D::ChainedTransform);
        ResRegTrans->AppendTransform(boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::MultiSectionTransform(3 * ngrid, ngrid, 2 * ngrid,
                    Copier)));

        //and then susceptibility, susceptibility also has a LogTransform in addition
        //to reduce the range of inversion parameters
        boost::shared_ptr<jif3D::ChainedTransform> SusceptibilityTransform(
            new jif3D::ChainedTransform);
        jif3D::rvec SusRefModel(ngrid);
        std::fill(SusRefModel.begin(), SusRefModel.end(), 1.0);
        SusceptibilityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::TanhTransform(std::log(minsuscep), std::log(maxsuscep))));
        SusceptibilityTransform->AppendTransform(
            boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::LogTransform(SusRefModel)));
        SuscepCrossTrans = boost::shared_ptr<jif3D::GeneralModelTransform>(
        		SusceptibilityTransform->clone());
        MagTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
            new jif3D::MultiSectionTransform(3 * ngrid, 2 * ngrid, 3 * ngrid,
            		SusceptibilityTransform));
        //we regularize on the raw model parameters as these are more evenly spaced than conductivities
        SuscepRegTrans = boost::shared_ptr<jif3D::ChainedTransform>(
            new jif3D::ChainedTransform);
        SuscepRegTrans->AppendTransform(boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::MultiSectionTransform(3 * ngrid, 2 * ngrid, 3 * ngrid,
                    Copier)));

        //if we want to regularize in the wavelet  domain
        //we need to add a wavelet transform to the regularization
        if (Wavelet)
          {
            SlowRegTrans->AppendTransform(WaveletTrans);
            ResRegTrans->AppendTransform(WaveletTrans);
            SuscepRegTrans->AppendTransform(WaveletTrans);
          }
        //finished setting up cross-gradient coupling

        SlowTrans = TomoTransform;
        ResTrans = DCTransform;
        SuscepTrans = MagTransform;
      }

    void SetupTomoDCMagCoupling::SetupCrossGradModel(jif3D::rvec &InvModel,
            const jif3D::ThreeDModelBase &ModelGeometry,
            const jif3D::ThreeDSeismicModel &SeisMod,
            const jif3D::ThreeDDCResistivityModel &DCMod, const jif3D::ThreeDMagneticModel &MagMod,
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

        jif3D::rvec DCModel(ngrid, 1.0);
        if (DCMod.GetNModelElements() == ngrid)
          {
            std::copy(DCMod.GetResistivities().origin(),
            		DCMod.GetResistivities().origin() + ngrid, DCModel.begin());
            std::cout << "Transforming resistivity model. " << std::endl;
            ublas::subrange(InvModel, ngrid, 2 * ngrid) = ublas::subrange(
            		ResTrans->PhysicalToGeneralized(DCModel), ngrid, 2 * ngrid);
          }

        jif3D::rvec MagModel(ngrid, 1.0);
        if (MagMod.GetNModelElements() == ngrid)
          {

            std::copy(MagMod.GetSusceptibilities().origin(),
            		MagMod.GetSusceptibilities().origin() + ngrid, MagModel.begin());
            std::cout << "Transforming susceptibility model. " << std::endl;
            ublas::subrange(InvModel, 2 * ngrid, 3 * ngrid) = ublas::subrange(
            		SuscepTrans->PhysicalToGeneralized(MagModel), 2 * ngrid, 3 * ngrid);
          }
        //then we construct the three cross gradient terms
        //the double section transform takes two sections of the model
        //and feeds them to the objective function
        boost::shared_ptr<jif3D::CrossGradient> SeisDCCross(
            new jif3D::CrossGradient(ModelGeometry));
        boost::shared_ptr<jif3D::MultiSectionTransform> SeisDCTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        SeisDCTrans->AddSection(0, ngrid, SlowCrossTrans);
        SeisDCTrans->AddSection(ngrid, 2 * ngrid, ResCrossTrans);

        //for each cross-gradient term we ask for a weight
        double seisdclambda = 1.0;
        std::cout << "Weight for seismic-dc cross-gradient term: ";
        std::cin >> seisdclambda;
        if (seisdclambda > 0.0)
          {
            Objective.AddObjective(SeisDCCross, SeisDCTrans, seisdclambda,
                "SeisDC",JointObjective::coupling);
          }

        boost::shared_ptr<jif3D::CrossGradient> SeisMagCross(
            new jif3D::CrossGradient(ModelGeometry));
        boost::shared_ptr<jif3D::MultiSectionTransform> SeisMagTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        SeisMagTrans->AddSection(0, ngrid, SlowCrossTrans);
        SeisMagTrans->AddSection(2 * ngrid, 3 * ngrid, SuscepCrossTrans);

        double seismaglambda = 1.0;
        std::cout << "Weight for seismic-mag cross-gradient term: ";
        std::cin >> seismaglambda;
        if (seismaglambda > 0.0)
          {
            Objective.AddObjective(SeisMagCross, SeisMagTrans, seismaglambda, "SeisMag",JointObjective::coupling);
          }

        boost::shared_ptr<jif3D::CrossGradient> DCMagCross(
            new jif3D::CrossGradient(ModelGeometry));
        boost::shared_ptr<jif3D::MultiSectionTransform> DCMagTrans(
            new jif3D::MultiSectionTransform(3 * ngrid));
        DCMagTrans->AddSection(ngrid, 2 * ngrid, ResCrossTrans);
        DCMagTrans->AddSection(2 * ngrid, 3 * ngrid, SuscepCrossTrans);

        double dcmaglambda = 1.0;
        std::cout << "Weight for dc-mag cross-gradient term: ";
        std::cin >> dcmaglambda;
        if (dcmaglambda > 0.0)
          {
            Objective.AddObjective(DCMagCross, DCMagTrans, dcmaglambda, "DCMag",JointObjective::coupling);
          }
        //finally we construct the regularization terms
        //we ask for a weight and construct a regularization object
        //for each type of physical parameter separately
        //first we set up seismic tomography
        jif3D::rvec Ones(DCModel.size(), 1.0);
        double seisreglambda = 1.0;
        std::cout << " Weight for seismic regularization: ";
        std::cin >> seisreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> SeisReg(Regularization->clone());
        jif3D::rvec TomoCovar(3 * ngrid);
        SetupTomoDCMagModelCovar(TomoCovar, Ones, SeisReg->GetDataError(), ngrid);
        SeisReg->SetDataError(TomoCovar);

        //then the regularization of resistivity
        double dcreglambda = 1.0;
        std::cout << " Weight for dcresistivity regularization: ";
        std::cin >> dcreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> DCReg(Regularization->clone());
        jif3D::rvec DCCovar(3 * ngrid);
        SetupTomoDCMagModelCovar(DCCovar, Ones, DCReg->GetDataError(), ngrid);
        DCReg->SetDataError(DCCovar);

        //and finally susceptibility
        double magreglambda = 1.0;
        std::cout << " Weight for magnetics regularization: ";
        std::cin >> magreglambda;
        boost::shared_ptr<jif3D::RegularizationFunction> MagReg(Regularization->clone());
        jif3D::rvec MagCovar(3 * ngrid);
        SetupTomoDCMagModelCovar(MagCovar, Ones, MagReg->GetDataError(), ngrid);
        MagReg->SetDataError(MagCovar);
        //if we specify on the command line that we want to subtract the
        //starting model, we set the corresponding reference model
        //in the regularization object
        if (substart)
          {
            SeisReg->SetReferenceModel(SlowRegTrans->GeneralizedToPhysical(InvModel));
            DCReg->SetReferenceModel(ResRegTrans->GeneralizedToPhysical(InvModel));
            MagReg->SetReferenceModel(SuscepRegTrans->GeneralizedToPhysical(InvModel));
          }
        if (seisreglambda > 0.0)
          {
            Objective.AddObjective(SeisReg, SlowRegTrans, seisreglambda, "SeisReg",JointObjective::regularization);
          }
        if (dcreglambda > 0.0)
          {
            Objective.AddObjective(DCReg, ResRegTrans, dcreglambda, "DCReg",JointObjective::regularization);
          }
        if (magreglambda > 0.0)
          {
            Objective.AddObjective(MagReg, SuscepRegTrans, magreglambda, "MagReg",JointObjective::regularization);
          }
      }

    void SetupTomoDCMagCoupling::SetupModelVector(const po::variables_map &vm, jif3D::rvec &InvModel,
            const jif3D::ThreeDModelBase &ModelGeometry,
            const jif3D::ThreeDSeismicModel &SeisMod,
            const jif3D::ThreeDDCResistivityModel &DCMod, const jif3D::ThreeDMagneticModel &MagMod,
            jif3D::JointObjective &Objective,
            boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {

    	SetupCrossGradModel(InvModel, ModelGeometry, SeisMod, DCMod, MagMod,
            Objective, Regularization, substart);
      }
  }
