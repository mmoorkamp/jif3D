//============================================================================
// Name        : SetupDCCoupling.h
// Author      : July 3, 2014
// Version     : 
// Copyright   : 2014, zhanjieshi
//============================================================================

#include "SetupDCCoupling.h"
#include "../Global/FileUtil.h"
#include "../Inversion/ModelTransforms.h"


namespace ublas = boost::numeric::ublas;

namespace jif3D
  {
    void SetupDCModelCovar(jif3D::rvec &Covar, const jif3D::rvec &InvModel,
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

    SetupDCCoupling::SetupDCCoupling() :
        mincond(0.0), maxcond(0.0), minslow(0.0), maxslow(0.0), minres(0.0), maxres(
            0.0), res_a(0.0), res_b(0.0), res_c(0.0), cond_a(0.0), cond_b(0.0), cond_c(0.0), ResReplace(
            0.0), CondReplace(3.0)
      {
      }

    SetupDCCoupling::~SetupDCCoupling()
      {
      }

    po::options_description SetupDCCoupling::SetupOptions()
      {
        po::options_description desc("Coupling options");
        desc.add_options()("minslow",
            po::value(&minslow)->default_value(5e-4))("maxslow",
            po::value(&maxslow)->default_value(0.01))("mincond",
            po::value(&mincond)->default_value(2e-3))("maxcond",
            po::value(&maxcond)->default_value(10))("minres",
            po::value(&minres)->default_value(0.1))("maxres",
            po::value(&maxres)->default_value(5e2))("res_a",
            po::value(&res_a)->default_value(-9.23e-7),
            "The quadratic term of the velocity-conductivity relationship")("res_b",
            po::value(&res_b)->default_value(4.2e-3),
            "The linear term of the velocity-conductivity relationship")("res_c",
            po::value(&res_c)->default_value(1.5),
            "The constant term of the velocity-conductivity relationship")("cond_a",
            po::value(&cond_a)->default_value(-9.23e-7),
            "The quadratic term of the velocity-conductivity relationship")("cond_b",
            po::value(&cond_b)->default_value(4.2e-3),
            "The linear term of the velocity-conductivity relationship")("cond_c",
            po::value(&cond_c)->default_value(1.5),
            "The constant term of the velocity-conductivity relationship")("relmodel",
            po::value(&RelModelName),
            "Name of a model file that specifies where to apply the parameter relationship")(
            "ResReplace", po::value(&ResReplace)->default_value(200.0),
            "Value to use for Density where relationship is not valid")("CondReplace",
            po::value(&CondReplace)->default_value(3.3),
            "Value to use for Conductivity where relationship is not valid");

        return desc;
      }

    void SetupDCCoupling::SetupTransforms(const po::variables_map &vm,
        ThreeDSeismicModel &GeometryModel,
        boost::shared_ptr<jif3D::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jif3D::GeneralModelTransform> &DCTransform,
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
                    GeometryModel.GetData().shape()[1], GeometryModel.GetData().shape()[2]));
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

            DCTransform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::ResistivityTransform(TomoTransform, RelModel, ResReplace,
                    res_a, res_b, res_c));
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
        ResTrans = DCTransform;
      }

    void SetupDCCoupling::SetupFixedCouplingModel(jif3D::rvec &InvModel,
        const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDSeismicModel &SeisMod, const jif3D::ThreeDDCResistivityModel &DCMod,
        const jif3D::ThreeDMTModel &MTMod, jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {
        const size_t ngrid = ModelGeometry.GetNModelElements();

        //if we want direct parameter coupling
        //we just need a single regularization term
        //as we invert for slowness, the regularization also works on slowness
        if (SeisMod.GetNModelElements() != ngrid)
          {
            throw jif3D::FatalException(
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
            jif3D::rvec TomoCovar(3 * ngrid);
            SetupDCModelCovar(TomoCovar, InvModel, Regularization->GetDataError(), ngrid);
            Regularization->SetDataError(TomoCovar);
          }
        if (substart)
          {

            jif3D::rvec RefModel = SlowTrans->PhysicalToGeneralized(InvModel);
            RefModel = SlowRegTrans->GeneralizedToPhysical(RefModel);
            Regularization->SetReferenceModel(RefModel);
          }
        Objective.AddObjective(Regularization, SlowRegTrans, reglambda, "Reg",JointObjective::regularization);
        InvModel = SlowTrans->PhysicalToGeneralized(InvModel);
      }

    void SetupDCCoupling::SetupModelVector(const po::variables_map &vm,
        jif3D::rvec &InvModel, const jif3D::ThreeDModelBase &ModelGeometry,
        const jif3D::ThreeDSeismicModel &SeisMod, const jif3D::ThreeDDCResistivityModel &DCMod,
        const jif3D::ThreeDMTModel &MTMod, jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart)
      {
         //depending on the type of coupling the model vector looks quite different
        //and we have to setup a number of different things
        //so we separate it into different functions even though it means
        //passing around many parameters
          SetupFixedCouplingModel(InvModel, ModelGeometry, SeisMod, DCMod, MTMod,
              Objective, Regularization, substart);
      }
  }
