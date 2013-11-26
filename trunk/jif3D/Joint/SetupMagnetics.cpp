//============================================================================
// Name        : SetupMagnetics.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupMagnetics.h"
#include "../Inversion/ModelTransforms.h"
#include "../Magnetics/ReadWriteMagneticData.h"
#include "../Magnetics/OMPMagneticImp.h"
#include "../Magnetics/MagneticTransforms.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"

namespace jif3D
  {

    SetupMagnetics::SetupMagnetics() :
        inclination(0.0), declination(0.0), fieldstrength(1.0), relerr(0.02), minerr(0.0)
      {

      }

    SetupMagnetics::~SetupMagnetics()
      {

      }

    po::options_description SetupMagnetics::SetupOptions()
      {
        po::options_description desc("Magnetics options");

        desc.add_options()("gpu", "Perform Magnetics calculation on GPU")("magrelerr",
            po::value(&relerr)->default_value(0.02),
            "The relative error for the Magnetics data")("magminerr",
            po::value(&minerr)->default_value(0.0),
            "The minimum absolute error for the Magnetics data")("inclination",
            po::value(&inclination)->default_value(0.0),
            "The inclination of the magnetic field in radian")("declination",
            po::value(&declination)->default_value(0.0),
            "The declination of the magnetic field in radian")("fieldstrength",
            po::value(&fieldstrength)->default_value(1.0),
            "The strength of the inducing magnetic field in nT");

        return desc;
      }

    bool SetupMagnetics::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::GeneralModelTransform> &Transform, double xorigin,
        double yorigin, boost::filesystem::path TempDir)
      {
        //if we want to use CUDA for forward modeling
        //we set a variable for easier access later and print
        //that information to the screen
        bool wantcuda = false;
        if (vm.count("gpu"))
          {
            std::cout << "Using GPU" << "\n";
            wantcuda = true;
          }

        jif3D::rvec MagData, MagError;
        double maglambda = 1.0;

        //we first ask for the weights for scalar and tensor Magnetics
        //as we might not have to do very much if the weights are zero
        std::cout << "Magnetics Lambda: ";
        std::cin >> maglambda;

        jif3D::ThreeDMagneticModel::tMeasPosVec PosX, PosY, PosZ;

        //if the weight is different from zero
        //we have to read in scalar Magnetics data
        if (maglambda > 0.0)
          {
            std::string magdatafilename = jif3D::AskFilename(
                "Total field magnetic Data Filename: ");
            jif3D::ReadTotalFieldMagneticMeasurements(magdatafilename, MagData, PosX,
                PosY, PosZ, MagError);
            std::string magmodelfilename = jif3D::AskFilename(
                "Magnetics Model Filename: ");
            Model.ReadNetCDF(magmodelfilename);
            Model.ClearMeasurementPoints();

            for (size_t i = 0; i < PosX.size(); ++i)
              {
                Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
              }

            Model.SetOrigin(xorigin, yorigin, 0.0);

            if (Transform.get() == NULL)
              {
                jif3D::rvec RefVec(Model.GetSusceptibilities().num_elements(), 1.0);
                Transform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                    new jif3D::NormalizeTransform(RefVec));
              }

            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<
                jif3D::ThreeDGravMagImplementation<jif3D::ThreeDMagneticModel> > Implementation;
            if (wantcuda)
              {
                throw jif3D::FatalException("No GPU support, yet !");
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::ThreeDMagneticModel> >(
                    new jif3D::OMPMagneticImp(inclination, declination, fieldstrength));
              }
            boost::shared_ptr<CalculatorType> Calculator(
                new CalculatorType(Implementation, TempDir));

            Calculator->SetDataTransform(
                boost::shared_ptr<jif3D::TotalField>(new jif3D::TotalField));
            MagObjective =
                boost::shared_ptr<jif3D::ThreeDModelObjective<CalculatorType> >(
                    new jif3D::ThreeDModelObjective<CalculatorType>(*Calculator));
            MagObjective->SetObservedData(MagData);
            MagObjective->SetCoarseModelGeometry(Model);
            jif3D::rvec Error(jif3D::ConstructError(MagData, MagError, relerr, minerr));
            std::cout << " MagData: " << MagData << std::endl;
            std::cout << " MagError: " << Error << std::endl;
            MagObjective->SetDataError(Error);

            Objective.AddObjective(MagObjective, Transform, maglambda, "Magnetics",
                JointObjective::datafit);
            std::cout << " Magnetics ndata: " << MagData.size() << std::endl;
            std::cout << " Magnetics lambda: " << maglambda << std::endl;
          }

        //indicate whether we added a Magnetics objective function
        //this way the caller can do additional consistency checks
        return (maglambda > 0.0);
      }
  }
