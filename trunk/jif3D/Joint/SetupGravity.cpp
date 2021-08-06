//============================================================================
// Name        : SetupGravity.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupGravity.h"
#include "../Inversion/ModelTransforms.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Gravity/ThreeDGravityFactory.h"
#include "../Gravity/ScalarGravityData.h"
#include "../Gravity/TensorGravityData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"

namespace jif3D
  {

    SetupGravity::SetupGravity() :
        scalrelerr(0.02), ftgrelerr(0.02), scalminerr(0.0), ftgminerr(1e-9), HaveScal(
            false), HaveFTG(false)
      {

      }

    SetupGravity::~SetupGravity()
      {

      }

    po::options_description SetupGravity::SetupOptions()
      {
        po::options_description desc("Gravity options");

        desc.add_options()("gpu", "Perform gravity calculation on GPU")("scalrelerr",
            po::value(&scalrelerr)->default_value(0.02),
            "The relative error for the scalar gravity data")("scalminerr",
            po::value(&scalminerr)->default_value(0.0),
            "The minimum absolute error for the scalar gravity data")("ftgrelerr",
            po::value(&ftgrelerr)->default_value(0.02),
            "The relative error for the FTG gravity data")("ftgminerr",
            po::value(&ftgminerr)->default_value(1e-9),
            "The minimum absolute error for the FTG gravity data");

        return desc;
      }

    bool SetupGravity::SetupObjective(const po::variables_map &vm,
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

        jif3D::ScalarGravityData ScalGravData;
        jif3D::TensorGravityData FTGData;
        double scalgravlambda = 1.0;
        double ftglambda = 1.0;
        //we first ask for the weights for scalar and tensor gravity
        //as we might not have to do very much if the weights are zero
        std::cout << "Scalar Gravimetry Lambda: ";
        std::cin >> scalgravlambda;
        std::cout << "FTG Lambda: ";
        std::cin >> ftglambda;

        //if the weight is different from zero
        //we have to read in scalar gravity data
        if (scalgravlambda > JointObjective::MinWeight)
          {
            std::string scalgravdatafilename = jif3D::AskFilename(
                "Scalar Gravity Data Filename: ");

            ScalGravData.ReadNetCDF(scalgravdatafilename);
            ScalGravData.WriteVTK(scalgravdatafilename);
            HaveScal = true;
          }

        //if the weight is different from zero
        //we have to read in ftg data
        if (ftglambda > JointObjective::MinWeight)
          {
            std::string ftgdatafilename = jif3D::AskFilename("FTG Data Filename: ");
            FTGData.ReadNetCDF(ftgdatafilename);
            HaveFTG = true;
          }
        //if the inversion includes any type of gravity data
        //we need the model geometry
        if (scalgravlambda > 0.0 || ftglambda > 0.0)
          {
            std::string gravmodelfilename = jif3D::AskFilename(
                "Gravity Model Filename: ");
            ScalGravModel.ReadNetCDF(gravmodelfilename);
            FTGGravModel.ReadNetCDF(gravmodelfilename);
            if (xorigin != 0.0 || yorigin != 0.0)
              {
                ScalGravModel.SetOrigin(xorigin, yorigin, 0.0);
                FTGGravModel.SetOrigin(xorigin, yorigin, 0.0);
              }
          }
        if (Transform.get() == NULL)
          {
            jif3D::rvec RefVec(ScalGravModel.GetDensities().num_elements());
            std::fill(RefVec.begin(), RefVec.end(), 1.0);
            Transform = boost::shared_ptr<jif3D::GeneralModelTransform>(
                new jif3D::NormalizeTransform(RefVec));
          }
        //now we setup the objective functions for each type of gravity data
        //the steps are the same for both, create new objective function,
        //set the observed data, model geometry and data error
        //finally output some basic information about the data to the screen
        //to signal the user that something happened.
        if (scalgravlambda > JointObjective::MinWeight)
          {
            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ScalarGravityData> > Implementation;
            if (wantcuda)
              {
#ifdef HAVEGPU
                Implementation = boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ScalarGravityData> >(
                    new jif3D::ScalarCudaGravityImp);
#else
                throw jif3D::FatalException(
                    "Code has been compiled without GPU support !", __FILE__, __LINE__);
#endif
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::ScalarGravityData> >(
                    new jif3D::ScalarOMPGravityImp);
              }
            boost::shared_ptr<ScalarCalculatorType> ScalarCalculator(
                new ScalarCalculatorType(Implementation, TempDir));

            ScalGravObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<ScalarCalculatorType> >(*ScalarCalculator);
            ScalGravObjective->SetObservedData(ScalGravData);
            ScalGravObjective->SetCoarseModelGeometry(ScalGravModel);
            ScalGravObjective->SetDataError(
                jif3D::ConstructError(ScalGravData.GetData(), ScalGravData.GetErrors(),
                    scalrelerr, scalminerr));

            Objective.AddObjective(ScalGravObjective, Transform, scalgravlambda,
                "ScalGrav", JointObjective::datafit);
            std::cout << "Scalar Gravity ndata: " << ScalGravData.GetData().size() << std::endl;
            std::cout << "Scalar Gravity lambda: " << scalgravlambda << std::endl;
          }
        if (ftglambda > JointObjective::MinWeight)
          {
            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::TensorGravityData> > Implementation;
            if (wantcuda)
              {
#ifdef HAVEGPU
                Implementation = boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::TensorGravityData> >(
                    new jif3D::TensorCudaGravityImp);
#else
                throw jif3D::FatalException(
                    "Code has been compiled without GPU support !", __FILE__, __LINE__);
#endif
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::TensorGravityData> >(
                    new jif3D::TensorOMPGravityImp);
              }
            boost::shared_ptr<TensorCalculatorType> TensorCalculator(
                new TensorCalculatorType(Implementation, TempDir));

            FTGObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<TensorCalculatorType> >(*TensorCalculator);
            FTGObjective->SetObservedData(FTGData);
            FTGObjective->SetCoarseModelGeometry(FTGGravModel);
            FTGObjective->SetDataError(
                jif3D::ConstructError(FTGData.GetData(), FTGData.GetErrors(), ftgrelerr,
                    ftgminerr));

            Objective.AddObjective(FTGObjective, Transform, ftglambda, "FTG",
                JointObjective::datafit);
            std::cout << "FTG ndata: " << FTGData.GetData().size() << std::endl;
            std::cout << "FTG lambda: " << ftglambda << std::endl;
          }
        //indicate whether we added a gravity objective function
        //this way the caller can do additional consistency checks
        return (ftglambda > JointObjective::MinWeight) || (scalgravlambda > JointObjective::MinWeight);
      }
  }
