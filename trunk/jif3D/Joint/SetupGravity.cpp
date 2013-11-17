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
            po::value(&scalrelerr)->default_value(0.02),
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

        jif3D::rvec ScalGravData, FTGData;
        double scalgravlambda = 1.0;
        double ftglambda = 1.0;
        //we first ask for the weights for scalar and tensor gravity
        //as we might not have to do very much if the weights are zero
        std::cout << "Scalar Gravimetry Lambda: ";
        std::cin >> scalgravlambda;
        std::cout << "FTG Lambda: ";
        std::cin >> ftglambda;

        jif3D::ThreeDGravityModel::tMeasPosVec ScalPosX, ScalPosY, ScalPosZ;
        jif3D::ThreeDGravityModel::tMeasPosVec FTGPosX, FTGPosY, FTGPosZ;
        //if the weight is different from zero
        //we have to read in scalar gravity data
        if (scalgravlambda > 0.0)
          {
            std::string scalgravdatafilename = jif3D::AskFilename(
                "Scalar Gravity Data Filename: ");

            jif3D::ReadScalarGravityMeasurements(scalgravdatafilename, ScalGravData,
                ScalPosX, ScalPosY, ScalPosZ);
            HaveScal = true;
          }

        //if the weight is different from zero
        //we have to read in ftg data
        if (ftglambda > 0.0)
          {
            std::string ftgdatafilename = jif3D::AskFilename("FTG Data Filename: ");
            jif3D::ReadTensorGravityMeasurements(ftgdatafilename, FTGData, FTGPosX,
                FTGPosY, FTGPosZ);
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

            ScalGravModel.ClearMeasurementPoints();
            FTGGravModel.ClearMeasurementPoints();
            if (ftglambda > 0.0)
              {
                for (size_t i = 0; i < FTGPosX.size(); ++i)
                  {
                    FTGGravModel.AddMeasurementPoint(FTGPosX.at(i), FTGPosY.at(i),
                        FTGPosZ.at(i));
                  }
              }
            if (scalgravlambda > 0.0)
              {
                for (size_t i = 0; i < ScalPosX.size(); ++i)
                  {
                    ScalGravModel.AddMeasurementPoint(ScalPosX.at(i), ScalPosY.at(i),
                        ScalPosZ.at(i));
                  }
              }
            ScalGravModel.SetOrigin(xorigin, yorigin, 0.0);
            FTGGravModel.SetOrigin(xorigin, yorigin, 0.0);
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
        if (scalgravlambda > 0.0)
          {
            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<
                jif3D::ThreeDGravMagImplementation<jif3D::ThreeDGravityModel> > Implementation;
            if (wantcuda)
              {
#ifdef HAVEGPU
                Implementation = boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ThreeDGravityModel> >(
                    new jif3D::ScalarCudaGravityImp);
#else
                throw jif3D::FatalException(
                    "Code has been compiled without GPU support !");
#endif
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::ThreeDGravityModel> >(
                    new jif3D::ScalarOMPGravityImp);
              }
            boost::shared_ptr<CalculatorType> ScalarCalculator(
                new CalculatorType(Implementation, TempDir));

            ScalGravObjective = boost::shared_ptr<
                jif3D::ThreeDModelObjective<CalculatorType> >(
                new jif3D::ThreeDModelObjective<CalculatorType>(*ScalarCalculator));
            ScalGravObjective->SetObservedData(ScalGravData);
            ScalGravObjective->SetCoarseModelGeometry(ScalGravModel);
            ScalGravObjective->SetDataError(
                jif3D::ConstructError(ScalGravData, scalrelerr, scalminerr));

            Objective.AddObjective(ScalGravObjective, Transform, scalgravlambda,
                "ScalGrav",JointObjective::datafit);
            std::cout << "Scalar Gravity ndata: " << ScalGravData.size() << std::endl;
            std::cout << "Scalar Gravity lambda: " << scalgravlambda << std::endl;
          }
        if (ftglambda > 0.0)
          {
            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<
                jif3D::ThreeDGravMagImplementation<jif3D::ThreeDGravityModel> > Implementation;
            if (wantcuda)
              {
#ifdef HAVEGPU
                Implementation = boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::ThreeDGravityModel> >(
                    new jif3D::TensorCudaGravityImp);
#else
                throw jif3D::FatalException(
                    "Code has been compiled without GPU support !");
#endif
              }
            else
              {
                Implementation = boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::ThreeDGravityModel> >(
                    new jif3D::TensorOMPGravityImp);
              }
            boost::shared_ptr<CalculatorType> TensorCalculator(
                new CalculatorType(Implementation, TempDir));

            FTGObjective =
                boost::shared_ptr<jif3D::ThreeDModelObjective<CalculatorType> >(
                    new jif3D::ThreeDModelObjective<CalculatorType>(*TensorCalculator));
            FTGObjective->SetObservedData(FTGData);
            FTGObjective->SetCoarseModelGeometry(FTGGravModel);
            FTGObjective->SetDataError(
                jif3D::ConstructError(FTGData, ftgrelerr, ftgminerr));

            Objective.AddObjective(FTGObjective, Transform, ftglambda, "FTG",JointObjective::datafit);
            std::cout << "FTG ndata: " << FTGData.size() << std::endl;
            std::cout << "FTG lambda: " << ftglambda << std::endl;
          }
        //indicate whether we added a gravity objective function
        //this way the caller can do additional consistency checks
        return (ftglambda > 0.0) || (scalgravlambda > 0.0);
      }
  }
