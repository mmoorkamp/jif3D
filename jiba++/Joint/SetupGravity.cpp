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

namespace jiba
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
        jiba::JointObjective &Objective,
        boost::shared_ptr<jiba::GeneralModelTransform> &Transform, double xorigin,
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

        jiba::rvec ScalGravData, FTGData;
        double scalgravlambda = 1.0;
        double ftglambda = 1.0;
        //we first ask for the weights for scalar and tensor gravity
        //as we might not have to do very much if the weights are zero
        std::cout << "Scalar Gravimetry Lambda: ";
        std::cin >> scalgravlambda;
        std::cout << "FTG Lambda: ";
        std::cin >> ftglambda;

        jiba::ThreeDGravityModel::tMeasPosVec ScalPosX, ScalPosY, ScalPosZ;
        jiba::ThreeDGravityModel::tMeasPosVec FTGPosX, FTGPosY, FTGPosZ;
        //if the weight is different from zero
        //we have to read in scalar gravity data
        if (scalgravlambda > 0.0)
          {
            std::string scalgravdatafilename = jiba::AskFilename(
                "Scalar Gravity Data Filename: ");

            jiba::ReadScalarGravityMeasurements(scalgravdatafilename, ScalGravData,
                ScalPosX, ScalPosY, ScalPosZ);
            HaveScal = true;
          }

        //if the weight is different from zero
        //we have to read in ftg data
        if (ftglambda > 0.0)
          {
            std::string ftgdatafilename = jiba::AskFilename("FTG Data Filename: ");
            jiba::ReadTensorGravityMeasurements(ftgdatafilename, FTGData, FTGPosX,
                FTGPosY, FTGPosZ);
            HaveFTG = true;
          }
        //if the inversion includes any type of gravity data
        //we need the model geometry
        if (scalgravlambda > 0.0 || ftglambda > 0.0)
          {
            std::string gravmodelfilename = jiba::AskFilename("Gravity Model Filename: ");
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
            jiba::rvec RefVec(ScalGravModel.GetDensities().num_elements());
            std::fill(RefVec.begin(), RefVec.end(), 1.0);
            Transform = boost::shared_ptr<jiba::GeneralModelTransform>(
                new jiba::NormalizeTransform(RefVec));
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
            boost::shared_ptr<jiba::ThreeDGravityImplementation> Implementation;
            if (wantcuda)
              {
                Implementation = boost::shared_ptr<jiba::ThreeDGravityImplementation>(
                    new jiba::ScalarCudaGravityImp);
              }
            else
              {
                Implementation = boost::shared_ptr<jiba::ThreeDGravityImplementation>(
                    new jiba::ScalarOMPGravityImp);
              }
            boost::shared_ptr<jiba::DiskGravityCalculator> ScalarCalculator(
                new jiba::DiskGravityCalculator(Implementation, TempDir));

            ScalGravObjective = boost::shared_ptr<
                jiba::ThreeDModelObjective<DiskGravityCalculator> >(
                new jiba::ThreeDModelObjective<DiskGravityCalculator>(*ScalarCalculator));
            ScalGravObjective->SetObservedData(ScalGravData);
            ScalGravObjective->SetCoarseModelGeometry(ScalGravModel);
            ScalGravObjective->SetDataError(
                jiba::ConstructError(ScalGravData, scalrelerr, scalminerr));

            Objective.AddObjective(ScalGravObjective, Transform, scalgravlambda,
                "ScalGrav");
            std::cout << "Scalar Gravity ndata: " << ScalGravData.size() << std::endl;
            std::cout << "Scalar Gravity lambda: " << scalgravlambda << std::endl;
          }
        if (ftglambda > 0.0)
          {
            //we want to set the path for temporary file storage
            //the factory function cannot perform this, so we
            //have to assemble the calculator object ourselves
            boost::shared_ptr<jiba::ThreeDGravityImplementation> Implementation;
            if (wantcuda)
              {
                Implementation = boost::shared_ptr<jiba::ThreeDGravityImplementation>(
                    new jiba::TensorCudaGravityImp);
              }
            else
              {
                Implementation = boost::shared_ptr<jiba::ThreeDGravityImplementation>(
                    new jiba::TensorOMPGravityImp);
              }
            boost::shared_ptr<jiba::DiskGravityCalculator> TensorCalculator(
                new jiba::DiskGravityCalculator(Implementation, TempDir));

            FTGObjective = boost::shared_ptr<
                jiba::ThreeDModelObjective<DiskGravityCalculator> >(
                new jiba::ThreeDModelObjective<DiskGravityCalculator>(*TensorCalculator));
            FTGObjective->SetObservedData(FTGData);
            FTGObjective->SetCoarseModelGeometry(FTGGravModel);
            FTGObjective->SetDataError(
                jiba::ConstructError(FTGData, ftgrelerr, ftgminerr));

            Objective.AddObjective(FTGObjective, Transform, ftglambda, "FTG");
            std::cout << "FTG ndata: " << FTGData.size() << std::endl;
            std::cout << "FTG lambda: " << ftglambda << std::endl;
          }
        //indicate whether we added a gravity objective function
        //this way the caller can do additional consistency checks
        return (ftglambda > 0.0) || (scalgravlambda > 0.0);
      }
  }
