//============================================================================
// Name        : SetupGravity.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include "SetupGravity.h"
#include "../Gravity/GravityObjective.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"

namespace jiba
  {

    SetupGravity::SetupGravity() :
      scalrelerr(0.02), ftgrelerr(0.02), scalminerr(0.0), ftgminerr(1e-9)
      {

      }

    SetupGravity::~SetupGravity()
      {

      }

    po::options_description SetupGravity::SetupOptions()
      {
        po::options_description desc("Gravity options");

        desc.add_options()("gpu", "Perform gravity calculation on GPU")(
            "scalrelerr", po::value(&scalrelerr)->default_value(0.02),
            "The relative error for the scalar gravity data")("scalminerr",
            po::value(&scalminerr)->default_value(0.0),
            "The minimum absolute error for the scalar gravity data")(
            "ftgrelerr", po::value(&scalrelerr)->default_value(0.02),
            "The relative error for the FTG gravity data")("ftgminerr",
            po::value(&ftgminerr)->default_value(1e-9),
            "The minimum absolute error for the FTG gravity data");

        return desc;
      }

    void SetupGravity::SetupObjective(const po::variables_map &vm,
        jiba::JointObjective &Objective, boost::shared_ptr<jiba::GeneralModelTransform> Transform)
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

        jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;
        //if the weight is different from zero
        //we have to read in scalar gravity data
        if (scalgravlambda > 0.0)
          {
            std::string scalgravdatafilename = jiba::AskFilename(
                "Scalar Gravity Data Filename: ");

            jiba::ReadScalarGravityMeasurements(scalgravdatafilename,
                ScalGravData, PosX, PosY, PosZ);
          }

        //if the weight is different from zero
        //we have to read in ftg data
        if (ftglambda > 0.0)
          {
            std::string ftgdatafilename = jiba::AskFilename(
                "FTG Data Filename: ");
            jiba::ReadTensorGravityMeasurements(ftgdatafilename, FTGData, PosX,
                PosY, PosZ);
          }
        //if the inversion includes any type of gravity data
        //we need the model geometry
        if (scalgravlambda > 0.0 || ftglambda > 0.0)
          {
            std::string gravmodelfilename = jiba::AskFilename(
                "Gravity Model Filename: ");
            GravModel.ReadNetCDF(gravmodelfilename);

            GravModel.ClearMeasurementPoints();
            for (size_t i = 0; i < PosX.size(); ++i)
              {
                GravModel.AddMeasurementPoint(PosX.at(i), PosY.at(i),
                    PosZ.at(i));
              }
          }
        //now we setup the objective functions for each type of gravity data
        //the steps are the same for both, create new objective function,
        //set the observed data, model geometry and data error
        //finally output some basic information about the data to the screen
        //to signal the user that something happened.
        if (scalgravlambda > 0.0)
          {
            boost::shared_ptr<jiba::GravityObjective> ScalGravObjective(
                new jiba::GravityObjective(false, wantcuda));
            ScalGravObjective->SetObservedData(ScalGravData);
            ScalGravObjective->SetModelGeometry(GravModel);
            ScalGravObjective->SetDataCovar(jiba::ConstructError(ScalGravData,
                scalrelerr, scalminerr));

            Objective.AddObjective(ScalGravObjective, Transform,
                scalgravlambda, "ScalGrav");
            std::cout << "Scalar Gravity ndata: " << ScalGravData.size()
                << std::endl;
            std::cout << "Scalar Gravity lambda: " << scalgravlambda
                << std::endl;
          }
        if (ftglambda > 0.0)
          {
            boost::shared_ptr<jiba::GravityObjective> FTGObjective(
                new jiba::GravityObjective(true, wantcuda));
            FTGObjective->SetObservedData(FTGData);
            FTGObjective->SetModelGeometry(GravModel);
            FTGObjective->SetDataCovar(jiba::ConstructError(FTGData, ftgrelerr,
                ftgminerr));

            Objective.AddObjective(FTGObjective, Transform, ftglambda, "FTG");
            std::cout << "FTG ndata: " << FTGData.size() << std::endl;
            std::cout << "FTG lambda: " << ftglambda << std::endl;
          }
      }
  }
