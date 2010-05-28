//============================================================================
// Name        : SetupGravity.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include "SetupGravity.h"
#include "../ModelBase/EqualGeometry.h"
#include "../Gravity/GravityObjective.h"
#include "../Gravity/ReadWriteGravityData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"

namespace jiba
  {

    SetupGravity::SetupGravity()
      {

      }

    SetupGravity::~SetupGravity()
      {

      }

    po::options_description SetupGravity::SetupOptions()
      {
        po::options_description desc("Gravity options");

        desc.add_options()("gpu", "Perform gravity calculation on GPU");

        return desc;
      }

    void SetupGravity::SetupObjective(const po::variables_map &vm,
        jiba::JointObjective &Objective, const ThreeDSeismicModel &StartModel,
        boost::shared_ptr<jiba::GeneralModelTransform> Transform)
      {
        bool wantcuda = false;
        if (vm.count("gpu"))
          {
            std::cout << "Using GPU" << "\n";
            wantcuda = true;
          }

        jiba::rvec ScalGravData, FTGData;
        double scalgravlambda = 1.0;
        double ftglambda = 1.0;
        std::cout << "Scalar Gravimetry Lambda: ";
        std::cin >> scalgravlambda;
        std::cout << "FTG Lambda: ";
        std::cin >> ftglambda;

        jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;

       // if (scalgravlambda > 0.0)
          {
            std::string scalgravdatafilename = jiba::AskFilename(
                "Scalar Gravity Data Filename: ");

            jiba::ReadScalarGravityMeasurements(scalgravdatafilename,
                ScalGravData, PosX, PosY, PosZ);
          }

        //if (ftglambda > 0.0)
          {
            std::string ftgdatafilename = jiba::AskFilename(
                "FTG Data Filename: ");
            jiba::ReadTensorGravityMeasurements(ftgdatafilename, FTGData, PosX,
                PosY, PosZ);
          }

        //if (scalgravlambda > 0.0 || ftglambda > 0.0)
          {
            std::string gravmodelfilename = jiba::AskFilename(
                "Gravity Model Filename: ");
            GravModel.ReadNetCDF(gravmodelfilename);
            if (!EqualGridGeometry(GravModel, StartModel))
              {
                throw jiba::FatalException(
                    "Gravity model does not have the same geometry as starting model");
              }
            GravModel.ClearMeasurementPoints();
            for (size_t i = 0; i < PosX.size(); ++i)
              {
                GravModel.AddMeasurementPoint(PosX.at(i), PosY.at(i),
                    PosZ.at(i));
              }
          }
        //if (scalgravlambda > 0.0)
          {
            boost::shared_ptr<jiba::GravityObjective> ScalGravObjective(
                new jiba::GravityObjective(false, wantcuda));
            ScalGravObjective->SetObservedData(ScalGravData);
            ScalGravObjective->SetModelGeometry(GravModel);
            ScalGravObjective->SetDataCovar(jiba::ConstructError(ScalGravData,
                0.0, 5e-7));

            Objective.AddObjective(ScalGravObjective, Transform, scalgravlambda,"ScalGrav");
            std::cout << "Scalar Gravity ndata: " << ScalGravData.size()
                << std::endl;
            std::cout << "Scalar Gravity lambda: " << scalgravlambda
                << std::endl;
          }
        //if (ftglambda > 0.0)
          {
            boost::shared_ptr<jiba::GravityObjective> FTGObjective(
                new jiba::GravityObjective(true, wantcuda));
            FTGObjective->SetObservedData(FTGData);
            FTGObjective->SetModelGeometry(GravModel);
            FTGObjective->SetDataCovar(
                jiba::ConstructError(FTGData, 0.02, 1e-9));

            Objective.AddObjective(FTGObjective, Transform, ftglambda,"FTG");
            std::cout << "FTG ndata: " << FTGData.size() << std::endl;
            std::cout << "FTG lambda: " << ftglambda << std::endl;
          }
      }
  }
