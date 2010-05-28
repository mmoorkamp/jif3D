//============================================================================
// Name        : SetupMT.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../MT/X3DObjective.h"
#include "../MT/ReadWriteImpedances.h"
#include "../ModelBase/EqualGeometry.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "SetupMT.h"

namespace jiba
  {

    SetupMT::SetupMT()
      {
      }

    SetupMT::~SetupMT()
      {
      }

    po::options_description SetupMT::SetupOptions()
      {
        po::options_description desc("MT options");
        return desc;
      }

    void SetupMT::SetupObjective(const po::variables_map &vm,
        jiba::JointObjective &Objective, const ThreeDSeismicModel &StartModel,
        boost::shared_ptr<jiba::GeneralModelTransform> Transform)
      {
        double mtlambda = 1.0;
        std::cout << "MT Lambda: ";
        std::cin >> mtlambda;

        std::string mtdatafilename = jiba::AskFilename("MT data filename: ");
        std::string mtmodelfilename = jiba::AskFilename("MT Model Filename: ");

        //read in MT data
        std::vector<double> MTXPos, MTYPos, MTZPos, Frequencies;
        jiba::rvec MTData;
        jiba::ReadImpedancesFromNetCDF(mtdatafilename, Frequencies, MTXPos,
            MTYPos, MTZPos, MTData);

        MTModel.ReadNetCDF(mtmodelfilename);
        if (!EqualGridGeometry(MTModel, StartModel))
          {
            throw jiba::FatalException(
                "MT model does not have the same geometry as starting model");
          }
        MTModel.ClearMeasurementPoints();
        std::copy(Frequencies.begin(), Frequencies.end(), std::back_inserter(
            MTModel.SetFrequencies()));
        for (size_t i = 0; i < MTXPos.size(); ++i)
          {
            MTModel.AddMeasurementPoint(MTXPos[i], MTYPos[i], MTZPos[i]);
          }
        if (mtlambda > 0.0)
          {
            boost::shared_ptr<jiba::X3DObjective> MTObjective(
                new jiba::X3DObjective());
            MTObjective->SetModelGeometry(MTModel);
            MTObjective->SetObservedData(MTData);
            MTObjective->SetDataCovar(jiba::ConstructMTError(MTData, 0.02));

            Objective.AddObjective(MTObjective, Transform, mtlambda, "MT");
            std::cout << "MT ndata: " << MTData.size() << std::endl;
            std::cout << "MT lambda: " << mtlambda << std::endl;
          }
      }
  }
