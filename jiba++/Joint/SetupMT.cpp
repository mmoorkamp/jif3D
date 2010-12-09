//============================================================================
// Name        : SetupMT.cpp
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include "../MT/ReadWriteImpedances.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "SetupMT.h"
#include <algorithm>

namespace jiba
  {

    SetupMT::SetupMT() :
      relerr(0.02)
      {
      }

    SetupMT::~SetupMT()
      {
      }

    po::options_description SetupMT::SetupOptions()
      {
        //set the program option description for the MT part of the inversion
        po::options_description desc("MT options");
        desc.add_options()("mtrelerr", po::value(&relerr)->default_value(0.02),
            "The relative error for the MT data")("mtfine", po::value(
            &FineModelName),
            "The name for the model with the MT forward geometry");
        return desc;
      }

    bool SetupMT::SetupObjective(const po::variables_map &vm,
        jiba::JointObjective &Objective, boost::shared_ptr<
            jiba::GeneralModelTransform> Transform)
      {
        //first we ask the user a few questions
        //these are all values that are likely to change from run to run
        // so we do not want them as command line arguments or in the
        //configuration file
        double mtlambda = 1.0;
        std::cout << "MT Lambda: ";
        std::cin >> mtlambda;
        //if lambda is negative we assume that no MT inversion is wanted
        //and there is nothing to do here
        if (mtlambda > 0.0)
          {
            //for inversion we need some data, so we ask for the filename
            std::string mtdatafilename =
                jiba::AskFilename("MT data filename: ");

            std::string mtmodelfilename = jiba::AskFilename(
                "MT Model Filename: ");
            //read in the model and check whether the geometry matches the one
            //of the tomography starting model
            MTModel.ReadNetCDF(mtmodelfilename);

            //read in MT data, the position of the measurement sites, frequencies and impedances
            std::vector<double> MTXPos, MTYPos, MTZPos, Frequencies;
            jiba::rvec MTData, MTError;
            jiba::ReadImpedancesFromNetCDF(mtdatafilename, Frequencies, MTXPos,
                MTYPos, MTZPos, MTData, MTError);

            //set the model object so that we can use it to calculate synthetic data
            // for each observation
            MTModel.ClearMeasurementPoints();
            std::copy(Frequencies.begin(), Frequencies.end(),
                std::back_inserter(MTModel.SetFrequencies()));
            for (size_t i = 0; i < MTXPos.size(); ++i)
              {
                MTModel.AddMeasurementPoint(MTXPos[i], MTYPos[i], MTZPos[i]);
              }
            //setup the objective function for the MT data
            jiba::X3DMTCalculator Calculator;

            MTObjective = boost::shared_ptr<jiba::ThreeDModelObjective<
                jiba::X3DMTCalculator> >(new jiba::ThreeDModelObjective<
                jiba::X3DMTCalculator>(Calculator));
            //if we specified the name for a refined model for forward calculations
            //we read in that model, set the measurement configuration for the observed
            //data and pass it to the objective function
            if (vm.count("mtfine"))
              {
                jiba::X3DModel FineModel;
                FineModel.ReadNetCDF(FineModelName);
                FineModel.CopyMeasurementConfigurations(MTModel);
                MTObjective->SetFineModelGeometry(FineModel);
              }
            MTObjective->SetCoarseModelGeometry(MTModel);
            MTObjective->SetObservedData(MTData);
            jiba::rvec MinErr(jiba::ConstructMTError(MTData, relerr));
            std::transform(MTError.begin(), MTError.end(), MinErr.begin(),
                MTError.begin(), std::max<double>);
            MTObjective->SetDataCovar(MTError);
            //add the MT part to the JointObjective that will be used
            //for the inversion
            Objective.AddObjective(MTObjective, Transform, mtlambda, "MT");
            //output some information to the screen
            //to signal that we added the MT data
            std::cout << "MT ndata: " << MTData.size() << std::endl;
            std::cout << "MT lambda: " << mtlambda << std::endl;
          }
        //return true if we added an MT objective function
        return (mtlambda > 0.0);
      }
  }
