//============================================================================
// Name        : SetupTomo.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include "../Tomo/TomographyObjective.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "SetupTomo.h"
#include <iostream>

namespace jiba
  {

    SetupTomo::SetupTomo()
      {
      }

    SetupTomo::~SetupTomo()
      {
      }

    po::options_description SetupTomo::SetupOptions()
      {
        po::options_description desc("Tomography options");
        return desc;
      }

    void SetupTomo::SetupObjective(const po::variables_map &vm,
        jiba::JointObjective &Objective, ThreeDSeismicModel &StartModel, boost::shared_ptr<jiba::GeneralModelTransform> Transform)
      {
        jiba::rvec TomoData;

        //first we read in the starting model and the measured data
        std::string modelfilename = jiba::AskFilename(
            "Tomography inversion model Filename: ");
        //we read in the starting modelfile
        jiba::ThreeDSeismicModel TomoFineGeometry;
        StartModel.ReadNetCDF(modelfilename);
        StartModel.WriteVTK(modelfilename + ".vtk");
        std::string tomogeometryfilename = jiba::AskFilename(
            "Tomography forward geometry filename: ");
        TomoFineGeometry.ReadNetCDF(tomogeometryfilename);

        //get the name of the file containing the data and read it in
        std::string tomodatafilename = jiba::AskFilename(
            "Tomography Data Filename: ");

        //read in data
        jiba::ReadTraveltimes(tomodatafilename, TomoData, StartModel);
        //copy measurement configuration to refined model
        TomoFineGeometry.CopyMeasurementConfigurations(StartModel);

        boost::shared_ptr<jiba::TomographyObjective> TomoObjective(
            new jiba::TomographyObjective());
        TomoObjective->SetObservedData(TomoData);
        TomoObjective->SetFineModelGeometry(TomoFineGeometry);
        TomoObjective->SetCoarseModelGeometry(StartModel);
        jiba::rvec TomoCovar(TomoData.size());
        //we assume a general error of 5 ms for the seismic data
        std::fill(TomoCovar.begin(), TomoCovar.end(), 5.0e-3);
        TomoObjective->SetDataCovar(TomoCovar);

        double tomolambda = 1.0;
        std::cout << "Tomography Lambda: ";
        std::cin >> tomolambda;

        if (tomolambda > 0.0)
          {
            Objective.AddObjective(TomoObjective, Transform,
                tomolambda);
            std::cout << "Tomo ndata: " << TomoData.size() << std::endl;
            std::cout << "Tomo lambda: " << tomolambda << std::endl;
          }
      }

  }
