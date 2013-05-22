//============================================================================
// Name        : SetupTomo.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Tomo/ReadWriteTomographyData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "SetupTomo.h"
#include <iostream>

namespace jif3D
  {

    SetupTomo::SetupTomo() :
        pickerr(5.0e-3)
      {
      }

    SetupTomo::~SetupTomo()
      {
      }

    po::options_description SetupTomo::SetupOptions()
      {
        po::options_description desc("Tomography options");
        desc.add_options()("pickerr", po::value(&pickerr)->default_value(5e-3),
            "The picking error for the travel time data")("tomofine",
            po::value(&FineModelName),
            "The name for the model with the MT forward geometry")("writerays",
            "Write out the rays for each seismic forward modelling");
        return desc;
      }

    bool SetupTomo::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::GeneralModelTransform> Transform, double xorigin,
        double yorigin)
      {
        jif3D::rvec TomoData;

        double tomolambda = 1.0;
        std::cout << "Tomography Lambda: ";
        std::cin >> tomolambda;

        if (tomolambda > 0.0)
          {
            //first we read in the starting model and the measured data
            std::string modelfilename = jif3D::AskFilename(
                "Tomography inversion model Filename: ");
            //we read in the starting modelfile
            //the starting model does not necessarily obey the gridding rules for seismic data
            //we can fix this with a grid refinement model
            TomoModel.ReadNetCDF(modelfilename, false);
            //write out the starting model as a .vtk file for plotting
            TomoModel.WriteVTK(modelfilename + ".vtk");

            //get the name of the file containing the data and read it in
            std::string tomodatafilename = jif3D::AskFilename(
                "Tomography Data Filename: ");

            //read in data
            jif3D::ReadTraveltimes(tomodatafilename, TomoData, TomoModel);
            TomoModel.SetOrigin(xorigin, yorigin, 0.0);
            bool writerays = false;
            if (vm.count("writerays"))
              {
                writerays = true;
              }
            jif3D::TomographyCalculator Calculator(writerays);

            TomoObjective = boost::shared_ptr<
                jif3D::ThreeDModelObjective<jif3D::TomographyCalculator> >(
                new jif3D::ThreeDModelObjective<jif3D::TomographyCalculator>(Calculator));
            TomoObjective->SetObservedData(TomoData);
            TomoObjective->SetCoarseModelGeometry(TomoModel);
            //we assume the same error for all measurements
            //this is either the default value set in the constructor
            //or set by the user
            jif3D::rvec TomoError(TomoData.size());
            std::fill(TomoError.begin(), TomoError.end(), pickerr);
            TomoObjective->SetDataError(TomoError);

            if (vm.count("tomofine"))
              {
                jif3D::ThreeDSeismicModel TomoFineGeometry;
                TomoFineGeometry.ReadNetCDF(FineModelName);
                //copy measurement configuration to refined model
                TomoFineGeometry.CopyMeasurementConfigurations(TomoModel);
                TomoObjective->SetFineModelGeometry(TomoFineGeometry);
              }
            Objective.AddObjective(TomoObjective, Transform, tomolambda, "Tomo");
            std::cout << "Tomo ndata: " << TomoData.size() << std::endl;
            std::cout << "Tomo lambda: " << tomolambda << std::endl;
          }
        //indicate whether we added a tomography objective
        return (tomolambda > 0.0);
      }

  }
