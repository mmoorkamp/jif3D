/*
 * SetupSW.cpp
 *
 *  Created on: 29 Oct 2019
 *      Author: bweise
 */

#include "../SurfaceWaves/SurfaceWaveData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../ModelTransforms/DensPVelTransform.h"
#include "SetupSW.h"
#include <iostream>
#include <vector>
#include <boost/make_shared.hpp>

namespace jif3D
  {

    SetupSW::SetupSW() :
        relerr(2.0e-2), minerr(1), swlambda(0.0)
      {
      }

    SetupSW::~SetupSW()
      {
      }

    po::options_description SetupSW::SetupOptions()
      {
        po::options_description desc("Surface wave tomography options");
        desc.add_options()("swmodel", po::value(&modelfilename),
            "The name of the starting model for the surface wave tomography")("swdata",
            po::value(&datafilename),
            "The name of the data for the surface wave tomography")("swlambda",
            po::value(&swlambda), "The weight for the surface wave tomography data")(
            "relerr", po::value(&relerr)->default_value(2.0e-2),
            "The relative error for the phase delay data")("minerr",
            po::value(&minerr)->default_value(1),
            "The minimal error for the phase delay data");
        return desc;
      }

    bool SetupSW::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::GeneralModelTransform> Transform, double xorigin,
        double yorigin)
      {
        jif3D::SurfaceWaveData SurfaceWaveData;

        swlambda = 1.0;
        if (!vm.count("swlamda"))
          {
            std::cout << "Surface wave tomography Lambda: ";
            std::cin >> swlambda;
          }
        if (swlambda > 0.0)
          {
            if (!vm.count("modelfilename"))
              {
                //first we read in the starting model and the measured data
                modelfilename = jif3D::AskFilename(
                    "Surface wave tomography inversion model Filename: ");
              }
            //we read in the starting modelfile
            SWModel.ReadNetCDF(modelfilename);
            //write out the starting model as a .vtk file for plotting
            SWModel.WriteVTK(modelfilename + ".vtk");

            if (!vm.count("swdata"))
              {
                //get the name of the file containing the data and read it in
                datafilename = jif3D::AskFilename(
                    "Surface wave tomography Data Filename: ");
              }
            //read in data
            SurfaceWaveData.ReadNetCDF(datafilename);

            if (xorigin != 0.0 || yorigin != 0.0)
              {
                SWModel.SetOrigin(xorigin, yorigin, 0.0);
              }
            //dynamic_cast<jif3D::DensPVelTransform*>(Transform.get())->SetBGDens(
            //    SWModel.GetBGDens());
            std::vector<double> SurfaceWaveError = ConstructError(
                SurfaceWaveData.GetData(), SurfaceWaveData.GetErrors(), relerr, minerr);
            jif3D::SurfaceWaveCalculator Calculator;
            Calculator.set_distance_tolerance(1.0);
            Calculator.set_vel_tolerance(0.001);
            Calculator.set_false_east(500000.0);
            Calculator.set_root_search_iterations(50);
            Calculator.set_mode_skip_iterations(2);
            SurfaceWaveData.SetDataAndErrors(SurfaceWaveData.GetData(), SurfaceWaveError);

            SurfaceWaveObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<jif3D::SurfaceWaveCalculator>>(Calculator);

            SurfaceWaveObjective->SetObservedData(SurfaceWaveData);
            SurfaceWaveObjective->SetCoarseModelGeometry(SWModel);
            //we assume the same error for all measurements
            //this is either the default value set in the constructor
            //or set by the user
            SurfaceWaveObjective->SetDataError(SurfaceWaveError);

            Objective.AddObjective(SurfaceWaveObjective, Transform, swlambda, "SWTomo",
                JointObjective::datafit);
            std::cout << "SurfaceWave tomo ndata: " << SurfaceWaveData.GetData().size()
                << std::endl;
            std::cout << "SurfaceWave tomo lambda: " << swlambda << std::endl;
          }
        //indicate whether we added a tomography objective
        return (swlambda > 0.0);
      }

  }

