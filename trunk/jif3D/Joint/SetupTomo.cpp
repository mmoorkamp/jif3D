//============================================================================
// Name        : SetupTomo.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "SetupTomo.h"
#include <iostream>
#include <boost/make_shared.hpp>

namespace jif3D
  {

    SetupTomo::SetupTomo() :
        pickerr(5.0e-3), CellSize(), tomolambda(0.0)
      {
      }

    SetupTomo::~SetupTomo()
      {
      }

    po::options_description SetupTomo::SetupOptions()
      {
        po::options_description desc("Tomography options");
        desc.add_options()("tomomodel", po::value(&modelfilename),
            "The name of the starting model for the seismic tomography")("tomodata",
            po::value(&datafilename), "The name of the data for the seismic tomography")(
            "tomolambda", po::value(&tomolambda),
            "The weight for the seismic tomography data")("pickerr",
            po::value(&pickerr)->default_value(5e-3),
            "The picking error for the travel time data")("tomofine",
            po::value(&CellSize), "The cell size in m for the refined tomography model")(
            "writerays", "Write out the rays for each seismic forward modelling");
        return desc;
      }

    bool SetupTomo::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::GeneralModelTransform> Transform, double xorigin,
        double yorigin)
      {
        jif3D::TomographyData TomoData;

        tomolambda = 1.0;
        if (!vm.count("tomolamda"))
          {
            std::cout << "Tomography Lambda: ";
            std::cin >> tomolambda;
          }
        if (tomolambda > 0.0)
          {
            if (!vm.count("modelfilename"))
              {
                //first we read in the starting model and the measured data
                modelfilename = jif3D::AskFilename(
                    "Tomography inversion model Filename: ");
              }
            //we read in the starting modelfile
            //the starting model does not necessarily obey the gridding rules for seismic data
            //we can fix this with a grid refinement model
            TomoModel.ReadNetCDF(modelfilename, false);
            //write out the starting model as a .vtk file for plotting
            TomoModel.WriteVTK(modelfilename + ".vtk");

            if (!vm.count("tomodata"))
              {
                //get the name of the file containing the data and read it in
                datafilename = jif3D::AskFilename("Tomography Data Filename: ");
              }
            //read in data
            TomoData.ReadNetCDF(datafilename);

            if (xorigin != 0.0 || yorigin != 0.0)
              {
                TomoModel.SetOrigin(xorigin, yorigin, 0.0);
              }

            jif3D::TomographyCalculator Calculator;
            TomoObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<jif3D::TomographyCalculator>>(Calculator);

            TomoObjective->SetObservedData(TomoData);
            TomoObjective->SetCoarseModelGeometry(TomoModel);
            //we assume the same error for all measurements
            //this is either the default value set in the constructor
            //or set by the user
            std::vector<double> TomoError = ConstructError(TomoData.GetData(), TomoData.GetErrors(), 0.0, pickerr);
            TomoObjective->SetDataError(TomoError);

            if (vm.count("tomofine") && CellSize > 0.0)
              {
                const double xmax = TomoModel.GetXCoordinates().back();
                const double ymax = TomoModel.GetYCoordinates().back();
                const double zmax = TomoModel.GetZCoordinates().back();
                const double xmin = TomoModel.GetXCoordinates().front();
                const double ymin = TomoModel.GetYCoordinates().front();
                const double zmin = TomoModel.GetZCoordinates().front();
                const double xextent = xmax - TomoModel.GetXCoordinates().front();
                const double yextent = ymax - TomoModel.GetYCoordinates().front();
                const double zextent = zmax - TomoModel.GetZCoordinates().front();
                const int nx = round(xextent / CellSize);
                const int ny = round(yextent / CellSize);
                const int nz = round(zextent / CellSize);
                //if the finely discretized grid does not fit into the inversion grid
                //with a tolerance of more than 10cm
                if (std::abs(xmin + nx * CellSize - xmax) > 0.1)
                  {
                    throw jif3D::FatalException(
                        "Refined grid does not fit in x-direction", __FILE__, __LINE__);
                  }
                if (std::abs(ymin + ny * CellSize - ymax) > 0.1)
                  {
                    throw jif3D::FatalException(
                        "Refined grid does not fit in y-direction", __FILE__, __LINE__);
                  }
                if (std::abs(zmin + nz * CellSize - zmax) > 0.1)
                  {
                    throw jif3D::FatalException(
                        "Refined grid does not fit in x-direction", __FILE__, __LINE__);
                  }

                jif3D::ThreeDSeismicModel TomoFineGeometry;
                TomoFineGeometry.SetCellSize(CellSize, nx, ny, nz);
                TomoFineGeometry.SetOrigin(xmin,ymin,zmin);
                std::cout << "Refined Model has " << nx << " * " << ny << " * " << nz
                    << "cells\n";
                //copy measurement configuration to refined model
                TomoObjective->SetFineModelGeometry(TomoFineGeometry);
              }
            Objective.AddObjective(TomoObjective, Transform, tomolambda, "Tomo",
                JointObjective::datafit);
            std::cout << "Tomo ndata: " << TomoData.GetData().size() << std::endl;
            std::cout << "Tomo lambda: " << tomolambda << std::endl;
          }
        //indicate whether we added a tomography objective
        return (tomolambda > 0.0);
      }

  }
