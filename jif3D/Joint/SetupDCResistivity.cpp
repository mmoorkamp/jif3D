//============================================================================
// Name        : SetupDCResistivity.h
// Author      : 12 Wen 2021
// Version     :
// Copyright   : 2021, zhanjie and mm489
//============================================================================

#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include <iostream>
#include <boost/make_shared.hpp>
#include "SetupDCResistivity.h"
#include "../DCResistivity/DCResistivityData.h"

namespace jif3D
  {
    SetupDCResistivity::SetupDCResistivity() :
		relerr(0.02), minerr(1e-3), CellSize(), dclambda(0.0)
      {
      }

    SetupDCResistivity::~SetupDCResistivity()
      {
      }



    po::options_description SetupDCResistivity::SetupOptions()
      {
        po::options_description desc("DC resistivity options");

        desc.add_options()("dcmodel", po::value(&dcmodelfilename),
                "The name of the starting model for the DCResistivity")("dcdata",
                po::value(&dcdatafilename), "The name of the data for the DCResistivity")
				("dcrelerr", po::value(&relerr)->default_value(0.02),
            "The relative error for the DC resistivity data")("dcminerr",
            po::value(&minerr)->default_value(1e-3))("dclambda", po::value(&dclambda),
                    "The weight for the DCResistivity data")("dcfine",
            po::value(&CellSize), "The cell size in m for the refined DC model");

        return desc;
      }

    bool SetupDCResistivity::SetupObjective(const po::variables_map &vm,
       jif3D::JointObjective &Objective,
       boost::shared_ptr<jif3D::GeneralModelTransform> Transform, double xorigin,
	   double yorigin)
      {
    	jif3D::DCResistivityData DCData;

        dclambda = 1.0;
        if (!vm.count("dclamda"))
          {
            std::cout << "DCResistivity Lambda: ";
            std::cin >> dclambda;
          }
        if (dclambda > 0.0)
          {
            if (!vm.count("dcmodelfilename"))
              {
                //first we read in the starting model of DCResisitivity
            	dcmodelfilename = jif3D::AskFilename(
                    "DC Resistivity Model Filename: ");
              }
            //we read in the starting modelfile
            //the starting model does not necessarily obey the gridding rules for seismic data
            //we can fix this with a grid refinement model
            DCModel.ReadNetCDF(dcmodelfilename);
            //write out the starting model as a .vtk file for plotting
            DCModel.WriteVTK(dcmodelfilename + ".vtk");

            if (!vm.count("dcdata"))
              {
                //get the name of the file containing the data and read it in
            	dcdatafilename = jif3D::AskFilename("DC Resistivity Data Filename: ");
              }
            //read in data
            DCData.ReadNetCDF(dcdatafilename);

            if (xorigin != 0.0 || yorigin != 0.0)
              {
            	DCModel.SetOrigin(xorigin, yorigin, 0.0);
              }

            jif3D::DCResistivityCalculator DCCalculator;
            DCObjective = boost::make_shared<
            		jif3D::ThreeDModelObjective<jif3D::DCResistivityCalculator>>(DCCalculator);

            DCObjective->SetObservedData(DCData);
            DCObjective->SetCoarseModelGeometry(DCModel);
            //we assume the same error for all measurements
            //this is either the default value set in the constructor
            //or set by the user
            std::vector<double> DCError = ConstructError(DCData.GetData(), DCData.GetErrors(), relerr, minerr);
            DCObjective->SetDataError(DCError);

            if (vm.count("dcfine") && CellSize > 0.0)
              {
            	const double xmax = DCModel.GetXCoordinates().back();
                const double ymax = DCModel.GetYCoordinates().back();
                const double zmax = DCModel.GetZCoordinates().back();
                const double xextent = xmax - DCModel.GetXCoordinates().front();
                const double yextent = ymax - DCModel.GetYCoordinates().front();
                const double zextent = zmax - DCModel.GetZCoordinates().front();
                const int nx = round(xextent / CellSize);
                const int ny = round(yextent / CellSize);
                const int nz = round(zextent / CellSize);
                //if the finely discretized grid does not fit into the inversion grid
                //with a tolerance of more than 10cm
                if (std::abs(nx * CellSize - xmax) > 0.1)
                  {
                    throw jif3D::FatalException(
                        "Refined grid does not fit in x-direction", __FILE__, __LINE__);
                  }
                if (std::abs(ny * CellSize - ymax) > 0.1)
                  {
                    throw jif3D::FatalException(
                        "Refined grid does not fit in y-direction", __FILE__, __LINE__);
                  }
                if (std::abs(nz * CellSize - zmax) > 0.1)
                  {
                    throw jif3D::FatalException(
                        "Refined grid does not fit in x-direction", __FILE__, __LINE__);
                  }

                jif3D::ThreeDDCResistivityModel DCFineGeometry;
                DCFineGeometry.SetCellSize(CellSize, nx, ny, nz);
                std::cout << "Refined Model has " << nx << " * " << ny << " * " << nz
                		<< "cells\n";
                //copy measurement configuration to refined model
                DCObjective->SetFineModelGeometry(DCFineGeometry);
              }
            Objective.AddObjective(DCObjective, Transform, dclambda, "DC Resistivity",
            		JointObjective::datafit);
            std::cout << "Resistivity ndata: " << DCData.GetData().size() << std::endl;
            std::cout << "Resistivity lambda: " << dclambda << std::endl;
          }
        //indicate whether we added a tomography objective
        return (dclambda > 0.0);
      }

  }
