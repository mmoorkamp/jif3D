//============================================================================
// Name        : SetupTomo.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupTomo.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyData.h"
#include "../Inversion/ModelTransforms.h"

#include "../Global/FileUtil.h"
#include "../Global/Noise.h"

#include <boost/make_shared.hpp>
#include <iostream>

namespace jif3D
  {

    const std::string Name = "Slowness";

    SetupTomo::SetupTomo() :
        GeneralDataSetup(Name), pickerr(5.0e-3), CellSize(), tomolambda(0.0), WriteRays(
            false)
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
            "writerays", po::value(&WriteRays)->default_value(false),
            "Write out the rays for each seismic forward modelling")("minslow",
            po::value(&minslow)->default_value(1e-4))("maxslow",
            po::value(&maxslow)->default_value(0.005));
        return desc;
      }

    bool SetupTomo::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective, jif3D::ThreeDModelBase &InversionMesh,
        jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
        std::vector<std::string> &SegmentNames, std::vector<parametertype> &SegmentTypes,
        boost::filesystem::path TempDir)
      {
        const size_t ngrid = InversionMesh.GetData().num_elements();

        tomolambda = 1.0;
        if (!vm.count("tomolamda"))
          {
            std::cout << "Tomography Lambda: ";
            std::cin >> tomolambda;
          }
        if (tomolambda > 0.0)
          {
            if (!vm.count("tomomodel"))
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

            StartingParameters.resize(ngrid);
            std::copy(TomoModel.GetSlownesses().origin(),
                TomoModel.GetSlownesses().origin() + ngrid,
                StartingParameters.begin());

            if (!vm.count("tomodata"))
              {
                //get the name of the file containing the data and read it in
                datafilename = jif3D::AskFilename("Tomography Data Filename: ");
              }
            //read in data
            TomoData.ReadNetCDF(datafilename);
            TomoData.WriteMeasurementPoints(modelfilename + ".rec.vtk");
            TomoData.WriteSourcePoints(modelfilename + ".sor.vtk");

            jif3D::TomographyCalculator Calculator;
            TomoObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<jif3D::TomographyCalculator>>(Calculator);

            TomoObjective->SetObservedData(TomoData);
            TomoObjective->SetCoarseModelGeometry(TomoModel);
            //we assume the same error for all measurements
            //this is either the default value set in the constructor
            //or set by the user
            std::vector<double> TomoError = ConstructError(TomoData.GetData(),
                TomoData.GetErrors(), 0.0, pickerr);
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
                TomoFineGeometry.SetOrigin(xmin, ymin, zmin);
                std::cout << "Refined Model has " << nx << " * " << ny << " * " << nz
                    << "cells\n";
                //copy measurement configuration to refined model
                TomoObjective->SetFineModelGeometry(TomoFineGeometry);
              }
            size_t start = startindices.back();
            size_t end = start + ngrid;
            boost::shared_ptr<jif3D::GeneralModelTransform> TomoTrans =
                boost::make_shared<jif3D::TanhTransform>(minslow, maxslow);
            Transform = boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid, start,
                end, TomoTrans);
            startindices.push_back(end);
            SegmentNames.push_back(Name);
            SegmentTypes.push_back(GeneralDataSetup::gridparameter);

            Objective.AddObjective(TomoObjective, Transform, tomolambda, "Tomo",
                JointObjective::datafit);
            std::cout << "Tomo ndata: " << TomoData.GetData().size() << std::endl;
            std::cout << "Tomo lambda: " << tomolambda << std::endl;
          }
        //indicate whether we added a tomography objective
        return (tomolambda > 0.0);
      }

    void SetupTomo::IterationOutput(const std::string &filename,
        const jif3D::rvec &ModelVector)
      {
        if (tomolambda > JointObjective::MinWeight)
          {
            jif3D::rvec TomoInvModel = Transform->GeneralizedToPhysical(ModelVector);
            std::copy(TomoInvModel.begin(), TomoInvModel.end(),
                TomoModel.SetSlownesses().origin());

            TomoModel.WriteVTK(filename + ".tomo.inv.vtk");
            TomoModel.WriteNetCDF(filename + ".tomo.inv.nc");
            if (WriteRays)
              {
                TomoObjective->GetCalculator().WriteRays(filename + "_rays.vtk");
              }
          }
      }

    void SetupTomo::FinalOutput(const std::string &filename, const jif3D::rvec &FinalModelVector)
      {
        if (tomolambda > JointObjective::MinWeight)
          {
            std::cout << "Writing final slowness models " << std::endl;
            jif3D::rvec TomoInvModel = Transform->GeneralizedToPhysical(FinalModelVector);

            std::copy(TomoInvModel.begin(), TomoInvModel.end(),
                TomoModel.SetSlownesses().origin());

            auto TomoDataVec = TomoObjective->GetSyntheticData();
            auto TomoErrVec = TomoObjective->GetDataError();
            if (TomoDataVec.size() > 0)
              {
                TomoData.SetDataAndErrors(
                    std::vector<double>(TomoDataVec.begin(), TomoDataVec.end()),
                    TomoErrVec);
                TomoData.WriteNetCDF(filename + ".inv_tt.nc");
              }
            TomoModel.WriteVTK(filename + ".tomo.inv.vtk");
            TomoModel.WriteNetCDF(filename + ".tomo.inv.nc");
          }
      }

  }
