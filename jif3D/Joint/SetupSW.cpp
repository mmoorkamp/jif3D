//
// Created by wangchao on 2022/4/14.
//

#include "../SurfaceWaves/SurfaceWaveData.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../ModelTransforms/DensPVelTransform.h"
#include "../Inversion/ModelTransforms.h"

#include "SetupSW.h"
#include <iostream>
#include <vector>
#include <boost/make_shared.hpp>

namespace jif3D
  {

    const std::string Name = "S-velocity";

    SetupSW::SetupSW() :
        GeneralDataSetup(Name), relerr(2.0e-2), minerr(0.0), dswlambda(0.0)
      {

      }

    SetupSW::~SetupSW()
      {

      }

    po::options_description SetupSW::SetupOptions()
      {
        po::options_description desc("Surface wave tomography options");
        desc.add_options()("swmodel", po::value(&modelfilename),
            "The name of the starting model for the surface wave tomography")("dswdata",
            po::value(&datafilename),
            "The name of the data for the surface wave tomography")("dswlambda",
            po::value(&dswlambda), "The weight for the surface wave tomography data")(
            "minvs", po::value(&minvs)->default_value(1.5),
            "The minimum shear wave velocivty")("maxvs",
            po::value(&maxvs)->default_value(6.5), "The maximum shear wave velocivty")(
            "relerr", po::value(&relerr)->default_value(2.0e-2),
            "The relative error for the phase delay data")("minerr",
            po::value(&minerr)->default_value(0.0),
            "The picker error for the phase delay data");
        return desc;
      }

    void SetupSW::IterationOutput(const std::string &filename,
        const jif3D::rvec &ModelVector)
      {

      }

    void SetupSW::FinalOutput(const jif3D::rvec &FinalModelVector)
      {

      }

    bool SetupSW::SetupObjective(const boost::program_options::variables_map &vm,
        jif3D::JointObjective &Objective, jif3D::ThreeDModelBase &InversionMesh,
        jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
        std::vector<std::string> &SegmentNames, std::vector<parametertype> &SegmentTypes,
        boost::filesystem::path TempDir)
      {

        jif3D::SurfaceWaveData DSurfaceWaveData;
        const size_t ngrid = InversionMesh.GetData().num_elements();

        dswlambda = 0.0;
        if (!vm.count("dswlambda"))
          {
            std::cout << "Surface wave tomography Lambda: ";
            std::cin >> dswlambda;
          }

        if (dswlambda > 0.0)
          {
            if (!vm.count("modelfilename"))
              {
                modelfilename = jif3D::AskFilename(
                    "Surface wave tomography inversion model Filename: ");
              }
            DSurfaceWaveModel.ReadNetCDF(modelfilename);
            DSurfaceWaveModel.WriteVTK(modelfilename + ".vtk");

            StartingParameters.resize(ngrid);
            std::copy(DSurfaceWaveModel.GetData().origin(),
                DSurfaceWaveModel.GetData().origin() + ngrid, StartingParameters.begin());

            if (!vm.count("swdata"))
              {
                datafilename = jif3D::AskFilename(
                    "Surface wave tomography Data Filename: ");
              }
            DSurfaceWaveData.ReadNetCDF(datafilename);

            std::vector<double> DSurfaceWaveError = ConstructError(
                DSurfaceWaveData.GetData(), DSurfaceWaveData.GetErrors(), relerr, minerr);
            DSurfaceWaveData.SetDataAndErrors(DSurfaceWaveData.GetData(),
                DSurfaceWaveError);

            jif3D::SurfaceWaveCalculator Calculator;
            DSurfaceWaveObjective = boost::make_shared<
                jif3D::ThreeDModelObjective<jif3D::SurfaceWaveCalculator>>(Calculator);

            size_t start = startindices.back();
            size_t end = start + ngrid;
            boost::shared_ptr<jif3D::GeneralModelTransform> SWTrans = boost::make_shared<
                jif3D::TanhTransform>(minvs, maxvs);
            Transform = boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid, start,
                end, SWTrans);
            startindices.push_back(end);
            SegmentNames.push_back(Name);
            SegmentTypes.push_back(GeneralDataSetup::gridparameter);

            DSurfaceWaveObjective->SetObservedData(DSurfaceWaveData);
            DSurfaceWaveObjective->SetCoarseModelGeometry(DSurfaceWaveModel);
            DSurfaceWaveObjective->SetDataError(DSurfaceWaveError);

            Objective.AddObjective(DSurfaceWaveObjective, Transform, dswlambda, "DSWTomo",
                JointObjective::datafit);
            std::cout << "DSurfaceWave tomo ndata: " << DSurfaceWaveData.GetData().size()
                << std::endl;
            std::cout << "DSurfaceWave tomo lambda: " << dswlambda << std::endl;

          }
        return (dswlambda > 0.0);
      }

  }
