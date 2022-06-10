/*
 * SetupMagnetization.cpp
 *
 *  Created on: Jun 2, 2022
 *      Author: max
 */

#include "SetupMagnetization.h"
#include "../Inversion/ModelTransforms.h"
#include "../Magnetics/ReadWriteMagneticData.h"
#include "../Magnetics/OMPMagnetizationImp.h"
#include "../Magnetics/ThreeComponentMagneticData.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/ReadAnyModel.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"

#include <boost/make_shared.hpp>

namespace jif3D
  {

    namespace po = boost::program_options;

    const std::string Name = "Magnetization";

    boost::program_options::options_description SetupMagnetization::SetupOptions()
      {
        po::options_description desc("Magnetization options");

        desc.add_options()("totalmagrelerr", po::value(&relerr)->default_value(0.02),
            "The relative error for the Magnetics data")("totalmagminerr",
            po::value(&minerr)->default_value(0.0),
            "The minimum absolute error for the Magnetics data")("magnetization_covmod",
            po::value<std::string>(),
            "A file containing the model covariance for susceptibilities")("magdepth",
            "Counteract the decay in sensitivities of magnetic data with depth")("minmag",
            po::value(&minmag)->default_value(-1.0))("maxmag",
            po::value(&maxmag)->default_value(1.0));

        return desc;
      }

    bool SetupMagnetization::SetupObjective(
        const boost::program_options::variables_map &vm, jif3D::JointObjective &Objective,
        jif3D::ThreeDModelBase &InversionMesh, jif3D::rvec &CovModVec,
        std::vector<size_t> &startindices, std::vector<std::string> &SegmentNames,
        std::vector<parametertype> &SegmentTypes, boost::filesystem::path TempDir)
      {

        const size_t ngrid = InversionMesh.GetData().num_elements();
        std::cout << "Magnetization Lambda: ";
        std::cin >> maglambda;
        if (maglambda > 0.0)
          {
            jif3D::ThreeComponentMagneticData MagData;

            //if the weight is really small we use this as a "secret setting"
            //for constrained inversion, so we only read in data when the
            //weight is larger than 1e-32, otherwise we create a single point
            //dummy dataset, the joint objective class is set that
            //objective functions with such a small weight are not evaluated anyway
            if (maglambda > JointObjective::MinWeight)
              {
                //read in data file
                std::string magdatafilename = jif3D::AskFilename(
                    "Three component magnetic Data Filename: ");
                MagData.ReadNetCDF(magdatafilename);
                MagData.WriteVTK(magdatafilename);
              }

            std::string magmodelfilename = jif3D::AskFilename(
                "Magnetization Model Filename: ");
            Model.ReadNetCDF(magmodelfilename);
            jif3D::rvec ModParm(Model.GetModelParameters());
            if (ModParm.size() > 0)
              {
                if (ngrid != Model.GetMagnetization_X().num_elements())
                  {
                    std::string err = "Magnetization model does not match grid size ! "
                        + std::to_string(Model.GetMagnetization_X().num_elements()) + " "
                        + std::to_string(ngrid);
                    throw jif3D::FatalException(err, __FILE__, __LINE__);
                  }
                StartingParameters.resize(3 * ngrid);
                std::copy(ModParm.begin(), ModParm.end(), StartingParameters.begin());
              }

            boost::shared_ptr<jif3D::GeneralModelTransform> MagTrans;
            MagTrans = boost::make_shared<jif3D::TanhTransform>(minmag, maxmag);

            size_t start = startindices.back();
            size_t end = start + ngrid;
            startindices.push_back(end);
            SegmentNames.push_back(Name + "_X");
            SegmentTypes.push_back(GeneralDataSetup::gridparameter);

            end += ngrid;
            startindices.push_back(end);
            SegmentNames.push_back(Name + "_Y");
            SegmentTypes.push_back(GeneralDataSetup::gridparameter);

            end += ngrid;
            startindices.push_back(end);
            SegmentNames.push_back(Name + "_Z");
            SegmentTypes.push_back(GeneralDataSetup::gridparameter);

            Transform = boost::make_shared<jif3D::MultiSectionTransform>(2 * ngrid, start,
                end, MagTrans);

            CovModVec.resize(3 * ngrid);
            std::fill(CovModVec.begin(), CovModVec.end(), 1.0);
            if (vm.count("mag_covmod"))
              {
                boost::shared_ptr<jif3D::ThreeDModelBase> CovModel = jif3D::ReadAnyModel(
                    vm["mag_covmod"].as<std::string>());

                const size_t ncovmod = CovModel->GetData().num_elements();
                if (ncovmod != ngrid)
                  {
                    std::cerr
                        << " Magnetization covariance does not have the same number of parameters "
                        << ncovmod << " as inversion model " << ngrid << std::endl;
                    return 100;
                  }
                else
                  {
                    std::cout << " Setting covariance for susceptibility from file. "
                        << std::endl;
                  }
                std::copy(CovModel->GetData().origin(),
                    CovModel->GetData().origin() + ncovmod, CovModVec.begin());
                std::copy(CovModel->GetData().origin(),
                    CovModel->GetData().origin() + ncovmod, CovModVec.begin() + ngrid);
                std::copy(CovModel->GetData().origin(),
                    CovModel->GetData().origin() + ncovmod,
                    CovModVec.begin() + 2 * ngrid);

              }
            if (maglambda > JointObjective::MinWeight)
              {
                //we want to set the path for temporary file storage
                //the factory function cannot perform this, so we
                //have to assemble the calculator object ourselves
                boost::shared_ptr<
                    jif3D::ThreeDGravMagImplementation<jif3D::ThreeComponentMagneticData> > Implementation =
                    boost::make_shared<jif3D::OMPMagnetizationImp>();
#ifdef MAGDISK
                  Calculator = boost::make_shared<MagCalculatorType>(Implementation, TempDir);
#else
                Calculator = boost::make_shared<MagCalculatorType>(Implementation);
#endif

                MagObjective = boost::make_shared<
                    jif3D::ThreeDModelObjective<MagCalculatorType> >(*Calculator);
                MagObjective->SetObservedData(MagData);
                MagObjective->SetCoarseModelGeometry(Model);
                std::vector<double> Error(
                    jif3D::ConstructError(MagData.GetData(), MagData.GetErrors(), relerr,
                        minerr));
                MagObjective->SetDataError(Error);

                Objective.AddObjective(MagObjective, Transform, maglambda, "Magnetics",
                    JointObjective::datafit);
                std::cout << " Magnetics ndata: " << MagData.GetData().size()
                    << std::endl;
                std::cout << " Magnetics lambda: " << maglambda << std::endl;
              }

          }

        //indicate whether we added a Magnetics objective function
        //this way the caller can do additional consistency checks
        return (maglambda > JointObjective::MinWeight);

      }

    void SetupMagnetization::IterationOutput(const std::string &filename,
        const jif3D::rvec &ModelVector)
      {
        if (maglambda > JointObjective::MinWeight)
          {
            jif3D::rvec MagInvModel = Transform->GeneralizedToPhysical(ModelVector);
            Model.SetModelParameter(MagInvModel);
            Model.WriteVTK(filename + ".mag.inv.vtk");
            Model.WriteNetCDF(filename + ".mag.inv.nc");
          }
      }

    void SetupMagnetization::FinalOutput(const jif3D::rvec &FinalModelVector)
      {
        if (maglambda > JointObjective::MinWeight)
          {
            std::string modelfilename = "result";
            jif3D::rvec MagInvModel = Transform->GeneralizedToPhysical(FinalModelVector);
            std::cout << "Writing final magnetization models " << std::endl;
            Model.SetModelParameter(MagInvModel);

            auto ObservedData = MagObjective->GetObservedData();
            jif3D::rvec MagInvData = MagObjective->GetSyntheticData();
            if (MagInvData.size() > 0)
              {
                jif3D::SaveMagneticComponentMeasurements(modelfilename + ".inv_mag.nc",
                    std::vector<double>(MagInvData.begin(), MagInvData.end()),
                    ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
                    ObservedData.GetMeasPosZ(), MagObjective->GetDataError());
                jif3D::Write3DVectorDataToVTK(modelfilename + ".inv_mag.vtk", "T",
                    std::vector<double>(MagInvData.begin(), MagInvData.end()),
                    ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
                    ObservedData.GetMeasPosZ());
                jif3D::rvec MagDiff(MagObjective->GetIndividualMisfit());
                jif3D::SaveMagneticComponentMeasurements(modelfilename + ".diff_mag.nc",
                    std::vector<double>(MagDiff.begin(), MagDiff.end()),
                    ObservedData.GetMeasPosX(), ObservedData.GetMeasPosY(),
                    ObservedData.GetMeasPosZ(), MagObjective->GetDataError());
              }
            Model.WriteVTK(modelfilename + ".mag.inv.vtk");
            Model.WriteNetCDF(modelfilename + ".mag.inv.nc");
          }
      }

    SetupMagnetization::SetupMagnetization() :
        GeneralDataSetup(Name)
      {
        // TODO Auto-generated constructor stub

      }

    SetupMagnetization::~SetupMagnetization()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
