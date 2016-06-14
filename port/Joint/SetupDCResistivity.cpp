//============================================================================
// Name        : SetupDCResistivity.cpp
// Author      : 4 Jun 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#include "SetupDCResistivity.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include "../DCResistivity/ReadWriteDCResistivityData.h"

namespace jif3D
  {

    po::options_description SetupDCResistivity::SetupOptions()
      {
        po::options_description desc("DC resistivity options");

        desc.add_options()("dcrelerr", po::value(&relerr)->default_value(0.02),
            "The relative error for the DC resistivity data")("dcminerr",
            po::value(&minerr)->default_value(1e-3))("dcfine",
                    po::value(&DCFineModelName),
                    "The name for the model with the DC forward geometry");

        return desc;
      }

    bool SetupDCResistivity::SetupObjective(const po::variables_map &vm,
        jif3D::JointObjective &Objective,
        boost::shared_ptr<jif3D::GeneralModelTransform> &Transform, double xorigin,
        double yorigin)
      {

        jif3D::rvec DCData, DCError;
        double dclambda = 1.0;

        std::cout << "DC Resistivity Lambda: ";
        std::cin >> dclambda;

        //if the weight is different from zero
        //we have to read in scalar Magnetics data
        if (dclambda > 0.0)
          {
            std::string dcdatafilename = jif3D::AskFilename(
                "DC Resistivity Data Filename: ");

            jif3D::ReadApparentResistivity(dcdatafilename, DCData, DCError, Model);

            std::string dcmodelfilename = jif3D::AskFilename(
                "DC Resistivity Model Filename: ");
            Model.ReadNetCDF(dcmodelfilename);

            Model.SetOrigin(xorigin, yorigin, 0.0);
            Calculator = boost::shared_ptr<DCResistivityCalculator>(new DCResistivityCalculator);
            DCObjective = boost::shared_ptr<
                jif3D::ThreeDModelObjective<DCResistivityCalculator> >(
                new jif3D::ThreeDModelObjective<DCResistivityCalculator>(*Calculator));
            DCObjective->SetObservedData(DCData);
            DCObjective->SetCoarseModelGeometry(Model);
            jif3D::rvec Error(jif3D::ConstructError(DCData, DCError, relerr, minerr));
            DCObjective->SetDataError(Error);

            if (vm.count("dcfine"))
              {
                jif3D::ThreeDDCResistivityModel DCFineGeometry;
                DCFineGeometry.ReadNetCDF(DCFineModelName);
                //copy measurement configuration to refined model
                DCFineGeometry.CopyMeasurementConfigurations(Model);
                DCObjective->SetFineModelGeometry(DCFineGeometry);
              }
            Objective.AddObjective(DCObjective, Transform, dclambda, "DC Resistivity",
                JointObjective::datafit);
            std::cout << " Resistivity ndata: " << DCData.size() << std::endl;
            std::cout << " Resistivity lambda: " << dclambda << std::endl;
          }

        //indicate whether we added a Magnetics objective function
        //this way the caller can do additional consistency checks
        return (dclambda > 0.0);
      }

    SetupDCResistivity::SetupDCResistivity()
      {
        // TODO Auto-generated constructor stub

      }

    SetupDCResistivity::~SetupDCResistivity()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
