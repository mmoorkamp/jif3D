//============================================================================
// Name        : SetupRegularization.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "SetupRegularization.h"
#include "../Global/convert.h"
#include "../Regularization/GradientRegularization.h"
#include "../Regularization/HOGradientRegularization.h"
#include "../Regularization/CurvatureRegularization.h"
#include "../Regularization/MinDiffRegularization.h"
#include "../Regularization/MinimumSupport.h"
#include "../Regularization/EntropyRegularization.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../ModelBase/ReadAnyModel.h"
#include <boost/make_shared.hpp>
#include <boost/filesystem.hpp>
#include <algorithm>

namespace jif3D
  {

    void SetTearModel(const po::variables_map &vm, const std::string &OptionName,
        const ThreeDModelBase &StartModel, jif3D::ThreeDModelBase &TearModel)
      {
        if (vm.count(OptionName))
          {
            //we use a seismic model file as a container for the tear information
            //but this model does not have to obey the gridding rules
            std::string Filename(vm[OptionName].as<std::string>());
            if (!boost::filesystem::exists(Filename))
              {
                throw jif3D::FatalException(
                    "Tear model has been specified but file: " + Filename
                        + " does not exist ", __FILE__, __LINE__);

              }
            TearModel = *ReadAnyModel(Filename).get();
            if (TearModel.GetModelShape()[0] != StartModel.GetModelShape()[0])
              {
                throw jif3D::FatalException(
                    Filename
                        + ": X-dimensions of TearModel do not match x-dimensions of starting model ",
                    __FILE__, __LINE__);
              }
            if (TearModel.GetModelShape()[1] != StartModel.GetModelShape()[1])
              {
                throw jif3D::FatalException(
                    Filename
                        + ": Y-dimensions of TearModel do not match y-dimensions of starting model ",
                    __FILE__, __LINE__);
              }
            if (TearModel.GetModelShape()[2] != StartModel.GetModelShape()[2])
              {
                throw jif3D::FatalException(
                    Filename
                        + ": Z-dimensions of TearModel do not match z-dimensions of starting model ",
                    __FILE__, __LINE__);
              }
          }
        else
          {
            TearModel.SetMeshSize(StartModel.GetXCellSizes().size(),
                StartModel.GetYCellSizes().size(), StartModel.GetZCellSizes().size());
            std::fill_n(TearModel.SetData().origin(), TearModel.GetData().num_elements(),
                1.0);
          }

      }

    SetupRegularization::SetupRegularization() :
        beta(0.0), substart(false), xweight(1.0), yweight(1.0), zweight(1.0), minsuppb(
            1.0), nbins(0)
      {
      }

    SetupRegularization::~SetupRegularization()
      {
      }

    po::options_description SetupRegularization::SetupOptions()
      {
        //add all possible options for regularization to an
        //options_description object that can be used for
        //parsing
        po::options_description desc("Regularization options");

        desc.add_options()("xreg", po::value(&xweight)->default_value(1.0),
            "The weight for the regularization in x-direction")("yreg",
            po::value(&yweight)->default_value(1.0),
            "The weight for the regularization in y-direction")("zreg",
            po::value(&zweight)->default_value(1.0),
            "The weight for the regularization in z-direction")("curvreg",
            "Use model curvature for regularization. If not set use gradient.")(
            "hogradient", "Use 4th order approximation for gradient")("mindiff",
            "Minize the model vector (or difference to starting model if substart is set")(
            "tearmodx", po::value<std::string>(),
            "Filename for a model containing information about tear zones in x-direction.")(
            "tearmody", po::value<std::string>(),
            "Filename for a model containing information about tear zones in y-direction.")(
            "tearmodz", po::value<std::string>(),
            "Filename for a model containing information about tear zones in z-direction.")(
            "beta", po::value(&beta)->default_value(0.0),
            "The weight for the model parameter minimization in the regularization")(
            "substart", po::value(&substart)->default_value(false),
            "Substract the starting model when calculating the roughness")("minsupp",
            po::value(&minsuppb), "Use minimum support regularization (EXPERIMENTAL)")(
            "mingradsupp", po::value(&minsuppb),
            "Use minimum gradient support regularization (EXPERIMENTAL)")("entropy",
            po::value(&nbins)->default_value(0),
            "Use entropy regularization (EXPERIMENTAL)")("considersize",
            "Consider the size of the cells when calculating regularization");

        return desc;
      }

    boost::shared_ptr<jif3D::RegularizationFunction> SetupRegularization::SetupObjective(
        const po::variables_map &vm, const ThreeDModelBase &StartModel)
      {
        //if we only want to use a minimum model, we do not
        //have to worry about tear models etc. so we create the object
        //and return
        if (vm.count("mindiff"))
          {
            return boost::shared_ptr<jif3D::MinDiffRegularization>(
                new jif3D::MinDiffRegularization(StartModel));
          }
        if (vm.count("minsupp"))
          {
            return boost::shared_ptr<jif3D::MinimumSupport>(
                new jif3D::MinimumSupport(
                    boost::shared_ptr<jif3D::MinDiffRegularization>(
                        new jif3D::MinDiffRegularization(StartModel)), minsuppb));
          }
        if (vm.count("entropy") && nbins > 0.0)
          {
            return boost::make_shared<jif3D::EntropyRegularization>(-2.0, 2.0, nbins);
          }

        //setup possible tearing for the regularization for the three directions
        jif3D::ThreeDSeismicModel TearModX, TearModY, TearModZ;
        SetTearModel(vm, "tearmodx", StartModel, TearModX);
        SetTearModel(vm, "tearmody", StartModel, TearModY);
        SetTearModel(vm, "tearmodz", StartModel, TearModZ);
        TearModX.WriteVTK("tearmodx.vtk");
        TearModY.WriteVTK("tearmody.vtk");
        TearModZ.WriteVTK("tearmodz.vtk");

        if (vm.count("mingradsupp"))
          {
            return boost::shared_ptr<jif3D::MinimumSupport>(
                new jif3D::MinimumSupport(
                    boost::shared_ptr<jif3D::GradientRegularization>(
                        new jif3D::GradientRegularization(StartModel, TearModX, TearModY,
                            TearModZ, beta)), minsuppb));
          }

        //decide whether we want to use gradient base regularization
        //or curvature based regularization, the default is gradient
        boost::shared_ptr<jif3D::MatOpRegularization> Regularization;
        if (vm.count("curvreg"))
          {
            Regularization = boost::shared_ptr<jif3D::MatOpRegularization>(
                new jif3D::CurvatureRegularization(StartModel, TearModX, TearModY,
                    TearModZ, beta));

          }
        else
          {
            if (vm.count("hogradient"))
              {
                Regularization = boost::shared_ptr<jif3D::MatOpRegularization>(
                    new jif3D::HOGradientRegularization(StartModel, TearModX, TearModY,
                        TearModZ, beta));
              }
            else
              {
                Regularization = boost::shared_ptr<jif3D::MatOpRegularization>(
                    new jif3D::GradientRegularization(StartModel, TearModX, TearModY,
                        TearModZ, beta, vm.count("considersize")));
              }
          }

        //we can directly use the values for the weights without checking
        //the options because we set the default value to 1
        Regularization->SetXWeight(xweight);
        Regularization->SetYWeight(yweight);
        Regularization->SetZWeight(zweight);

        return Regularization;
      }

    boost::shared_ptr<jif3D::RegularizationFunction> SetupRegularization::SetupObjective(
        const po::variables_map &vm, const ThreeDModelBase &StartModel,
        const jif3D::rvec &CovModVec, jif3D::ThreeDModelBase &TearModelX,
        jif3D::ThreeDModelBase &TearModelY, jif3D::ThreeDModelBase &TearModelZ)
      {
        SetTearModel(vm, "tearmodx", StartModel, TearModelX);
        SetTearModel(vm, "tearmody", StartModel, TearModelY);
        SetTearModel(vm, "tearmodz", StartModel, TearModelZ);
        return SetupObjective(vm, StartModel);

      }

  }
