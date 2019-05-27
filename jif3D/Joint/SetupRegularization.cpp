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
#include "../Tomo/ThreeDSeismicModel.h"
#include "../ModelBase/ReadAnyModel.h"

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
            TearModel = *ReadAnyModel(Filename).get();
            if (TearModel.GetModelShape()[0] != StartModel.GetModelShape()[0])
              {
                throw jif3D::FatalException(
                    "X-dimensions of TearModel do not match x-dimensions of starting model ",
                    __FILE__, __LINE__);
              }
            if (TearModel.GetModelShape()[1] != StartModel.GetModelShape()[1])
              {
                throw jif3D::FatalException(
                    "Y-dimensions of TearModel do not match x-dimensions of starting model ",
                    __FILE__, __LINE__);
              }
            if (TearModel.GetModelShape()[2] != StartModel.GetModelShape()[2])
              {
                throw jif3D::FatalException(
                    "Z-dimensions of TearModel do not match x-dimensions of starting model ",
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
            1.0)
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
            "Use minimum gradient support regularization (EXPERIMENTAL)")("considersize",
            "Consider the size of the cells when calculating regularization");

        return desc;
      }

    boost::shared_ptr<jif3D::RegularizationFunction> SetupRegularization::SetupObjective(
        const po::variables_map &vm, const ThreeDModelBase &StartModel,
        const jif3D::rvec &CovModVec)
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

        //setup possible tearing for the regularization for the three directions
        jif3D::ThreeDSeismicModel TearModX, TearModY, TearModZ;
        SetTearModel(vm, "tearmodx", StartModel, TearModX);
        SetTearModel(vm, "tearmody", StartModel, TearModY);
        SetTearModel(vm, "tearmodz", StartModel, TearModZ);

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
        //We either pass an empty covariance vector then the regularization class
        //takes care of setting the covariance to 1
        //or we set a covariance vector that has one value for each model cell
        if (!CovModVec.empty())
          {
            //as we treat each direction separately, the covariance
            //vector has to have a length 3 times the length of the model vector
            //here we copy the covariances to the appropriate position
            const size_t ngrid = StartModel.GetNModelElements();
            if (CovModVec.size() != ngrid)
              {
                throw jif3D::FatalException(
                    "Size of model covariance: " + jif3D::stringify(CovModVec.size())
                        + " does not match model size: " + jif3D::stringify(ngrid),
                    __FILE__, __LINE__);
              }
            jif3D::rvec Cov(ngrid * 3);
            ublas::subrange(Cov, 0, ngrid) = CovModVec;
            ublas::subrange(Cov, ngrid, 2 * ngrid) = CovModVec;
            ublas::subrange(Cov, 2 * ngrid, 3 * ngrid) = CovModVec;
            Regularization->SetDataError(std::vector<double>(Cov.begin(),Cov.end()));
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
        return SetupObjective(vm, StartModel, CovModVec);

      }

  }
