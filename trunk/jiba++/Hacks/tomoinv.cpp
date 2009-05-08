//============================================================================
// Name        : gravinv.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================


#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <boost/bind.hpp>
#include "../Global/convert.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/VectorTransform.h"
#include "../ModelBase/VTKTools.h"
#include "../ModelBase/NetCDFTools.h"
#include "../Inversion/LimitedMemoryQuasiNewton.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/MinDiffRegularization.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../Tomo/TomographyObjective.h"
#include "../Tomo/TomographyCalculator.h"

namespace ublas = boost::numeric::ublas;

class LogTransform: public jiba::GeneralModelTransform
  {
public:
  virtual jiba::rvec Transform(const jiba::rvec &FullModel)
    {
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = pow(10, FullModel(i));
      return Output;
    }

  virtual jiba::rvec Derivative(const jiba::rvec &FullModel,
      const jiba::rvec &Derivative)
    {
      const double factor = std::log(10.0);
      jiba::rvec Output(FullModel.size());
      for (size_t i = 0; i < FullModel.size(); ++i)
        Output( i) = factor * pow(10, FullModel(i)) * Derivative(i);

      return Output;
    }
  LogTransform()
    {
    }
  virtual ~LogTransform()
    {
    }
  };

int main(int argc, char *argv[])
  {
    //these objects hold information about the measurements and their geometry
    jiba::rvec Data;

    //first we read in the starting model and the measured data
    std::string modelfilename, datafilename;
    std::cout << "Starting model Filename: ";
    std::cin >> modelfilename;
    //we read in the starting modelfile
    jiba::ThreeDSeismicModel Model;
    Model.ReadNetCDF(modelfilename);
    //get the name of the file containing the data and read it in
    std::cout << "Data Filename: ";
    std::cin >> datafilename;

    //read in data
    jiba::ReadTraveltimes(datafilename, Data, Model);
    //if we don't have data inversion doesn't make sense;
    if (Data.empty())
      {
        std::cerr << "No measurements defined" << std::endl;
        exit(100);
      }

    //we define a few constants that are used throughout the inversion

    const size_t ndata = Data.size();
    const double errorlevel = 0.02;
    const double maxdata = *std::max_element(Data.begin(), Data.end(),
        jiba::absLess<double, double>());
    //create objects for the misfit and a very basic error estimate
        jiba::rvec DataError(ndata);
    for (size_t i = 0; i < ndata; ++i)
      {
        DataError( i) = std::max(std::abs(Data(i) * errorlevel), 1e-2 * maxdata
            * errorlevel);
      }

    jiba::rvec InvModel(Model.GetSlownesses().num_elements());
    std::copy(Model.GetSlownesses().origin(), Model.GetSlownesses().origin()
        + Model.GetSlownesses().num_elements(), InvModel.begin());

    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel( i) = std::log10(InvModel(i));
    boost::shared_ptr<jiba::TomographyObjective> TomoObjective(
        new jiba::TomographyObjective());
    TomoObjective->SetObservedData(Data);
    TomoObjective->SetModelGeometry(Model);
    TomoObjective->SetDataCovar(DataError);

    boost::shared_ptr<jiba::JointObjective> Objective(
        new jiba::JointObjective());
    boost::shared_ptr<jiba::MinDiffRegularization> Regularization(
        new jiba::MinDiffRegularization());
    jiba::rvec RefModel(InvModel);
    Regularization->SetReferenceModel(RefModel);

    double lambda = 1.0;
    std::cout << "Lambda: ";
    std::cin >> lambda;
    Objective->AddObjective(TomoObjective, boost::shared_ptr<
        jiba::GeneralModelTransform>(new LogTransform()));
    Objective->AddObjective(Regularization, boost::shared_ptr<
        jiba::GeneralModelTransform>(new LogTransform()), lambda);

    std::cout << "Performing inversion." << std::endl;

    jiba::LimitedMemoryQuasiNewton LBFGS(Objective, 5);

    for (size_t i = 0; i < 5; ++i)
      {
        LBFGS.MakeStep(InvModel);
        std::cout << std::endl;

      }

    for (size_t i = 0; i < InvModel.size(); ++i)
      InvModel( i) = pow(10,InvModel(i));
    std::copy(InvModel.begin(), InvModel.begin()
        + Model.GetSlownesses().num_elements(), Model.SetSlownesses().origin());

    //calculate the predicted data
    std::cout << "Calculating response of inversion model." << std::endl;
    jiba::rvec InvData(jiba::TomographyCalculator().Calculate(Model));
    jiba::SaveTraveltimes(modelfilename + ".inv_tt.nc", InvData, Model);
    //and write out the data and model
    //here we have to distinguish again between scalar and ftg data
    std::cout << "Writing out inversion results." << std::endl;

    Model.WriteVTK(modelfilename + ".inv.vtk");
    Model.WriteNetCDF(modelfilename + ".inv.nc");
    std::cout << std::endl;
  }
