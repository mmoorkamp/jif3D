//============================================================================
// Name        : gravgrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file gravgrid.cpp
 * A simple program to calculate evenly spaced gravity data from a model stored in a netcdf file.
 * The spacing and region for the measurements is specified interactively. The model is taken from
 * a netcdf file.
 */

#include <iostream>
#include <string>
#include "ThreeDGravityModel.h"
#include "ScalarGravityData.h"
#include "TensorGravityData.h"
#include "ReadWriteGravityData.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "ThreeDGravityFactory.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include <boost/cast.hpp>
#include <boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
  {
    bool wantcuda = false;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("cuda",
        po::value(&wantcuda)->default_value(false), "Use cuda for forward calculations.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }

    jif3D::ThreeDGravityModel GravModel;
    jif3D::ScalarGravityData ScalData;
    jif3D::TensorGravityData TensData;
    double minx, miny, maxx, maxy, deltax, deltay, z;
    //ask for the measurement grid specifications
    //first x-direction
    cout << "Minimum x-position: ";
    cin >> minx;
    cout << "Maximum x-position: ";
    cin >> maxx;
    cout << "Delta x: ";
    cin >> deltax;
    //then y-direction
    cout << "Minimum y-position: ";
    cin >> miny;
    cout << "Maximum y-position: ";
    cin >> maxy;
    cout << "Delta y: ";
    cin >> deltay;
    //all measurements have to be at the same height.
    cout << "Z-level: ";
    cin >> z;
    //setup the measurements in the forward modelling code
    const size_t nmeasx = boost::numeric_cast<size_t>((maxx - minx) / deltax);
    const size_t nmeasy = boost::numeric_cast<size_t>((maxy - miny) / deltay);
    for (size_t i = 0; i <= nmeasx; ++i)
      {
        for (size_t j = 0; j <= nmeasy; ++j)
          {
            ScalData.AddMeasurementPoint(minx + i * deltax, miny + j * deltay, z);
            TensData.AddMeasurementPoint(minx + i * deltax, miny + j * deltay, z);
          }
      }
    //ask for the name of the netcdf file containing the model
    std::string ModelFilename = jif3D::AskFilename("Model Filename: ");

    //determine the extension to find out the type
    std::string extension = jif3D::GetFileExtension(ModelFilename);
    //read in the file
    if (extension == ".nc")
      {
        GravModel.ReadNetCDF(ModelFilename);
      }
    else
      {
        GravModel.ReadIgmas(ModelFilename);
        GravModel.WriteNetCDF(ModelFilename + ".nc");
      }
    typedef typename jif3D::MinMemGravMagCalculator<jif3D::ScalarGravityData> ScalCalculatorType;
    typedef typename jif3D::MinMemGravMagCalculator<jif3D::TensorGravityData> TensCalculatorType;
    //save the measurements and some plots
    boost::shared_ptr<TensCalculatorType> TensorCalculator(
        jif3D::CreateGravityCalculator<TensCalculatorType>::MakeTensor(wantcuda));
    boost::shared_ptr<ScalCalculatorType> ScalarCalculator(
        jif3D::CreateGravityCalculator<ScalCalculatorType>::MakeScalar(wantcuda));
    jif3D::rvec ScalarResults(ScalarCalculator->Calculate(GravModel, ScalData));
    jif3D::rvec TensorResults(TensorCalculator->Calculate(GravModel, TensData));

    jif3D::rvec ScalErr(ScalarResults.size());
    jif3D::rvec TensErr(TensorResults.size());
    double scalrelnoise = 0.0;
    double ftgrelnoise = 0.0;
    double scalabsnoise = 0.0;
    double ftgabsnoise = 0.0;
    cout << "Relative noise level (scalar data): ";
    cin >> scalrelnoise;
    cout << "Relative noise level (ftg data): ";
    cin >> ftgrelnoise;
    cout << "Absolute noise level (scalar data): ";
    cin >> scalabsnoise;
    cout << "Absolute noise level (ftg data): ";
    cin >> ftgabsnoise;
    jif3D::AddNoise(ScalarResults, scalrelnoise, scalabsnoise);
    std::transform(ScalarResults.begin(), ScalarResults.end(), ScalErr.begin(),
        [scalrelnoise, scalabsnoise] (const double data)
          { return std::max(scalabsnoise,data*scalrelnoise);});
    jif3D::AddNoise(TensorResults, ftgrelnoise, ftgabsnoise);
    std::transform(TensorResults.begin(), TensorResults.end(), TensErr.begin(),
        [ftgrelnoise, ftgabsnoise] (const double data)
          { return std::max(ftgabsnoise,data*ftgrelnoise);});
    ScalData.SetDataAndErrors(
        std::vector<double>(ScalarResults.begin(), ScalarResults.end()),
        std::vector<double>(ScalErr.begin(), ScalErr.end()));
    TensData.SetDataAndErrors(
        std::vector<double>(TensorResults.begin(), TensorResults.end()),
        std::vector<double>(TensErr.begin(), TensErr.end()));

    ScalData.WriteNetCDF(ModelFilename + ".sgd.nc");
    TensData.WriteNetCDF(ModelFilename + ".ftg.nc");

    //write the model in .vtk format, at the moment the best plotting option
    GravModel.WriteVTK(ModelFilename + ".vtk");

    jif3D::Write3DDataToVTK(ModelFilename + ".data.vtk", "grav_accel",
        std::vector<double>(ScalarResults.begin(), ScalarResults.end()),
        ScalData.GetMeasPosX(), ScalData.GetMeasPosY(), ScalData.GetMeasPosZ());
    jif3D::Write3DTensorDataToVTK(ModelFilename + ".ftgdata.vtk", "U",
        std::vector<double>(TensorResults.begin(), TensorResults.end()),
        TensData.GetMeasPosX(), TensData.GetMeasPosY(), TensData.GetMeasPosZ());

  }
