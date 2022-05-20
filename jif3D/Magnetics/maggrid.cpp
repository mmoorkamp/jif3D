//============================================================================
// Name        : maggrid.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2013, MM
//============================================================================

/*! \file maggrid.cpp
 * A simple program to calculate evenly spaced magnetic data from a model stored in a netcdf file.
 * The spacing and region for the measurements is specified interactively. The model is taken from
 * a netcdf file.
 */

#include <iostream>
#include <string>
#include "ReadWriteMagneticData.h"
#include "MagneticTransforms.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/ThreeDGravMagImplementation.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"
#include "../Global/Noise.h"
#include <boost/cast.hpp>
#include <boost/program_options.hpp>

#include "OMPMagneticSusceptibilityImp.h"
#include "ThreeDSusceptibilityModel.h"
#include "TotalFieldMagneticData.h"

using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
  {
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }

    jif3D::ThreeDSusceptibilityModel MagModel;
    jif3D::TotalFieldMagneticData Data;
    double minx, miny, maxx, maxy, deltax, deltay, z;
    double inclination, declination, fieldstrength;
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
    cout << "Inclination: ";
    cin >> inclination;
    cout << "Declination: ";
    cin >> declination;
    cout << "Field Strength: ";
    cin >> fieldstrength;
    //setup the measurements in the forward modelling code
    const size_t nmeasx = boost::numeric_cast<size_t>((maxx - minx) / deltax);
    const size_t nmeasy = boost::numeric_cast<size_t>((maxy - miny) / deltay);
    for (size_t i = 0; i <= nmeasx; ++i)
      {
        for (size_t j = 0; j <= nmeasy; ++j)
          {
            Data.AddMeasurementPoint(minx + i * deltax, miny + j * deltay, z);
          }
      }
    //ask for the name of the netcdf file containing the model
    std::string ModelFilename = jif3D::AskFilename("Model Filename: ");

    MagModel.ReadNetCDF(ModelFilename);

    typedef typename jif3D::MinMemGravMagCalculator<jif3D::TotalFieldMagneticData> CalculatorType;
    boost::shared_ptr<jif3D::ThreeDGravMagImplementation<jif3D::TotalFieldMagneticData> > Implementation(
        new jif3D::OMPMagneticSusceptibilityImp(inclination, declination, fieldstrength));

    boost::shared_ptr<CalculatorType> Calculator(new CalculatorType(Implementation));

    jif3D::rvec VecResults(Calculator->Calculate(MagModel, Data));

    jif3D::Write3DVectorDataToVTK(ModelFilename + ".mag.vtk", "B",
            std::vector<double>(VecResults.begin(), VecResults.end()), Data.GetMeasPosX(),
            Data.GetMeasPosY(), Data.GetMeasPosZ());


    Calculator->SetDataTransform(
        boost::shared_ptr<jif3D::TotalFieldAnomaly>(
            new jif3D::TotalFieldAnomaly(inclination, declination)));

    jif3D::rvec Results(Calculator->Calculate(MagModel, Data));
    std::vector<double> Err(Results.size());

    double relnoise = 0.0;
    double absnoise = 0.0;

    cout << "Relative noise level : ";
    cin >> relnoise;
    cout << "Absolute noise level : ";
    cin >> absnoise;

    jif3D::AddNoise(Results, relnoise, absnoise);
    std::transform(Results.begin(), Results.end(), Err.begin(),
        [relnoise, absnoise] (const double data)
          { return std::max(absnoise,data*relnoise);});

    jif3D::SaveTotalFieldMagneticMeasurements(ModelFilename + ".mag.nc",
        std::vector<double>(Results.begin(), Results.end()), Data.GetMeasPosX(),
        Data.GetMeasPosY(), Data.GetMeasPosZ(), Err);

    //write the model in .vtk format, at the moment the best plotting option
    MagModel.WriteVTK(ModelFilename + ".vtk");

    jif3D::Write3DDataToVTK(ModelFilename + ".data.vtk", "T",
        std::vector<double>(Results.begin(), Results.end()), Data.GetMeasPosX(),
        Data.GetMeasPosY(), Data.GetMeasPosZ());

  }
