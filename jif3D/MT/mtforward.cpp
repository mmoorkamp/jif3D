//============================================================================
// Name        : mtforward.cpp
// Author      : Jul 14, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include "../Global/Jif3DGlobal.h"
#include "../Global/FileUtil.h"
#include "../Global/convert.h"
#include "../Global/Noise.h"
#include "../ModelBase/VTKTools.h"
#include "X3DModel.h"
#include "X3DMTCalculator.h"
#include "X3DTipperCalculator.h"
#include "ReadWriteImpedances.h"
#include "MTData.h"
#include "TipperData.h"

#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>

#ifdef HAVEOPENMP
#include <omp.h>
#endif

namespace po = boost::program_options;

int main(int argc, char *argv[])
  {

    std::string X3DName = "x3d";
    double DistDeviation = 0;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("threads", po::value<int>(),
        "The number of openmp threads")("x3dname",
        po::value(&X3DName)->default_value("x3d"), "The name of the executable for x3d")(
        "debug", "Show debugging output.")("tempdir", po::value<std::string>(),
        "The name of the directory to store temporary files in")("opt",
        "Use opt for Green's function calculation in x3d.")("dist",
        po::value(&DistDeviation)->default_value(0.0),
        "Standard deviation for random distortion, 0 means no distortion");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }
#ifdef HAVEOPENMP
    if (vm.count("threads"))
      {
        omp_set_num_threads(vm["threads"].as<int>());
      }
#endif

    boost::filesystem::path TempDir = boost::filesystem::current_path();
    if (vm.count("tempdir"))
      {
        TempDir = vm["tempdir"].as<std::string>();
        if (!boost::filesystem::is_directory(TempDir))
          {
            std::cerr << TempDir.string() << " is not a directory or does not exist ! \n";
            return 500;
          }
      }

    std::string modelfilename = jif3D::AskFilename("Filename: ", true);
    jif3D::X3DModel MTModel;
    jif3D::MTData DataMT;
    jif3D::TipperData DataTipper;
    MTModel.ReadNetCDF(modelfilename);
    std::cout << "Model size: " << MTModel.GetXCellSizes().size() << " "
        << MTModel.GetYCellSizes().size() << " " << MTModel.GetZCellSizes().size()
        << std::endl;
    double minx, miny, maxx, maxy, deltax, deltay, z;
    //ask for the measurement grid specifications
    //first x-direction
    std::cout << "Minimum x-position: ";
    std::cin >> minx;
    std::cout << "Maximum x-position: ";
    std::cin >> maxx;
    std::cout << "Delta x: ";
    std::cin >> deltax;
    //then y-direction
    std::cout << "Minimum y-position: ";
    std::cin >> miny;
    std::cout << "Maximum y-position: ";
    std::cin >> maxy;
    std::cout << "Delta y: ";
    std::cin >> deltay;
    //all measurements have to be at the same height.
    std::cout << "Z-level: ";
    std::cin >> z;
    //setup the measurements in the forward modelling code
    const size_t nmeasx = boost::numeric_cast<size_t>((maxx - minx) / deltax);
    const size_t nmeasy = boost::numeric_cast<size_t>((maxy - miny) / deltay);
    for (size_t i = 0; i <= nmeasx; ++i)
      {
        for (size_t j = 0; j <= nmeasy; ++j)
          {
            DataMT.AddMeasurementPoint(minx + i * deltax, miny + j * deltay, z);
            DataTipper.AddMeasurementPoint(minx + i * deltax, miny + j * deltay, z);
          }
      }

    std::string outfilename = jif3D::AskFilename("Output filename: ", false);
    std::cout << "Frequencies: ";
    double currfreq = 1.0;
    std::vector<double> frequencies;
    try
      {
        while (true)
          {
            std::string input;
            std::cin >> input;
            jif3D::convert(input, currfreq);
            frequencies.push_back(currfreq);
          }
      } catch (jif3D::BadConversion &e)
      {

      }
    std::sort(frequencies.begin(), frequencies.end());
    DataMT.SetFrequencies(frequencies);
    DataTipper.SetFrequencies(frequencies);

    std::vector<double> StatNum(DataMT.GetMeasPosX().size());
    std::iota(StatNum.begin(), StatNum.end(), 0);
    //! Write scalar data with 3D coordinate information into a .vtk file for plotting
    jif3D::Write3DDataToVTK(outfilename + ".vtk", "MTStats", StatNum,
        DataMT.GetMeasPosX(), DataMT.GetMeasPosY(), DataMT.GetMeasPosZ());
    DataMT.CompleteObject();
    std::cout << "Calculating forward response " << std::endl;
    jif3D::X3DMTCalculator Calculator(TempDir, X3DName);
    if (vm.count("opt"))
      {
        Calculator.SetGreenType1(jif3D::GreenCalcType::opt);
        Calculator.SetGreenType4(jif3D::GreenCalcType::opt);
      }
    jif3D::rvec Impedances(Calculator.Calculate(MTModel, DataMT));

    double relnoise = 0.0;
    double absnoise = 0.0;
    double range = 0.0;
    std::cout << "Relative noise level: ";
    std::cin >> relnoise;
    std::cout << "Absolute noise level: ";
    std::cin >> absnoise;
    std::cout << "Dynamic range: ";
    std::cin >> range;

    const size_t nimp = Impedances.size();
    jif3D::rvec Errors(nimp, 0.0);

    for (size_t i = 0; i < nimp; i += 8)
      {
        double maximp = *std::max_element(Impedances.begin() + i,
            Impedances.begin() + i + 8, [](double a, double b)
              { return std::abs(a) < std::abs(b);});
        double threshold = std::max(absnoise, std::abs(range * maximp));
        std::transform(Impedances.begin() + i, Impedances.begin() + i + 8,
            Errors.begin() + i, [&](double d) -> double
              { return std::max(std::abs(d * relnoise),threshold);});
      }

    jif3D::AddNoise(Impedances, relnoise, Errors);
    std::vector<double> C, NoC;
    std::vector<std::string> Names(DataMT.GetMeasPosX().size(), "");
    if (DistDeviation > 0.0)
      {
        boost::lagged_fibonacci607 generator(static_cast<unsigned int>(std::time(0)));

        boost::normal_distribution<> diag_dist(1.0, DistDeviation);
        boost::normal_distribution<> offdiag_dist(0.0, DistDeviation);
        boost::variate_generator<boost::lagged_fibonacci607&, boost::normal_distribution<> > diag_noise(
            generator, diag_dist);
        boost::variate_generator<boost::lagged_fibonacci607&, boost::normal_distribution<> > offdiag_noise(
            generator, offdiag_dist);
        const size_t nstats = StatNum.size();
        const size_t nfreq = frequencies.size();

        for (size_t i = 0; i < nstats; ++i)
          {
            double Cxx = diag_noise();
            double Cxy = offdiag_noise();
            double Cyx = offdiag_noise();
            double Cyy = diag_noise();
            for (size_t j = 0; j < nfreq; ++j)
              {
                size_t offset = (j * nstats + i) * 8;

                std::vector<double> NewZ;
                NewZ.push_back(Cxx * Impedances(offset) + Cxy * Impedances(offset + 4));
                NewZ.push_back(
                    Cxx * Impedances(offset + 1) + Cxy * Impedances(offset + 5));
                NewZ.push_back(
                    Cxx * Impedances(offset + 2) + Cxy * Impedances(offset + 6));
                NewZ.push_back(
                    Cxx * Impedances(offset + 3) + Cxy * Impedances(offset + 7));
                NewZ.push_back(Cyx * Impedances(offset) + Cyy * Impedances(offset + 4));
                NewZ.push_back(
                    Cyx * Impedances(offset + 1) + Cyy * Impedances(offset + 5));
                NewZ.push_back(
                    Cyx * Impedances(offset + 2) + Cyy * Impedances(offset + 6));
                NewZ.push_back(
                    Cyx * Impedances(offset + 3) + Cyy * Impedances(offset + 7));
                std::copy(NewZ.begin(), NewZ.end(), Impedances.begin() + offset);
              }
            C.push_back(Cxx);
            C.push_back(Cxy);
            C.push_back(Cyx);
            C.push_back(Cyy);
            NoC.push_back(1.0);
            NoC.push_back(0.0);
            NoC.push_back(0.0);
            NoC.push_back(1.0);
          }
        DataMT.WriteNetCDF(outfilename + "dist.nc");

      }
    DataMT.WriteNetCDF(outfilename);
    DataMT.WriteModEM(outfilename + ".dat");

    jif3D::X3DTipperCalculator TipCalc(TempDir, X3DName);
    jif3D::rvec Tipper = TipCalc.Calculate(MTModel, DataTipper);

    const size_t ntip = Tipper.size();
    jif3D::rvec TipErrors(ntip, 0.0);

    for (size_t i = 0; i < ntip; i += 4)
      {
        double maxtip = *std::max_element(Tipper.begin() + i, Tipper.begin() + i + 4,
            [](double a, double b)
              { return std::abs(a) < std::abs(b);});
        double threshold = std::max(absnoise, std::abs(range * maxtip));
        std::transform(Tipper.begin() + i, Tipper.begin() + i + 4, TipErrors.begin() + i,
            [&](double d) -> double
              { return std::max(std::abs(d * relnoise),threshold);});
      }

    jif3D::AddNoise(Tipper, relnoise, TipErrors);

    DataTipper.SetDataAndErrors(std::vector<double>(Tipper.begin(), Tipper.end()),
        std::vector<double>(TipErrors.begin(), TipErrors.end()));
    DataTipper.WriteNetCDF(outfilename + ".tip.nc");

    MTModel.WriteVTK(modelfilename + ".vtk");
    MTModel.WriteModEM(modelfilename + ".dat");
  }

