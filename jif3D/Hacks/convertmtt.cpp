//============================================================================
// Name        : convertmtt.cpp
// Author      : Sep 21, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/math/special_functions/relative_difference.hpp>
#include <GeographicLib/UTMUPS.hpp>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTData.h"
#include "../MT/TipperData.h"
#include "../MT/MTEquations.h"

using namespace GeographicLib;
namespace po = boost::program_options;

void Interpolate(const std::vector<double> &xvalues, const std::vector<double> &yvalues,
    const std::vector<double> &err, const std::vector<double> &InterpX,
    std::vector<double> &InterpY)
  {
    const double invalid = 1e10;
    const size_t ninter = InterpX.size();
    size_t currindex = 0;
    for (size_t i = 0; i < ninter; ++i)
      {
        while (currindex < xvalues.size() && xvalues[currindex] >= InterpX[i])
          ++currindex;
        if (currindex == 0 && InterpX[i] == xvalues[currindex])
          {
            InterpY.push_back(yvalues[currindex]);
          }
        else
          {
            if (currindex == 0 && InterpX[i] != xvalues[currindex])
              {
                InterpY.push_back(0.0);
              }
            else
              {
                if (currindex == xvalues.size() && InterpX[i] != xvalues[currindex - 1])
                  {
                    InterpY.push_back(0.0);
                  }
                else
                  {
                    if (InterpX[i] == xvalues[currindex - 1])
                      {
                        InterpY.push_back(yvalues[currindex - 1]);

                      }
                    else
                      {
                        if (std::abs(yvalues[currindex - 1]) > invalid
                            || std::abs(yvalues[currindex]) > invalid
                            || err[currindex - 1] >= std::abs(yvalues[currindex - 1])
                            || err[currindex] >= std::abs(yvalues[currindex]))
                          {
                            InterpY.push_back(0.0);
                          }
                        else
                          {
                            double curry = yvalues[currindex - 1]
                                + (yvalues[currindex] - yvalues[currindex - 1])
                                    / (xvalues[currindex] - xvalues[currindex - 1])
                                    * (InterpX[i] - xvalues[currindex - 1]);
                            if (std::isnan(curry))
                              {
                                std::cout << " Encountered nan, ignoring " << std::endl;
                                InterpY.push_back(0.0);
                              }
                            else
                              {
                                InterpY.push_back(curry);
                              }
                          }
                      }
                  }
              }
          }

      }
  }

void InterpolateError(const std::vector<double> &xvalues,
    const std::vector<double> &yvalues, const std::vector<double> &InterpX,
    std::vector<double> &InterpY)
  {
    const size_t ninter = InterpX.size();
    size_t currindex = 0;
    const double largeerr = 1e10;

    for (size_t i = 0; i < ninter; ++i)
      {
        double error = largeerr;
        while (xvalues[currindex] >= InterpX[i] && currindex < xvalues.size())
          ++currindex;

        if (currindex == 0 && InterpX[i] == xvalues[currindex])
          {
            error = yvalues[currindex];
          }
        else
          {
            if (currindex == 0 && InterpX[i] != xvalues[currindex])
              {
                error = largeerr;
              }
            else
              {
                if (currindex == xvalues.size() && InterpX[i] != xvalues[currindex - 1])
                  {
                    error = largeerr;
                  }
                else
                  {
                    if (InterpX[i] == xvalues[currindex - 1])
                      {
                        error = yvalues[currindex - 1];
                      }
                    else
                      {

                        double curry = yvalues[currindex - 1]
                            + (yvalues[currindex] - yvalues[currindex - 1])
                                / (xvalues[currindex] - xvalues[currindex - 1])
                                * (InterpX[i] - xvalues[currindex - 1]);
                        if (std::isnan(curry))
                          {
                            std::cout << " Encountered nan, setting large error "
                                << std::endl;
                            error = largeerr;
                          }
                        else
                          {
                            error = curry;
                          }

                      }
                  }
              }
          }

        if (error <= 0.0)
          {
            std::cerr << "Zero error" << std::endl;
            error = largeerr;
          }
        InterpY.push_back(error);
      }

  }

void SetMissing(std::complex<double> &value, double &error)
  {
    value = std::complex<double>(1.0, 1.0);
    error = 1000.0;
  }

void SetMissingTip(std::complex<double> &value, double &error)
  {
    value = std::complex<double>(0.0, 0.0);
    error = 1.0;
  }

int main(int argc, char *argv[])
  {

    std::string DistFilename, FreqFilename;
    int utmzone = -1;
    double minfreq, maxfreq;
    double depth;
    double weightlevel1, weightlevel2;
    bool rotate;
    bool rotangle;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("distfile",
        po::value(&DistFilename),
        "Read in ascii file with distortion estimates for each station (4 values per line)")(
        "toutm", po::value(&utmzone)->default_value(-1),
        "Project coordinates to utm zone, -1 means no projection, 0 automatic zone detection")(
        "freqfile", po::value(&FreqFilename), "File with frequencies for interpolation")(
        "minfreq", po::value(&minfreq)->default_value(-1),
        "Minimum frequency kept for interpolation")("maxfreq",
        po::value(&maxfreq)->default_value(1e32), "Max frequency kept for interpolation")(
        "depth", po::value(&depth)->default_value(0.0),
        "Set depth in m (negative of elevation) to a fixed value. "
            "If >= 100,000 keep elevations from input files")("weightlevel1",
        po::value(&weightlevel1)->default_value(0.3),
        "The threshold for the weight in the j-file, error will be multiplied by inverse")(
        "weightlevel2", po::value(&weightlevel2)->default_value(0.2),
        "The threshold for the weight in the j-file, error will be set extremely high")(
        "rotate", po::value(&rotate)->default_value(true),
        "Rotate all data to a consistent angle")("rotangle",
        po::value(&rotangle)->default_value(0.0), "The angle to rotate to");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
    if (vm.count("help"))
      {
        std::cout << desc << "\n";
        return 1;
      }

    std::string filename = jif3D::AskFilename("File with station names/positions: ");
    std::ifstream StationFile;

    StationFile.open(filename.c_str());
    char dummy[2048];
    size_t nstats = 0;
    while (StationFile.good())
      {
        StationFile.getline(dummy, 2048);
        if (StationFile.good())
          {
            ++nstats;
          }
      }
    StationFile.close();

    std::vector<double> C;
    if (vm.count("distfile"))
      {
        std::ifstream distfile(DistFilename.c_str());
        std::copy(std::istream_iterator<double>(distfile),
            std::istream_iterator<double>(), std::back_inserter(C));
        if (C.size() != nstats * 4)
          {
            std::cerr << "Number of elements in Distortion file " << C.size()
                << " does not match 4* number of stations " << nstats << std::endl;
            return 100;
          }
      }
    else
      {
        C.resize(nstats * 4, 0);
        for (size_t i = 0; i < nstats; ++i)
          {
            C.at(i * 4) = 1.0;
            C.at(i * 4 + 3) = 1.0;
          }
      }

    std::cout << "Number of stations: " << nstats << std::endl;
    if (nstats > 0)
      {
        std::vector<double> MasterFreqs;
        std::vector<double> StatXCoord(nstats), StatYCoord(nstats), StatZCoord(nstats),
            Angles(nstats);
        std::vector<std::string> Names;
        StationFile.open(filename.c_str());
        std::string StationName;
        double XC, YC, ZC;

        size_t stationindex = 0;
        //read in all the data
        std::vector<std::vector<double>> AllFreqs, AllTipFreqs, AllImpedances,
            AllImpErrors, AllTipper, AllTipErr;
        while (StationFile.good() && stationindex < nstats)
          {
            std::vector<double> CurrFrequencies, TipCurrFrequencies;
            std::vector<double> CurrImpedances, CurrErrors, CurrTip, CurrTipErr,
                CurrWeight;
            double CurrAngle;
            StationFile >> StationName;

            if (StationFile.good())
              {
                std::cout << stationindex << " " << StationName << " "
                    << CurrFrequencies.size() << " " << CurrAngle << std::endl;
                Names.push_back(StationName);
                std::string extension = jif3D::GetFileExtension(StationName);
                if (extension == ".mtt")
                  {
                    jif3D::ReadImpedancesFromMTT(StationName, CurrFrequencies,
                        CurrImpedances, CurrErrors, CurrTip, CurrTipErr);
                  }
                else
                  {
                    jif3D::ReadImpedancesFromJ(StationName, CurrFrequencies, XC, YC, ZC,
                        CurrImpedances, CurrErrors, CurrWeight, CurrAngle);
                    for (size_t i = 0; i < CurrErrors.size(); ++i)
                      {
                        if (CurrWeight.at(i) > 0.0)
                          {
                            if (CurrWeight.at(i) < weightlevel1)
                              {
                                std::cout << " Reweighting level 1: " << CurrWeight.at(i)
                                    << std::endl;
                                CurrErrors.at(i) /= CurrWeight.at(i);
                              }
                            if (CurrWeight.at(i) < weightlevel2)
                              {
                                std::cout << " Reweighting level 2 " << CurrWeight.at(i)
                                    << std::endl;
                                CurrErrors.at(i) *= 1e6;
                              }
                          }
                      }
                    jif3D::ReadTipperFromJ(StationName, TipCurrFrequencies, XC, YC, ZC,
                        CurrTip, CurrTipErr, CurrAngle);
                  }
                if (rotate)
                  {
                    double rangle = rotangle - CurrAngle;
                    std::cout << "Rotating by " << rangle << "  degrees." << std::endl;
                    rangle = rangle / 180.0 * boost::math::constants::pi<double>();
                    CurrImpedances = jif3D::RotateImpedanceVector(rangle, CurrImpedances);
                    CurrTip = jif3D::RotateTipperVector(rangle, CurrTip);
                    CurrAngle = rotangle;
                  }
                StatXCoord.at(stationindex) = XC;
                StatYCoord.at(stationindex) = YC;
                StatZCoord.at(stationindex) = ZC;
                Angles.at(stationindex) = CurrAngle;
                std::copy(CurrFrequencies.begin(), CurrFrequencies.end(),
                    std::back_inserter(MasterFreqs));
                AllFreqs.push_back(CurrFrequencies);
                AllTipFreqs.push_back(TipCurrFrequencies);
                AllImpedances.push_back(CurrImpedances);
                AllImpErrors.push_back(CurrErrors);
                AllTipper.push_back(CurrTip);
                AllTipErr.push_back(CurrTipErr);
              }
            std::cout << std::endl;
            ++stationindex;
          }

        if (vm.count("freqfile"))
          {
            std::vector<double> filefreqs;
            std::ifstream freqfile(FreqFilename.c_str());
            std::copy(std::istream_iterator<double>(freqfile),
                std::istream_iterator<double>(), std::back_inserter(filefreqs));
            MasterFreqs = filefreqs;
          }
        else
          {
            //sort the frequencies across all files and keep only the ones that differ by more than 1%
            std::sort(MasterFreqs.begin(), MasterFreqs.end(), std::greater<double>());
            MasterFreqs.erase(
                std::unique(MasterFreqs.begin(), MasterFreqs.end(), [](double a, double b)
                  { return boost::math::relative_difference(a,b) < 0.01;}),
                MasterFreqs.end());
            auto end = std::upper_bound(MasterFreqs.begin(), MasterFreqs.end(), minfreq,
                std::greater<double>());
            auto start = std::lower_bound(MasterFreqs.begin(), MasterFreqs.end(), maxfreq,
                std::greater<double>());
            std::vector<double> tmp(start, end);
            MasterFreqs = tmp;
          }
        std::cout << "Master frequencies \n" << std::endl;
        std::cout << std::setw(12) << "Filename ";
        for (double freq : MasterFreqs)
          {
            std::cout << std::setw(10) << freq << std::endl;
          }
        std::cout << std::endl;

        std::vector<double> Impedances, Errors, Tipper, TipErr;
        const size_t nfreq = MasterFreqs.size();
        Impedances.resize(nstats * nfreq * 8);
        Errors.resize(nstats * nfreq * 8);
        Tipper.resize(nstats * nfreq * 4, 0.0);
        TipErr.resize(nstats * nfreq * 4, 1.0);

        for (size_t i = 0; i < nstats; ++i)
          {
            std::vector<double> SiteFreqs(AllFreqs.at(i));
            const size_t ncurrfreq = SiteFreqs.size();
            std::vector<double> ZxxR, ZxxI, ZxyR, ZxyI, ZyxR, ZyxI, ZyyR, ZyyI, dZxx,
                dZxy, dZyx, dZyy, TxR, TxI, TyR, TyI, dTx, dTy;
            for (size_t j = 0; j < ncurrfreq; ++j)
              {
                ZxxR.push_back(AllImpedances.at(i).at(j * 8));
                ZxxI.push_back(AllImpedances.at(i).at(j * 8 + 1));
                dZxx.push_back(AllImpErrors.at(i).at(j * 8));

                ZxyR.push_back(AllImpedances.at(i).at(j * 8 + 2));
                ZxyI.push_back(AllImpedances.at(i).at(j * 8 + 3));
                dZxy.push_back(AllImpErrors.at(i).at(j * 8 + 2));

                ZyxR.push_back(AllImpedances.at(i).at(j * 8 + 4));
                ZyxI.push_back(AllImpedances.at(i).at(j * 8 + 5));
                dZyx.push_back(AllImpErrors.at(i).at(j * 8 + 4));

                ZyyR.push_back(AllImpedances.at(i).at(j * 8 + 6));
                ZyyI.push_back(AllImpedances.at(i).at(j * 8 + 7));
                dZyy.push_back(AllImpErrors.at(i).at(j * 8 + 6));
              }
            std::vector<double> IZxxR, IZxxI, IZxyR, IZxyI, IZyxR, IZyxI, IZyyR, IZyyI,
                IdZxx, IdZxy, IdZyx, IdZyy, ITxR, ITxI, ITyR, ITyI, IdTx, IdTy;

            std::vector<double> TipFreqs(AllTipFreqs.at(i));
            size_t ntipfreq = TipFreqs.size();
            for (size_t j = 0; j < ntipfreq; ++j)
              {
                TxR.push_back(AllTipper.at(i).at(j * 4));
                TxI.push_back(AllTipper.at(i).at(j * 4 + 1));
                dTx.push_back(AllTipErr.at(i).at(j * 4));

                TyR.push_back(AllTipper.at(i).at(j * 4 + 2));
                TyI.push_back(AllTipper.at(i).at(j * 4 + 3));
                dTy.push_back(AllTipErr.at(i).at(j * 4 + 2));
              }

            Interpolate(SiteFreqs, ZxxR, dZxx, MasterFreqs, IZxxR);
            Interpolate(SiteFreqs, ZxxI, dZxx, MasterFreqs, IZxxI);
            Interpolate(SiteFreqs, ZxyR, dZxy, MasterFreqs, IZxyR);
            Interpolate(SiteFreqs, ZxyI, dZxy, MasterFreqs, IZxyI);
            Interpolate(SiteFreqs, ZyxR, dZyx, MasterFreqs, IZyxR);
            Interpolate(SiteFreqs, ZyxI, dZyx, MasterFreqs, IZyxI);
            Interpolate(SiteFreqs, ZyyR, dZyy, MasterFreqs, IZyyR);
            Interpolate(SiteFreqs, ZyyI, dZyy, MasterFreqs, IZyyI);
            InterpolateError(SiteFreqs, dZxx, MasterFreqs, IdZxx);
            InterpolateError(SiteFreqs, dZxy, MasterFreqs, IdZxy);
            InterpolateError(SiteFreqs, dZyx, MasterFreqs, IdZyx);
            InterpolateError(SiteFreqs, dZyy, MasterFreqs, IdZyy);
            if (TipFreqs.size() > 0)
              {
                Interpolate(TipFreqs, TxR, dTx, MasterFreqs, ITxR);
                Interpolate(TipFreqs, TxI, dTx, MasterFreqs, ITxI);
                Interpolate(TipFreqs, TyR, dTy, MasterFreqs, ITyR);
                Interpolate(TipFreqs, TyI, dTy, MasterFreqs, ITyI);
                InterpolateError(TipFreqs, dTx, MasterFreqs, IdTx);
                InterpolateError(TipFreqs, dTy, MasterFreqs, IdTy);
              }
            else
              {
                ITxR.resize(nfreq, 0);
                ITxI.resize(nfreq, 0);
                ITyR.resize(nfreq, 0);
                ITyI.resize(nfreq, 0);
                IdTx.resize(nfreq, 100);
                IdTy.resize(nfreq, 100);
              }

            for (size_t j = 0; j < nfreq; ++j)
              {
                size_t index = j * nstats * 8 + i * 8;
                Impedances.at(index) = IZxxR.at(j);
                Impedances.at(index + 1) = IZxxI.at(j);
                Impedances.at(index + 2) = IZxyR.at(j);
                Impedances.at(index + 3) = IZxyI.at(j);
                Impedances.at(index + 4) = IZyxR.at(j);
                Impedances.at(index + 5) = IZyxI.at(j);
                Impedances.at(index + 6) = IZyyR.at(j);
                Impedances.at(index + 7) = IZyyI.at(j);
                Errors.at(index) = IdZxx.at(j);
                Errors.at(index + 1) = IdZxx.at(j);
                Errors.at(index + 2) = IdZxy.at(j);
                Errors.at(index + 3) = IdZxy.at(j);
                Errors.at(index + 4) = IdZyx.at(j);
                Errors.at(index + 5) = IdZyx.at(j);
                Errors.at(index + 6) = IdZyy.at(j);
                Errors.at(index + 7) = IdZyy.at(j);
                index = j * nstats * 4 + i * 4;
                Tipper.at(index) = ITxR.at(j);
                Tipper.at(index + 1) = ITxI.at(j);
                Tipper.at(index + 2) = ITyR.at(j);
                Tipper.at(index + 3) = ITyI.at(j);
                TipErr.at(index) = IdTx.at(j);
                TipErr.at(index + 1) = IdTx.at(j);
                TipErr.at(index + 2) = IdTy.at(j);
                TipErr.at(index + 3) = IdTy.at(j);
              }
          }

        if (utmzone >= 0)
          {
            int zone;
            bool northp;
            double x, y, gamma, k;
            std::cout << "Converting coordinates " << std::endl;
            for (size_t i = 0; i < StatXCoord.size(); ++i)
              {
                if (utmzone == 0)
                  {
                    //automatic zone determination
                    UTMUPS::Forward(StatXCoord.at(i), StatYCoord.at(i), zone, northp, x,
                        y);
                    std::string zonestr = UTMUPS::EncodeZone(zone, northp);
                    std::cout << std::fixed << std::setprecision(2) << zonestr << " " << x
                        << " " << y << "\n";
                  }
                else
                  {
                    zone = utmzone;
                    UTMUPS::Forward(StatXCoord.at(i), StatYCoord.at(i), zone, northp, x,
                        y, gamma, k, utmzone);
                  }
                StatXCoord.at(i) = y;
                StatYCoord.at(i) = x;
              }

          }
        if (depth < 1e5)
          {
            std::fill(StatZCoord.begin(), StatZCoord.end(), depth);
          }

        std::string outfilename = jif3D::AskFilename("Output file: ", false);
        jif3D::MTData MTData;
        MTData.SetMeasurementPoints(StatXCoord, StatYCoord, StatZCoord);
        MTData.SetDataAndErrors(Impedances, Errors);
        MTData.SetFrequencies(MasterFreqs);
        MTData.SetRotAngles(Angles);
        MTData.SetNames(Names);
        MTData.CompleteObject();
        MTData.WriteNetCDF(outfilename);
        MTData.WriteMeasurementPoints(outfilename + ".vtk");

        jif3D::TipperData TipperData;
        TipperData.SetMeasurementPoints(StatXCoord, StatYCoord, StatZCoord);
        TipperData.SetDataAndErrors(Tipper, TipErr);
        TipperData.SetFrequencies(MasterFreqs);
        TipperData.CompleteObject();
        TipperData.WriteNetCDF(outfilename + ".tip.nc");

      }
    else
      {
        std::cerr << "No usable station information in file !" << std::endl;
        return 100;
      }
  }
