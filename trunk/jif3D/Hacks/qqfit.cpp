//============================================================================
// Name        : qqfit.cpp
// Author      : Dec 6, 2013
// Version     :
// Copyright   : 2013, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTEquations.h"

void QQPlot(std::ofstream &outfile, size_t index, std::vector<double> &Misfit)
  {
    const size_t nimp = Misfit.size();
    const size_t nval = nimp / 8;
    std::vector<double> sorted(nval * 2);
    std::cout << "Processing " << nval << " impedance elements" << std::endl;
    for (size_t i = 0; i < nval; ++i)
      {
        sorted.at(2 * i) = Misfit.at(8 * i + index);
        sorted.at(2 * i + 1) = Misfit.at(8 * i + index + 1);
      }
    std::sort(sorted.begin(), sorted.end());
    std::cout << "Writing out " << sorted.size() << " values for qq-plot" << std::endl;
    for (double q : sorted)
      {
        outfile << q << std::endl;
      }
    //std::copy(sorted.begin(), sorted.end(),
    //		std::ostream_iterator<double>(outfile, "\n"));
  }

int main()
  {
    std::string misfitfilename = jif3D::AskFilename("Name of misfit file: ");

    std::vector<double> Impedances, Errors, Misfit;
    std::vector<double> Frequencies, StatX, StatY, StatZ, C;
    jif3D::ReadImpedancesFromNetCDF(misfitfilename, Frequencies, StatX, StatY, StatZ,
        Misfit, Errors, C);

    std::cout << " Read " << Misfit.size() << " misfit values " << std::endl;
    std::ofstream Zxxqq("Zxx.qq");
    std::ofstream Zxyqq("Zxy.qq");
    std::ofstream Zyxqq("Zyx.qq");
    std::ofstream Zyyqq("Zyy.qq");

    QQPlot(Zxxqq, 0, Misfit);
    QQPlot(Zxyqq, 2, Misfit);
    QQPlot(Zyxqq, 4, Misfit);
    QQPlot(Zyyqq, 6, Misfit);

    double minthresh = -1e6;
    double maxthresh = 1e6;

    std::cout << "Minimum threshold: ";
    std::cin >> minthresh;
    std::cout << "Maximum threshold: ";
    std::cin >> maxthresh;
    std::vector<size_t> Indices;
    for (size_t i = 0; i < Misfit.size(); ++i)
      {
        if (Misfit.at(i) < minthresh || Misfit.at(i) > maxthresh)
          {
            Indices.push_back(i);
          }
      }
    std::string impfilename = jif3D::AskFilename("Name of data file: ");
    jif3D::ReadImpedancesFromNetCDF(impfilename, Frequencies, StatX, StatY, StatZ,
        Impedances, Errors, C);
    std::ofstream indexfile("ind.out");
    const size_t nstat = StatX.size();
    for (size_t ind : Indices)
      {
        Errors.at(ind) = std::abs(Impedances.at(ind));
        size_t stati = ind % (nstat * 8) / 8;
        size_t freqi = ind / (nstat * 8);
        indexfile << ind << " " << stati << " " << " " << freqi << std::endl;
      }

    for (size_t i = 0; i < Errors.size() - 1; i += 2)
      {
        Errors.at(i) = std::max(Errors.at(i),Errors.at(i+1));
      }
    std::cout << "Modified: " << Indices.size() << " out of " << Impedances.size()
        << " data " << std::endl;
    jif3D::WriteImpedancesToNetCDF(impfilename + ".qq.nc", Frequencies, StatX, StatY,
        StatZ, Impedances, Errors, C);
  }
