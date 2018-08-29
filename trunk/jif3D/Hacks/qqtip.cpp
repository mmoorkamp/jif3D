//============================================================================
// Name        : qqtip.cpp
// Author      : Jun 21, 2018
// Version     :
// Copyright   : 2018, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTEquations.h"

void QQPlot(std::ofstream &outfile, size_t index, jif3D::rvec &Misfit)
  {
    const size_t nimp = Misfit.size();
    const size_t nval = nimp / 4;
    jif3D::rvec sorted(nval * 2);
    std::cout << "Processing " << nval << " tipper elements" << std::endl;
    for (size_t i = 0; i < nval; ++i)
      {
        sorted(2 * i) = Misfit(4 * i + index);
        sorted(2 * i + 1) = Misfit(4 * i + index + 1);
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

    jif3D::rvec Tipper, Errors, Misfit;
    std::vector<double> Frequencies, StatX, StatY, StatZ;
    jif3D::ReadTipperFromNetCDF(misfitfilename, Frequencies, StatX, StatY, StatZ,
        Misfit, Errors);

    std::cout << " Read " << Misfit.size() << " misfit values " << std::endl;
    std::ofstream Txqq("Tx.qq");
    std::ofstream Tyqq("Ty.qq");


    QQPlot(Txqq, 0, Misfit);
    QQPlot(Tyqq, 2, Misfit);

    double minthresh = -1e6;
    double maxthresh = 1e6;

    std::cout << "Minimum threshold: ";
    std::cin >> minthresh;
    std::cout << "Maximum threshold: ";
    std::cin >> maxthresh;
    std::vector<size_t> Indices;
    for (size_t i = 0; i < Misfit.size(); ++i)
      {
        if (Misfit(i) < minthresh || Misfit(i) > maxthresh)
          {
            Indices.push_back(i);
          }
      }
    std::string tipfilename = jif3D::AskFilename("Name of data file: ");
    jif3D::ReadTipperFromNetCDF(tipfilename, Frequencies, StatX, StatY, StatZ,
        Tipper, Errors);
    std::ofstream indexfile("ind.out");

    const size_t nstat = StatX.size();
    for (size_t ind : Indices)
      {
        Errors(ind) = 10 * std::abs(Tipper(ind));
        size_t stati = ind % (nstat * 4) / 4;
        size_t freqi = ind / (nstat * 4);
        indexfile << ind << " " << stati << " " << " " << freqi << std::endl;
      }

    for (size_t i = 0; i < Errors.size() - 1; i += 2)
      {
        Errors(i) = std::max(Errors(i),Errors(i+1));
      }
    std::cout << "Modified: " << Indices.size() << " out of " << Tipper.size()
        << " data " << std::endl;
    jif3D::WriteTipperToNetCDF(tipfilename + ".qq.nc", Frequencies, StatX, StatY,
        StatZ, Tipper, Errors);
  }
