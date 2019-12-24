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
#include "../MT/TipperData.h"

void QQPlot(std::ofstream &outfile, size_t index, std::vector<double> &Misfit)
  {
    const size_t nimp = Misfit.size();
    const size_t nval = nimp / 4;
    std::vector<double> sorted(nval * 2);
    std::cout << "Processing " << nval << " tipper elements" << std::endl;
    for (size_t i = 0; i < nval; ++i)
      {
        sorted.at(2 * i) = Misfit.at(4 * i + index);
        sorted.at(2 * i + 1) = Misfit.at(4 * i + index + 1);
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
    jif3D::TipperData TipperData;
    std::vector<double> Tipper, Errors, Misfit;
    TipperData.ReadNetCDF(misfitfilename);
    Misfit = TipperData.GetData();
    Errors = TipperData.GetErrors();

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
        if (Misfit.at(i) < minthresh || Misfit.at(i) > maxthresh)
          {
            Indices.push_back(i);
          }
      }
    std::string tipfilename = jif3D::AskFilename("Name of data file: ");
    TipperData.ReadNetCDF(tipfilename);
    Tipper = TipperData.GetData();
    std::ofstream indexfile("ind.out");

    const size_t nstat = TipperData.GetMeasPosX().size();
    for (size_t ind : Indices)
      {
        Errors.at(ind) = 10 * std::abs(Tipper.at(ind));
        size_t stati = ind % (nstat * 4) / 4;
        size_t freqi = ind / (nstat * 4);
        indexfile << ind << " " << stati << " " << " " << freqi << std::endl;
      }

    for (size_t i = 0; i < Errors.size() - 1; i += 2)
      {
        Errors.at(i) = std::max(Errors.at(i), Errors.at(i + 1));
        Errors.at(i + 1) = Errors.at(i);
      }
    std::cout << "Modified: " << Indices.size() << " out of " << Tipper.size() << " data "
        << std::endl;
    TipperData.SetErrors(Errors);
    TipperData.WriteNetCDF(tipfilename + ".qq.nc");
  }
