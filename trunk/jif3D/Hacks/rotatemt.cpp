//============================================================================
// Name        : writemtt.cpp
// Author      : Feb 21, 2011
// Version     :
// Copyright   : 2011, mmoorkamp
//============================================================================

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <boost/math/constants/constants.hpp>
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../Global/VecMat.h"
#include "../Global/ReadWriteSparseMatrix.h"
#include "../Global/convert.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTEquations.h"
#include "../MT/MTCovar.h"
#include "../Inversion/MatrixTools.h"
#include "../ModelBase/VTKTools.h"


namespace po = boost::program_options;

int main(int argc, char *argv[])
  {
    bool RotCenter = false;
    po::options_description desc("General options");
    desc.add_options()("help", "produce help message")("center",
        po::value(&RotCenter)->default_value(false),
        "Rotate around the center of all the sites, otherwise rotates around zero.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    std::string ncfilename = jif3D::AskFilename("Name of netcdf file: ");

    jif3D::rvec Impedances, Errors;
    std::vector<double> Frequencies, StatX, StatY, StatZ, C;
    jif3D::ReadImpedancesFromNetCDF(ncfilename, Frequencies, StatX, StatY, StatZ,
        Impedances, Errors, C);

    std::cout << "Rotation angle [degree]: ";
    double dangle = 0.0;
    std::cin >> dangle;
    double rangle = dangle / 180.0 * boost::math::constants::pi<double>();
    jif3D::rvec RotImp = jif3D::RotateImpedanceVector(rangle, Impedances);
    const size_t nelem = 8;
    const size_t ntensor = Impedances.size() / nelem;
    jif3D::comp_mat InvCov(Errors.size(), Errors.size());
    for (size_t i = 0; i < ntensor; ++i)
      {
        jif3D::rmat InvCovOrig(4, 4, 0.0), CovOrig(4, 4, 0.0);
        CovOrig(0, 0) = pow(Errors(i * nelem), 2);
        CovOrig(1, 1) = pow(Errors(i * nelem + 2), 2);
        CovOrig(2, 2) = pow(Errors(i * nelem + 4), 2);
        CovOrig(3, 3) = pow(Errors(i * nelem + 6), 2);

        InvCovOrig(0, 0) = 1.0 / CovOrig(0, 0);
        InvCovOrig(1, 1) = 1.0 / CovOrig(1, 1);
        InvCovOrig(2, 2) = 1.0 / CovOrig(2, 2);
        InvCovOrig(3, 3) = 1.0 / CovOrig(3, 3);

        jif3D::rmat RotInvCov = jif3D::RotateMTCovar(rangle, InvCovOrig, true);
        jif3D::rmat RotCov = jif3D::RotateMTCovar(rangle, CovOrig, false);
        for (size_t j = 0; j < 4; ++j)
          {
            for (size_t k = 0; k < 4; ++k)
              {
                InvCov(i * nelem + j * 2, i * nelem + k * 2) = RotInvCov(j, k);
                InvCov(i * nelem + j * 2 + 1, i * nelem + k * 2 + 1) = RotInvCov(j, k);
              }
          }
        Errors(i * nelem) = sqrt(RotCov(0, 0));
        Errors(i * nelem + 2) = sqrt(RotCov(1, 1));
        Errors(i * nelem + 4) = sqrt(RotCov(2, 2));
        Errors(i * nelem + 6) = sqrt(RotCov(3, 3));

      }

    double MeanX = 0.0;
    double MeanY = 0.0;

    if (RotCenter)
      {
        MeanX = std::accumulate(StatX.begin(), StatX.end(), 0.0) / StatX.size();
        MeanY = std::accumulate(StatY.begin(), StatY.end(), 0.0) / StatY.size();
      }
    //we rotate the coordinates clockwise
    for (size_t i = 0; i < StatX.size(); ++i)
      {

        double newx = (StatX.at(i) - MeanX) * cos(rangle)
            + (StatY.at(i) - MeanY) * sin(rangle);
        double newy = (StatY.at(i) - MeanY) * cos(rangle)
            - (StatX.at(i) - MeanX) * sin(rangle);
        StatX.at(i) = MeanX + newx;
        StatY.at(i) = MeanY + newy;
      }

    jif3D::WriteSparseMatrixToNetcdf(ncfilename + "_invcov.nc", InvCov, "InvCovariance");
    jif3D::WriteImpedancesToNetCDF(ncfilename + "_rot" + jif3D::stringify(dangle) + ".nc",
        Frequencies, StatX, StatY, StatZ, RotImp, Errors);

    jif3D::Write3DDataToVTK(ncfilename + "_rot.statpos.vtk", "Station",
        jif3D::rvec(StatX.size(), 1.0), StatX, StatY, StatZ);
  }
