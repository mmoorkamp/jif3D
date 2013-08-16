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
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../Global/VecMat.h"
#include "../Global/ReadWriteSparseMatrix.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTEquations.h"
#include "../MT/MTCovar.h"
#include "../Inversion/MatrixTools.h"
#include "../ModelBase/VTKTools.h"

int main()
  {
    std::string ncfilename = jif3D::AskFilename("Name of netcdf file: ");

    jif3D::rvec Impedances, Errors;
    std::vector<double> Frequencies, StatX, StatY, StatZ;
    jif3D::ReadImpedancesFromNetCDF(ncfilename, Frequencies, StatX, StatY, StatZ,
        Impedances, Errors);

    std::cout << "Rotation angle [degree]: ";
    double dangle = 0.0;
    std::cin >> dangle;
    double rangle = dangle / 180.0 * M_PI;
    jif3D::rvec RotImp = jif3D::RotateImpedanceVector(rangle, Impedances);
    const size_t nelem = 8;
    const size_t ntensor = Impedances.size() / nelem;
    jif3D::comp_mat InvCov(Errors.size(), Errors.size());
    for (size_t i = 0; i < ntensor; ++i)
      {
        jif3D::rmat InvCovOrig(4, 4);
        InvCovOrig(0, 0) = pow(1.0 / Errors(i * nelem), 2);
        InvCovOrig(1, 1) = pow(1.0 / Errors(i * nelem + 2), 2);
        InvCovOrig(2, 2) = pow(1.0 / Errors(i * nelem + 4), 2);
        InvCovOrig(3, 3) = pow(1.0 / Errors(i * nelem + 6), 2);
        jif3D::rmat RotInvCov = jif3D::RotateMTInvCovar(rangle, InvCovOrig);

        for (size_t j = 0; j < 4; ++j)
          {
            for (size_t k = 0; k < 4; ++k)
              {
                InvCov(i * nelem + j * 2, i * nelem + k * 2) = RotInvCov(j, k);
                InvCov(i * nelem + j * 2 + 1, i * nelem + k * 2 + 1) = RotInvCov(j, k);
              }
          }
      }

    double MeanX = std::accumulate(StatX.begin(), StatX.end(), 0.0) / StatX.size();
    double MeanY = std::accumulate(StatY.begin(), StatY.end(), 0.0) / StatY.size();
    for (size_t i = 0; i < StatX.size(); ++i)
      {
        double newx = MeanX + (StatX.at(i) - MeanX) * cos(rangle)
            - (StatY.at(i) - MeanY) * sin(rangle);
        double newy = MeanY + (StatX.at(i) - MeanX) * sin(rangle)
            + (StatY.at(i) - MeanY) * cos(rangle);
        StatX.at(i) = newx;
        StatY.at(i) = newy;
      }

    jif3D::WriteSparseMatrixToNetcdf(ncfilename + "_invcov.nc", InvCov, "InvCovariance");
    jif3D::WriteImpedancesToNetCDF(ncfilename + "_rot.nc", Frequencies, StatX, StatY,
        StatZ, RotImp, Errors);

    jif3D::Write3DDataToVTK(ncfilename + "_rot.statpos.vtk", "Station",
        jif3D::rvec(ntensor, 1.0), StatX, StatY, StatZ);
  }
