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
        jif3D::covmat CovOrig(4, 4);
        CovOrig(0, 0) = Errors(i * nelem);
        CovOrig(1, 1) = Errors(i * nelem + 2);
        CovOrig(2, 2) = Errors(i * nelem + 4);
        CovOrig(3, 3) = Errors(i * nelem + 6);
        jif3D::covmat RotCov = jif3D::RotateMTCovar(rangle, CovOrig);
        jif3D::rmat InvCovLocal(RotCov.size1(), RotCov.size2());
        jif3D::InvertMatrix(RotCov, InvCovLocal);
        for (size_t j = 0; j < 4; ++j)
          {
            for (size_t k = 0; k < 4; ++k)
              {
                InvCov(i * nelem + j * 2, i * nelem + k * 2) = InvCovLocal(j, k);
                InvCov(i * nelem + j * 2 + 1, i * nelem + k * 2 + 1) = InvCovLocal(j, k);
              }
          }

      }
    jif3D::WriteSparseMatrixToNetcdf(ncfilename + "_invcov.nc", InvCov, "InvCovariance");
    jif3D::WriteImpedancesToNetCDF(ncfilename + "_rot.nc", Frequencies, StatX, StatY,
        StatZ, RotImp, Errors);

  }
