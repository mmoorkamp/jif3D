//============================================================================
// Name        : X3DMTCalculator.cpp
// Author      : Jul 7, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <complex>
#include "../Global/convert.h"
#include "X3DMTCalculator.h"
#include "ReadWriteX3D.h"
#include "MTEquations.h"

namespace jiba
  {

    X3DMTCalculator::X3DMTCalculator()
      {

      }

    X3DMTCalculator::~X3DMTCalculator()
      {

      }

    rvec X3DMTCalculator::Calculate(const X3DModel &Model)
      {
        const std::string modelfilename("x3d.model");
        const std::string resultfilename("x3d.result");
        WriteProjectFile(Model.GetFrequencies(), X3DModel::MT, resultfilename,
            modelfilename);
        Write3DModelForX3D(modelfilename, Model.GetXCellSizes(),
            Model.GetYCellSizes(), Model.GetZCellSizes(),
            Model.GetConductivities(), Model.GetBackgroundConductivities(),
            Model.GetBackgroundThicknesses());
        system("x3d");
        jiba::rvec result(Model.GetXCellSizes().size()
            * Model.GetYCellSizes().size() * Model.GetFrequencies().size() * 8);
        result.clear();
        const size_t nfreq = Model.GetFrequencies().size();
        for (size_t i = 0; i < nfreq; ++i)
          {
            std::string emoAname = resultfilename + jiba::stringify(i)
                + "a.emo";
            std::string emoBname = resultfilename + jiba::stringify(i)
                + "b.emo";
            std::vector<std::complex<double> > Ex1, Ex2, Ey1, Ey2, Hx1, Hx2,
                Hy1, Hy2;
            std::complex<double> Zxx, Zxy, Zyx, Zyy;
            ReadEMO(emoAname, Ex1, Ey1, Hx1, Hy1);
            ReadEMO(emoBname, Ex2, Ey2, Hx2, Hy2);
            const size_t nobs = Ex1.size();
            const size_t freq_index = nobs * i;
            for (size_t j = 0; j < nobs; ++j)
              {
                const size_t obs_index = freq_index + j * 8;
                FieldsToImpedance(Ex1[i], Ex2[i], Ey1[i], Ey2[i], Hx1[i],
                    Hx2[i], Hy1[i], Hy2[i], Zxx, Zxy, Zyx, Zyy);
                result(obs_index) = Zxx.real();
                result(obs_index + 1) = Zxx.imag();
                result(obs_index + 2) = Zxy.real();
                result(obs_index + 3) = Zxy.imag();
                result(obs_index + 4) = Zyx.real();
                result(obs_index + 5) = Zyx.imag();
                result(obs_index + 6) = Zyy.real();
                result(obs_index + 7) = Zyy.imag();
              }
            //finished with one frequency
          }
        return result;
      }

    rvec X3DMTCalculator::LQDerivative(const X3DModel &Model,
        const rvec &Misfit)
      {

      }
  }
