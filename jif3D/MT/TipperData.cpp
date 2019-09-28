/*
 * TipperData.cpp
 *
 *  Created on: May 22, 2019
 *      Author: max
 */

#include "TipperData.h"
#include "ReadWriteImpedances.h"
#include "../ModelBase/VTKTools.h"


namespace jif3D
  {

    TipperData::TipperData()
      {
        // TODO Auto-generated constructor stub

      }

    TipperData::~TipperData()
      {
        // TODO Auto-generated destructor stub
      }

    //! Write all information to a netcdf file
    void TipperData::ReadNetCDF(const std::string &filename)
      {
        std::vector<double> PosX, PosY, PosZ, D, E, Freqs;
        ReadTipperFromNetCDF(filename, Freqs, PosX, PosY, PosZ, HxIndices, HyIndices,
            HzIndices, D, E);
        SetFrequencies(Freqs);

        SetMeasurementPoints(PosX, PosY, PosZ);
        SetDataAndErrors(D, E);
      }
    void TipperData::WriteNetCDF(const std::string &filename)
      {
        WriteTipperToNetCDF(filename, GetFrequencies(), GetMeasPosX(), GetMeasPosY(),
            GetMeasPosZ(), HxIndices, HyIndices, HzIndices, GetData(), GetErrors());
      }

    void TipperData::CompleteObject()
      {
        const size_t nmeas = GetMeasPosX().size();
        if (nmeas == 0)
          {
            throw jif3D::FatalException(
                "No measurements specified, cannot complete MT Object", __FILE__,
                __LINE__);
          }
        const size_t nfreq = GetFrequencies().size();
        if (nfreq == 0)
          {
            throw jif3D::FatalException(
                "No frequencies specified, cannot complete MT Object", __FILE__,
                __LINE__);
          }

        if (HxIndices.size() != nmeas * nfreq)
          {
            HxIndices.resize(nmeas * nfreq);

            for (size_t ifr = 0; ifr < nfreq; ++ifr)
              {
                const size_t ind_shift = nmeas * ifr;
                for (size_t i = 0; i < nmeas; ++i)
                  {
                    HxIndices[i + ind_shift] = i;
                  }
              }
            HyIndices = HxIndices;
            HzIndices = HxIndices;
          }
        if (GetRotAngles().size() != nmeas)
          {
            std::vector<double> RA(nmeas, 0.0);
            SetRotAngles(RA);
          }
        if (GetData().size() != nmeas)
          {
            std::vector<double> dummy(nmeas, 0.0);
            SetDataAndErrors(dummy, dummy);
          }
      }

    void TipperData::ReadModEM(const std::string &filename)
      {

      }
    void TipperData::WriteModEM(const std::string &filename)
      {

      }
    void TipperData::PlotMeasurementConfiguration(const std::string &filename)
      {
        std::vector<double> FirstFreq(GetMeasPosX().size());

        const size_t nstats = HxIndices.size() / GetFrequencies().size();
        std::vector<double> tmpx, tmpy, tmpz;

        for (unsigned int i = 0; i < nstats; ++i)
          {
            tmpx.push_back(GetMeasPosX()[HxIndices[i]]);
            tmpy.push_back(GetMeasPosY()[HxIndices[i]]);
            tmpz.push_back(GetMeasPosZ()[HxIndices[i]]);
          }

        FirstFreq.resize(tmpx.size());
        std::iota(FirstFreq.begin(), FirstFreq.end(), 1);
        jif3D::Write3DDataToVTK(filename + "_Hx.vtk", "HxSites", FirstFreq, tmpx, tmpy,
            tmpz);

        tmpx.clear();
        tmpy.clear();
        tmpz.clear();

        for (unsigned int i = 0; i < nstats; ++i)
          {
            tmpx.push_back(GetMeasPosX()[HyIndices[i]]);
            tmpy.push_back(GetMeasPosY()[HyIndices[i]]);
            tmpz.push_back(GetMeasPosZ()[HyIndices[i]]);
          }
        jif3D::Write3DDataToVTK(filename + "_Hy.vtk", "HySites", FirstFreq, tmpx, tmpy,
            tmpz);

        tmpx.clear();
        tmpy.clear();
        tmpz.clear();

        for (unsigned int i = 0; i < nstats; ++i)
          {
            tmpx.push_back(GetMeasPosX()[HzIndices[i]]);
            tmpy.push_back(GetMeasPosY()[HzIndices[i]]);
            tmpz.push_back(GetMeasPosZ()[HzIndices[i]]);
          }
        jif3D::Write3DDataToVTK(filename + "_Hz.vtk", "HzSites", FirstFreq, tmpx, tmpy,
            tmpz);
      }

  } /* namespace jif3D */
