/*
 * MTData.cpp
 *
 *  Created on: May 21, 2019
 *      Author: max
 */

#include "MTData.h"
#include "ReadWriteTitanData.h"
#include "../ModelBase/VTKTools.h"
namespace jif3D
  {

    //! Write all information to a netcdf file
    void MTData::ReadNetCDF(const std::string &filename)
      {
        std::vector<double> PosX, PosY, PosZ, Freqs, Imp, Err, RA;
        ReadTitanDataFromNetCDF(filename, Freqs, PosX, PosY, PosZ, ExIndices, EyIndices,
            HxIndices, Imp, Err, Distortion, RA);
        HyIndices = HxIndices;
        SetFrequencies(Freqs);
        SetRotAngles(RA);
        SetMeasurementPoints(PosX, PosY, PosZ);
        SetDataAndErrors(Imp, Err);
      }
    void MTData::WriteNetCDF(const std::string &filename)
      {
        WriteTitanDataToNetCDF(filename, GetFrequencies(), GetMeasPosX(), GetMeasPosY(),
            GetMeasPosZ(), ExIndices, EyIndices, HxIndices, GetData(), GetErrors(),
            Distortion, GetRotAngles());
      }

    void MTData::SetDefaultDistortion()
      {
        const size_t nstats =
            ExIndices.empty() ?
                GetMeasPosX().size() : ExIndices.size() / GetFrequencies().size();
        Distortion.resize(nstats *4);
        for (size_t i = 0; i < nstats; ++i)
          {
            Distortion[i * 4] = 1.0;
            Distortion[i * 4 + 1] = 0.0;
            Distortion[i * 4 + 2] = 0.0;
            Distortion[i * 4 + 3] = 1.0;
          }
      }

    void MTData::CompleteObject()
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
        const size_t ndist = nmeas * 4;
        if (GetDistortion().size() != ndist)
          {
            SetDefaultDistortion();
          }
        if (ExIndices.size() != nmeas * nfreq)
          {
            ExIndices.resize(nmeas * nfreq);
            for (size_t ifr = 0; ifr < nfreq; ++ifr)
              {
                const size_t ind_shift = nmeas * ifr;
                for (size_t i = 0; i < nmeas; ++i)
                  {
                    ExIndices[i + ind_shift] = i;
                  }
              }
            EyIndices = ExIndices;
            HxIndices = ExIndices;
            HyIndices = ExIndices;
          }
        if (GetRotAngles().size() != nmeas)
          {
            std::vector<double> RA(nmeas, 0.0);
            SetRotAngles(RA);
          }
        if (GetData().size() != nmeas)
          {
            std::vector<double> dummy(nmeas,0.0);
            SetDataAndErrors(dummy,dummy);
          }
      }

    void MTData::ReadModEM(const std::string &filename)
      {

      }

    void MTData::WriteModEM(const std::string &filename)
      {

      }

    void MTData::PlotMeasurementConfiguration(const std::string &filename)
      {
        std::vector<double> FirstFreq(GetMeasPosX().size());

        const size_t nstats = ExIndices.size() / GetFrequencies().size();
        std::vector<double> tmpx, tmpy, tmpz;

        for (unsigned int i = 0; i < nstats; ++i)
          {
            tmpx.push_back(GetMeasPosX()[ExIndices[i]]);
            tmpy.push_back(GetMeasPosY()[ExIndices[i]]);
            tmpz.push_back(GetMeasPosZ()[ExIndices[i]]);
          }

        FirstFreq.resize(tmpx.size());
        std::iota(FirstFreq.begin(), FirstFreq.end(), 1);
        jif3D::Write3DDataToVTK(filename + "_Ex.vtk", "ExSites", FirstFreq, tmpx, tmpy,
            tmpz);

        tmpx.clear();
        tmpy.clear();
        tmpz.clear();

        for (unsigned int i = 0; i < nstats; ++i)
          {
            tmpx.push_back(GetMeasPosX()[EyIndices[i]]);
            tmpy.push_back(GetMeasPosY()[EyIndices[i]]);
            tmpz.push_back(GetMeasPosZ()[EyIndices[i]]);
          }

        FirstFreq.resize(tmpx.size());
        std::iota(FirstFreq.begin(), FirstFreq.end(), 1);
        jif3D::Write3DDataToVTK(filename + "_Ey.vtk", "EySites", FirstFreq, tmpx, tmpy,
            tmpz);

        tmpx.clear();
        tmpy.clear();
        tmpz.clear();

        for (unsigned int i = 0; i < nstats; ++i)
          {
            tmpx.push_back(GetMeasPosX()[HxIndices[i]]);
            tmpy.push_back(GetMeasPosY()[HxIndices[i]]);
            tmpz.push_back(GetMeasPosZ()[HxIndices[i]]);
          }

        FirstFreq.resize(tmpx.size());
        std::iota(FirstFreq.begin(), FirstFreq.end(), 1);
        jif3D::Write3DDataToVTK(filename + "_Hx.vtk", "HSites", FirstFreq, tmpx, tmpy,
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

        FirstFreq.resize(tmpx.size());
        std::iota(FirstFreq.begin(), FirstFreq.end(), 1);
        jif3D::Write3DDataToVTK(filename + "_Hy.vtk", "HSites", FirstFreq, tmpx, tmpy,
            tmpz);

      }

    MTData::MTData()
      {
        // TODO Auto-generated constructor stub

      }

    MTData::~MTData()
      {
        // TODO Auto-generated destructor stub
      }

  } /* namespace jif3D */
