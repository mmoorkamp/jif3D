/*
 * TipperData.cpp
 *
 *  Created on: May 22, 2019
 *      Author: max
 */

#include "TipperData.h"
#include "ReadWriteImpedances.h"
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
            GetMeasPosZ(), GetData(), GetErrors());
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
        const size_t ndist = nmeas * 4;

        if (HxIndices.size() != nmeas * nfreq)
          {
            HxIndices.resize(nmeas * nfreq);

            for (int ifr = 0; ifr < nfreq; ++ifr)
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
      }

    void TipperData::ReadModEM(const std::string &filename)
      {

      }
    void TipperData::WriteModEM(const std::string &filename)
      {

      }
    void TipperData::PlotMeasurementConfiguration(const std::string &filename)
      {

      }

  } /* namespace jif3D */
