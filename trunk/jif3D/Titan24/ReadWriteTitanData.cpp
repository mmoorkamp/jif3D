//============================================================================
// Name        : ReadWriteTitanData.cpp
// Author      : Feb 08, 2017
// Version     :
// Copyright   : 2017, aavdeeva
//============================================================================
#include "ReadWriteTitanData.h"
#include "../MT/ReadWriteImpedances.h"
#include "../MT/MTEquations.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/convert.h"
#include "../ModelBase/NetCDFModelTools.h"
#include "../Global/NetCDFTools.h"
#include "../Global/NetCDFPortHelper.h"

#include <fstream>
#include <iomanip>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <boost/cast.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using netCDF::NcVar;
using netCDF::NcFile;
using netCDF::NcDim;

namespace jif3D
  {
    using namespace std;

    //! The name used for the index of Stations with Ex measurement used to compute Titan TF data in netcdf files
    static const std::string ExIndicesName = "ExIndices";
    //! The name used for the index of Stations with Ey measurement used to compute Titan TF data in netcdf files
    static const std::string EyIndicesName = "EyIndices";
    //! The name used for the index of Stations with H measurement used to compute Titan TF data in netcdf files
    static const std::string HIndicesName = "HIndices";

    static const std::string FreqDimName = "Frequency";
    static const std::string DistortionName = "C";
    static const std::string AngleName = "RotationAngle";

    //write the Titan TF indices matrix to a netcdf file
    //this is an internal helper function
    void WriteTitanTFIndices(NcFile &NetCDFFile, NcDim &TitanTFNumDim, NcDim &FreqDim,
        const std::vector<int> &Indices, const std::string &IndexName)
      {
        std::vector<NcDim> dimVec;
        dimVec.push_back(FreqDim);
        dimVec.push_back(TitanTFNumDim);

        NcVar IndexVar = NetCDFFile.addVar(IndexName, netCDF::ncInt, dimVec);

        //        CompVar.put(&Component[0], FreqDim.getSize(), StatNumDim.getSize());
        cxxport::put_legacy_ncvar(IndexVar, &Indices[0], FreqDim.getSize(),
            TitanTFNumDim.getSize());
      }

    //read the Titan TF indices from a netcdf file
    //indices are matrices of size TitanTFNumDim x FreqDim
    //this is an internal helper function
    void ReadTitanTFIndices(NcFile &NetCDFFile, std::vector<int> &Indices,
        const std::string &IndexName, const bool MustExist = true)
      {
        try
          {
            NcVar SizeVar = NetCDFFile.getVar(IndexName);
            if (!SizeVar.isNull())
              {

                const std::vector<long> edges = cxxport::get_legacy_var_edges(SizeVar);
                Indices.resize(edges[0] * edges[1]);
//              SizeVar.get(&Temp[0], SizeVar->edges()[0], SizeVar->edges()[1]);
                cxxport::get_legacy_ncvar(SizeVar, &Indices[0], edges[0], edges[1]);

              }
          } catch (const netCDF::exceptions::NcException &ex)
          {
            if (MustExist)
              {
                throw std::runtime_error(
                    "Call to ReadTitanTFIndices with MustExist failed.");
              }
          }
      }

    void WriteTitanDataToNetCDF(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &MeasXCoord,
        const std::vector<double> &MeasYCoord, const std::vector<double> &MeasZCoord,
        const std::vector<int> &ExIndices, const std::vector<int> &EyIndices,
        const std::vector<int> &HIndices, const jif3D::rvec &Impedances,
        const jif3D::rvec &Errors, const std::vector<double> &Distortion,
        const std::vector<double> &RotAngles)
      {
        const size_t nmeas = MeasXCoord.size();
        const size_t nTitanTFs = ExIndices.size();
        const size_t nfreqs = Frequencies.size();
        const size_t nstats = nTitanTFs / nfreqs;

        const size_t nimp = nTitanTFs * 8;
        assert(nmeas == MeasYCoord.size());
        assert(nmeas == MeasZCoord.size());
        assert(nTitanTFs == EyIndices.size());
        assert(nTitanTFs == HIndices.size());
        assert(Impedances.size() == nimp);
        //create a netcdf file
        NcFile DataFile(filename, NcFile::replace);
        //Create the dimensions for the measurements
        NcDim MeasNumDim = DataFile.addDim(MeasNumberName, nmeas);

        //write out the measurement coordinates
        WriteVec(DataFile, MeasPosXName, MeasXCoord, MeasNumDim, "m");
        WriteVec(DataFile, MeasPosYName, MeasYCoord, MeasNumDim, "m");
        WriteVec(DataFile, MeasPosZName, MeasZCoord, MeasNumDim, "m");

        //write out the frequencies that we store
        NcDim FreqDim = DataFile.addDim(FreqDimName, Frequencies.size());
        NcVar FreqVar = DataFile.addVar(FreqDimName, netCDF::ncDouble, FreqDim);

//        FreqVar->put(&Frequencies[0], nfreqs);
        cxxport::put_legacy_ncvar(FreqVar, Frequencies.data(), nfreqs);

        //Create the dimensions for the stations
        NcDim StatNumDim = DataFile.addDim(StationNumberName, nstats);

        // write indices for Titan TF

        WriteTitanTFIndices(DataFile, StatNumDim, FreqDim, ExIndices, ExIndicesName);
        WriteTitanTFIndices(DataFile, StatNumDim, FreqDim, EyIndices, EyIndicesName);
        WriteTitanTFIndices(DataFile, StatNumDim, FreqDim, HIndices, HIndicesName);

        //and now we can write all the impedance components
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxx_re", 0);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxx_im", 1);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxy_re", 2);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxy_im", 3);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyx_re", 4);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyx_im", 5);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyy_re", 6);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyy_im", 7);
        //now we deal with the errors, if no parameter has been explicitly passed
        //Errors is empty, so we just fill the vector with zeros
        jif3D::rvec ZErr(Errors);
        if (ZErr.empty())
          {
            ZErr.resize(nimp);
            std::fill_n(ZErr.begin(), nimp, 0.0);
          }
        //in any case we write some error information even if it is all zero
        //as the error is identical for real and imaginary part
        //we only have to write the elements corresponding to the real part
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, ZErr, "dZxx", 0);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, ZErr, "dZxy", 2);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, ZErr, "dZyx", 4);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, ZErr, "dZyy", 6);
        if (!Distortion.empty())
          {
            NcDim CDim = DataFile.addDim("Celem", 4);

            std::vector<NcDim> dimVec;
            dimVec.push_back(StatNumDim);
            dimVec.push_back(CDim);

            NcVar CVar = DataFile.addVar(DistortionName, netCDF::ncDouble, dimVec);

//            CVar.put(&Distortion[0], nstats, 4);
            cxxport::put_legacy_ncvar(CVar, Distortion.data(), nstats, 4);
          }
        if (!RotAngles.empty())
          {
            std::vector<NcDim> dimVec;
            dimVec.push_back(StatNumDim);
            NcVar AngleVar = DataFile.addVar(AngleName, netCDF::ncDouble, dimVec);

//            CVar.put(&Distortion[0], nstats, 4);
            cxxport::put_legacy_ncvar(AngleVar, RotAngles.data(), nstats);
          }

      }

    void ReadTitanDataFromNetCDF(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &MeasXCoord,
        std::vector<double> &MeasYCoord, std::vector<double> &MeasZCoord,
        std::vector<int> &ExIndices, std::vector<int> &EyIndices,
        std::vector<int> &HIndices, jif3D::rvec &Impedances, jif3D::rvec &ImpError,
        std::vector<double> &Distortion, std::vector<double> &RotAngles)
      {
        //open the netcdf file readonly
        NcFile DataFile(filename, NcFile::read);
        //read in the station coordinates in the three directions
        ReadVec(DataFile, MeasPosXName, MeasXCoord);
        ReadVec(DataFile, MeasPosYName, MeasYCoord);
        ReadVec(DataFile, MeasPosZName, MeasZCoord);
        //read in the frequencies
        ReadVec(DataFile, FreqDimName, Frequencies);

        ReadTitanTFIndices(DataFile, ExIndices, ExIndicesName, false);
        ReadTitanTFIndices(DataFile, EyIndices, EyIndicesName, false);
        ReadTitanTFIndices(DataFile, HIndices, HIndicesName, false);

        //for each frequency and each data index we have 4 titan transfer function elements
        //currently it is not possible to have a different number of frequencies
        //at each site
        const size_t nimp = ExIndices.size() * 8;

        Impedances.resize(nimp);
        ImpError.resize(nimp);
        //read the impedances
        ReadImpedanceComp(DataFile, Impedances, "Zxx_re", 0);
        ReadImpedanceComp(DataFile, Impedances, "Zxx_im", 1);
        ReadImpedanceComp(DataFile, Impedances, "Zxy_re", 2);
        ReadImpedanceComp(DataFile, Impedances, "Zxy_im", 3);
        ReadImpedanceComp(DataFile, Impedances, "Zyx_re", 4);
        ReadImpedanceComp(DataFile, Impedances, "Zyx_im", 5);
        ReadImpedanceComp(DataFile, Impedances, "Zyy_re", 6);
        ReadImpedanceComp(DataFile, Impedances, "Zyy_im", 7);
        //now read the errors, we make their presence in the file optional
        ReadImpedanceComp(DataFile, ImpError, "dZxx", 0, false);
        ReadImpedanceComp(DataFile, ImpError, "dZxy", 2, false);
        ReadImpedanceComp(DataFile, ImpError, "dZyx", 4, false);
        ReadImpedanceComp(DataFile, ImpError, "dZyy", 6, false);
        //we only read in the errors for the real part
        //and then assign the same error to the imaginary part
        for (size_t i = 0; i < nimp - 1; i += 2)
          {
            ImpError(i + 1) = ImpError(i);
          }

        const int nstats = ExIndices.size() / Frequencies.size();
        // Distortion values are optional, so we have to catch exceptions thrown by the netcdf API
        //and ignore them
        try
          {
            NcVar DistVar = DataFile.getVar(DistortionName);
            const int nvalues = nstats * 4;
            Distortion.resize(nvalues);
            if (!DistVar.isNull())
              {
                const std::vector<long> edges = cxxport::get_legacy_var_edges(DistVar);
                if (nvalues != edges[0] * edges[1])
                  {
                    throw jif3D::FatalException(
                        "Number of distortion parameters does not match number of stations !",
                        __FILE__, __LINE__);
                  }

                //              DistVar.get(&Distortion[0], edges[0], edges[1]);
                cxxport::get_legacy_ncvar(DistVar, Distortion.data(), edges[0], edges[1]);
              }
            else
              {
                const size_t nstats = ExIndices.size() / Frequencies.size();
                for (size_t i = 0; i < nstats; ++i)
                  {
                    Distortion[i * 4] = 1.0;
                    Distortion[i * 4 + 1] = 0.0;
                    Distortion[i * 4 + 2] = 0.0;
                    Distortion[i * 4 + 3] = 1.0;
                  }
              }
          } catch (const netCDF::exceptions::NcException &ex)
          {
            // ignore
          }

        // rotation angles are optional, so we have to catch exceptions thrown by the netcdf API
        //and ignore them
        try
          {
            NcVar AngleVar = DataFile.getVar(AngleName);
            RotAngles.resize(nstats);
            if (!AngleVar.isNull())
              {
                const std::vector<long> edges = cxxport::get_legacy_var_edges(AngleVar);
                if (nstats != edges[0])
                  {
                    throw jif3D::FatalException(
                        "Number of rotation angle values does not match number of stations !",
                        __FILE__, __LINE__);
                  }
                cxxport::get_legacy_ncvar(AngleVar, RotAngles.data(), edges[0]);
              }
            else
              {
                std::fill(RotAngles.begin(), RotAngles.end(), 0.0);
              }

          } catch (const netCDF::exceptions::NcException &ex)
          {
            // ignore
          }

      }

    void WriteTitanDataToModEM(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<int> &ExIndices, const jif3D::rvec &Imp, const jif3D::rvec &Err)
      {
        const double convfactor = 4.0 * 1e-4 * acos(-1.0);
        std::ofstream outfile(filename.c_str());
        outfile << "# Description: \n";
        outfile
            << "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error \n";
        outfile << "> Full_Impedance\n";
        outfile << "> exp(+i\\omega t)\n ";
        outfile << "> [mV/km]/[nT]\n";
        outfile << "> 0.00\n";
        outfile << "> 0.0 0.0\n";

        size_t nfreqs = Frequencies.size();
        size_t nsites = ExIndices.size() / nfreqs;
        outfile << "> " << nfreqs << " " << nsites << "\n";
        for (size_t i = 0; i < nsites; ++i)
          {
            std::string SiteName = "Site" + std::to_string(i);
            std::ostringstream SiteLine;
            SiteLine << " " << SiteName << "  0.0  0.0 " << StatXCoord[ExIndices[i]]
                << " " << StatYCoord[ExIndices[i]] << " " << StatZCoord[ExIndices[i]]
                << " ";
            for (size_t j = 0; j < nfreqs; ++j)
              {
                size_t index = 8 * (nsites * j + i);
                double period = 1.0 / Frequencies.at(j);
                outfile << period << SiteLine.str();
                outfile << " ZXX " << Imp(index) / convfactor << " "
                    << Imp(index + 1) / convfactor << " " << Err(index) / convfactor
                    << "\n";

                outfile << period << SiteLine.str();
                outfile << " ZXY " << Imp(index + 2) / convfactor << " "
                    << Imp(index + 3) / convfactor << " " << Err(index + 2) / convfactor
                    << "\n";

                outfile << period << SiteLine.str();
                outfile << " ZYX " << Imp(index + 4) / convfactor << " "
                    << Imp(index + 5) / convfactor << " " << Err(index + 4) / convfactor
                    << "\n";

                outfile << period << SiteLine.str();
                outfile << " ZYY " << Imp(index + 6) / convfactor << " "
                    << Imp(index + 6) / convfactor << " " << Err(index + 6) / convfactor
                    << "\n";
              }
          }
      }
  }

