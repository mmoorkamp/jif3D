//============================================================================
// Name        : ReadWriteImpedance.cpp
// Author      : Jul 13, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================
#include "ReadWriteImpedances.h"
#include "MTEquations.h"
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../Global/convert.h"
#include "../Global/FileUtil.h"
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
    static const std::string FreqDimName = "Frequency";
    static const std::string DistortionName = "C";

    //! The name used for the index of Stations with Ex measurement used to compute Titan TF data in netcdf files
    static const std::string HxIndicesName = "HxIndices";
    //! The name used for the index of Stations with Ey measurement used to compute Titan TF data in netcdf files
    static const std::string HyIndicesName = "HyIndices";
    //! The name used for the index of Stations with H measurement used to compute Titan TF data in netcdf files
    static const std::string HzIndicesName = "HzIndices";
//write one compoment of the impedance tensor to a netcdf file
//this is an internal helper function
    void WriteImpedanceComp(NcFile &NetCDFFile, NcDim &StatNumDim, NcDim &FreqDim,
        const std::vector<double> &Impedances, const std::string &CompName,
        const size_t compindex)
      {
        std::vector<NcDim> dimVec;
        dimVec.push_back(FreqDim);
        dimVec.push_back(StatNumDim);

        NcVar CompVar = NetCDFFile.addVar(CompName, netCDF::ncDouble, dimVec);

        std::vector<double> Component(FreqDim.getSize() * StatNumDim.getSize());

        for (size_t i = 0; i < Component.size(); ++i)
          {
            Component.at(i) = Impedances.at(i * 8 + compindex);
          }

//        CompVar.put(&Component[0], FreqDim.getSize(), StatNumDim.getSize());
        cxxport::put_legacy_ncvar(CompVar, &Component[0], FreqDim.getSize(),
            StatNumDim.getSize());
        //we write the impedance in units of Ohm
        CompVar.putAtt("units", "Ohm");
      }

    void WriteTipperComp(NcFile &NetCDFFile, NcDim &StatNumDim, NcDim &FreqDim,
        const std::vector<double> &Tipper, const std::string &CompName,
        const size_t compindex)
      {
        std::vector<NcDim> dimVec;
        dimVec.push_back(FreqDim);
        dimVec.push_back(StatNumDim);

        NcVar CompVar = NetCDFFile.addVar(CompName, netCDF::ncDouble, dimVec);

        jif3D::rvec Component(FreqDim.getSize() * StatNumDim.getSize());

        for (size_t i = 0; i < Component.size(); ++i)
          {
            Component(i) = Tipper.at(i * 4 + compindex);
          }

//        CompVar.put(&Component[0], FreqDim.getSize(), StatNumDim.getSize());
        cxxport::put_legacy_ncvar(CompVar, &Component[0], FreqDim.getSize(),
            StatNumDim.getSize());
        //we write the impedance in units of Ohm
        CompVar.putAtt("units", "");
      }

    void ReadTFIndices(NcFile &NetCDFFile, std::vector<int> &Indices,
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

    //write the Titan TF indices matrix to a netcdf file
    //this is an internal helper function
    void WriteTFIndices(NcFile &NetCDFFile, NcDim &TitanTFNumDim, NcDim &FreqDim,
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

//read one component of the impedance tensor from a netcdf file
//this is an internal helper function
    void ReadImpedanceComp(NcFile &NetCDFFile, std::vector<double> &Impedances,
        const std::string &CompName, const size_t compindex, const bool MustExist)
      {
        try
          {
            NcVar SizeVar = NetCDFFile.getVar(CompName);
            if (!SizeVar.isNull())
              {
                const size_t nvalues = Impedances.size() / 8;
                const std::vector<long> edges = cxxport::get_legacy_var_edges(SizeVar);

                assert(nvalues == boost::numeric_cast<size_t>(edges[0] * edges[1]));
                std::vector<double> Temp(nvalues);

//              SizeVar.get(&Temp[0], SizeVar->edges()[0], SizeVar->edges()[1]);
                cxxport::get_legacy_ncvar(SizeVar, &Temp[0], edges[0], edges[1]);

                for (size_t i = 0; i < nvalues; ++i)
                  {
                    Impedances.at(i * 8 + compindex) = Temp.at(i);
                  }
              }
          } catch (const netCDF::exceptions::NcException &ex)
          {
            if (MustExist)
              {
                throw std::runtime_error(
                    "Call to ReadImpedanceComp with MustExist failed.");
              }
          }
      }

    void ReadTipperComp(NcFile &NetCDFFile, std::vector<double> &Tipper,
        const std::string &CompName, const size_t compindex, const bool MustExist = true)
      {
        try
          {
            NcVar SizeVar = NetCDFFile.getVar(CompName);
            if (!SizeVar.isNull())
              {
                const size_t nvalues = Tipper.size() / 4;
                const std::vector<long> edges = cxxport::get_legacy_var_edges(SizeVar);

                assert(nvalues == boost::numeric_cast<size_t>(edges[0] * edges[1]));
                jif3D::rvec Temp(nvalues);

//              SizeVar.get(&Temp[0], SizeVar->edges()[0], SizeVar->edges()[1]);
                cxxport::get_legacy_ncvar(SizeVar, &Temp[0], edges[0], edges[1]);

                for (size_t i = 0; i < nvalues; ++i)
                  {
                    Tipper.at(i * 4 + compindex) = Temp(i);
                  }
              }
          } catch (const netCDF::exceptions::NcException &ex)
          {
            if (MustExist)
              {
                throw std::runtime_error(
                    "Call to ReadImpedanceComp with MustExist failed.");
              }
          }
      }

    void WriteImpedancesToNetCDF(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<double> &Impedances, const std::vector<double> &Errors,
        const std::vector<double> &Distortion)
      {
        const size_t nstats = StatXCoord.size();
        const size_t nfreqs = Frequencies.size();
        const size_t nimp = nstats * nfreqs * 8;
        assert(nstats == StatYCoord.size());
        assert(nstats == StatYCoord.size());
        assert(Impedances.size() == nimp);
        //create a netcdf file
        NcFile DataFile(filename, NcFile::replace);
        //Create the dimensions for the stations
        NcDim StatNumDim = DataFile.addDim(StationNumberName, nstats);

        //write out the measurement coordinates
        WriteVec(DataFile, MeasPosXName, StatXCoord, StatNumDim, "m");
        WriteVec(DataFile, MeasPosYName, StatYCoord, StatNumDim, "m");
        WriteVec(DataFile, MeasPosZName, StatZCoord, StatNumDim, "m");
        //write out the frequencies that we store
        NcDim FreqDim = DataFile.addDim(FreqDimName, Frequencies.size());
        NcVar FreqVar = DataFile.addVar(FreqDimName, netCDF::ncDouble, FreqDim);

//        FreqVar->put(&Frequencies[0], nfreqs);
        cxxport::put_legacy_ncvar(FreqVar, Frequencies.data(), nfreqs);
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
        std::vector<double> ZErr(Errors);
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
      }

    void ReadImpedancesFromNetCDF(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        std::vector<double> &Impedances, std::vector<double> &ImpError,
        std::vector<double> &Distortion)
      {
        //open the netcdf file readonly
        NcFile DataFile(filename, NcFile::read);
        //read in the station coordinates in the three directions
        ReadVec(DataFile, MeasPosXName, StatXCoord);
        ReadVec(DataFile, MeasPosYName, StatYCoord);
        ReadVec(DataFile, MeasPosZName, StatZCoord);
        //read in the frequencies
        ReadVec(DataFile, FreqDimName, Frequencies);
        //for each frequency and each station we have 8 impedances
        //currently it is not possible to have a different number of frequencies
        //at each site
        const size_t nimp = Frequencies.size() * StatXCoord.size() * 8;
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
            ImpError.at(i + 1) = ImpError.at(i);
          }

        try
          {
            NcVar DistVar = DataFile.getVar(DistortionName);

            const int nvalues = StatXCoord.size() * 4;
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
                for (size_t i = 0; i < StatXCoord.size(); ++i)
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
      }

    void WriteTipperToNetCDF(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<int> &HxIndices, const std::vector<int> &HyIndices,
        const std::vector<int> &HzIndices, const std::vector<double> &Tipper,
        const std::vector<double> &Errors)
      {
        const size_t nstats = StatXCoord.size();
        const size_t nfreqs = Frequencies.size();
        const size_t nimp = nstats * nfreqs * 4;
        assert(nstats == StatYCoord.size());
        assert(nstats == StatYCoord.size());
        assert(Tipper.size() == nimp);
        //create a netcdf file
        NcFile DataFile(filename, NcFile::replace);
        //Create the dimensions for the stations
        NcDim StatNumDim = DataFile.addDim(StationNumberName, nstats);

        //write out the measurement coordinates
        WriteVec(DataFile, MeasPosXName, StatXCoord, StatNumDim, "m");
        WriteVec(DataFile, MeasPosYName, StatYCoord, StatNumDim, "m");
        WriteVec(DataFile, MeasPosZName, StatZCoord, StatNumDim, "m");

        //write out the frequencies that we store
        NcDim FreqDim = DataFile.addDim(FreqDimName, Frequencies.size());
        NcVar FreqVar = DataFile.addVar(FreqDimName, netCDF::ncDouble, FreqDim);

        //        FreqVar->put(&Frequencies[0], nfreqs);
        cxxport::put_legacy_ncvar(FreqVar, Frequencies.data(), nfreqs);
        //and now we can write all the impedance components

        WriteTFIndices(DataFile, StatNumDim, FreqDim, HxIndices, HxIndicesName);
        WriteTFIndices(DataFile, StatNumDim, FreqDim, HyIndices, HyIndicesName);
        WriteTFIndices(DataFile, StatNumDim, FreqDim, HzIndices, HzIndicesName);
        WriteTipperComp(DataFile, StatNumDim, FreqDim, Tipper, "Tx_re", 0);
        WriteTipperComp(DataFile, StatNumDim, FreqDim, Tipper, "Tx_im", 1);
        WriteTipperComp(DataFile, StatNumDim, FreqDim, Tipper, "Ty_re", 2);
        WriteTipperComp(DataFile, StatNumDim, FreqDim, Tipper, "Ty_im", 3);

        //in any case we write some error information even if it is all zero
        //as the error is identical for real and imaginary part
        //we only have to write the elements corresponding to the real part
        WriteTipperComp(DataFile, StatNumDim, FreqDim, Errors, "dTx", 0);
        WriteTipperComp(DataFile, StatNumDim, FreqDim, Errors, "dTy", 2);

      }

    void ReadTipperFromNetCDF(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        std::vector<int> &HxIndices, std::vector<int> &HyIndices,
        std::vector<int> &HzIndices, std::vector<double> &Tipper,
        std::vector<double> &Error)
      {
        //open the netcdf file readonly
        NcFile DataFile(filename, NcFile::read);
        //read in the station coordinates in the three directions
        ReadVec(DataFile, MeasPosXName, StatXCoord);
        ReadVec(DataFile, MeasPosYName, StatYCoord);
        ReadVec(DataFile, MeasPosZName, StatZCoord);
        //read in the frequencies
        ReadVec(DataFile, FreqDimName, Frequencies);
        //for each frequency and each station we have 4 tipper values
        //currently it is not possible to have a different number of frequencies
        //at each site
        ReadTFIndices(DataFile, HxIndices, HxIndicesName, false);
        ReadTFIndices(DataFile, HyIndices, HyIndicesName, false);
        ReadTFIndices(DataFile, HzIndices, HzIndicesName, false);

        if (HyIndices.size() != HxIndices.size() || HzIndices.size() != HxIndices.size())
          {
            throw jif3D::FatalException("Inconsistent index information for Tipper Data",
            __FILE__, __LINE__);
          }

        const size_t nmeas = StatXCoord.size();
        const size_t nfreq = Frequencies.size();
        size_t ind_shift = 0;
        if (HxIndices.empty())
          {
            HxIndices.resize(nmeas * nfreq);
            HyIndices.resize(nmeas * nfreq);
            HzIndices.resize(nmeas * nfreq);
            for (size_t ifr = 0; ifr < nfreq; ++ifr)
              {
                ind_shift = nmeas * ifr;
                for (size_t i = 0; i < nmeas; ++i)
                  {
                    HxIndices[i + ind_shift] = i;
                  }
              }
            HyIndices = HxIndices;
            HzIndices = HxIndices;
          }

        const size_t nimp = Frequencies.size() * StatXCoord.size() * 4;
        Tipper.resize(nimp);
        Error.resize(nimp);
        //read the impedances
        ReadTipperComp(DataFile, Tipper, "Tx_re", 0);
        ReadTipperComp(DataFile, Tipper, "Tx_im", 1);
        ReadTipperComp(DataFile, Tipper, "Ty_re", 2);
        ReadTipperComp(DataFile, Tipper, "Ty_im", 3);

        //now read the errors, we make their presence in the file optional
        ReadTipperComp(DataFile, Error, "dTx", 0, false);
        ReadTipperComp(DataFile, Error, "dTy", 2, false);
        //we only read in the errors for the real part
        //and then assign the same error to the imaginary part
        for (size_t i = 0; i < nimp - 1; i += 2)
          {
            Error.at(i + 1) = Error.at(i);
          }
      }

    void ReadImpedancesFromMTT(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &Impedances,
        std::vector<double> &Errors, std::vector<double> &Tipper,
        std::vector<double> &TippErr)
      {
        std::ifstream infile;
        double currentreal;
        int nentries = 0;
        infile.open(filename.c_str());
        while (infile.good())
          {
            infile >> currentreal;
            ++nentries;
          }
        infile.close();
        if (nentries == 0)
          throw FatalException("No data in file: " + filename,
          __FILE__, __LINE__);
        if (((nentries - 1) % 23) != 0)
          throw FatalException("Number of records does not match expected: " + filename,
          __FILE__, __LINE__);
        const int nrecords = (nentries - 1) / 23;
        Impedances.resize(nrecords * 8);
        Errors.resize(nrecords * 8);
        Tipper.resize(nrecords * 4);
        TippErr.resize(nrecords * 4);
        Frequencies.resize(nrecords);
        infile.open(filename.c_str());
        int currentrecord = 0;
        int currentip = 0;
        if (infile)
          {
            while (infile.good()) // read in inputfile
              {
                infile >> Frequencies[currentrecord / 8];
                if (infile.good())
                  {
                    //read number of degrees of freedom in file an throw away
                    infile >> currentreal;

                    infile >> Impedances.at(currentrecord)
                        >> Impedances.at(currentrecord + 1)
                        >> Impedances.at(currentrecord + 2)
                        >> Impedances.at(currentrecord + 3)
                        >> Impedances.at(currentrecord + 4)
                        >> Impedances.at(currentrecord + 5)
                        >> Impedances.at(currentrecord + 6)
                        >> Impedances.at(currentrecord + 7);
                    // read in the impedance errors
                    infile >> Errors.at(currentrecord) >> Errors.at(currentrecord + 2)
                        >> Errors.at(currentrecord + 4) >> Errors.at(currentrecord + 6);
                    //fpr the moment we ignore these values in the .mtt file
                    //Tx
                    infile >> Tipper.at(currentip) >> Tipper.at(currentip + 1);
                    //Ty
                    infile >> Tipper.at(currentip + 2) >> Tipper.at(currentip + 3);
                    //dTx
                    infile >> TippErr.at(currentip);
                    //dTy
                    infile >> TippErr.at(currentip + 2);
                    //Coherence Rx
                    infile >> currentreal;
                    //Coherence Ry
                    infile >> currentreal;
                    //Coherence Rz
                    infile >> currentreal;

                    currentrecord += 8;
                    currentip += 4;
                  }
              }
            infile.close();

            //we read the error information into the vector components corresponding to the real part
            //we have to copy that error to the components corresponding to the imaginary part
            for (size_t i = 0; i < Errors.size() - 1; i += 2)
              {
                Errors.at(i + 1) = Errors.at(i);
              }
            for (size_t i = 0; i < TippErr.size() - 1; i += 2)
              {
                TippErr.at(i + 1) = TippErr.at(i);
              }

            //convert the units in the .mtt file (km/s) into S.I units (Ohm)
            const double convfactor = 4.0 * 1e-4 * acos(-1.0);
            auto ConvLambda = [convfactor](double val)
              { return val * convfactor;};
            std::transform(Impedances.begin(), Impedances.end(), Impedances.begin(),
                ConvLambda);
            std::transform(Errors.begin(), Errors.end(), Errors.begin(), ConvLambda);

          }
        else
          {
            throw jif3D::FatalException("File not found: " + filename, __FILE__,
            __LINE__);
          }
      }

    void WriteImpedancesToMtt(const std::string &filenamebase,
        const std::vector<double> &Frequencies, const std::vector<double> &Imp,
        const std::vector<double> &Err, const std::vector<double> &Tipper,
        const std::vector<double> &TipErr)
      {
        const double convfactor = 1.0 / (4.0 * 1e-4 * acos(-1.0));
        std::vector<double> Impedances(Imp.size());
        auto ConvLambda = [convfactor](double val)
          { return convfactor * val;};
        std::transform(Imp.begin(), Imp.end(), Impedances.begin(), ConvLambda);
        std::vector<double> Errors(Err.size());
        std::transform(Err.begin(), Err.end(), Errors.begin(), ConvLambda);

        const size_t nfreq = Frequencies.size();
        const size_t nimp = Impedances.size();
        const size_t ndatapersite = nfreq * 8;
        const size_t nsites = nimp / ndatapersite;
        assert(nimp % ndatapersite == 0);
        for (size_t i = 0; i < nsites; ++i)
          {

            ofstream outfile;
            std::string currfilename = filenamebase + jif3D::stringify(i) + ".mtt";
            outfile.open(currfilename.c_str());

            for (unsigned int j = 0; j < nfreq; ++j) //write mtt-file
              {
                const size_t startindex = (j * nsites + i) * 8;
                const size_t tipindex = (j * nsites + i) * 4;
                outfile << Frequencies.at(j);
                outfile << "   1 \n";

                outfile << Impedances.at(startindex) << " ";
                outfile << Impedances.at(startindex + 1) << " ";
                outfile << Impedances.at(startindex + 2) << " ";
                outfile << Impedances.at(startindex + 3) << " ";
                outfile << Impedances.at(startindex + 4) << " ";
                outfile << Impedances.at(startindex + 5) << " ";
                outfile << Impedances.at(startindex + 6) << " ";
                outfile << Impedances.at(startindex + 7) << " ";
                outfile << "\n";
                outfile << Errors.at(startindex) << " ";
                outfile << Errors.at(startindex + 2) << " ";
                outfile << Errors.at(startindex + 4) << " ";
                outfile << Errors.at(startindex + 6) << " ";
                outfile << Tipper.at(tipindex) << " ";
                outfile << Tipper.at(tipindex + 1) << " ";
                outfile << Tipper.at(tipindex + 2) << " ";
                outfile << Tipper.at(tipindex + 3) << " ";
                outfile << "\n";
                outfile << TipErr.at(tipindex) << " ";
                outfile << TipErr.at(tipindex + 1) << " ";
                outfile << 0.0 << " ";
                outfile << 0.0 << " ";
                outfile << 0.0 << " ";
                outfile << "\n";
                //write out mtt file entries

              }
            outfile.close();
          }
      }

    void WriteImpedancesToMtt(const std::string &filenamebase,
        const std::vector<double> &Frequencies, const std::vector<double> &Imp,
        const std::vector<double> &Err)
      {
        std::vector<double> Dummy(Imp.size() / 2, 0.0);
        WriteImpedancesToMtt(filenamebase, Frequencies, Imp, Err, Dummy, Dummy);
      }

    void WriteJBlock(const std::vector<double> &Frequencies,
        const std::vector<double> &Imp, const std::vector<double> &Err, ofstream &outfile,
        const double convfactor, int CompIndex, int SiteIndex)
      {
        const size_t nfreq = Frequencies.size();
        const size_t nimp = Imp.size();
        const size_t ndatapersite = nfreq * 8;
        const size_t nsites = nimp / ndatapersite;
        for (unsigned int i = 0; i < Frequencies.size(); ++i)
          {
            const size_t startindex = (i * nsites + SiteIndex) * 8 + CompIndex;
            outfile << setfill(' ') << setw(15) << resetiosflags(ios::fixed)
                << 1. / Frequencies.at(i) << " ";
            outfile << setfill(' ') << setw(15) << convfactor * Imp.at(startindex) << " ";
            outfile << setfill(' ') << setw(15) << convfactor * Imp.at(startindex + 1)
                << " ";
            outfile << setfill(' ') << setw(15) << convfactor * Err.at(startindex) << " ";
            outfile << setfill(' ') << setw(15) << setiosflags(ios::fixed) << 1.000 << " "
                << endl;

          }
      }

    void WriteImpedancesToJ(const std::string &filenamebase,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<double> &Imp, const std::vector<double> &Err)
      {
        //const double convfactor = 4.0 * 1e-4 * acos(-1.0);
        const double convfactor = 1.0;
        const size_t nfreq = Frequencies.size();
        const size_t nimp = Imp.size();
        const size_t ndatapersite = nfreq * 8;
        const size_t nsites = nimp / ndatapersite;
        assert(nimp % ndatapersite == 0);
        for (size_t i = 0; i < nsites; ++i)
          {

            ofstream outfile;
            std::string currfilename = filenamebase + jif3D::stringify(i) + ".dat";
            outfile.open(currfilename.c_str());

            outfile << "# J-File Produced by jif3D" << endl;
            outfile << ">LATITUDE  = " << StatXCoord.at(i) << endl;
            outfile << ">LONGITUDE = " << StatYCoord.at(i) << endl;
            outfile << ">ELEVATION = " << StatZCoord.at(i) << endl;
            outfile << ">AZIMUTH   = " << 0.0 << endl;
            outfile << currfilename << endl;

            outfile << "ZXX S.I." << endl;
            outfile << nfreq << endl;
            //WriteJBlock(DataXX,outfile,convfactor);
            WriteJBlock(Frequencies, Imp, Err, outfile, convfactor, 0, i);
            outfile << "ZXY S.I." << endl;
            outfile << nfreq << endl;
            WriteJBlock(Frequencies, Imp, Err, outfile, convfactor, 2, i);

            outfile << "ZYX S.I." << endl;
            outfile << nfreq << endl;
            WriteJBlock(Frequencies, Imp, Err, outfile, convfactor, 4, i);

            outfile << "ZYY S.I." << endl;
            outfile << nfreq << endl;
            WriteJBlock(Frequencies, Imp, Err, outfile, convfactor, 6, i);

            outfile.close();
          }
      }

    inline double rad(double d)
      {
        return d / 180.0 * boost::math::constants::pi<double>();
      }

    void ReadAppResFromAscii(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        std::vector<double> &Imp, std::vector<double> &Err)
      {
        std::ifstream infile;
        infile.open(filename.c_str());
        char dummy[5000];
        //skip header line
        infile.getline(dummy, 5000);
        std::string major, link, site;
        double x, y, elevation, fre, rxy, ryx, pxy, pyx, erxy, eryx, epxy, epyx;

        std::vector<double> ImpTemp, ErrTemp, tmpx, tmpy, tmpz;
        while (infile.good())
          {
            infile >> major >> link >> site >> x >> y >> elevation >> fre >> rxy >> ryx
                >> pxy >> pyx >> erxy >> eryx >> epxy >> epyx;
            if (infile.good())
              {
                std::cout << fre << " " << rxy << " " << pxy << std::endl;
                ImpTemp.push_back(0.0);
                ImpTemp.push_back(0.0);
                fre = std::pow(10, fre);
                const double absZxy = sqrt(twopimu * fre * rxy);
                const double cpxy = cos(rad(pxy));
                const double spxy = sin(rad(pxy));
                double cpyx = cos(rad(pyx));
                double spyx = sin(rad(pyx));
                cpyx = cpyx > 0.0 ? -cpyx : cpyx;
                spyx = spyx > 0.0 ? -spyx : spyx;
                ImpTemp.push_back(absZxy * cpxy);
                ImpTemp.push_back(absZxy * spxy);
                std::cout << "Fre: " << fre << " Rho: "
                    << AppRes(std::complex<double>(absZxy * cpxy, absZxy * spxy), fre)
                    << std::endl;
                std::cout << "Zxy_re: " << absZxy * cpxy << " Zxy_im " << absZxy * spxy
                    << std::endl;
                const double absZyx = sqrt(twopimu * fre * ryx);
                ImpTemp.push_back(absZyx * cpyx);
                ImpTemp.push_back(absZyx * spyx);
                ImpTemp.push_back(0.0);
                ImpTemp.push_back(0.0);
                ErrTemp.push_back(0.0);
                ErrTemp.push_back(0.0);
                double dR = sqrt(
                    twopimu * fre / (4 * rxy) * std::pow(erxy * cpxy, 2)
                        + std::pow(absZxy * spxy * rad(epxy), 2));
                double dI = sqrt(
                    twopimu * fre / (4 * rxy) * std::pow(erxy * spxy, 2)
                        + std::pow(absZxy * cpxy * rad(epxy), 2));
                ErrTemp.push_back(std::max(dR, dI));
                ErrTemp.push_back(std::max(dR, dI));
                dR = sqrt(
                    twopimu * fre / (4 * ryx) * std::pow(eryx * cpxy, 2)
                        + std::pow(absZyx * spyx * rad(epyx), 2));
                dI = sqrt(
                    twopimu * fre / (4 * ryx) * std::pow(eryx * spxy, 2)
                        + std::pow(absZyx * cpyx * rad(epyx), 2));
                ErrTemp.push_back(std::max(dR, dI));
                ErrTemp.push_back(std::max(dR, dI));
                ErrTemp.push_back(0.0);
                ErrTemp.push_back(0.0);
                Frequencies.push_back(fre);
                tmpx.push_back(x);
                tmpy.push_back(y);
                tmpz.push_back(elevation);
              }

          }

        Imp.resize(ImpTemp.size());
        Err.resize(ErrTemp.size());
        //std::copy(ImpTemp.begin(), ImpTemp.end(), Imp.begin());
        //std::copy(ErrTemp.begin(), ErrTemp.end(), Err.begin());
        bool Ascending = true;
        if (Frequencies.at(0) > Frequencies.at(1))
          {
            Ascending = false;
          }
        if (Ascending)
          {
            std::sort(Frequencies.begin(), Frequencies.end(), std::less<double>());
          }
        else
          {
            std::sort(Frequencies.begin(), Frequencies.end(), std::greater<double>());
          }
        Frequencies.erase(
            std::unique(Frequencies.begin(), Frequencies.end(), [](double a, double b)
              { return std::abs((a-b)/std::max(a,b)) < 0.001;}), Frequencies.end());
        //StatXCoord.erase(std::unique(StatXCoord.begin(), StatXCoord.end()),
        //    StatXCoord.end());
        //StatYCoord.erase(std::unique(StatYCoord.begin(), StatYCoord.end()),
        //    StatYCoord.end());
        //StatZCoord.erase(std::unique(StatZCoord.begin(), StatZCoord.end()),
        //    StatZCoord.end());
        const size_t nfreq = Frequencies.size();
        const size_t nstat = tmpx.size() / nfreq;
        for (size_t i = 0; i < nstat; ++i)
          {
            StatXCoord.push_back(tmpx.at(i * nfreq));
            StatYCoord.push_back(tmpy.at(i * nfreq));
            StatZCoord.push_back(tmpz.at(i * nfreq));
            for (size_t j = 0; j < nfreq; ++j)
              {
                size_t Impindex = (j * nstat + i) * 8;
                size_t TempIndex = (i * nfreq + j) * 8;
                std::copy(ImpTemp.begin() + TempIndex, ImpTemp.begin() + TempIndex + 8,
                    Imp.begin() + Impindex);
                std::copy(ErrTemp.begin() + TempIndex, ErrTemp.begin() + TempIndex + 8,
                    Err.begin() + Impindex);
              }
          }
        std::cout << "Stations: " << StatXCoord.size() << std::endl;
        std::cout << "Frequencies: " << Frequencies.size() << std::endl;
        std::copy(Frequencies.begin(), Frequencies.end(),
            std::ostream_iterator<double>(std::cout, " "));
        std::cout << std::endl;
        std::cout << "Impendances: " << Imp.size() << std::endl;
      }

    void WriteAppResToAscii(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<double> &Imp, const std::vector<double> &Err)
      {
        std::ofstream outfile(filename.c_str());
        outfile.precision(3);
        outfile
            << " Major    Link    Site           x            y          elevation       Fre         rxy         ryx";
        outfile
            << "          pxy          pyx            erxy         eryx         epxy         epyx\n";
        const size_t nfreq = Frequencies.size();
        const size_t nstat = StatXCoord.size();

        for (size_t j = 0; j < nstat; ++j)
          {
            for (size_t i = 0; i < nfreq; ++i)
              {
                const size_t currindex = 8 * (i * nstat + j);
                const double CurrFreq = Frequencies.at(i);
                double rhoxy = jif3D::AppRes(
                    std::complex<double>(Imp.at(currindex + 2), Imp.at(currindex + 3)),
                    CurrFreq);
                double pxy = jif3D::ImpedancePhase(
                    std::complex<double>(Imp.at(currindex + 2), Imp.at(currindex + 3)));
                double rhoyx = jif3D::AppRes(
                    std::complex<double>(Imp.at(currindex + 4), Imp.at(currindex + 5)),
                    CurrFreq);
                double pyx = jif3D::ImpedancePhase(
                    std::complex<double>(Imp.at(currindex + 4), Imp.at(currindex + 5)));
                outfile << std::fixed;
                outfile << std::setw(10) << j << std::setw(10) << j << std::setw(10) << j;
                outfile << std::setw(15) << StatXCoord[j];
                outfile << std::setw(15) << StatYCoord[j];
                outfile << std::setw(15) << StatZCoord[j];
                outfile << std::setw(15) << log10(CurrFreq);
                outfile << std::setw(15) << rhoxy;
                outfile << std::setw(15) << rhoyx;
                outfile << std::setw(15) << pxy;
                outfile << std::setw(15) << pyx;
                outfile << " 0      0     0     0\n";
              }
          }
      }

    void ReadImpedancesFromModEM(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        std::vector<double> &Imp, std::vector<double> &Err)
      {
        const double convfactor = 4.0 * 1e-4 * acos(-1.0);
        std::ifstream infile(filename.c_str());
        char dummy[1024];
        //swallow up the first 7 lines with header information
        for (size_t i = 0; i < 7; ++i)
          infile.getline(dummy, 1024);
        //ignore the fist character (>) of the last header line
        infile.ignore(1);
        size_t nfreq = 0;
        size_t nsites = 0;
        infile >> nfreq >> nsites;
        //swallow up the rest of the line
        infile.getline(dummy, 1024);
        Imp.resize(nfreq * nsites * 8);
        std::fill(Imp.begin(), Imp.end(), 1.0);
        Err.resize(Imp.size());
        std::fill(Err.begin(), Err.end(), 100.0);
        Frequencies.resize(nfreq);
        StatXCoord.resize(nsites);
        StatYCoord.resize(nsites);
        StatZCoord.resize(nsites);
        bool CanConvert = true;

        std::vector<double> AllFreq, AllX, AllY, AllZ, AllImpReal, AllImpImag, AllErr;
        std::vector<std::string> StatName, CompName;
        while (CanConvert && infile.good())
          {
            try
              {
                string line;
                std::getline(infile, line);
                typedef std::vector<std::string> split_vector_type;

                split_vector_type SplitVec; // #2: Search for tokens
                split(SplitVec, line, boost::is_any_of(" \n\r"),
                    boost::token_compress_on);

                double p;
                jif3D::convert(SplitVec[0], p);
                AllFreq.push_back(1.0 / p);
                StatName.push_back(SplitVec[1]);
                jif3D::convert(SplitVec[4], p);
                AllX.push_back(p);
                jif3D::convert(SplitVec[5], p);
                AllY.push_back(p);
                jif3D::convert(SplitVec[6], p);
                AllZ.push_back(p);
                CompName.push_back(SplitVec[7]);
                jif3D::convert(SplitVec[8], p);
                AllImpReal.push_back(p);
                jif3D::convert(SplitVec[9], p);
                AllImpImag.push_back(p);
                jif3D::convert(SplitVec[10], p);
                AllErr.push_back(p);
              } catch (std::exception &e)
              {
                CanConvert = false;
              }
          }
        std::vector<std::string> UniqStats(StatName);
        UniqStats.erase(std::unique(UniqStats.begin(), UniqStats.end()), UniqStats.end());

        Frequencies = AllFreq;
        auto FuzzComp = [](double a, double b)
          { return std::abs((a-b)/std::max(a,b)) < 0.001;};
        std::sort(Frequencies.begin(), Frequencies.end(), std::greater<double>());
        Frequencies.erase(std::unique(Frequencies.begin(), Frequencies.end(), FuzzComp),
            Frequencies.end());
        std::vector<std::string> Comps =
          { "ZXX", "ZXY", "ZYX", "ZYY" };
        size_t ndata = AllFreq.size();
        for (size_t i = 0; i < ndata; ++i)
          {
            auto FreqIter = std::find(Frequencies.begin(), Frequencies.end(), AllFreq[i]);
            auto StatIter = std::find(UniqStats.begin(), UniqStats.end(), StatName[i]);
            auto CompIter = std::find(Comps.begin(), Comps.end(), CompName[i]);
            size_t FreqIndex = std::distance(Frequencies.begin(), FreqIter);
            size_t StatIndex = std::distance(UniqStats.begin(), StatIter);
            size_t CompIndex = std::distance(Comps.begin(), CompIter);
            StatXCoord[StatIndex] = AllX[i];
            StatYCoord[StatIndex] = AllY[i];
            StatZCoord[StatIndex] = AllZ[i];
            size_t ImpIndex = 2 * CompIndex + 8 * StatIndex + 8 * nsites * FreqIndex;
            Imp[ImpIndex] = convfactor * AllImpReal[i];
            Imp[ImpIndex + 1] = convfactor * AllImpImag[i];
            Err[ImpIndex] = convfactor * AllErr[i];
            Err[ImpIndex + 1] = convfactor * AllErr[i];
          }

      }

    void ReadTipperFromModEM(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        std::vector<double> &Tip, std::vector<double> &Err)
      {
        std::ifstream infile(filename.c_str());
        char dummy[1024];
        //swallow up the first 7 lines with header information
        for (size_t i = 0; i < 7; ++i)
          infile.getline(dummy, 1024);
        //ignore the fist character (>) of the last header line
        infile.ignore(1);
        size_t nfreq = 0;
        size_t nsites = 0;
        infile >> nfreq >> nsites;
        //swallow up the rest of the line
        infile.getline(dummy, 1024);
        Tip.resize(nfreq * nsites * 4);
        std::fill(Tip.begin(), Tip.end(), 1.0);
        Err.resize(Tip.size());
        std::fill(Err.begin(), Err.end(), 100.0);
        Frequencies.resize(nfreq);
        StatXCoord.resize(nsites);
        StatYCoord.resize(nsites);
        StatZCoord.resize(nsites);
        bool CanConvert = true;

        std::vector<double> AllFreq, AllX, AllY, AllZ, AllTipReal, AllTipImag, AllErr;
        std::vector<std::string> StatName, CompName;
        while (CanConvert && infile.good())
          {
            try
              {
                string line;
                std::getline(infile, line);
                typedef std::vector<std::string> split_vector_type;

                split_vector_type SplitVec; // #2: Search for tokens
                split(SplitVec, line, boost::is_any_of(" \n\r"),
                    boost::token_compress_on);

                double p;
                jif3D::convert(SplitVec[0], p);
                AllFreq.push_back(1.0 / p);
                StatName.push_back(SplitVec[1]);
                jif3D::convert(SplitVec[4], p);
                AllX.push_back(p);
                jif3D::convert(SplitVec[5], p);
                AllY.push_back(p);
                jif3D::convert(SplitVec[6], p);
                AllZ.push_back(p);
                CompName.push_back(SplitVec[7]);
                jif3D::convert(SplitVec[8], p);
                AllTipReal.push_back(p);
                jif3D::convert(SplitVec[9], p);
                AllTipImag.push_back(p);
                jif3D::convert(SplitVec[10], p);
                AllErr.push_back(p);
              } catch (std::exception &e)
              {
                CanConvert = false;
              }
          }
        std::vector<std::string> UniqStats(StatName);
        UniqStats.erase(std::unique(UniqStats.begin(), UniqStats.end()), UniqStats.end());

        Frequencies = AllFreq;
        auto FuzzComp = [](double a, double b)
          { return std::abs((a-b)/std::max(a,b)) < 0.001;};
        std::sort(Frequencies.begin(), Frequencies.end(), std::greater<double>());
        Frequencies.erase(std::unique(Frequencies.begin(), Frequencies.end(), FuzzComp),
            Frequencies.end());
        std::vector<std::string> Comps =
          { "TX", "TY" };
        size_t ndata = AllFreq.size();
        for (size_t i = 0; i < ndata; ++i)
          {
            auto FreqIter = std::find(Frequencies.begin(), Frequencies.end(), AllFreq[i]);
            auto StatIter = std::find(UniqStats.begin(), UniqStats.end(), StatName[i]);
            auto CompIter = std::find(Comps.begin(), Comps.end(), CompName[i]);
            size_t FreqIndex = std::distance(Frequencies.begin(), FreqIter);
            size_t StatIndex = std::distance(UniqStats.begin(), StatIter);
            size_t CompIndex = std::distance(Comps.begin(), CompIter);
            StatXCoord[StatIndex] = AllX[i];
            StatYCoord[StatIndex] = AllY[i];
            StatZCoord[StatIndex] = AllZ[i];
            size_t TipIndex = 2 * CompIndex + 4 * StatIndex + 4 * nsites * FreqIndex;
            Tip[TipIndex] = AllTipReal[i];
            Tip[TipIndex + 1] = AllTipImag[i];
            Err[TipIndex] = AllErr[i];
            Err[TipIndex + 1] = AllErr[i];
          }
      }

    void WriteModEMLine(std::ofstream &outfile, double Zr, double Zi, double Err,
        double convfactor = 1.0)
      {
        outfile << std::setw(15) << Zr / convfactor << " " << std::setw(15)
            << Zi / convfactor << " " << std::setw(15) << Err / convfactor << "\n";
      }

    void WriteImpedancesToModEM(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<double> &Imp, const std::vector<double> &Err)
      {
        std::ofstream outfile(filename.c_str());
        outfile.precision(6);

        outfile << "# Description: \n";
        outfile
            << "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error \n";
        outfile << "> Full_Impedance\n";
        outfile << "> exp(+i\\omega t)\n";
        outfile << "> [mV/km]/[nT]\n";
        outfile << "> 0.00\n";
        outfile << "> 0.0 0.0\n";
        size_t nsites = StatXCoord.size();
        size_t nfreqs = Frequencies.size();
        outfile << "> " << nfreqs << " " << nsites << "\n";
        outfile.setf(std::ios::scientific);
        const double convfactor = 4.0 * 1e-4 * acos(-1.0);
        for (size_t i = 0; i < nsites; ++i)
          {
            std::string SiteName = "Site" + std::to_string(i);
            std::ostringstream SiteLine;
            SiteLine << " " << SiteName << "  0.0  0.0 " << StatXCoord.at(i) << " "
                << StatYCoord.at(i) << " " << StatZCoord.at(i) << " ";
            for (size_t j = 0; j < nfreqs; ++j)
              {
                size_t index = 8 * (nsites * j + i);
                double period = 1.0 / Frequencies.at(j);
                outfile << period << SiteLine.str();
                outfile << " ZXX ";
                WriteModEMLine(outfile, Imp.at(index), Imp.at(index + 1), Err.at(index),
                    convfactor);

                outfile << period << SiteLine.str();
                outfile << " ZXY ";
                WriteModEMLine(outfile, Imp.at(index + 2), Imp.at(index + 3),
                    Err.at(index + 2), convfactor);

                outfile << period << SiteLine.str();
                outfile << " ZYX ";
                WriteModEMLine(outfile, Imp.at(index + 4), Imp.at(index + 5),
                    Err.at(index + 4), convfactor);

                outfile << period << SiteLine.str();
                outfile << " ZYY ";
                WriteModEMLine(outfile, Imp.at(index + 6), Imp.at(index + 7),
                    Err.at(index + 6), convfactor);
              }
          }
      }

    J3DEXPORT void WriteTipperToModEM(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const jif3D::rvec &Tip, const jif3D::rvec &Err)
      {
        std::ofstream outfile(filename.c_str());
        outfile.precision(6);

        outfile << "# Description: \n";
        outfile
            << "# Period(s) Code GG_Lat GG_Lon X(m) Y(m) Z(m) Component Real Imag Error \n";
        outfile << "> Full_Vertical_Components\n";
        outfile << "> exp(+i\\omega t)\n";
        outfile << "> []\n";
        outfile << "> 0.00\n";
        outfile << "> 0.0 0.0\n";
        size_t nsites = StatXCoord.size();
        size_t nfreqs = Frequencies.size();
        outfile << "> " << nfreqs << " " << nsites << "\n";
        outfile.setf(std::ios::scientific);
        for (size_t i = 0; i < nsites; ++i)
          {
            std::string SiteName = "Site" + std::to_string(i);
            std::ostringstream SiteLine;
            SiteLine << " " << SiteName << "  0.0  0.0 " << StatXCoord.at(i) << " "
                << StatYCoord.at(i) << " " << StatZCoord.at(i) << " ";
            for (size_t j = 0; j < nfreqs; ++j)
              {
                size_t index = 4 * (nsites * j + i);
                double period = 1.0 / Frequencies.at(j);
                outfile << period << SiteLine.str();
                outfile << " TX ";
                WriteModEMLine(outfile, Tip(index), Tip(index + 1), Err(index));

                outfile << period << SiteLine.str();
                outfile << " TY ";
                WriteModEMLine(outfile, Tip(index + 2), Tip(index + 3), Err(index + 2));

              }
          }

      }

    J3DEXPORT void ReadImpedancesFromJ(const std::string &filename,
        std::vector<double> &Frequencies, double &StatXCoord, double &StatYCoord,
        double &StatZCoord, std::vector<double> &Imp, std::vector<double> &Err)
      {
        std::ifstream infile;
        infile.open(filename.c_str());
        std::string line;
        std::vector<std::string> SplitVec; // #2: Search for tokens

        line = jif3D::FindToken(infile, ">LATITUDE");
        split(SplitVec, line, boost::is_any_of("="), boost::token_compress_on);
        StatXCoord = std::stod(SplitVec.at(1));
        infile.seekg(0);
        line = jif3D::FindToken(infile, ">LONGITUDE");
        split(SplitVec, line, boost::is_any_of("="), boost::token_compress_on);
        StatYCoord = std::stod(SplitVec.at(1));
        infile.seekg(0);
        line = jif3D::FindToken(infile, ">LATITUDE");
        split(SplitVec, line, boost::is_any_of("="), boost::token_compress_on);
        StatZCoord = std::stod(SplitVec.at(1));
        infile.seekg(0);
        bool HasZ = false;
        bool HasRhoPhi = false;
        try
          {
            jif3D::FindToken(infile, "RXX");
            HasRhoPhi = true;
          } catch (jif3D::FatalException &e)
          {

          }
        infile.clear();
        infile.seekg(0);
        try
          {
            jif3D::FindToken(infile, "ZXX");
            HasZ = true;
          } catch (jif3D::FatalException &e)
          {

          }
        infile.clear();
        infile.seekg(0);
        if (HasZ)
          {
            size_t nfreq = 0;
            jif3D::FindToken(infile, "ZXX");
            infile >> nfreq;
            Imp.resize(nfreq * 8);
            Err.resize(nfreq * 8);
            double dummy;
            Frequencies.resize(nfreq);
            for (size_t i = 0; i < nfreq; ++i)
              {
                infile >> Frequencies.at(i);
                Frequencies.at(i) =
                    Frequencies.at(i) > 0.0 ?
                        1.0 / Frequencies.at(i) : std::abs(Frequencies.at(i));
                infile >> Imp.at(i * 8);
                infile >> Imp.at(i * 8 + 1);
                infile >> Err.at(i * 8);
                Err.at(i * 8 + 1) = Err.at(i * 8);
                infile >> dummy;
              }
            jif3D::FindToken(infile, "ZXY");
            infile >> dummy;
            for (size_t i = 0; i < nfreq; ++i)
              {
                infile >> dummy;
                infile >> Imp.at(i * 8 + 2);
                infile >> Imp.at(i * 8 + 3);
                infile >> Err.at(i * 8 + 2);
                Err.at(i * 8 + 3) = Err.at(i * 8 + 2);
                infile >> dummy;
              }
            jif3D::FindToken(infile, "ZYX");
            infile >> dummy;
            for (size_t i = 0; i < nfreq; ++i)
              {
                infile >> dummy;
                infile >> Imp.at(i * 8 + 4);
                infile >> Imp.at(i * 8 + 5);
                infile >> Err.at(i * 8 + 4);
                Err.at(i * 8 + 5) = Err.at(i * 8 + 4);
                infile >> dummy;
              }
            jif3D::FindToken(infile, "ZYY");
            infile >> dummy;
            for (size_t i = 0; i < nfreq; ++i)
              {
                infile >> dummy;
                infile >> Imp.at(i * 8 + 6);
                infile >> Imp.at(i * 8 + 7);
                infile >> Err.at(i * 8 + 6);
                Err.at(i * 8 + 7) = Err.at(i * 8 + 6);
                infile >> dummy;
              }
            //std::transform(Imp.begin(), Imp.end(), Imp.begin(), [convfactor](double val)
            //  { return val/convfactor;});
            //std::transform(Err.begin(), Err.end(), Err.begin(), [convfactor](double val)
            //              { return val/convfactor;});
          }
        else
          {
            if (HasRhoPhi)
              {

              }
            else
              {
                throw jif3D::FatalException("No valid MT data in file", __FILE__,
                __LINE__);
              }
          }
      }
  }

