//============================================================================
// Name        : ReadWriteImpedance.cpp
// Author      : Jul 13, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <fstream>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "../Global/FatalException.h"
#include "../Global/NumUtil.h"
#include "../ModelBase/NetCDFTools.h"
#include "ReadWriteImpedances.h"

namespace jiba
  {

    const std::string FreqDimName = "Frequency";
    //write one compoment of the impedance tensor to a netcdf file
    //this is an internal helper function
    void WriteImpedanceComp(NcFile &NetCDFFile, NcDim *StatNumDim,
        NcDim *FreqDim, const jiba::rvec &Impedances,
        const std::string &CompName, const size_t compindex)
      {
        NcVar *CompVar = NetCDFFile.add_var(CompName.c_str(), ncDouble,
            FreqDim, StatNumDim);
        jiba::rvec Component(FreqDim->size() * StatNumDim->size());
        for (size_t i = 0; i < Component.size(); ++i)
          {
            Component(i) = Impedances(i * 8 + compindex);
          }
        CompVar->put(&Component[0], FreqDim->size(), StatNumDim->size());
      }
    //read one compoment of the impedance tensor from a netcdf file
    //this is an internal helper function
    void ReadImpedanceComp(NcFile &NetCDFFile, jiba::rvec &Impedances,
        const std::string &CompName, const size_t compindex)
      {
        NcVar *SizeVar = NetCDFFile.get_var(CompName.c_str());
        const size_t nvalues = Impedances.size() / 8;
        jiba::rvec Temp(nvalues);
        SizeVar->get(&Temp[0], SizeVar->edges()[0], SizeVar->edges()[1]);
        for (size_t i = 0; i < nvalues; ++i)
          {
            Impedances(i * 8 + compindex) = Temp(i);
          }

      }

    void WriteImpedancesToNetCDF(const std::string &filename,
        const std::vector<double> &Frequencies,
        const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord,
        const std::vector<double> &StatZCoord, const jiba::rvec &Impedances)
      {
        const size_t nstats = StatXCoord.size();
        const size_t nfreqs = Frequencies.size();
        assert(nstats == StatYCoord.size());
        assert(nstats == StatYCoord.size());
        assert(Impedances.size() == nstats * nfreqs * 8);
        //create a netcdf file
        NcFile DataFile(filename.c_str(), NcFile::Replace);
        //each station gets an index number that we use as a dimension
        //in the netcdf file
        NcDim *StatNumDim = DataFile.add_dim(StationNumberName.c_str(), nstats);
        std::vector<int> StationNumber;
        std::generate_n(back_inserter(StationNumber), nstats, IntSequence(0));
        NcVar *StatNumVar = DataFile.add_var(StationNumberName.c_str(), ncInt,
            StatNumDim);
        StatNumVar->put(&StationNumber[0], nstats);
        //write out the measurement coordinates
        WriteVec(DataFile, MeasPosXName, StatXCoord, StatNumDim, "m");
        WriteVec(DataFile, MeasPosYName, StatYCoord, StatNumDim, "m");
        WriteVec(DataFile, MeasPosZName, StatZCoord, StatNumDim, "m");
        //write out the frequencies that we store
        NcDim *FreqDim = DataFile.add_dim(FreqDimName.c_str(),
            Frequencies.size());
        NcVar *FreqVar = DataFile.add_var(FreqDimName.c_str(), ncDouble,
            FreqDim);
        FreqVar->put(&Frequencies[0], nfreqs);
        //and now we can write all the impedance components
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxx_re",
            0);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxx_im",
            1);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxy_re",
            2);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zxy_im",
            3);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyx_re",
            4);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyx_im",
            5);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyy_re",
            6);
        WriteImpedanceComp(DataFile, StatNumDim, FreqDim, Impedances, "Zyy_im",
            7);
      }

    void ReadImpedancesFromNetCDF(const std::string &filename, std::vector<
        double> &Frequencies, std::vector<double> &StatXCoord, std::vector<
        double> &StatYCoord, std::vector<double> &StatZCoord,
        jiba::rvec &Impedances)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        ReadVec(DataFile, MeasPosXName, StationNumberName, StatXCoord);
        ReadVec(DataFile, MeasPosYName, StationNumberName, StatYCoord);
        ReadVec(DataFile, MeasPosZName, StationNumberName, StatZCoord);
        ReadVec(DataFile, FreqDimName, FreqDimName, Frequencies);
        Impedances.resize(Frequencies.size() * StatXCoord.size() * 8);
        ReadImpedanceComp(DataFile, Impedances, "Zxx_re", 0);
        ReadImpedanceComp(DataFile, Impedances, "Zxx_im", 1);
        ReadImpedanceComp(DataFile, Impedances, "Zxy_re", 2);
        ReadImpedanceComp(DataFile, Impedances, "Zxy_im", 3);
        ReadImpedanceComp(DataFile, Impedances, "Zyx_re", 4);
        ReadImpedanceComp(DataFile, Impedances, "Zyx_im", 5);
        ReadImpedanceComp(DataFile, Impedances, "Zyy_re", 6);
        ReadImpedanceComp(DataFile, Impedances, "Zyy_im", 7);

      }

    void ReadImpendacesFromMTT(const std::string &filename,
        std::vector<double> &Frequencies, jiba::rvec &Impedances)
      {
        std::ifstream infile;
        double currentreal, currentimag;
        int nentries = 0;
        infile.open(filename.c_str());
        while (infile.good())
          {
            infile >> currentreal;
            ++nentries;
          }
        infile.close();
        if (((nentries - 1) % 23) != 0)
          throw FatalException("Number of records does not match expected: "
              + filename);
        const int nrecords = (nentries - 1) / 23;
        Impedances.resize(nrecords * 8);
        Frequencies.resize(nrecords);
        infile.open(filename.c_str());
        int currentrecord = 0;
        if (infile)
          {
            while (infile.good()) // read in inputfile
              {
                infile >> Frequencies[currentrecord];
                if (infile.good())
                  {
                    //read number of degrees of freedom in file an throw away
                    infile >> currentreal;

                    infile >> Impedances(currentrecord*8) >> Impedances(currentrecord*8+1)
                         >> Impedances(currentrecord*8+2) >> Impedances(currentrecord*8+3)
                         >> Impedances(currentrecord*8+4) >> Impedances(currentrecord*8+5)
                         >> Impedances(currentrecord*8+6) >> Impedances(currentrecord*8+7);
                    //we skip the rest of the information contained in the file for now
                    //dZxx
                    infile >> currentreal;
                    //dZxy
                    infile >> currentreal;
                    //dZyx
                    infile >> currentreal;
                    //dZyy
                    infile >> currentreal;
                    //Tx
                    infile >> currentreal >> currentimag;
                    //Ty
                    infile >> currentreal >> currentimag;
                    //dTx
                    infile >> currentreal;
                    //dTy
                    infile >> currentreal;
                    //Coherence Rx
                    infile >> currentreal;
                    //Coherence Ry
                    infile >> currentreal;
                    //Coherence Rz
                    infile >> currentreal;

                    ++currentrecord;
                  }
              }
            infile.close();

          }
        else
          {
            throw jiba::FatalException("File not found: " + filename);
          }
      }
  }
