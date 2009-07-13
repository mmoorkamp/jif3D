//============================================================================
// Name        : ReadWriteImpedance.cpp
// Author      : Jul 13, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "../Global/NumUtil.h"
#include "../ModelBase/NetCDFTools.h"
#include "ReadWriteImpedances.h"

namespace jiba
  {

    const std::string FreqDimName = "Frequency";

    void WriteImpedanceComp(NcFile &NetCDFFile, NcDim *StatNumDim,
        NcDim *FreqDim, const jiba::rvec &Impedances,
        const std::string &CompName, const size_t compindex)
      {
        NcVar *CompVar = NetCDFFile.add_var(CompName.c_str(), ncDouble,
            StatNumDim, FreqDim);
        jiba::rvec Component(StatNumDim->size() * FreqDim->size());
        for (size_t i = 0; i < Component.size(); ++i)
          {
            Component(i) = Impedances(i * 8 + compindex);
          }
        CompVar->put(&Component[0], Component.size());
      }

    void ReadImpedanceComp(NcFile &NetCDFFile, jiba::rvec &Impedances,
        const std::string &CompName, const size_t compindex)
      {
        NcVar *SizeVar = NetCDFFile.get_var(CompName.c_str());
        const size_t nvalues = Impedances.size() / 8;
        jiba::rvec Temp(nvalues);
        SizeVar->get(&Temp[0], nvalues);
        for (size_t i = 0; i < nvalues; ++i)
          {
            Impedances(i * 8 + compindex) = Temp(i);
          }

      }

    void WriteImpedancesToNetCDF(const std::string &filename,
        const jiba::rvec &Frequencies, const jiba::rvec &StatXCoord,
        const jiba::rvec &StatYCoord, const jiba::rvec &StatZCoord,
        const jiba::rvec &Impedances)
      {
        const size_t nstats = StatXCoord.size();
        const size_t nfreqs = Frequencies.size();
        assert(nstats == StatYCoord.size());
        assert(nstats == StatYCoord.size());
        assert(Impedances.size() == nstats * nfreqs * 8);
        //create a netcdf file
        NcFile DataFile(filename.c_str(), NcFile::Replace);

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

        NcDim *FreqDim = DataFile.add_dim(FreqDimName.c_str(),
            Frequencies.size());
        NcVar *FreqVar = DataFile.add_var(FreqDimName.c_str(), ncDouble,
            FreqDim);
        FreqVar->put(&Frequencies[0], nfreqs);
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

    void ReadImpedancesFromNetCDF(const std::string &filename,
        jiba::rvec &Frequencies, jiba::rvec &StatXCoord,
        jiba::rvec &StatYCoord, jiba::rvec &StatZCoord, jiba::rvec &Impedances)
      {
        NcFile DataFile(filename.c_str(), NcFile::ReadOnly);
        ReadVec(DataFile, MeasPosXName, StationNumberName, StatXCoord);
        ReadVec(DataFile, MeasPosXName, StationNumberName, StatYCoord);
        ReadVec(DataFile, MeasPosXName, StationNumberName, StatZCoord);
        ReadVec(DataFile, FreqDimName, FreqDimName, Frequencies);

      }
  }
