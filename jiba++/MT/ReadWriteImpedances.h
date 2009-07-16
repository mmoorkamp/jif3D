//============================================================================
// Name        : ReadWriteImpedances.h
// Author      : Jul 13, 2009
// Version     : 
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef READWRITEIMPEDANCES_H_
#define READWRITEIMPEDANCES_H_

#include "../Global/VecMat.h"

namespace jiba
  {

    void WriteImpedancesToNetCDF(const std::string &filename,
        const jiba::rvec &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord,
        const std::vector<double> &StatZCoord, const jiba::rvec &Impedances);

    void ReadImpedancesFromNetCDF(const std::string &filename,
        jiba::rvec &Frequencies, std::vector<double> &StatXCoord, std::vector<
            double> &StatYCoord, std::vector<double> &StatZCoord,
        jiba::rvec &Impedances);
  }

#endif /* READWRITEIMPEDANCES_H_ */
