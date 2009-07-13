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
        const jiba::rvec &Frequencies, const jiba::rvec &StatXCoord,
        const jiba::rvec &StatYCoord, const jiba::rvec &StatZCoord,
        const jiba::rvec &Impedances);

    void ReadImpedancesFromNetCDF(const std::string &filename,
        jiba::rvec &Frequencies, jiba::rvec &StatXCoord,
        jiba::rvec &StatYCoord, jiba::rvec &StatZCoord, jiba::rvec &Impedances);

  }

#endif /* READWRITEIMPEDANCES_H_ */
