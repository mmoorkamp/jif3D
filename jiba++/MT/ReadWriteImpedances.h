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
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! Write magnetotelluric impedances to a netcdf file
    /*! We can save MT impedances for several stations in a netcdf file for storage
     * and analysis with external programs.
     * @param filename The name for the netcdf file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param Impedances The impedances as a vector of real numbers. 8 consecutive elements form the impedance matrix for one frequency and site, all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     */
    void WriteImpedancesToNetCDF(const std::string &filename,
        const std::vector<double> &Frequencies,
        const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord,
        const std::vector<double> &StatZCoord, const jiba::rvec &Impedances);

    //! Read magnetotelluric impedances from a netcdf file
    /*! Read MT impedances for several stations from a netcdf file.
     * @param filename The name for the netcdf file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param Impedances The impedances as a vector of real numbers. 8 consecutive elements form the impedance matrix for one frequency and site, all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     */
    void ReadImpedancesFromNetCDF(const std::string &filename, std::vector<
        double> &Frequencies, std::vector<double> &StatXCoord, std::vector<
        double> &StatYCoord, std::vector<double> &StatZCoord,
        jiba::rvec &Impedances);
  /* @} */
  }

#endif /* READWRITEIMPEDANCES_H_ */
