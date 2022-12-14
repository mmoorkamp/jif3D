//============================================================================
// Name        : ReadWriteImpedances.h
// Author      : Jul 13, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef READWRITEIMPEDANCES_H_
#define READWRITEIMPEDANCES_H_

#include "../Global/Jif3DGlobal.h"
#include "../Global/Serialization.h"
#include "../Global/VecMat.h"

namespace jif3D
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
     * @param Impedances The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Errors Optional parameter containing the impedance errors with the same number of Elements with Impedances.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     */
    J3DEXPORT void WriteImpedancesToNetCDF(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const jif3D::rvec &Impedances, const jif3D::rvec &Errors = jif3D::rvec(),
        const std::vector<double> &Distortion = std::vector<double>());

    //! Read magnetotelluric impedances from a netcdf file
    /*! Read MT impedances for several stations from a netcdf file.
     * @param filename The name for the netcdf file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param Impedances The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param ImpError The error of the impedances, has the same number of elements as Impedances
     *         but contains all zeros if no error information present in the file.
     *         Also the error for the real and imaginary parts are always identical.
     */
    J3DEXPORT void ReadImpedancesFromNetCDF(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        jif3D::rvec &Impedances, jif3D::rvec &ImpError, std::vector<double> &Distortion);

    //! A very basic routine to read impedances at a single site from a .mtt file in the format used by University of Goettingen
    /*! A very basic routine to read impedances at a single site from a .mtt file in the
     * format used by University of Goettingen.
     * @param filename The name of the .mtt file
     * @param Frequencies The frequencies contained in the file
     * @param Impedances The impedances in the same convention as above
     * @param Errors The impedance errors recorded in the file, will have the same number of elements as Impedances
     *        but the values for the real and imaginary part of each element are always identical
     */
    J3DEXPORT void ReadImpedancesFromMTT(const std::string &filename,
        std::vector<double> &Frequencies, jif3D::rvec &Impedances, jif3D::rvec &Errors);

    //! A very basic routine to write impedances for several sites to  .mtt files in the format used by University of Goettingen.
    /*! A very basic routine to write impedances for several sites to  .mtt files in the format used by University of Goettingen.
     * @param filenamebase The start of the name of the .mtt file, each file will get a number appended
     * @param Frequencies The frequencies contained in the file
     * @param Imp The impedances in the same convention as above
     * @param Err The impedance errors recorded in the file, will have the same number of elements as Impedances
     *        but the values for the real and imaginary part of each element are always identical
     */
    J3DEXPORT void WriteImpedancesToMtt(const std::string &filenamebase,
        const std::vector<double> &Frequencies, const jif3D::rvec &Imp,
        const jif3D::rvec &Err);

    //! Reads apparent resistivity and phase information from an ascii file with all stations joined together
    /*!
     * @param filename The name of the ascii file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param Imp The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Err Impedance errors with the same number of Elements with Impedances.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     */
    J3DEXPORT void ReadAppResFromAscii(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        jif3D::rvec &Imp, jif3D::rvec &Err);

    //! Write apparent resistivity and phase information to an ascii file with all stations joined together
    /*!
     * @param filename The name of the ascii file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param Imp The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Err Impedance errors with the same number of Elements with Impedances.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     */
    J3DEXPORT void WriteAppResToAscii(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const jif3D::rvec &Imp, const jif3D::rvec &Err);

    //! Reads impedances from an ascii file as written by ModEM
    /*!
     * @param filename The name of the ascii file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param Imp The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Err Impedance errors with the same number of Elements with Impedances.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     */
    J3DEXPORT void ReadImpedancesFromModEM(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        jif3D::rvec &Imp, jif3D::rvec &Err);

    //! Write impedances to an ascii file as written by ModEM
    /*!
     * @param filename The name of the ascii file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param Imp The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Err Impedance errors with the same number of Elements with Impedances.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     */
    J3DEXPORT void WriteImpedancesToModEM(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const jif3D::rvec &Imp, const jif3D::rvec &Err);
  /* @} */
  }

#endif /* READWRITEIMPEDANCES_H_ */
