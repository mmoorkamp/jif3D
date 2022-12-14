//============================================================================
// Name        : ReadWriteImpedances.h
// Author      : Jul 13, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef READWRITEIMPEDANCES_H_
#define READWRITEIMPEDANCES_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"

#include "../Global/NetCDFPortHelper.h"

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */

    //!write one component of the impedance tensor to a netcdf file
    /*!this is an internal helper function
     * @param NetCDFFile The name for the netcdf file
     * @param StatNumDim The number of sites
     * @param FreqDim The number of frequencies
     * @param Impedances The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param CompName Name of the component to be written to netcdf file ('Zxx_re' ... 'Zyy_im')
     * @param compindex component index (0-7)
     */
    J3DEXPORT void WriteImpedanceComp(netCDF::NcFile &NetCDFFile,
        netCDF::NcDim &StatNumDim, netCDF::NcDim &FreqDim,
        const std::vector<double> &Impedances, const std::string &CompName,
        const size_t compindex);

    //!read one component of the impedance tensor from a netcdf file
    /*!this is an internal helper function
     * @param NetCDFFile The name for the netcdf file
     * @param Impedances The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param CompName Name of the component to be written to netcdf file ('Zxx_re' ... 'Zyy_im')
     * @param compindex component index (0-7)
     */
    J3DEXPORT void ReadImpedanceComp(netCDF::NcFile &NetCDFFile,
        std::vector<double> &Impedances, const std::string &CompName,
        const size_t compindex, const bool MustExist = true);

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
        const std::vector<double> &Impedances, const std::vector<double> &Errors,
        const std::vector<double> &Distortion, const std::vector<std::string> &Names,
        const std::vector<double> &RotationAngles);

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
        std::vector<double> &Impedances, std::vector<double> &ImpError,
        std::vector<double> &Distortion, std::vector<std::string> &Names,
        std::vector<double> &RotationAngles);

    //! Write tipper to a netcdf file
    /*! We can save tipper for several stations in a netcdf file for storage
     * and analysis with external programs.
     * @param filename The name for the netcdf file
     * @param Frequencies The vector of frequencies in Hz for the tipper
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the tipper in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the tipper in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the tipper in m
     * @param Tipper The tipper (dimensionless) as a vector of real numbers.
     *        4 consecutive elements form the tipper vectore for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Errors Optional parameter containing the errors with the same number of Elements as tipper.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     */
    J3DEXPORT void WriteTipperToNetCDF(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<int> &HxIndices, const std::vector<int> &HyIndices,
        const std::vector<int> &HzIndices, const std::vector<double> &Tipper,
        const std::vector<double> &Errors, const std::vector<std::string> &Names,
        const std::vector<double> &RotationAngles);

    //! Read magnetotelluric tipper from a netcdf file
    /*! Read MT tipper for several stations from a netcdf file.
     * @param filename The name for the netcdf file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param Tipper The tipper (dimensionless) as a vector of real numbers.
     *        4 consecutive elements form the tipper vector for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Error The error  has the same number of elements as tipper
     *         but contains all zeros if no error information present in the file.
     *         Also the error for the real and imaginary parts are always identical.
     */
    J3DEXPORT void ReadTipperFromNetCDF(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        std::vector<int> &HxIndices, std::vector<int> &HyIndices,
        std::vector<int> &HzIndices, std::vector<double> &Tipper,
        std::vector<double> &Error, std::vector<std::string> &Names,
        std::vector<double> &RotationAngles);

    //! A very basic routine to read impedances at a single site from a .mtt file in the format used by University of Goettingen
    /*! A very basic routine to read impedances at a single site from a .mtt file in the
     * format used by University of Goettingen.
     * @param filename The name of the .mtt file
     * @param Frequencies The frequencies contained in the file
     * @param Impedances The impedances in the same convention as above
     * @param Errors The impedance errors recorded in the file, will have the same number of elements as Impedances
     *        but the values for the real and imaginary part of each element are always identical
     * @param Tipper The tipper for the same frequencies as the impedances, convention Txr,Txi,Tyr,Tyi,... with same convention as impedances, i.e. frequency varies slowest
     * @param TippErr The error for the tipper
     */
    J3DEXPORT void ReadImpedancesFromMTT(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &Impedances,
        std::vector<double> &Errors, std::vector<double> &Tipper,
        std::vector<double> &TippErr);

    //! A very basic routine to write impedances for several sites to  .mtt files in the format used by University of Goettingen.
    /*! A very basic routine to write impedances for several sites to  .mtt files in the format used by University of Goettingen.
     * @param filenamebase The start of the name of the .mtt file, each file will get a number appended
     * @param Frequencies The frequencies contained in the file
     * @param Imp The impedances in the same convention as above
     * @param Err The impedance errors recorded in the file, will have the same number of elements as Impedances
     *        but the values for the real and imaginary part of each element are always identical
     * @param Tipper The tipper for the same frequencies as the impedances, convention Txr,Txi,Tyr,Tyi,... with same convention as impedances, i.e. frequency varies slowest
     * @param TippErr The error for the tipper
     */
    J3DEXPORT void WriteImpedancesToMtt(const std::string &filenamebase,
        const std::vector<double> &Frequencies, const std::vector<double> &Imp,
        const std::vector<double> &Err, const std::vector<double> &Tipper,
        const std::vector<double> &TipErr);

    J3DEXPORT void WriteImpedancesToMtt(const std::string &filenamebase,
        const std::vector<double> &Frequencies, const std::vector<double> &Imp,
        const std::vector<double> &Err);

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
        std::vector<double> &Imp, std::vector<double> &Err);

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
    void WriteAppResToAscii(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<double> &Imp, const std::vector<double> &Err);

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
     * @param StatNames The name of each station
     */
    J3DEXPORT void ReadImpedancesFromModEM(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        std::vector<double> &Imp, std::vector<double> &Err,
        std::vector<std::string> &StatNames);

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
        const std::vector<double> &Imp, const std::vector<double> &Err,
        const std::vector<std::string> &Names);

    J3DEXPORT void ReadTipperFromModEM(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &StatXCoord,
        std::vector<double> &StatYCoord, std::vector<double> &StatZCoord,
        std::vector<double> &Tip, std::vector<double> &Err);

    J3DEXPORT void WriteTipperToModEM(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<double> &Tip, const std::vector<double> &Err,
        const std::vector<std::string> &Names);

    //! Read magnetotelluric impedances from a .j file
    /*!
     * @param filename The name of the .j file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinate (North) of the measurement station for the impedances
     * @param StatYCoord The y-coordinate (East) of the measurement stations for the impedances
     * @param StatZCoord The z-coordinate (Down) of the measurement stations for the impedances
     * @param Imp The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Err Impedance errors with the same number of Elements with Impedances.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     *
     * @param Additional weighting information that can be used to adjust errors
     * @param Tipper The tipper for the same frequencies as the impedances, convention Txr,Txi,Tyr,Tyi,... with same convention as impedances, i.e. frequency varies slowest
     * @param TippErr The error for the tipper
     */
    J3DEXPORT void ReadImpedancesFromJ(const std::string &filename,
        std::vector<double> &Frequencies, double &StatXCoord, double &StatYCoord,
        double &StatZCoord, std::vector<double> &Imp, std::vector<double> &Err, std::vector<double> &Weight,  double &RotAngle);
    J3DEXPORT void ReadTipperFromJ(const std::string &filename,
        std::vector<double> &Frequencies, double &StatXCoord, double &StatYCoord,
        double &StatZCoord, std::vector<double> &Tipper, std::vector<double> &TippErr, double &RotAngle);
    void WriteImpedancesToJ(const std::string &filenamebase,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<double> &Imp, const std::vector<double> &Err);

    J3DEXPORT void ReadImpedancesFromXML(const std::string &filename,
        std::vector<double> &Frequencies, double &StatXCoord, double &StatYCoord,
        double &StatZCoord, std::vector<double> &Imp, std::vector<double> &Err, std::vector<double> &Weight, std::vector<double> &Tipper, std::vector<double> &TippErr,  double &RotAngle);

  /* @} */
  }

#endif /* READWRITEIMPEDANCES_H_ */
