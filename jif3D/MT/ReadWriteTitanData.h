//============================================================================
// Name        : ReadWriteTitanData.h
// Author      : Feb 9, 2017
// Version     :
// Copyright   : 2017, aavdeeva
//============================================================================

#ifndef READWRITETITANDATA_H_
#define READWRITETITANDATA_H_

#include "../Global/Serialization.h"
#include "../Global/Jif3DGlobal.h"
#include "../Global/VecMat.h"
#include <vector>
namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */

    //! Write Titan24 Data transfer functions to a netcdf file
    /*! We can save Titan24 Data for several stations in a netcdf file for storage
     * and analysis with external programs.
     * @param filename The name for the netcdf file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param MeasXCoord The x-coordinates (North) of the measurement stations for all Ex, Ey or H fields in m
     * @param MeasYCoord The y-coordinates (East) of the measurement stations for all Ex, Ey or H fields in m
     * @param MeasZCoord The z-coordinates (Down) of the measurement stations for all Ex, Ey or H fields in m
     * @param ExIndices The indices of the measurements corresponding to Ex fields in involved in impedance computation,
     * these are integer vec of size nFreqs x nSites, the frequencies vary slowest.
     * @param EyIndices The indices of the measurements corresponding to Ey fields in involved in impedance computation,
     * these are integer vec of size nFreqs x nSites, the frequencies vary slowest.
     * @param HIndices The indices of the measurements corresponding to H fields in involved in impedance computation,
     * these are integer vec of size nFreqs x nSites, the frequencies vary slowest.
     * @param Impedances The Titan24 transfer functions (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Errors Optional parameter containing errors of the Titan24 transfer functions with the same number of Elements with Impedances.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     */
    J3DEXPORT void WriteTitanDataToNetCDF(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &MeasXCoord,
        const std::vector<double> &MeasYCoord, const std::vector<double> &MeasZCoord,
        const std::vector<int> &ExIndices, const std::vector<int> &EyIndices,
        const std::vector<int> &HIndices, const std::vector<double> &Impedances,
        const std::vector<double> &Errors, const std::vector<double> &Distortion,
        const std::vector<double> &RotAngles);

    //! Read Titan24 Data transfer functions from a netcdf file
    /*! Read Titan24 Data for several stations from a netcdf file.
     * @param filename The name for the netcdf file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param MeasXCoord The x-coordinates (North) of the measurement stations for all Ex, Ey or H fields in m
     * @param MeasYCoord The y-coordinates (East) of the measurement stations for all Ex, Ey or H fields in m
     * @param MeasZCoord The z-coordinates (Down) of the measurement stations for all Ex, Ey or H fields in m
     * @param ExIndices The indices of the measurements corresponding to Ex fields in involved in impedance computation,
     * these are integer vec of size nFreqs x nSites, the frequencies vary slowest.
     * @param EyIndices The indices of the measurements corresponding to Ey fields in involved in impedance computation,
     * these are integer vec of size nFreqs x nSites, the frequencies vary slowest.
     * @param HIndices The indices of the measurements corresponding to H fields in involved in impedance computation,
     * these are integer vec of size nFreqs x nSites, the frequencies vary slowest.
     * @param Impedances The Titan24 transfer functions (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param ImpError The error of the Titan24 transfer functions, has the same number of elements as Impedances
     *         but contains all zeros if no error information present in the file.
     *         Also the error for the real and imaginary parts are always identical.
     */
    J3DEXPORT void ReadTitanDataFromNetCDF(const std::string &filename,
        std::vector<double> &Frequencies, std::vector<double> &MeasXCoord,
        std::vector<double> &MeasYCoord, std::vector<double> &MeasZCoord,
        std::vector<int> &ExIndices, std::vector<int> &EyIndices,
        std::vector<int> &HIndices, std::vector<double>  &Impedances, std::vector<double>  &ImpError,
        std::vector<double> &Distortion, std::vector<double> &RotAngles);

    //! Write Titan24 transfer functions into an ascii file as written by ModEM
    /*! Write Titan24 transfer functions into an ascii file as written by ModEM,
     * the coordinates of the sites are the coordinates of Ex measurements
     * @param filename The name of the ascii file
     * @param Frequencies The vector of frequencies in Hz for the impedances in the vector Impedances
     * @param StatXCoord The x-coordinates (North) of the measurement stations for the impedances in m
     * @param StatYCoord The y-coordinates (East) of the measurement stations for the impedances in m
     * @param StatZCoord The z-coordinates (Down) of the measurement stations for the impedances in m
     * @param ExIndices The indices of the measurements corresponding to Ex fields in involved in impedance computation,
     * these are integer vec of size nFreqs x nSites, the frequencies vary slowest.
     * @param Imp The impedances (in Ohm, i.e. E/H) as a vector of real numbers.
     *        8 consecutive elements form the impedance matrix for one frequency and site,
     *        all impedances for one frequency and all stations form a contiguous block, the frequencies vary slowest.
     * @param Err Impedance errors with the same number of Elements with Impedances.
     *        As we only have one error estimate per element we write only the components corresponding to the real parts.
     */
    J3DEXPORT void WriteTitanDataToModEM(const std::string &filename,
        const std::vector<double> &Frequencies, const std::vector<double> &StatXCoord,
        const std::vector<double> &StatYCoord, const std::vector<double> &StatZCoord,
        const std::vector<int> &ExIndices, const jif3D::rvec &Imp,
        const jif3D::rvec &Err);

  /* @} */
  }

#endif /* READWRITETITANDATA_H_ */
