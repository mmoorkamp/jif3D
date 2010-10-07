//============================================================================
// Name        : ReadWriteTomographyData.h
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef READWRITETOMOGRAPHYDATA_H_
#define READWRITETOMOGRAPHYDATA_H_

#include <fstream>
#include "../Global/VecMat.h"
#include "ThreeDSeismicModel.h"
#include "modeling_seismic.h"

namespace jiba
  {
    /** \addtogroup tomo Seismic tomography classes and functions */
    /* @{ */
    //!Save a collection of travel times and associated positions to a netcdf file
    /*! Saves traveltimes together with the source receiver configuration to a netcdf file.
     * @param filename The name of the file to store the data
     * @param Data The traveltimes in s, the length of this vector matches the number of
     *             source-receiver combinations in the Model object
     * @param Model The model object containing information about the source and receiver setup
     */
    void SaveTraveltimes(const std::string &filename, const jiba::rvec &Data,
        const jiba::ThreeDSeismicModel &Model);
    //!Read a collection of travel times and associated positions from a netcdf file
    /*! Read traveltimes together with the source receiver configuration from a netcdf file.
     * @param filename The name of the file to store the data
     * @param Data The traveltimes in s, the length of this vector matches the number of
     *         source-receiver combinations in the Model object
     * @param Model The model object containing information about the source and receiver setup
     */
    void ReadTraveltimes(const std::string &filename, jiba::rvec &Data,
        jiba::ThreeDSeismicModel &Model);

    //! Plot the travel time field, i.e. the travel time at each grid cell in a vtk file
    /*! Write the travel time field to a .vtk file. This is a template so that we can
     * use it with both the original Podvin code as well as our own class.
     * @param filename The filename to write the time field to
     * @param Times The travel time for each model cell as a 1D vector has to have nx*ny*nz elements, nz varies fastest
     * @param gridspacing The spacing in m between adjacent grid cells, has to be constant throughout the grid
     * @param nx The number of cells in x-direction
     * @param ny The number of cells in y-direction
     * @param nz The number of cells in z-direction
     */
    template<typename floattype>
    void PlotTimeField(const std::string &filename, const floattype *Times,
        const double gridspacing, const size_t nx, const size_t ny,
        const size_t nz)
      {
        std::ofstream outfile(filename.c_str());
        outfile << "# vtk DataFile Version 2.0\n";
        outfile << "Traveltimes\n";
        outfile << "ASCII\n";
        outfile << "DATASET STRUCTURED_POINTS\n";
        outfile << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
        outfile << "ORIGIN 0 0 0\n";
        outfile << "SPACING " << gridspacing << " " << gridspacing << " "
            << gridspacing << std::endl;
        outfile << "POINT_DATA " << nx * ny * nz << std::endl;
        outfile << "SCALARS traveltimes float\n";
        outfile << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < nz; ++i)
          for (size_t j = 0; j < ny; ++j)
            {
              for (size_t k = 0; k < nx; ++k)
                {
                  outfile << Times[k * (nz * ny) + j * nz + i] << " ";
                }
              outfile << std::endl;
            }
      }
    //! Plot the raypaths for a 3D forward calculation in a .vtk file
    /*! Plot the rays as calculated by the tomographic forward code. This uses
     * the low level information and internal structures of the TomographyCalculator class.
     * @param filename The name of the file to store the raypath information
     * @param raypath The structure containing the information about the raypaths
     * @param nmeas The number of measurements
     * @param gridspacing The spacing of the grid cells
     * @param nairlayers The number of airlayers
     */
    void PlotRaypath(const std::string &filename, jiba::RP_STRUCT *raypath,
        const size_t nmeas, const double gridspacing, const size_t nairlayers);
  /* @} */
  }
#endif /* READWRITETOMOGRAPHYDATA_H_ */
