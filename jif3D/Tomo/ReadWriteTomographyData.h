//============================================================================
// Name        : ReadWriteTomographyData.h
// Author      : May 5, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef READWRITETOMOGRAPHYDATA_H_
#define READWRITETOMOGRAPHYDATA_H_


#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"
#include <fstream>
#include "../Tomo/tomo_types.h"
#include "../Tomo/ThreeDSeismicModel.h"

namespace jif3D
  {


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
    J3DEXPORT void PlotTimeField(const std::string &filename, const floattype *Times,
        const double gridspacing, const size_t nx, const size_t ny, const size_t nz)
      {
        std::ofstream outfile(filename.c_str());
        outfile << "# vtk DataFile Version 2.0\n";
        outfile << "Traveltimes\n";
        outfile << "ASCII\n";
        outfile << "DATASET STRUCTURED_POINTS\n";
        outfile << "DIMENSIONS " << nx << " " << ny << " " << nz << std::endl;
        outfile << "ORIGIN 0 0 0\n";
        outfile << "SPACING " << gridspacing << " " << gridspacing << " " << gridspacing
            << std::endl;
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
     * @param minxindex The shift in number of grid cells between the actual model and the forward grid in x-direction
     * @param minyindex The shift in number of grid cells between the actual model and the forward grid in y-direction
     */
    J3DEXPORT void PlotRaypath(const std::string &filename, const std::vector<RP_STRUCT> &raypath,
        const size_t nmeas, const double gridspacing, const size_t nairlayers, int minxindex, int minyindex);
  /* @} */
  }
#endif /* READWRITETOMOGRAPHYDATA_H_ */
