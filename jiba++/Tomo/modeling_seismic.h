#ifndef MODELING_SEISMIC_H
#define MODELING_SEISMIC_H

#include <vector>

/*! \file
 * This file contains some data structures and function definitions for tomography modeling.
 * These are based on the code by B. Heincke with some modifications to avoid c-style allocations
 * and screen output.
 *  */
/*--------------------------------------------------------------*/
/* Define function:*/

namespace jiba
  {
    /** \addtogroup tomo Seismic tomography classes and functions */
    /* @{ */

    //! Parameters of the grid structure:
    class GRID_STRUCT
      {
    public:
      size_t nx; /*!< Number of grid cell center in x-direction*/
      size_t ny; /*!< Number of grid cell center in y-direction*/
      size_t nz; /*!< Number of grid cell center in z-direction*/

      float h; /*!< Size of the grid cells in m (!The same in all three directions!)*/

      /*!< Seismic Velocity parameters*/
      std::vector<float> slow; /*!< Slowness model used for the forward model (normalized by the grid cell size)*/
      GRID_STRUCT() :
        nx(0), ny(0), nz(0), h(0), slow()
        {
        }
      virtual ~GRID_STRUCT()
        {
        }
      };

    //! Geometry of the  the seismic shots and receivers
    class GEOMETRY
      {
    public:
      std::vector<float> x; /*!< x-coordinates of the shot/receiver locations  in m*/
      std::vector<float> y; /*!< y-coordinates of the shot/receiver locations  in m*/
      std::vector<float> z; /*!< z-coordinates of the shot/receiver locations in m*/
      size_t nshot; /*!< Number of shot positions*/
      size_t nrec; /*!< Number of receiver positions*/
      GEOMETRY() :
        x(NULL), y(NULL), z(NULL), nshot(0), nrec(0)
        {
        }
      virtual ~GEOMETRY()
        {
        }
      };

    //! Structure including informations about the seismic data
    class DATA_STRUCT
      {
    public:
      /*!< Seismic parameters*/
      size_t ndata_seis; /*!< Number of picked first-breaks*/
      size_t ndata_seis_act; /*!< Number of rays that could be traced back*/

      std::vector<size_t> sno; /*!< List of the shot position numbers of the traces for which the first breaks were picked*/
      std::vector<size_t> rno; /*!< List of the receiver position numbers of the traces for which the first breaks were picked*/

      std::vector<double> tcalc; /*!< Calculated travel times for the different shot-receiver combinations in ms*/

      std::vector<int> lshots; /*!< Number of active receivers for the corresponding shot position number (related to shots)*/
      DATA_STRUCT() :
        ndata_seis(0), ndata_seis_act(0), sno(NULL), rno(NULL), tcalc(NULL),
            lshots(NULL)
        {
        }
      virtual ~DATA_STRUCT()
        {
        }
      };

    //! Structure for the ray paths
    class RP_STRUCT
      {
    public:
      size_t n; /*!< Number specifying the row in the matrix for the shot-receiver combination*/
      size_t nray; /*!< Number of segments for the ray*/
      std::vector<double> x; /*!<               Forward grid: x-position of the begin of the ray-path segment*/
      /*!< BUT: Inversion grid: x-component of the ray in the cell*/
      std::vector<double> y; /*!<               Forward grid: y-position of the begin of the ray-path segment*/
      /*!< BUT: Inversion grid: y-component of the ray in the cell*/
      std::vector<double> z; /*!<               Forward grid: z-position of the begin of the ray-path segment*/
      /*!< BUT: Inversion grid: x-component of the ray in the cell*/
      std::vector<double> len; /*!< Length of the ray-path segment*/
      std::vector<long> ele; /*!< Grid cell position in the grid*/
      RP_STRUCT() :
        n(0), nray(0), x(), y(), z(), len(), ele()
        {
        }
      virtual ~RP_STRUCT()
        {
        }
      };
    //! The basic forward modeling function, uses the Podvin-Lecomte algorithm to calculate traveltimes and then constructs rays from these times
    int ForwardModRay(const GEOMETRY &geo, const GRID_STRUCT &grid,
        DATA_STRUCT *data, RP_STRUCT *raypath);
    /*in podvin-lecomte-3D.c*/
    float interpolate(float x, float y, float z, const GRID_STRUCT &grid,
        float *data);
  /* @} */
  }
#endif