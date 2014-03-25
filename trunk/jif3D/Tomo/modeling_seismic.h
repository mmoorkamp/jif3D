#ifndef MODELING_SEISMIC_H
#define MODELING_SEISMIC_H

#include <vector>
#include <boost/serialization/serialization.hpp>

/*! \file
 * This file contains some data structures and function definitions for tomography modeling.
 * These are based on the code by B. Heincke with some modifications to avoid c-style allocations
 * and screen output.
 *  */
/*--------------------------------------------------------------*/
/* Define function:*/

namespace jif3D
  {
    /** \addtogroup tomo Seismic tomography classes and functions */
    /* @{ */

    //! Parameters of the grid structure:
    class GRID_STRUCT
      {
    private:
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & nx;
          ar & ny;
          ar & nz;
          ar & h;
        }
    public:
      std::size_t nx; /*!< Number of grid cell center in x-direction*/
      std::size_t ny; /*!< Number of grid cell center in y-direction*/
      std::size_t nz; /*!< Number of grid cell center in z-direction*/

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
      std::size_t nshot; /*!< Number of shot positions*/
      std::size_t nrec; /*!< Number of receiver positions*/
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & x;
          ar & y;
          ar & z;
          ar & nshot;
          ar & nrec;
        }
      GEOMETRY() :
          x(), y(), z(), nshot(0), nrec(0)
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
      std::size_t ndata_seis; /*!< Number of picked first-breaks*/
      std::size_t ndata_seis_act; /*!< Number of rays that could be traced back*/

      std::vector<std::size_t> sno; /*!< List of the shot position numbers of the traces for which the first breaks were picked*/
      std::vector<std::size_t> rno; /*!< List of the receiver position numbers of the traces for which the first breaks were picked*/

      std::vector<double> tcalc; /*!< Calculated travel times for the different shot-receiver combinations in ms*/

      std::vector<int> lshots; /*!< Number of active receivers for the corresponding shot position number (related to shots)*/
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & ndata_seis;
          ar & ndata_seis_act;
          ar & sno;
          ar & rno;
          ar & tcalc;
          ar & lshots;
        }
      DATA_STRUCT() :
          ndata_seis(0), ndata_seis_act(0), sno(), rno(), tcalc(), lshots()
        {
        }
      virtual ~DATA_STRUCT()
        {
        }
      };

    //! Structure for the ray paths
    class RP_STRUCT
      {
    private:
      friend class boost::serialization::access;
      //! Provide serialization to be able to store objects and, more importantly for simpler MPI parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & n;
          ar & nray;
          ar & x;
          ar & y;
          ar & z;
          ar & len;
          ar & ele;
        }
    public:
      std::size_t n; /*!< Number specifying the row in the matrix for the shot-receiver combination*/
      std::size_t nray; /*!< Number of segments for the ray*/
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
    int ForwardModRay(const GEOMETRY &geo, const GRID_STRUCT &grid, DATA_STRUCT &data,
        std::vector<RP_STRUCT> &raypath);
    int ForwardModShot(int i, const std::vector<size_t> &uniqueshots, const GEOMETRY &geo, const GRID_STRUCT &grid, DATA_STRUCT &data,
    		std::vector<RP_STRUCT> &raypath);
    /*in podvin-lecomte-3D.c*/
    float interpolate(float x, float y, float z, const GRID_STRUCT &grid, float *data);
  /* @} */
  }
#endif
