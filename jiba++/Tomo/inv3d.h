/* Local libary of "inv3d": */
#include <cstdlib>
#include <vector>
#ifndef _INV3D_
#define _INV3D_

/*-------------------------------------------------------------*/
/* Define "Structures": */
namespace jiba
  {

    /*! Parameters of the grid structure:*/
    class GRID_STRUCT
      {
    public:
      size_t nx; /*!< Number of grid cell center in x-direction*/
      size_t ny; /*!< Number of grid cell center in y-direction*/
      size_t nz; /*!< Number of grid cell center in z-direction*/

      float h; /*!< Size of the grid cells in m (!The same in all three directions!)*/

      float org[3]; /*!< Components of the origin (center of the first cell in the upper left corner) of the grid in m*/

      /*!< Seismic Velocity parameters*/
      size_t nborder; /*!< Number of grid cells for the boundary*/

      std::vector<float> slow; /*!< Slowness model used for the forward model (normalized by the grid cell size)*/
      GRID_STRUCT() :
        nx(0), ny(0), nz(0), h(0), org(), nborder(0), slow()
        {
        }
      virtual ~GRID_STRUCT()
        {
        }
      };

    /*! Geometry of the  the seismic shots and receivers */

    class GEOMETRY
      {
    public:
      std::vector<float> x, y, z; /*!< Positions of the shot/receiver locations and fixpoints in m*/
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

    /*! Structure including informations about the seismic data*/
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

    /*! Structure for the ray paths*/
    class RP_STRUCT
      {
    public:
      size_t n; /*!< Number specifying the row in the matrix for the shot-receiver combination*/
      size_t nray; /*!< Number of segments for the ray*/
      std::vector<double> x; /*!< 		Forward grid: x-position of the begin of the ray-path segment*/
      /*!< BUT:	Inversion grid: x-component of the ray in the cell*/
      std::vector<double> y; /*!< 		Forward grid: y-position of the begin of the ray-path segment*/
      /*!< BUT:	Inversion grid: y-component of the ray in the cell*/
      std::vector<double> z; /*!< 		Forward grid: z-position of the begin of the ray-path segment*/
      /*!< BUT:	Inversion grid: x-component of the ray in the cell*/
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

  }

#endif
