/* Local libary of "inv3d": */
#include <cstdlib>

#ifndef _INV3D_
#define _INV3D_

/*-------------------------------------------------------------*/
/* Define "Structures": */
namespace jiba
  {

    /*! Parameters of the grid structure:*/
    typedef struct _GRID_STRUCT_
      {
      int nx; /*!< Number of grid cell center in x-direction*/
      int ny; /*!< Number of grid cell center in y-direction*/
      int nz; /*!< Number of grid cell center in z-direction*/

      float h; /*!< Size of the grid cells in m (!The same in all three directions!)*/

      float org[3]; /*!< Components of the origin (center of the first cell in the upper left corner) of the grid in m*/

      /*!< Seismic Velocity parameters*/
      int nborder; /*!< Number of grid cells for the boundary*/

      double *slow; /*!< Slowness model used for the forward model (normalized by the grid cell size)*/

      } GRID_STRUCT;

    /*! Geometry of the  the seismic shots and receivers */

    typedef struct _GEOMETRY_
      {
      float *x, *y, *z; /*!< Positions of the shot/receiver locations and fixpoints in m*/
      int nshot; /*!< Number of shot positions*/
      int nrec; /*!< Number of receiver positions*/
      } GEOMETRY;

    /*! Structure including informations about the seismic data*/
    typedef struct _DATA_STRUCT_
      {
      /*!< Seismic parameters*/
      long ndata_seis; /*!< Number of picked first-breaks*/
      long ndata_seis_act; /*!< Number of rays that could be traced back*/

      int *sno; /*!< List of the shot position numbers of the traces for which the first breaks were picked*/
      int *rno; /*!< List of the receiver position numbers of the traces for which the first breaks were picked*/

      double *tcalc; /*!< Calculated travel times for the different shot-receiver combinations in ms*/

      int *lshots; /*!< Number of active receivers for the corresponding shot position number (related to shots)*/
      int *lrecs; /*!< Number of used shots for the corresponding receicer position number (related to recs)*/
      } DATA_STRUCT;

    /*! Structure for the ray paths*/
    typedef struct _RP_STRUCT_
      {
      long n; /*!< Number specifying the row in the matrix for the shot-receiver combination*/
      long nray; /*!< Number of segments for the ray*/
      double *x; /*!< 		Forward grid: x-position of the begin of the ray-path segment*/
      /*!< BUT:	Inversion grid: x-component of the ray in the cell*/
      double *y; /*!< 		Forward grid: y-position of the begin of the ray-path segment*/
      /*!< BUT:	Inversion grid: y-component of the ray in the cell*/
      double *z; /*!< 		Forward grid: z-position of the begin of the ray-path segment*/
      /*!< BUT:	Inversion grid: x-component of the ray in the cell*/
      double *len; /*!< Length of the ray-path segment*/
      long *ele; /*!< Grid cell position in the grid*/

      } RP_STRUCT;


  }

#endif
