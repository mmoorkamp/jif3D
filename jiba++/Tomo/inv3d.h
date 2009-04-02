/* Local libary of "inv3d": */
#include <cstdlib>

#ifndef _INV3D_
#define _INV3D_

/*-------------------------------------------------------------*/
/* Define "Structures": */

/*! Parameters of the grid structure:*/
#ifndef  _GRID_STRUCT_
typedef struct _GRID_STRUCT_
  {
  int nx; /*!< Number of grid cell center in x-direction*/
  int ny; /*!< Number of grid cell center in y-direction*/
  int nz; /*!< Number of grid cell center in z-direction*/

  float h; /*!< Size of the grid cells in m (!The same in all three directions!)*/

  float org[3]; /*!< Components of the origin (center of the first cell in the upper left corner) of the grid in m*/

  int *border_index; /*!< grid cell belongs to the border or to the "air"(yes=0,no=1)*/

  /*!< Seismic Velocity parameters*/
  int nborder; /*!< Number of grid cells for the boundary*/

  double *slow; /*!< Slowness model used for the forward model (normalized by the grid cell size)*/

  } GRID_STRUCT;
#define  _GRID_STRUCT_
#endif

/*! Geometry of the topography points, the seismic shots and receivers, the gravimetric and MT stations*/
#ifndef  _GEOMETRY_
typedef struct _GEOMETRY_
  {
  float *x, *y, *z; /*!< Positions of the shot/receiver locations and fixpoints in m*/
  int nshot; /*!< Number of shot positions*/
  int nrec; /*!< Number of receiver positions*/
  } GEOMETRY;
#define  _GEOMETRY_
#endif

/*! Structure including informations about the seismic data*/
#ifndef  _DATA_STRUCT_
typedef struct _DATA_STRUCT_
  {
  /*!< Seismic parameters*/
  long ndata_seis; /*!< Number of picked first-breaks*/
  long ndata_seis_act; /*!< Number of rays that could be traced back*/

  int *sno; /*!< List of the shot position numbers of the traces for which the first breaks were picked*/
  int *rno; /*!< List of the receiver position numbers of the traces for which the first breaks were picked*/

  double *tcalc; /*!< Calculated travel times for the different shot-receiver combinations in ms*/
  /*!< REMARK: sno, rno, tobs and tcalc have the same number of elements (ndata_seis) and are linked to each other*/
  float *xdist; /*!< Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/
  float *ydist; /*!< Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/
  float *zdist; /*!< Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/

  int *shots; /*!< List of all shot positions used (number of elements: nshot)*/
  int *recs; /*!< List of all receiver positions used (number of elements: nrec)*/

  int *lshots; /*!< Number of active receivers for the corresponding shot position number (related to shots)*/
  int *lrecs; /*!< Number of used shots for the corresponding receicer position number (related to recs)*/
  } DATA_STRUCT;
#define  _DATA_STRUCT_
#endif

/*! Structure for the ray paths*/
#ifndef _RP_STRUCT_
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
#define _RP_STRUCT_
#endif

/*! Structure for the fat-rays*/
#ifndef _F_RP_STRUCT_
typedef struct _F_RP_STRUCT_
  {
  long n; /*!< Number specifying the row in the matrix for the shot-receiver combination*/
  long ncell; /*!< Number of cells which contribute to the fat ray*/
  float fatthres; /*!< Threshold for the fat-ray in ms*/
  float fatb; /*!< Exponent governing the decrease of the weight with increasing "difference"-time in 1/ms*/
  double *weight; /*!< Normalized contribution of the cell to the fat-ray*/
  long *ele; /*!< Grid cell position in the grid*/

  } F_RP_STRUCT;
#define _F_RP_STRUCT_
#endif

/*! Structure to organize the cell parameters during back tracing the rays*/
#ifndef _CELL_STRUCT_
typedef struct _CELL_STRUCT_
  {
  int xno; /*!< position number in x-direction of the cell in the grid*/
  int yno; /*!< position number in y-direction of the cell in the grid*/
  int zno; /*!< position number in z-direction of the cell in the grid*/
  int dirx_i; /*!< The ray runs in negative x-direction into the cell = 1; the ray runs in positive x-direction into the cell = 1; else=0 */
  int diry_i; /*!< The ray runs in negative y-direction into the cell = 1; the ray runs in positive y-direction into the cell = 1; else=0 */
  int dirz_i; /*!< The ray runs in negative z-direction into the cell = 1; the ray runs in positive z-direction into the cell = 1; else=0 */
  double xpos; /*!< The position of the starting point of the ray in x-direction (normalized by the position in the grid)*/
  double ypos; /*!< The position of the starting point of the ray in y-direction (normalized by the position in the grid)*/
  double zpos; /*!< The position of the starting point of the ray in z-direction (normalized by the position in the grid)*/

  } CELL_STRUCT;
#define _CELL_STRUCT_
#endif

#endif
