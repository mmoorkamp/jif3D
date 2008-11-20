/* Local libary of "inv3d": */
#include <cstdlib>

#ifndef _INV3D_
#define _INV3D_

/*-------------------------------------------------------------*/
/* Define "Structures": */


/*! Parameters of the initial velocity/gravity or resistivity model:*/
#ifndef  _GRADIENT_
typedef struct _GRADIENT_
  {
  short index_topo; /*!< This index governs if the velocity/density/resistivity field is determined by the surface topography(1) or */
  /*!< by a defined gridpoint "org"*/

  float value_0; /*!< Starting value (velocity[m/s], density[g/cm^3] or resistivity[ohmm]) (at the defined gridpoint "org" or at the surface)*/
  float min, max; /*!< Minimum and maximum value (velocity[m/s], density[g/cm^3], resistivity[ohmm]) for the starting model*/
  float org[3]; /*!< Position, where "value_o" is defined*/
  float g[3]; /*!< Gradient for velocity[(m/s)/m],density[g/cm^3] or resistivity[ohmm/m]*/

  } GRADIENT;
#define  _GRADIENT_
#endif

/*! Parameters of the grid structure:*/
#ifndef  _GRID_STRUCT_
typedef struct _GRID_STRUCT_
  {
  int nx; /*!< Number of grid cell center in x-direction*/
  int ny; /*!< Number of grid cell center in y-direction*/
  int nz; /*!< Number of grid cell center in z-direction*/

  float h; /*!< Size of the grid cells in m (!The same in all three directions!)*/

  float org[3]; /*!< Components of the origin (center of the first cell in the upper left corner) of the grid in m*/

  short topo_index; /*!< The index governs if the topography will be ignored(0) or introduced by means of triangulation(1)*/

  int grid_cell_factor[3]; /*!< Sizes of "inversion" grid cell dimensions compared to the forward" grid cell dimensions*/

  int *border_index; /*!< grid cell belongs to the border or to the "air"(yes=0,no=1)*/

  /*!< Seismic Velocity parameters*/
  float vborder; /*!< velocity of the boundaries m/s*/
  int nborder; /*!< Number of grid cells for the boundary*/

  float v_air_water; /*!< velocity [m/s] in the air or water layer*/

  double *slow; /*!< Slowness model used for the forward model (normalized by the grid cell size)*/
  GRADIENT grad_vel; /*!< Structure including the parameters of the initial (gradient) velocity model*/

  /*!< Density parameters*/
  float dens_air_water; /*!< density [g/cm^3] of the air or the water layer*/

  } GRID_STRUCT;
#define  _GRID_STRUCT_
#endif

/*! Geometry of the topography points, the seismic shots and receivers, the gravimetric and MT stations*/
#ifndef  _GEOMETRY_
typedef struct _GEOMETRY_
  {
  int nstat; /*!< Number of receiver/shot stations and fixpoints*/
  int *coor_info; /*!< Coordinate information used for the topography (yes=1,no=0)*/
  float *x, *y, *z; /*!< Positions of the shot/receiver locations and fixpoints in m*/
  float eps; /*!< Maximum distance between two coordinates location in m so
   that they are assigned to the same shot resp. receiver location */
  /*!<  seismic parameters*/
  int nshot; /*!< Number of shot positions*/
  int nrec; /*!< Number of receiver positions*/
  } GEOMETRY;
#define  _GEOMETRY_
#endif

/*! Structure with the topography information*/
#ifndef _TOPO_STRUCT_
typedef struct _TOPO_STRUCT_
  {
  long nr_of_topo_points; /*!< Number of topography points*/
  long nr_of_triangles; /*!< Number of triangles after the delaunay triangulation*/

  double *x; /*!< x coordinates of the topography points*/
  double *y; /*!< y coordinates of the topography points*/
  double *z; /*!< z coordinates of the topography points*/

  int *index; /*!< Index of the topography points that specify the corners of triangles (first triangle = first three entries)*/

  } TOPO_STRUCT;
#define  _TOPO_STRUCT_
#endif

/*! Structure including informations about the seismic data*/
#ifndef  _DATA_STRUCT_
typedef struct _DATA_STRUCT_
  {
  /*!< Seismic parameters*/
  long ndata_seis; /*!< Number of picked first-breaks*/
  long ndata_seis_act; /*!< Number of rays that could be traced back*/

  double rms_seis; /*!< RMS value of the traveltime residuals in ms*/

  int *sno; /*!< List of the shot position numbers of the traces for which the first breaks were picked*/
  int *rno; /*!< List of the receiver position numbers of the traces for which the first breaks were picked*/
  double *tobs; /*!< Observed travel times for the different shot-receiver combinations in ms*/
  double *tcalc; /*!< Calculated travel times for the different shot-receiver combinations in ms*/
  /*!< REMARK: sno, rno, tobs and tcalc have the same number of elements (ndata_seis) and are linked to each other*/
  float *xdist; /*!< Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/
  float *ydist; /*!< Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/
  float *zdist; /*!< Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/

  int *shots; /*!< List of all shot positions used (number of elements: nshot)*/
  int *recs; /*!< List of all receiver positions used (number of elements: nrec)*/

  int *lshots; /*!< Number of active receivers for the corresponding shot position number (related to shots)*/
  int *lrecs; /*!< Number of used shots for the corresponding receicer position number (related to recs)*/

  double timeshift; /*!< Constant timeshifts of all traces in ms*/

  double *weigth_seis; /*!< Weighting of the rays for the inversion (1.0 = is the usually weigthing)*/
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


/*! Parameters for each inversion cell*/
#ifndef _BIDX_STRUCT_
typedef struct _BIDX_STRUCT_
  {
  double bdamp; /*!< Damping parameter*/
  double xo; /*!< x-coordinate of the center of the inversion cell in m*/
  double yo; /*!< y-coordinate of the center of the inversion cell in m*/
  double zo; /*!< x-coordinate of the center of the inversion cell in m*/
  double xdim; /*!< Dimension of the cell in x-direction in m*/
  double ydim; /*!< Dimension of the cell in y-direction in m*/
  double zdim; /*!< Dimension of the cell in z-direction in m*/

  int border_inv_cell_x_left; /*!< Specify, if the cell is located at the left x-border (yes==1;no==0)*/
  int border_inv_cell_x_right; /*!< Specify, if the cell is located at the right x-border (yes==1;no==0)*/

  int border_inv_cell_y_front; /*!< Specify, if the cell is located at the front y-border (yes==1;no==0)*/
  int border_inv_cell_y_back; /*!< Specify, if the cell is located at the back y-border (yes==1;no==0)*/

  int border_inv_cell_z; /*!< Specify, if the cell is located at the z-border (yes==1;no==0)*/

  double val_slow; /*!< slowness in the cell*/
  double val_dens; /*!< density in the cell*/
  double val_res; /*!< resistivity in the cell*/

  int use; /*!< inversion cell used for inversion (yes=1;no=0)*/
  long used_nr; /*!< number after renumbering the ACTIVE inversion cells*/

  int nele; /*!< number of forward-cells included in the inversion cell (ONLY the ones, which are activated for the inversion)*/
  long *ele; /*!< the positions of the forward cells in the grid*/

  int nele_inv; /*!< Number of former inversion cells in the "new" inversion cells*/
  long *ele_inv; /*!< number of the former inversion cells included in the "new" inversion cell*/

  } BIDX_STRUCT;
#define _BIDX_STRUCT_
#endif



#endif
