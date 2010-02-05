//============================================================================
// Name        : PodvinTime3D.h
// Author      : May 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef PODVINTIME3D_H_
#define PODVINTIME3D_H_

namespace jiba
  {
    /** \addtogroup tomo Seismic tomography classes and functions */
    /* @{ */
    //! This is the 3D forward algorithm by Podvin and Lecomte encapsulated in a C++ class.
    /*! We moved the Podvin and Lecomte code into a C++ class. This allows us to
     * replace global static variables with class variables. Now we can do several
     * calculations in parallel at the cost of separate memory for each class.
     *
     * Otherwise we made no changes to the functions or the data structures.
     */
    class PodvinTime3D
      {
    private:
      int pre_init(void), init_point(void), /* routine modified in Dec 2005 */
      recursive_init(void), propagate_point(int), x_side(int, int, int, int,
          int, int), y_side(int, int, int, int, int, int), z_side(int, int,
          int, int, int, int), scan_x_ff(int, int, int, int, int, int, int),
          scan_x_fb(int, int, int, int, int, int, int), scan_x_bf(int, int,
              int, int, int, int, int), scan_x_bb(int, int, int, int, int, int,
              int), scan_y_ff(int, int, int, int, int, int, int), scan_y_fb(
              int, int, int, int, int, int, int), scan_y_bf(int, int, int, int,
              int, int, int), scan_y_bb(int, int, int, int, int, int, int),
          scan_z_ff(int, int, int, int, int, int, int), scan_z_fb(int, int,
              int, int, int, int, int), scan_z_bf(int, int, int, int, int, int,
              int), scan_z_bb(int, int, int, int, int, int, int);
      /* the only fully commented "side" functions are x_side() and scan_x_ff() */

      void error(int), init_nearest(void), /* routine modified in Dec 2005 */
      init_cell(float, float, float, int, int, int), /* routine modified in Dec 2005 */
      free_ptrs(int);

      float
      /* new function init_cellh(): see patches 271205[1] and 271205[2] */
      init_cellh(float vh, float vv, float hsc, float hsn), exact_delay(float,
          float, float, int, int, int);

      int t_1d(int, int, int, float, float, float, float, float), t_2d(int,
          int, int, float, float, float, float), diff_2d(int, int, int, float,
          float, float), t_3d_(int, int, int, float, float, float, float,
          float, int), t_3d_part1(int, int, int, float, float, float, float),
          point_diff(int, int, int, float, float), edge_diff(int, int, int,
              float, float, float);
      int nmesh_x, nmesh_y, nmesh_z; /* Model dimensions (cells) */
      float ***hs, *hs_buf, /* 1D and 3D arrays */
      *hs_keep; /* to save boundary values */

      /* TIMEFIELD */

      int nx, ny, nz; /* Timefield dimensions (nodes) */
      float ***t, *t_buf; /* 1D and 3D arrays */
      float timeshift; /* required by "fuzzy tests" */
      /* for more comments, see init_point() */

      /* SOURCE */

      float fxs, fys, fzs; /* Point source coordinates */
      int xs, ys, zs; /* Nearest node */

      int messages, /* message flag (0:silent)              */
      source_at_node, /* are source coordinate int's ? (0/1)  */
      init_stage, /* level of recursivity during init.    */
      current_side_limit, /* actual boundary of computations      */
      X0, X1, Y0, Y1, Z0, Z1, /* inclusive boundaries of timed region */
      reverse_order, /* level of recursivity in FD scheme    */
      *longflags, /* local headwave flags.                */
      flag_fb, x_start_fb, y_start_fb, z_start_fb, flag_bf, x_start_bf,
          y_start_bf, z_start_bf, flag_ff, x_start_ff, y_start_ff, z_start_ff,
          flag_bb, x_start_bb, y_start_bb, z_start_bb;
      /* control current side scanning.       */

      float hs_eps_init; /* tolerance on homogeneity
       (fraction of slowness at source point) */

      //this class cannot be copied
      PodvinTime3D(const PodvinTime3D&);
      PodvinTime3D &operator=(const PodvinTime3D&);
    public:
      //! This is the definition of the original function by Podvin and Lecomte
      int time_3d(float *HS, float *T, int NX, int NY, int NZ, float XS,
          float YS, float ZS, float HS_EPS_INIT, int MSG);
      PodvinTime3D();
      virtual ~PodvinTime3D();
      };
  /* @} */
  }

#endif /* PODVINTIME3D_H_ */
