//============================================================================
// Name        : PodvinTime3D.h
// Author      : May 8, 2009
// Version     :
// Copyright   : 2009, mmoorkamp (modifications, otherwise see below)
//============================================================================

/*---------------------------------------------------------------------------*/
/*  TIME_3D: FINITE DIFFERENCE COMPUTATION OF FIRST ARRIVAL TIMES IN 3D.     */
/*---------------------------------------------------------------------------*/
/*  P.PODVIN, Geophysique, Ecole des Mines de Paris, Fontainebleau.          */
/*  e-mail: Pascal.Podvin@ensmp.fr          Tel: 33-(1) 64 69 49 25.         */
/*                                                                           */
/*  First release: April 1991; Previous revision: 7 May 1993                 */
/*  A version with improved portability was released 19 January 2004         */
/*  Previous update was dated 16 February 2006 (with corrections in the      */
/*  initialization procedure, with improved consistency when source point    */
/*  is at the immediate vicinity of a velocity heterogeneity)                */
/*  This version is dated 16 October 2006 (two patches, seek '16102006')     */
/*                                                                           */
/* PROTOTYPE : (see ../include/fdtimes.h)                                    */
/*                                                                           */
/* int time_3d(float *hs, float *t, int nx, int ny, int nz,                  */
/*             float xs, float ys, float zs, float eps_init, int messages);  */
/*                                                                           */
/* ARGUMENTS (C description; all FORTRAN arguments are pointers)             */
/*                                                                           */
/*       (int)     nx,ny,nz        : dimensions of the time field (number of */
/*                                   grid points). Cells are cubic. No       */
/*                                   dimension may be lower than 2.          */
/*                                                                           */
/*       (float *) hs,t            : 1D arrays of nx*ny*nz elements.         */
/*                                   t will contain the computed time field  */
/*                                   arranged as a succession of planes x    */
/*                                   (resp. z)=Cst, each plane consisting of */
/*                                   a suite of y=Cst vectors when this      */
/*                                   routine is called from a program in C   */
/*                                   (resp. in FORTRAN). hs contains model   */
/*                                   data (slowness*mesh spacing) organized  */
/*                                   in the same manner. Within the routine, */
/*                                   values stored in hs are implicitly      */
/*                                   ascribed to cell centers, whereas times */
/*                                   computed and stored in t refer to cell  */
/*                                   corners (grid points). Cells located at */
/*                                   coordinates x=nx-1, y=ny-1 or z=nz-1    */
/*                                   (x=nx, y=ny, z=nz in FORTRAN dialect)   */
/*                                   are thus dummy cells (out of the model).*/
/*                                   The corresponding values will be        */
/*                                   ignored and treated as "infinite".      */
/*                              Note:                                        */
/*                                   Values stored in hs may not be negative */
/*                                   and must not exceed 0.499e+19 (see      */
/*                                   macro FD_HUGE defined below).           */
/*                                                                           */
/*       (float)   xs,ys,zs      :   Point source coordinates referred to    */
/*                                   "upper left" corner of model (grid      */
/*                                   point with lowest indices), expressed   */
/*                                   in mesh spacing (h) unit.               */
/*                                   Licit ranges: [0.0,nx-1.][0.0,ny-1.]    */
/*                                   [0.,nz-1.]                              */
/*                                   If source point is found to be out      */
/*                                   of the licit range, the timefield is    */
/*                                   treated as already initialized in the   */
/*                                   calling program (extended source, e.g., */
/*                                   exploding reflector). In such case, the */
/*                                   timefield must not be uniform, and      */
/*                                   nodes where times are computed must be  */
/*                                   initialized with a value certainly      */
/*                                   higher than their final arrival time.   */
/*                              Note: although this is not required when the */
/*                                    source is licit, you should always     */
/*                                    initialize array t with whatever       */
/*                                    constant value (why not 0?), so that   */
/*                                    you'll be warned if a wrong (illicit)  */
/*                                    source location is entered.            */
/*                                                                           */
/*       (float)   eps_init      :   tolerance on relative inhomogeneity     */
/*                                   used (only) during initialization.      */
/*                                   Although licit values are in [0,1],     */
/*                                   relevant values should be <0.01.        */
/*                                   0.001 is a reasonable choice.           */
/*                                   (arg. ignored if source is illicit)     */
/*                                                                           */
/*       (int)     messages        : 0 makes the routine silent (except on   */
/*                                   diagnosed error); 1: small talk mode;   */
/*                                   >1: debug mode (verbose).               */
/*                                   (a negative value is treated as 1).     */
/*                                                                           */
/* VALUE :                                                                   */
/*                                                                           */
/*       time_3d() returns a nonzero value if an error was detected.         */
/*       An explicit error message is printed on 'stderr'.                   */
/*                                                                           */
/* CALLING TIME_2D FROM A PROGRAM WRITTEN IN FORTRAN :                       */
/*                                                                           */
/*       The C-function time_3d_() is provided as the interface to Fortran.  */
/*       In this routine, all arguments are pointers as required by Fortran  */
/*       (where routine arguments are passed by address, not by value), and  */
/*       dimensions 'x','z' are swapped to mimic standard Fortran memory     */
/*       mapping (hs[i][j][k] in C "means" HS(K,J,I) in Fortran).            */
/*       Normally, calling TIME_3D (without the trailing underscore) in      */
/*       Fortran will automatically result in a call to time_3d_() in C.     */
/*    Compiling :                                                            */
/*       With Sun compilers, nothing particular is required to obtain this   */
/*       automatic behaviour.                                                */
/*       With the GNU compilers (gcc, g77), you will need to compile your    */
/*       Fortran program with option -fno-second-underscore to obtain this.  */
/*       Because the C/FORTRAN interface is not standardized, you might      */
/*       have to browse your documentation on other platforms.               */
/*       If you get into trouble :                                           */
/*       Compiling this program with the option -DDEBUG_ARGS allows to check */
/*       what C-function is actually called (time_3d or time_3d_) and how    */
/*       its arguments are interpreted.                                      */
/*     Declarations in the calling program :                                 */
/*       As seen from Fortran, TIME_3D is an INTEGER FUNCTION.               */
/*       It should thus be declared as such in the calling program.          */
/*       Please note that arguments passed to TIME_3D MUST be declared with  */
/*       their correct types (e.g., XS,YS,ZS MUST NOT be declared INTEGER,   */
/*       while NX,NY,NZ MUST NOT be declared REAL !)                         */
/*       Not conforming to this may result into incomprehensible situations  */
/*       (once again, compile option -DDEBUG_ARGS may help...)               */
/*       All scalar arguments are read-only and may thus be declared         */
/*       constant (using PARAMETER statements).                              */
/*     Program skeleton :                                                    */
/*       INTEGER NX,NY,NZ                                                    */
/*       REAL HS(NX,NY,NZ),T(NX,NY,NZ)                                       */
/* C or  REAL HS(NX*NY*NZ),T(NX*NY*NZ)                                       */
/*       REAL XS,YS,ZS,EPS_INIT                                              */
/*       INTEGER MESSAGES,TIME_3D,STATUS                                     */
/*       .....                                                               */
/*       STATUS=TIME_3D(HS,T,NX,NY,NZ,XS,YS,ZS,EPS_INIT,MESSAGES)            */
/*       IF(STATUS.NE.0)                                                     */
/*         STOP "time_3d diagnosed a (presumably fatal) error"               */
/*       ENDIF                                                               */
/*       .....                                                               */
/*       STOP                                                                */
/*       END                                                                 */
/*                                                                           */
/* RECENT UPDATES :                                                          */
/*                                                                           */
/*        Although time_3d has been used by dozens of people worldwide       */
/*        for many years, some changes were recently introduced (2003).      */
/*        These changes do not affect the routine's usage (except for        */
/*        some values returned on error, and the treatment of "infinite"     */
/*        slowness values, now forbidden in hs). Their only justification    */
/*        is improved portability.                                           */
/*        The changes are the following :                                    */
/*        - The program now conforms to ANSI-C (you must include fdtimes.h   */
/*          in your calling program, if it is written in C).                 */
/*        - I decided to drop Sun-specific calls to routines implementing    */
/*          the IEEE standards for infinity and related tests. This is non   */
/*          portable, and would even create problems on Sun platforms when   */
/*          the calling program is in Fortran !                              */
/*          As a consequence, a finite threshold is now used to test whether */
/*          a float is treated as infinite. No value in array hs is allowed  */
/*          to exceed this threshold (see macro FD_HUGE below).              */
/*        - Unpredictible crashes were seldom encountered on Intel based PCs */
/*          due to the fact that the routine's behaviour strongly depended   */
/*          on tests comparing floats in an "exact" manner (e.g., replacing  */
/*          a test x<y by x<=y altered significantly the code behaviour).    */
/*          These tests are now "fuzzified" with no significant consequence  */
/*          on precision (using macro EPS_FUZZY below).                      */
/*                                                                           */
/*        In 2005, Ari Trygvasson and Bjorn Bergman signaled an unexpected   */
/*        asymmetry of results in a particular case. The initialization      */
/*        procedure was found to provide wrong results when the source was   */
/*        located at the immediate vicinity of a velocity heterogeneity      */
/*        (with significant but not severe impact on precision). A quite     */
/*        significantly revised version was released on 2 February 2006.     */
/*                                                                           */
/* REFERENCE : Podvin & Lecomte, Geophys.J.Intern. 105, 271-284, 1991.       */
/*----------------------------------------------------------------------------*/
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
