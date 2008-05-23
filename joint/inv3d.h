/* Local libary of "inv3d": */
#include "meschachd_no_complex/sparse.h"

#ifndef _INV3D_
#define _INV3D_


/*-------------------------------------------------------------*/
/* Define "Structures": */

/*Filenames of geo/dat-files*/
#ifndef  _FILENAME_
typedef struct _FILENAME_ 
{ 
	  char name[128];  /* Input-filenames of geo-files */

} FILENAME;
#define _FILENAME_
#endif


/* Input parameters:*/
#ifndef  _PARAMETER_
typedef struct _PARAMETER_
{
	 short binary_index; /*1= A model from a binary file will be used; 0 = the parameters, set in the parameter file govern the model*/

	 int index_tseis;	    /*specifying if the seismic traveltime data will be inverted (yes=1/no=0)*/
	 int index_grav;	    /*specifying if the gravity data will be inverted (yes=1/no=0)*/
	 int index_mt;          /*specifying if the MT data will be inverted (yes=1/no=0)*/

         int nx;			    /*Number of grid cells in x-direction*/
	 int ny;			    /*Number of grid cells in y-direction*/
	 int nz;			    /*Number of grid cells in z-direction*/
		
	 float h;			    /*Size of the grid cells in m*/

	 float org[3];		    /*components of the origin of the grid in m*/

	 short topo_index;      /*The index governs if the topography will be ignored(0) or introduced by means of triangulation(1)*/

	 int grid_cell_factor[3]; /*Sizes of "inversion" grid cell dimensions compared to the forward" grid cell dimensions in all three directions*/
	 int flex_grid_size;	  /*spatially invariant grid size (yes=1/no=0)*/

	 /*Seismic parameters*/
	 short index_velocity_field; /*This index governs if the velocity field is determined by the surface topography(1)*/

	 float vo;			    /*Starting velocity in m/s (at a defined gridpoint or at the surface)*/
	 float vmin,vmax;       /*Minimum and maximum velocity for the starting model in m/s */
	 float org_vo[3];	    /*Position, where v0 is defined*/
	 float vg[3];		    /*velocity-gradient (m/s)/m */

	 int nborder;		    /*Number of grid cells for the boundary*/
	 float vborder;		    /*velocity of the boundaries m/s*/
	 float v_air_water;	    /*velocity in the air or water layer*/

	 
	 short ray_index;		/*Kind of rays used (normal rays = 1; fat rays = 2)*/
	 float fatthres;		/*Threshold for the fat-ray in ms*/
	 float fatb;			/*Exponent governing the decrease of the weight with increasing "difference"-time in 1/ms*/	

	 int ray_density_index;	 /*Determine the hit matrix and the derivative weighted sum(DWS) for (fat-)rays (yes=1/no=0)*/
	 int ray_density_index2; /*Determine the ray density tensor (yes=1/no=0)*/

	 float vmina;			/*min.allowed velocity in m/s*/
	 float vmaxa;			/*max.allowed velocity in m/s*/
	 float maxadj;			/*max.velocity change per iteration*/

	 int para_rays_index;	/*Parallel rays are downscaled (following the idea of Hansruedi) (yes=1/no=0)*/
	 float angle_para;		/*Angle between two rays (in degree), where they are considered as "parallel"*/

	 /*Gravity parameters*/
	 short index_density_field; /*This index governs if the density field is determined by the surface topography(=1)*/

	 float go;				/*Starting density in g/cm^3 (at a defined gridpoint or at the surface)*/
	 float gmin,gmax;       /*Minimum and maximum density for the starting model in g/cm^3*/
	 float org_go[3];	    /*Position, where g0 is defined*/
	 float gg[3];		    /*density-gradient (g/cm^3)/m*/

	 float g_air_water;		/*density in the air or in the water layer (g/cm^3)*/

	 float gmina;			/*min.allowed density in g/cm^3*/
	 float gmaxa;			/*max.allowed density in g/cm^3*/

	 /*extra cell*/
	 float extra_cell;      /*The weighting of the extra cell for the gravity that is used to remove the data shift in the gravity data */

	 /*MT parameters*/

     short dimension_mt;    /*Specify the dimension of the MT forward calculation (1=1D/2=2D/3=2D forward calculation and 1D inversion)*/

	 short direc_2D_mt;     /*Specify for the 2D measurements along which horizontal direction the modelling will be performed (x=1/y=2)*/
	 int  nr_cells_2d_mt;	/*Specify the thickness of the layers for the individual 2-D MT models in numbers of forward cells*/
	 short kind_of_data_mt; /*Specify the kind of input data used: 1= TE-Mode, 2= TM-Mode, 3= 1D: Berdichewsky average/2D: Both Modes*/

	 short index_resistivity_field; /*This index governs if the elec.resistivity field is determined by the surface topography(=1)*/

	 float ro;				/*Starting resistivity in ohmm (at a defined gridpoint or at the surface)*/
	 float rmin,rmax;       /*Minimum and maximum elec. resistivity for the starting model in ohmm*/
	 float org_ro[3];	    /*Position, where r0 is defined*/
	 float rg[3];		    /*elec.resistivity-gradient ohmm/m*/

	 float r_air_water;		/*elec. resistivity in the air or in the water layer (ohmm)*/

	 float rmina;			/*min.allowed resistivity in ohmm*/
	 float rmaxa;			/*max.allowed resistivity in ohmm*/

	 int index_mt_weighting; /*Using impedances that are individually weighted for different frequencies*/

	 /*Other parameters*/
	 float eps;             /*Maximum distance between two coordinates location in m so that they are assigned to the same shot resp. receiver location */

	 float damp;			/*Parameter that governing the damping*/
	 float smoothfac;		/*Relationship between smoothing and damping (0= only smoothing/ 1= only damping)*/
	 float smoothratio_y;	/*Smoothing ratio y/x*/
	 float smoothratio_z;	/*Smoothing ratio z/x*/

	 /*Jin Chen changed new*/
	 int extra_smooth;      /*Specify if the extra smooth of the border cells will be introduced. (yes=1/no=0)*/
	 float val_extra_smooth;   /*Specify the value of the extra smooth of the border cells.*/
	/*Jin Chen changed 02.2007*/
	 float dist_extra_smooth;  /*Special the distance, that declare which inner cells will be considered in extra smooth.*/

	 int inv_max;			/*Max. number of iterations*/

	 int write_sens_out;	/*Write out the sensitivities of the different methods in binary files (yes =1; no=0)*/
	 int sens_in_perc;		/*Write the summed sensitivities in percentages (YES=1/N0=0)*/
							/*(if NO, the real weighted entries of the matrix are used)*/

	 FILENAME *files_geo;			/*Pointer on the structure with the names of the geo-input*/
	 FILENAME *files_dat;			/*Pointer on the structure with the names of the dat-input*/
	 FILENAME *binary_mod_dat;		/*Name of a binary model file*/
	 FILENAME *ascii_data_dat;		/*Name of the ascii data file*/
	 FILENAME *binary_data_dat;		/*Name of a binary data file; NOT USED ANY MORE !!*/
	 FILENAME *file_grid_size;		/*Parameter-file including the spatially dependent grid cell sizes (for inversion)*/
	 FILENAME *file_seis_grav_rel;	/*Parameter-file including the relationship of velocity and density*/
	 FILENAME *file_seis_res_rel;	/*Parameter-file including the relationship of velocity and elec.resistivity*/

} PAR;
#define _PARAMETER_
#endif

/* Parameters of the initial velocity/gravity or resistivity model:*/
#ifndef  _GRADIENT_
typedef struct _GRADIENT_
{
	short index_topo;     /*This index governs if the velocity/density/resistivity field is determined by the surface topography(1) or */
						  /*by a defined gridpoint "org"*/

	float value_0;		  /*Starting value (velocity[m/s], density[g/cm^3] or resistivity[ohmm]) (at the defined gridpoint "org" or at the surface)*/
	float min, max;		  /*Minimum and maximum value (velocity[m/s], density[g/cm^3], resistivity[ohmm]) for the starting model*/
	float org[3];		  /*Position, where "value_o" is defined*/
	float g[3];			  /*Gradient for velocity[(m/s)/m],density[g/cm^3] or resistivity[ohmm/m]*/

} GRADIENT;
#define  _GRADIENT_
#endif

/* Parameters of the grid structure:*/
#ifndef  _GRID_STRUCT_
typedef struct _GRID_STRUCT_
{
	    int nx;				/*Number of grid cell center in x-direction*/
	    int ny;				/*Number of grid cell center in y-direction*/
	    int nz;				/*Number of grid cell center in z-direction*/

		float h;			/*Size of the grid cells in m (!The same in all three directions!)*/

	    float org[3];		/*Components of the origin (center of the first cell in the upper left corner) of the grid in m*/
						
		short topo_index;   /*The index governs if the topography will be ignored(0) or introduced by means of triangulation(1)*/

		int grid_cell_factor[3]; /*Sizes of "inversion" grid cell dimensions compared to the forward" grid cell dimensions*/

		int *border_index;	/*grid cell belongs to the border or to the "air"(yes=0,no=1)*/

		/*Seismic Velocity parameters*/
		float vborder;		/*velocity of the boundaries m/s*/
		int nborder;		/*Number of grid cells for the boundary*/

		float v_air_water;	/*velocity [m/s] in the air or water layer*/

		double *slow;		/*Slowness model used for the forward model (normalized by the grid cell size)*/
		GRADIENT grad_vel;  /*Structure including the parameters of the initial (gradient) velocity model*/

		/*Density parameters*/
		float dens_air_water; /*density [g/cm^3] of the air or the water layer*/

		double *dens;		/*density model used for the forward calculation*/
		GRADIENT grad_dens; /*Structure including the parameters of the initial (gradient) density model*/

		/*Electromagnetical parameters*/
		float res_air_water; /*resistivity [ohmm] of the air or the water layer*/
		
		double *res;			/*resistivity model used for the forward calculation*/
		GRADIENT grad_res;	/*Structure including the parameters of the initial (gradient) elc. resistivity model*/

} GRID_STRUCT;
#define  _GRID_STRUCT_
#endif



/* Geometry of the topography points, the seismic shots and receivers, the gravimetric and MT stations*/
#ifndef  _GEOMETRY_
typedef struct _GEOMETRY_
{
		int nstat;      /*Number of receiver/shot stations and fixpoints*/
		int *coor_info; /*Coordinate information used for the topography (yes=1,no=0)*/
	    float *x,*y,*z; /*Positions of the shot/receiver locations and fixpoints in m*/
		float eps;      /*Maximum distance between two coordinates location in m so 
						  that they are assigned to the same shot resp. receiver location */
		/* seismic parameters*/
		int nshot;      /*Number of shot positions*/
		int nrec;       /*Number of receiver positions*/
		
		/*gravimetric parameters*/
		int nstat_grav;	/*Number of gravimetric stations*/
		/*MT parameters*/
		int nstat_mt;   /*Number of MT stations*/

} GEOMETRY;
#define  _GEOMETRY_
#endif


/*Structure with the topography information*/
#ifndef _TOPO_STRUCT_
typedef struct _TOPO_STRUCT_
{
		long nr_of_topo_points; /*Number of topography points*/
		long nr_of_triangles;	/*Number of triangles after the delaunay triangulation*/

		double *x;  /*x coordinates of the topography points*/
		double *y;  /*y coordinates of the topography points*/
		double *z;  /*z coordinates of the topography points*/

		int *index; /*Index of the topography points that specify the corners of triangles (first triangle = first three entries)*/

} TOPO_STRUCT;
#define  _TOPO_STRUCT_
#endif


/* Structure including informations about the seismic, gravity and MT data*/
#ifndef  _DATA_STRUCT_
typedef struct _DATA_STRUCT_
{
	/*Seismic parameters*/
		long ndata_seis; /*Number of picked first-breaks*/
		long ndata_seis_act; /*Number of rays that could be traced back*/

		double rms_seis;	/*RMS value of the traveltime residuals in ms*/

		int *sno;		/*List of the shot position numbers of the traces for which the first breaks were picked*/
		int *rno;		/*List of the receiver position numbers of the traces for which the first breaks were picked*/
	 double *tobs;		/*Observed travel times for the different shot-receiver combinations in ms*/
	 double *tcalc;		/*Calculated travel times for the different shot-receiver combinations in ms*/
						/*REMARK: sno, rno, tobs and tcalc have the same number of elements (ndata_seis) and are linked to each other*/
	 float *xdist;		/*Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/
	 float *ydist;		/*Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/
	 float *zdist;		/*Distance (in x-direction) between the shots and receivers locations for all shot-receiver combinations in m*/

	    int *shots;     /*List of all shot positions used (number of elements: nshot)*/
		int *recs;		/*List of all receiver positions used (number of elements: nrec)*/
        
		int *lshots;    /*Number of active receivers for the corresponding shot position number (related to shots)*/
		int *lrecs;     /*Number of used shots for the corresponding receicer position number (related to recs)*/

	double timeshift;   /*Constant timeshifts of all traces in ms*/

	double *weigth_seis; /*Weighting of the rays for the inversion (1.0 = is the usually weigthing)*/

	/*Gravity parameters*/
	long ndata_grav;	 /*Number of gravity data measurements*/

	double rms_grav;	/*RMS value of the gravity residuals in mgal*/

	  int  *gno;		 /*List of gravity position numbers (in the data files)*/
	double *obs_grav; 	 /*Observed gravity data (in mgal) */ 
	double *calc_grav;	 /*Calculated gravity data (in mgal)*/
						 /*REMARK: gno, obs_grav and calc_grav have the same number of elements (ndata_grav) and are linked to each other*/

	  int  *gravs;		 /*List of all gravity stations used*/	

	double *weigth_grav; /*Weighting of the gravity measurements (1.0 = is the usually weighting)*/

	/*MT parameters*/
	long ndata_mt;		  /*Number of MT measurements*/
	long ndata_freq_mt;   /*Number of all data (sum of ALL frequencies)*/

	double rms_mt_real;	 /*RMS value of the real part of the MT residuals of Z*/
	double rms_mt_imag;	 /*RMS value of the imaginary part of the MT residuals of Z*/
    
     int  *mno;			 /*List of MT position numbers (in the data files)*/

	 int  *mts;		 /*List of all MT stations used*/

	 int *nfreq_mt;		 /*Number of MT frequency for each station*/

	double **freq_mt;	 /*Frequencies of the MT measurements (HZ)*/
    
			/*TE-response for 2D or common response for 1-D*/
	double **real_mt_TE;	 /*Obs.real part of the  MT impedance E_y/H_x (Ohmm)*/ 
	double **imag_mt_TE;	 /*Obs.imaginary part of the MT impedance E_y/H_x (Ohmm) */
	double **calc_real_mt_TE; /*Calc.real part of the MT impedance E_y/H_x (Ohmm)*/
	double **calc_imag_mt_TE; /*Calc.imag part of the MT impedance E_y/H_x (Ohmm)*/
			/*TM-response for 2D; NO influence on 1D modelling*/
	double **real_mt_TM;	 /*Obs.real part of the  MT impedance E_x/H_y (Ohmm)*/ 
	double **imag_mt_TM;	 /*Obs.imaginary part of the MT impedance E_x/H_y (Ohmm) */
	double **calc_real_mt_TM; /*Calc.real part of the MT impedance E_x/H_y (Ohmm)*/
	double **calc_imag_mt_TM; /*Calc.imag part of the MT impedance E_x/H_y (Ohmm)*/
    
	double **weigth_mt;	 /*Weighting of the MT measurements; every frequency for its own (1.0 = is the usually weighting)*/

} DATA_STRUCT; 
#define  _DATA_STRUCT_
#endif

/* Structure for the 2-D MT foward calculation*/
#ifndef _CALC_2D_MT_
typedef struct _CALC_2D_MT_
{
	    long nr_cells_2d_mt;	 /*Specify the thickness of the slices for the individual 2-D MT models in numbers of forward cells*/
		double min_y, max_y; /*Minimum and maximum y-coordinates of the slice (vertical to the 2-D plane)*/
	
		double *org[2];		/*Components of the origin (upper-left corner of the first cell(upper-left)) of the slice in m*/
		long *nx,*nz;			 /*Number of cells in horizontal and vertical direction within the 2D plane (without border cells and grid refinement)*/
		double **hx,**hz;			 /*Grid sizes in x- and z-direction m*/

		long *nz_atm,*nz_ionos;  /*Number of used cells in the atmosphere and in the ionosphere*/
		long *nz_basement;       /*Number of cells at the basement*/
		long *nx_right_border,*nx_left_border; /*Number of cells at the left and right border*/

		double *x,*z;        /*x and z coordinates of the stations in the slice*/
		long nstat_mt;		 /*Nr. of mt stations in a slice*/
		long *mts;			 /*List of station numbers in the slice*/

		long nfreq;			 /*Nr.of frequencies used for the stations in the considered slice*/
		double *freq;		 /*List of measured frequencies for the slice*/
		long *freq_ndata;     /*Nr of measurements using the frequency in the slice*/
		long **freq_data;     /*List of the indices of the data measurements for this frequency*/
		long **freq_stat;     /*List of the station numbers for this frequency*/

		double **Ex_r,**Ex_i;	/*The electric fields at the stations (x-component; real and imaginary part)*/
		double **Ey_r,**Ey_i;	/*The electric fields at the stations (y-component; real and imaginary part)*/
		double **Hx_r,**Hx_i;	/*The magnetic fields at the stations (x-component; real and imaginary part)*/
		double **Hy_r,**Hy_i;	/*The magnetic fields at the stations (y-component; real and imaginary part)*/
		
		int    **index_slice; /*Index that specify if the cell is active (==1) or NOT (==0)*/
		double **res_slice;	  /*Resistivity values within the cells of the slices*/

		int *cell_scal_factor; /*Factor for down-scaling the grid cell sizes; (controlled by the skin-depth)*/

} CALC_2D_MT;
#define  _CALC_2D_MT_
#endif


/* Structure for the relationship of velocity and density and velocity and resistivity*/
#ifndef _REL_VDR_STRUCT_
typedef struct _REL_VDR_STRUCT_
{
		int index;	/*Specify the kind of relationship*/
					/*		= 1: Linear relationship v[m/s] = a + rho*b
							= 2: v-res relationship (Dell' Aversana, 2001, First Break, 19,335-341)
								v[m/s] = a(ln(ln(res[ohmm]))) + b
							= 9: table of densities[mgal](resistivities[ohmm]) and velocities [m/s]   */             

		/*Linear relationship*/
		float lin_para_a;	/*First free parameter of the equation*/
		float lin_para_b;  /*Second free parameter of the equation*/

		/*Dell' Aversana*/
		float dell_para_a;	/*First free parameter of the equation*/
		float dell_para_b; /*Second free parameter of the equation*/

		/*Tabular link*/
		long tab_nr_pairs;	/*Number of pairs*/
		double *tab_v;		/*velocities in [m/s]*/
		double *tab_dr;		/*densities in [mgal] or resistivities [ohmm]*/

} REL_VDR_STRUCT;
#define _REL_VDR_STRUCT_
#endif

/* Structure for the ray paths*/
#ifndef _RP_STRUCT_
typedef struct _RP_STRUCT_
{
		long n;				/*Number specifying the row in the matrix for the shot-receiver combination*/
		long nray;			/*Number of segments for the ray*/
		double *x;			/*		Forward grid: x-position of the begin of the ray-path segment*/
							/*BUT:	Inversion grid: x-component of the ray in the cell*/
		double *y;			/*		Forward grid: y-position of the begin of the ray-path segment*/
							/*BUT:	Inversion grid: y-component of the ray in the cell*/
		double *z;			/*		Forward grid: z-position of the begin of the ray-path segment*/
							/*BUT:	Inversion grid: x-component of the ray in the cell*/
		double *len;		/*Length of the ray-path segment*/
		long *ele;			/*Grid cell position in the grid*/

} RP_STRUCT;
#define _RP_STRUCT_
#endif

/*Structure for the fat-rays*/
#ifndef _F_RP_STRUCT_
typedef struct _F_RP_STRUCT_
{
		long n;				/*Number specifying the row in the matrix for the shot-receiver combination*/
		long ncell;			/*Number of cells which contribute to the fat ray*/
		float fatthres;		/*Threshold for the fat-ray in ms*/
		float fatb;			/*Exponent governing the decrease of the weight with increasing "difference"-time in 1/ms*/	
		double *weight;		/*Normalized contribution of the cell to the fat-ray*/
		long *ele;			/*Grid cell position in the grid*/

} F_RP_STRUCT;
#define _F_RP_STRUCT_
#endif

/*Structure for the gravity derivatives*/
#ifndef _GRAV_STRUCT_
typedef struct _GRAV_STRUCT_
{
		long n;				/*Number that specifies the row in the matrix for the measurement*/
		long ncell;			/*Number of cells that contribute to the measurement*/
		double *deriv;		/*Value of the derivative dG/ds or dG/drho*/
		long *ele;			/*Grid cell position in the grid*/

} GRAV_STRUCT;
#define _GRAV_STRUCT_
#endif

/*Structure for the gravity derivatives*/
#ifndef _MT_STRUCT_
typedef struct _MT_STRUCT_
{
		long *n;			/*Numbers that specify the rows in the matrix for the measurement (size: size: 2*nfreq because of the real and the imaginary part of Z)*/
		int nfreq;  		/*Number of frequencies*/
		long ncell;			/*Number of cells that contribute to the measurement*/
		double **deriv;		/*Values of the derivatives dZ/ds or dZ/drhoR* (First index specify the cell; second index specify the frequency (size: 2*nfreq because of the real and the imaginary part of Z; even numbers = real part; odd numbers = imag. part)*/
		long *ele;			/*Grid cell position in the grid*/

} MT_STRUCT;
#define _MT_STRUCT_
#endif

/*Structure to organize the cell parameters during back tracing the rays*/
#ifndef _CELL_STRUCT_
typedef struct _CELL_STRUCT_
{
		int xno;			/*position number in x-direction of the cell in the grid*/
		int yno;			/*position number in y-direction of the cell in the grid*/
		int zno;			/*position number in z-direction of the cell in the grid*/
		int dirx_i;			/*The ray runs in negative x-direction into the cell = 1; the ray runs in positive x-direction into the cell = 1; else=0 */
		int diry_i;			/*The ray runs in negative y-direction into the cell = 1; the ray runs in positive y-direction into the cell = 1; else=0 */
		int dirz_i;			/*The ray runs in negative z-direction into the cell = 1; the ray runs in positive z-direction into the cell = 1; else=0 */
		double xpos;		/*The position of the starting point of the ray in x-direction (normalized by the position in the grid)*/ 
		double ypos;        /*The position of the starting point of the ray in y-direction (normalized by the position in the grid)*/
		double zpos;        /*The position of the starting point of the ray in z-direction (normalized by the position in the grid)*/
	
} CELL_STRUCT;
#define _CELL_STRUCT_
#endif

/* Inversion dependent parameters*/
#ifndef  _INV_
typedef struct _INV_
{
		int ninv;			/*Number of iterations*/
		int ninv_max;		/*max.number of iteration*/

		long nvel;			/*Number of ALL cells in the inversion grid*/
		long nvel_used;		/*Number of the active cells in the inversion grid*/

		long ndata;			 /*Number of all data measurements (corresponds to the number of "data" rows in the matrix)*/
		long seis_first_row; /*First row with seismic data*/
		long seis_nr_row;	 /*Number of rows with seismic data*/
		long grav_first_row; /*First row with gravity data*/
		long grav_nr_row;	 /*Number of rows with gravity data*/
		long mt_first_row;	 /*First row with MT data*/
		long mt_nr_row;		 /*Number of rows with MT data*/

		long nrow;			/*Number of rows in the matrix*/
		long ncol;			/*Number of columm in the matrix*/

							/*Scaling of the different data sets to balance the different methods to each other (default = 1)*/
		double rel_scal_seis; /*seismic*/
		double rel_scal_grav; /*gravity*/
		double rel_scal_mt;	  /*MT*/

		int nx;				/*Number of inversion cells in x-direction (for a constant grid)*/
		int ny;				/*Number of inversion cells in y-direction (for a constant grid)*/
		int nz;				/*Number of inversion cells in z-direction (for a constant grid)*/

		VEC *v_res;			/*(Residual) data vector*/
		double *para_mod;	/*Final model vector (For joint inversions and seismic traveltime inversion are the entries slownesses in [s/m])*/
							/*For simple gravity inversion are the entries densities in [g/cm^3]*/
							/*For simple MT inversion are the entries resistivities in [ohmm]*/
		double *para_mod2,*para_mod3;
							/*Final model vectors of the other parameters (gravity [mgal] and/or resistivity [ohmm]) that are determined from the velocities in the joint inversions*/
		double *mx;			/*Vector consisting of the original model parameter*/
		SPMAT *g;			/*Matrix*/

		double vmin;		/*min.allowed velocity in m/s*/
		double vmax;		/*max.allowed velocity in m/s*/
		double maxadj;		/*max.velocity change per iteration*/
		double gmin;		/*min.allowed density in g/cm^3*/
		double gmax;		/*max.allowed density in g/cm^3*/
		double rmin;		/*min.allowed resistivity in ohmm*/
		double rmax;		/*max.allowed resistivity in ohmm*/

		double angle_para;  /*Angle between two rays (in degree), where they are considered as "parallel"*/

		/*Jin Chen changed new*/
		long x_left_border_nr;  /*the number of the actived left x_border cells.*/
		long x_right_border_nr;  /*the number of the actived right x_border cells.*/

		long y_front_border_nr;  /*the number of the actived front y_border cells.*/
		long y_back_border_nr;  /*the number of the actived back y_border cells.*/
		
		/*extra cell (BJOERN_MOD)*/
		double dens_extra_cell; /*density [g/cm^3] of the extra cell below the model*/

		
} INV;
#define  _INV_
#endif

/* Parameters for each inversion cell*/
#ifndef _BIDX_STRUCT_
typedef struct _BIDX_STRUCT_
{
		double bdamp;		/*Damping parameter*/
		double xo;			/*x-coordinate of the center of the inversion cell in m*/
		double yo;			/*y-coordinate of the center of the inversion cell in m*/
		double zo;			/*x-coordinate of the center of the inversion cell in m*/
		double xdim;		/*Dimension of the cell in x-direction in m*/
		double ydim;		/*Dimension of the cell in y-direction in m*/
		double zdim;		/*Dimension of the cell in z-direction in m*/

		int border_inv_cell_x_left;	/*Specify, if the cell is located at the left x-border (yes==1;no==0)*/
		int border_inv_cell_x_right;	/*Specify, if the cell is located at the right x-border (yes==1;no==0)*/

		int border_inv_cell_y_front;	/*Specify, if the cell is located at the front y-border (yes==1;no==0)*/
		int border_inv_cell_y_back;	/*Specify, if the cell is located at the back y-border (yes==1;no==0)*/

		int border_inv_cell_z;	/*Specify, if the cell is located at the z-border (yes==1;no==0)*/

		double val_slow;	/*slowness in the cell*/
		double val_dens;	/*density in the cell*/
		double val_res;		/*resistivity in the cell*/

		int use;			/*inversion cell used for inversion (yes=1;no=0)*/
		long used_nr;		/*number after renumbering the ACTIVE inversion cells*/

		int nele;			/*number of forward-cells included in the inversion cell (ONLY the ones, which are activated for the inversion)*/
		long *ele;			/*the positions of the forward cells in the grid*/

		int nele_inv;		/*Number of former inversion cells in the "new" inversion cells*/
		long *ele_inv;			/*number of the former inversion cells included in the "new" inversion cell*/

} BIDX_STRUCT;
#define _BIDX_STRUCT_
#endif

/*Parameters governing the spatially variations of the inversion grid*/
 #ifndef _FLEX_STRUCT_
 typedef struct _FLEX_STRUCT_
 {
		int nx;				/*Number of inversion grid changes in x-direction*/
		int ny;				/*Number of inversion grid changes in y-direction*/
		int nz;				/*Number of inversion grid changes in z-direction*/
		double *xpos;		/*x-positions of the inversion grid, where the modified grid starts*/
		double *ypos;		/*y-positions of the inversion grid, where the modified grid starts*/
		double *zpos;		/*z-positions of the inversion grid, where the modified grid starts*/
		int *xratio[3];		/*new grid cell ratios*/
		int *yratio[3];		/*new grid cell ratios*/
		int *zratio[3];		/*new grid cell ratios*/
} FLEX_STRUCT;
#define _FLEX_STRUCT_
#endif

/* Parameters for the regularistion*/
#ifndef _REGU_STRUCT_
typedef struct _REGU_STRUCT_
{
		double damp;		/*Parameter that governing the damping*/
		/*Jin Chen changed new*/
		double smooth;      /*Parameter that governing the smoothing*/
		double val_extra_smooth; /*Parameter that governing the extra smoothing*/
		/*Jin Chen changed 02.2007*/
	    double dist_extra_smooth;  /*Special the distance of inner cells will be considered in extra smooth.*/
		/*extra cell*/
		double w_extra_cell; /*The weighting of the extra cell that is used to remove data shifts in the gravity data, if it is smaller than 0, then the extra cell will
							 not be intruduced.*/
		/*BJOERN_MOD5*/

		double damp2;		/*Scaling of the damping by the so-called slowness-regularization factor*/
							/*(Sum of all sensitivities/Sum of all inversion cells with non-zero entries)*/
		double smoothratio_y; /*Smoothing ratio y/x*/
		double smoothratio_z; /*Smoothing ratio z/x*/
} REGU_STRUCT;
#define _REGU_STRUCT_
#endif


/* Flags*/
#ifndef _FLAG_STRUCT_
typedef struct _FLAG_STRUCT_
{
		int index_tseis;		/*specifying if the seismic traveltime data will be inverted (yes=1/no=0)*/
		int index_grav;			/*specifying if the gravity data will be inverted (yes=1/no=0)*/
		int index_mt;			/*specifying if the MT data will be inverted (yes=1/no=0)*/

		int index_tseis_grav;	/*specifiying if the seismic traveltime AND gravity data will be jointed inverted (yes=1/no=0)*/
		int index_tseis_mt;		/*specifiying if the seismic traveltime AND MT data will be jointed inverted (yes=1/no=0)*/
		int index_grav_mt;		/*specifiying if the gravity AND MT data will be jointed inverted (yes=1/no=0)*/

		int index_tseis_grav_mt; /*specifiying if the seismic traveltime, gravity AND MT data will be jointed inverted (yes=1/no=0)*/

		int index_joint;		/*Joint inversion (more than one kind of data will be inverted) (yes=1/no=0)*/

		int link_seis_grav_mod;	/*Specify, if the density starting model is linked with seismic velocity model (no=0/yes=1)*/
		int link_seis_mt_mod;	/*Specify, if the elec.resistivity starting model is linked with seismic velocity/gravity model (no=0/yes(seismic velocity)=1/yes(gravity)=2)*/

		int kind_of_rays;		/*Kind of ray used (1=conventional rays/2=fat-rays)*/
		int balance_para_rays;  /*Parallel rays are downscaled (following the idea of Hansruedi) (yes=1/no=0)*/
		int ray_density;		/*Determine the hit-matrix and the DWS (yes=1/no=0)*/
		int ray_density2;		/*Determine the ray density tensor (yes=1/no=0)*/

		int dimension_mt;    /*Specify the dimension of the MT forward calculation (1=1D/2=2D)*/
		int kind_of_data_mt; /*Specify the kind of input data used: 1= TE-Mode, 2= TM-Mode, 3= 1D: Berdichewsky average/2D: Both Modes*/
		int nr_of_modes;	 /*Specify the number of modes that will be inverted: Only: TM or TE ==1 , Both ==2*/	
		int direc_2D_mt;     /*Specify for the 2D measurements along which horizontal direction the modelling will be performed (x=1/y=2)*/
		int nr_cells_2d_mt;	 /*Specify the thickness of the slices for the individual 2-D MT models in numbers of forward cells*/
	
		int index_mt_weighting; /*Using impedances that are individually weighted for different frequencies */

		int flex_grid_size;		/*spatially invariant grid size (yes=1/no=0)*/

		int do_smooth;			/*Perform smoothing (yes=1/no=0)*/
		int do_damp;			/*Perform damping (yes=1/no=0)*/

		int write_sens_out;		/*Write out the total (balanced) sensitivities of the different methods in binary files (yes =1; no=0)*/
		int sens_in_perc;		/*Write the summed sensitivities in percentages (YES=1/N0=0)*/
								/*(if NO, the real weighted entries of the matrix are used)*/

		/*Jin Chen changed new*/
		int do_extra_smooth;    /*Perform extra smooth to the border inversion cells. (yes=1/no=0)*/


} FLAG_STRUCT;
#define _FLAG_STRUCT_
#endif



#endif
