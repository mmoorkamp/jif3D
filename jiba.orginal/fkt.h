#include <time.h>
#include "meschachd_no_complex/sparse.h"
/*--------------------------------------------------------------*/
/* Define "makros": */
#define PI 4.*atan(1.e0)

/*--------------------------------------------------------------*/
/* Define function:*/
#ifndef _CALLING_

/*in file toolbox2.c:*/
char *memory (char *prev_addr,int n,int size,char *progname);
int GetNext(char *line,FILE *inf);
int ReadInput(char *fname,PAR *par,int *num_geo, int *num_mod);
int WriteModSeisOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid);
int WriteModGravOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid);
int WriteModResOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid);
int WriteSensSeisOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid, double *sens);
int WriteSensGravOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid, double *sens);
int WriteSensMTOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid, double *sens);
int WriteRayOut(RP_STRUCT *rp, DATA_STRUCT data, GRID_STRUCT grid);
int WriteFatRayOut(F_RP_STRUCT *fr,DATA_STRUCT data,GRID_STRUCT grid);
int WriteDatSeisOut(DATA_STRUCT data, int inv, GEOMETRY geo);
int WriteDatGravOut(DATA_STRUCT data, int ninv, GEOMETRY geo);
int WriteDatMTOut(DATA_STRUCT data, int ninv, GEOMETRY geo, int dimension_mt, int kind_of_data_mt);
int WriteTravelTimeModOut(float *tt, int loc,int nx,int ny,int nz);
int WriteMergedDatOut(DATA_STRUCT data, int mt_data_format);
float *ReadTravelTimeModIn(char *fname);
int WriteRayDenseOut(GRID_STRUCT grid, int ninv, long ninv_cell, long *nrays, double *length, BIDX_STRUCT *inv_cell);
int	WriteRayDensityTensorOut(GRID_STRUCT grid, int ninv, long ninv_cell, double *rel_minmax_eig, double *rel_medmax_eig, BIDX_STRUCT *inv_cell);
int WriteMTResistivityOut(long layer_i, double freq, long nsamples, double *parameter, int kind_of_parameter);
int ReadTResistivityIn(long layer_i, double freq, long nsamples, double *parameter, int kind_of_parameter);
int WriteMTResistivityIndexOut(long layer_i, double freq, long nsamples, int *parameter);

/*set_parameters.c*/
int set_flags(PAR par,FLAG_STRUCT *flag);
int set_regu(PAR par,REGU_STRUCT  *regu);
int set_inv(PAR par, INV *inv);

/*in file data_and_geometry1.c*/
int MakeGeoandDataStructure(GEOMETRY *geo,FILENAME *files_geo,GEOMETRY *geo_all,DATA_STRUCT *data,FILENAME *files_dat, DATA_STRUCT *data_all, int ngeo_file, float eps, int *kind_of_data_mt);
int ReadGeoFiles(char *fname , GEOMETRY *geo);
GEOMETRY MergeGeometry(GEOMETRY *geo, int num_geo);
int ReadDataFilesOLD(char *infile,DATA_STRUCT *ttdata,GEOMETRY *geo); /*For Hansruedis ".dat-Formats"*/
int ReadDataFiles(char *infile,DATA_STRUCT *ttdata,GEOMETRY *geo, long *count, int *kind_of_data_mt); /*For Bjoerns ".dat-Formats"*/
DATA_STRUCT MergeData(DATA_STRUCT *data, int num_dat, GEOMETRY *geo);
int OptiData(DATA_STRUCT *data,GEOMETRY *geo, float eps);
int CalcShotRecDist(DATA_STRUCT *data,GEOMETRY geo);
int StatData(DATA_STRUCT *data,GEOMETRY geo);

/*in file binary_data_and_geometry.c*/
int ReadFromBinary(int nmod_file, GRID_STRUCT *grid, PAR par, FLAG_STRUCT *flag, GEOMETRY geo, char *fname);
int ReadModIn(GRID_STRUCT *grid, char *fname);

/*in file grid_routines.c:*/
int MakeGridStruct(PAR par,GRID_STRUCT *grid_struct);
int MakeGradientStruct(short tmp_index_velocity_field, float tmp_vgx,float tmp_vgy,float tmp_vgz,float tmp_orgx_vo, float tmp_orgy_vo,float tmp_orgz_vo,float tmp_vo,float tmp_vmin,float tmp_vmax,GRADIENT *grad);

/*in file init_cube4.c:*/
int MakeInitModels(GRID_STRUCT *grid, GEOMETRY *geo, FLAG_STRUCT flag, char *fname_dens, char *fname_res, TOPO_STRUCT *topo);
long *TopoCalc(GEOMETRY *geo, GRID_STRUCT *grid, TOPO_STRUCT *topo);
int MakeIndividualModels(GRID_STRUCT *grid, int kind_of_model, long *sample_topo, int vel_dens_link, char *fname, char *fname2);
long *FitGridtoTopo(double *Z, GRID_STRUCT *grid, GEOMETRY geo);
int ModCubeCalc(double *vel ,int *border_index, GRID_STRUCT *grid, long *sample_topo, int kind_of_model);
int CalcVelDepMod(GRID_STRUCT *grid, double *mod, REL_VDR_STRUCT rel_vdr);
int CalcGravDepMod(GRID_STRUCT *grid, double *mod, REL_VDR_STRUCT rel_vdr, REL_VDR_STRUCT rel_vdr2);
int CalcSurfDepMod(long *sample_topo,GRID_STRUCT *grid,double *vel, int kind_of_model);
int	CalcOrgDepMod(GRID_STRUCT *grid, double *mod, int kind_of_mod);
int AddBoundaries(GRID_STRUCT *grid, double *vel1, int *border_index1, double *vel2, int *border_index2);

/*in file delaunay_triangulation.c*/
int delaunay (double *x_in, double *y_in, int n, int **link);
int TriangleInterpol(double *X,double *Y,double *Z, long numXYZ ,double *x ,double *y ,double *z ,int numxyz ,int *corner_triangles,int ntriangles);

/*in file vel_dens_res_rel.c*/
int ReadRelVelDensRes(char *fname, REL_VDR_STRUCT *rel_vdr);
int CheckMonotom(long nr_pairs,double *x,double *y);
double VeltoDensRes(double vel_value, REL_VDR_STRUCT rel_vdr);
double DensRestoVel(double dens_res_value, REL_VDR_STRUCT rel_vdr);

/* in file prepare_modeling.c*/
int CheckCoordinates(GEOMETRY geo, DATA_STRUCT data, GRID_STRUCT grid, FLAG_STRUCT flag);
int CheckCoordinatesTopo(GEOMETRY geo,DATA_STRUCT data,GRID_STRUCT grid);

/*in file modeling_seismic.c*/
int ForwardModRay(GEOMETRY geo,GRID_STRUCT grid,DATA_STRUCT *data, RP_STRUCT *raypath, time_t time_start);
int ForwardModFatRay(GEOMETRY geo,GRID_STRUCT grid, DATA_STRUCT *data, F_RP_STRUCT *fat_rays, time_t time_start);
int FatRayCalc(float *tt1,float *tt2,long nx,long ny,long nz, float tcalc,F_RP_STRUCT *fr,double *slow);
float interpolate(float x,float y,float z,GRID_STRUCT *grid,float *data);
int RayCalc(float *tt, int nx,int ny,int nz, float Xs, float Ys, float Zs,float *Xr,float *Yr, float *Zr,int nrec, RP_STRUCT *rp);
double *TimeGrad(int x, int y, int z, float *tt,int ny, int nz);
CELL_STRUCT RayBackTrace(double gradx, double grady, double gradz, CELL_STRUCT cell, float *tt,int ny,int nz);
int ResortRays(RP_STRUCT *raypath,DATA_STRUCT data,GRID_STRUCT grid);

/*in podvin-lecomte-3D.c*/
int time_3d(float *HS,float *T,int NX,int NY,int NZ,float XS,float YS,float ZS,float HS_EPS_INIT,int MSG);

/*in file modeling_gravity.c */
int ForwardModGrav(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, TOPO_STRUCT topo, INV inv);/*Jin,12_03*/
double Response_3D_Rectangle(double dens, double x1, double x2, double y1, double y2, double z1, double z2);
double Response_2D_Rectangle(double dens, double x1, double x2, double z1, double z2);
int RemoveMeanGrav(DATA_STRUCT *data);
int AddEffectExtraCell(double dens_extra_cell, DATA_STRUCT *data, GRID_STRUCT grid); /*BJOERN_MOD5*/
double Response_2D_lower_triangle(double dens, double x1, double x2, double z1, double z2);
double Response_2D_upper_triangle(double dens, double x1, double x2, double z1, double z2);
double Response_2D_lower_triangle1(double dens, double x1, double x2, double z1, double z2);
double Response_2D_upper_triangle1(double dens, double x1, double x2, double z1, double z2);
double Response_3D_Rectangle_1(double dens, double x1, double x2, double y1, double z1, double z2);
double Response_3D_Rectangle_2(double dens, double x1, double y1, double z1, double z2);


/*in file modeling_mt.c */
int ForwardModMT_1D(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, int kind_of_data);
int	MT_1D_CALC(double ang_freq, double *z_real_part, double *z_imag_part, double *dz, long nr_layer, double *resist, double p);
int ForwardModMT_2D(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, CALC_2D_MT *mt_2D_struct_out, FLAG_STRUCT flag, TOPO_STRUCT topo);
int DetThickness2DSlice(CALC_2D_MT *mt_2D_struct, GRID_STRUCT grid, long layer, int direc_2D_mt, long nr_cells_2D_mt);
int DetStation2DSlice(DATA_STRUCT *data,GEOMETRY geo, CALC_2D_MT *mt_2D_struct, int direc_2D_mt);
int CalcAveValues2DSlice(GRID_STRUCT grid, long layer, int direc_2D_mt, long nr_cells_2D_mt, long nx, long nz, double *slice, double *cube);
int DetTopography2DSlice(GRID_STRUCT grid, long first_sample, int direc_2D_mt, long nr_cells_2d_mt, long nx, double *hx, TOPO_STRUCT topo, double *grid_points);
int DetTopographyMTStat2DSlice(GEOMETRY geo, double *grid_points, int direc_2D_mt, long nx, double *hx, double orgx, long nr_stat_mt, long *mts);
int AddTopography2DSlice(long nx, long nz, double *hz, double *grid_points, double orgz, double *res_value, float res_water_air, int *index_topo);
int AdjustFreq2Dslice(CALC_2D_MT *mt_2D_struct, DATA_STRUCT *data);
int MakeRefinedGridMT(double grid_h, CALC_2D_MT *mt_2D_struct, double min_res, double *org_res_slice, double *org_index_slice, long freq_index);
int grid_refinement(long *nx, long *nz, double *hx, double *hz, double org_cell_size, long cell_scal_factor, double *org_slice1, double *new_slice1, double *org_slice2, double *new_slice2);
int MakeLocalRefinedGridMT(double grid_h, CALC_2D_MT *mt_2D_struct,double min_res, long nr_layer, long freq_index);
int AddCellsAtBoundaryMT(CALC_2D_MT *mt_2D_struct, double max_res_basement, double max_res_left, double max_res_right,long freq_index, long *nx_left_border, long *nx_right_border, long *nz_basement);
int AddCellsAtBoundaryAirMT(CALC_2D_MT *mt_2D_struct, long freq_index,  long *nz_air, long *nz_ionos);
int DetImpAtStations(CALC_2D_MT *mt_structure, DATA_STRUCT *data, double *H_r, double *H_i, double *E_r, double *E_i, int index_mode, long nz_atm, long frequency_nr);
int MakeGlobMTStruct(long layer_i, long freq_i, CALC_2D_MT *mt_2D_struct, CALC_2D_MT *mt_2D_struct_out);

/*in file balance_data.c*/
int	ScalParallelRays( RP_STRUCT *raypath, BIDX_STRUCT *inv_cell, long ncell, long nray, double angle_threshold, double *weighting_factor, int kind_of_weight);
int BalanceSparseMatrix(INV *inv, FLAG_STRUCT flag, SPMAT *A, double *res);
int BalanceMTfrequencies(INV *inv, SPMAT *A, double *res);
int ReweightSparseMatrix(FLAG_STRUCT flag, F_RP_STRUCT *fatray,RP_STRUCT *raypath, GRAV_STRUCT *grav, MT_STRUCT *mt, BIDX_STRUCT *inv_cell, DATA_STRUCT data, SPMAT *A);
int ReweightResdDataVector(FLAG_STRUCT flag,VEC *vec, DATA_STRUCT data, RP_STRUCT *raypath, F_RP_STRUCT *fatray, GRAV_STRUCT *grav, MT_STRUCT *mt);

/*in inversion_new27a.c*/
int InvRoutine(GEOMETRY geo, GRID_STRUCT *grid,RP_STRUCT *raypath, F_RP_STRUCT *fatray, GRAV_STRUCT *grav, MT_STRUCT *mt , INV *inv, DATA_STRUCT *data, FLAG_STRUCT *flag, REGU_STRUCT *regu, CALC_2D_MT *calc_2D_mt, long nr_of_slices_2D_mt, char *fname, char *fname_vel_dens, char *fname_vel_res);
int SetupInvMatrix(DATA_STRUCT *data, GRAV_STRUCT *grav, MT_STRUCT *mt, INV *inv, FLAG_STRUCT *flag);
int FillSparseMatrix(FLAG_STRUCT flag, F_RP_STRUCT *fatray,RP_STRUCT *raypath, GRAV_STRUCT *grav, MT_STRUCT *mt, BIDX_STRUCT *inv_cell, DATA_STRUCT data, SPMAT *A);
int FillResdDataVector(FLAG_STRUCT flag,VEC *vec, DATA_STRUCT data, RP_STRUCT *raypath, F_RP_STRUCT *fatray, GRAV_STRUCT *grav, MT_STRUCT *mt);
int SetRegu(INV *inv, BIDX_STRUCT *inv_cell, SPMAT *A, FLAG_STRUCT flag, REGU_STRUCT *regu, double *res);
int SetRegu1(INV *inv, BIDX_STRUCT *inv_cell, SPMAT *A, FLAG_STRUCT flag, REGU_STRUCT *regu);
double RegFactor(SPMAT *A, int first_row, int last_row, int first_col, int last_col, int kind_of_reg ,long *c);
int MinMaxConstrainsSeis(INV *inv);
int MinMaxConstrainsGrav(INV *inv);
int MinMaxConstrainsRes(INV *inv);
int CalcResDensfromVel(double *seis_mod, double *mod2, long nvel, double vmin, double vmax, REL_VDR_STRUCT rel_vdr);
int MakeNewSlowMod(GRID_STRUCT *grid,long nvel, double *mod, BIDX_STRUCT *inv_cell);
int MakeNewDensResMod(GRID_STRUCT *grid,long nvel, double *mod, BIDX_STRUCT *inv_cell, int kind_of_model);
int CopyResultsToInvParaMod(INV *inv, FLAG_STRUCT *flag, VEC *v_mod, char *fname_vel_dens, char *fname_vel_res);
int AssignValueToForward(FLAG_STRUCT *flag, GRID_STRUCT *grid, INV *inv, BIDX_STRUCT *inv_cell);
int CalcNumBorder(INV *inv, BIDX_STRUCT *inv_cell);
int InvExtraCell(GRID_STRUCT *grid, INV *inv, GRAV_STRUCT *grav, SPMAT *A, DATA_STRUCT *data, REGU_STRUCT *regu, double *res, FLAG_STRUCT  flag, double weight_grav); /*BJOERN_MOD5*/
int TransExtraCellValue(INV *inv, FLAG_STRUCT *flag, VEC *v_mod, int index); /*(BJOERN_MOD)*/

/*in file make_inv_struct.c*/
int DetNumInvCells(GRID_STRUCT grid, INV *inv);
long MakeRegularInvStruct(GRID_STRUCT grid,BIDX_STRUCT *inv_cell,INV inv, FLAG_STRUCT flag);
int InvRayStruct(GRID_STRUCT grid,RP_STRUCT *raypath, INV *inv, DATA_STRUCT data);
int InvFatRayStruct(GRID_STRUCT grid,F_RP_STRUCT *fatray, INV *inv, DATA_STRUCT data);
FLEX_STRUCT ReadFlexGridIn(char *fname);
int MakeIrregularStructuresI(FLEX_STRUCT flex,BIDX_STRUCT *inv_cell,GRID_STRUCT grid, INV *inv, FLAG_STRUCT flag);
int MakeIrregularRayStructures(RP_STRUCT *raypath,F_RP_STRUCT *fatray, BIDX_STRUCT *inv_cell, long nr_inv_cell, long nrays, int kind_of_rays);
int InvGravStruct(GEOMETRY geo,DATA_STRUCT *data, GRAV_STRUCT *grav, BIDX_STRUCT *inv_cell, long nr_inv_cell_used, long nr_inv_cell);
int InvMTStruct(DATA_STRUCT *data, MT_STRUCT *mt, MT_STRUCT *forward_mt, BIDX_STRUCT *inv_cell, long nr_inv_cell, long nr_of_forward_cells, int kind_of_data_mt, int dimension_mt);
int ReadinTmpMTfile(double *derivTE, double *derivTM, long i_stat, long j_freq, long nr_of_cells, int kind_of_data_mt);

/*in file resolution_tests.c*/
int RayDensity(GRID_STRUCT grid,RP_STRUCT *raypath, BIDX_STRUCT *inv_cell, long nray, long ncell, double *weight_ray, int ninv);
int FatRayDensity(GRID_STRUCT grid,F_RP_STRUCT *fatray, BIDX_STRUCT *inv_cell, long nray, long ncell, double *weight_ray, int ninv);
int RayDensityTensor(GRID_STRUCT grid, RP_STRUCT *raypath, BIDX_STRUCT *inv_cell, long nray, long ncell,double *weight_ray, int ninv);
int DetSensSeisRays(INV *inv, RP_STRUCT *raypath, F_RP_STRUCT *fatray, int kind_of_rays ,BIDX_STRUCT *inv_cell, GRID_STRUCT *grid, DATA_STRUCT *data, GEOMETRY geo, int sens_in_perc);
int DetSensGrav(INV *inv, GRAV_STRUCT *grav, BIDX_STRUCT *inv_cell, GRID_STRUCT *grid, DATA_STRUCT *data, GEOMETRY geo, int sens_in_perc, int kind_of_para);
int DetSensMT(INV *inv, MT_STRUCT *mt,BIDX_STRUCT *inv_cell, GRID_STRUCT *grid, DATA_STRUCT *data, GEOMETRY geo, int sens_in_perc, int kind_of_para, FLAG_STRUCT *flag);

/*in file derivatives.c*/
double CalcDerivDensSlow(GRID_STRUCT *grid,double *Y2, double *Y1, REL_VDR_STRUCT rel_vdr);
int FillwithValues(GRID_STRUCT *grid, double *values, double eps, int kind_of_model);
int DerivativesGrav(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, double *Rho2, double *Rho1, double deltaSlow);
int Derivatives1DMT(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, double *R2, double *R1, double deltaSlow, MT_STRUCT *mt);
int Derivavtives2DMTAnaly(CALC_2D_MT *calc_2D_mt, long nr_of_slices_2D_mt, GRID_STRUCT grid, DATA_STRUCT *data, double *R2, double *R1, double deltaSlow, MT_STRUCT *mt, FLAG_STRUCT *flag);
int Deriv_MT_analytical(CALC_2D_MT calc_2D_mt, double *perturb, long i_layer,long j_freq, int kind_of_mt, double **d_Er, double **d_Ei, double **d_Hr,double **d_Hi, double grid_h);
int AssignRefinedGridMT(CALC_2D_MT calc_2D_mt, double *perturb, double grid_h, long nx, long nz, long j);
int SumSensToForwardCells(CALC_2D_MT calc_2D_mt, double *dSens, double grid_h, long nx, long nz, long j_freq);
int makeMTstructure2D(CALC_2D_MT calc_2D_mt, double *d_Zr, double *d_Zi, double *b_index, GRID_STRUCT grid, MT_STRUCT *mt, long j_freq, long k_stat, int index_mode, int direc_2D_mt, DATA_STRUCT *data,  int *already_used, long first_cell_strike);
int DetWaterColumn(double *waterdepths, double *res_water, double res_air, long z_int_stat_index, long x_int_stat_index, double *res_slice, CALC_2D_MT calc_2D_mt, double first_cell_dz, long nz, long j_freq, int *index_no_water, int flag_ray_based);
int DetMeanResWithinLayer(CALC_2D_MT calc_2D_mt,double *res_slice, double *average_res, double *dz, long z_int_stat_index, double first_cell_dz, long iz, long nz, long j_freq, double *station_res, double x_index, double x_index_stat, int kind_of_calc);
int CalcBx_Ey(long nr_of_samples, double *p,  double *fft_input_Bx, double *fft_input_Ey, double *average_res, double *hz, long nhz, double freq, double waterdepths, double res_water_layer, int index_water_layer, double sample_int, double res_station, int index_imped);
int CalcBy_Ex_Ix(long nr_of_samples, long nr_of_samples_Ex, double *p, double *p_Ex, double *fft_input_Bx, double *fft_input_Ey, double *int_res, double *hz, long nhz, double freq, double waterdepths, double res_water_layer, double sample_int, double sample_int_Ex, double res_station, int index_imped);
int CalcBy_Ex_Iz(long nr_of_samples, long nr_of_samples_Ex, double *p, double *p_Ex, double *fft_input_By, double *fft_input_Ex, double *int_res, double *hz, long nhz, double freq, double waterdepths, double res_water_layer, double sample_int, double sample_int_Ex, double res_station, int index_imped);
int Reflex_TE(double waterdepths, double theta_r, double theta_i, double theta_0_r, double theta_0_i, double theta_air_r, double theta_air_i, double *R_TE_r, double *R_TE_i, int index_water_layer);
int Reflex_TM(double waterdepths, double sigma, double sigma_0, double theta_r, double theta_i, double theta_0_r, double theta_0_i, double *R_TM_r, double *R_TM_i);

/*Numerical recipes*/
/*in dsvdcmp.c*/
void dsvdcmp(double **a, int m, int n, double w[], double **v);
/*in dpythag.c*/
double dpythag(double a, double b);
/*in sort2.c*/
void sort2(unsigned long n, double arr[], double brr[]);
/*in dfour1.c*/
void four1(double *data,unsigned long nn, int isign);

/*Remark: - in dcomplex.c/d.complex.h are functions defined/declarated to work with complex numbers*/
/*        - in nrutil.c/nrutil.h are some functions defined/declarated that are required for the numerical recipes functions*/   


/************************************************/
/*Fortran Routinen (modified from Tartis)*/
/*in files B_hpol.for and B_epol.for*/
extern void HPOL (double* T, long* NX, long* NZ, double* DXJ, double* DZI, double* RO, double* HY_REAL, double* HY_IMAG, double* EX_REAL, double* EX_IMAG, double* EZ_REAL, double* EZ_IMAG);
extern void EPOL (double* T, long* NX, long* NZ, long* IONOS, long* IATM, double* DXJ, double* DZI, double* RO, double* RIONOS, double* HX_REAL, double* HY_IMAG, double* EY_REAL, double* EY_IMAG);

#define _CALLING_
#endif
