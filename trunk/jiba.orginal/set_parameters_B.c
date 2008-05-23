#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "inv3d.h"
#include "fkt.h"

/*---------------------------------------------------------------------------*/
/* Set flags from the input parameters*/
/* parameters:  par       := Structure of input parameters */                                     
/*				*flag	  := flag structure*/

int set_flags(PAR par,FLAG_STRUCT *flag)
{
   flag->index_tseis = par.index_tseis;
   flag->index_grav = par.index_grav;
   flag->index_mt = par.index_mt;

	/*Specifying the kind of joint inversions*/
   if(flag->index_grav != 0 && flag->index_tseis != 0 && flag->index_mt == 0)
	   flag->index_tseis_grav = 1;
   else
	   flag->index_tseis_grav = 0;

   if(flag->index_mt != 0 && flag->index_tseis != 0 && flag->index_grav == 0)
	   flag->index_tseis_mt = 1;
   else
	   flag->index_tseis_mt = 0;

   if(flag->index_grav != 0 && flag->index_mt != 0 && flag->index_tseis == 0)
	   flag->index_grav_mt = 1;
   else
	   flag->index_grav_mt = 0;

   if(flag->index_tseis != 0 && flag->index_grav != 0 && flag->index_mt != 0)
	   flag->index_tseis_grav_mt = 1;
   else
	   flag->index_tseis_grav_mt = 0;

   if(flag->index_tseis_grav != 0 || flag->index_tseis_mt != 0 || flag->index_grav_mt != 0 || flag->index_tseis_grav_mt != 0)
	   flag->index_joint = 1;
   else	
	   flag->index_joint = 0;

   if(flag->index_grav != 0 && flag->index_tseis != 0)
		flag->link_seis_grav_mod = 1; /*The gravity model will be calculated by means of the seismic model in the joint inversion*/
   else
	   flag->link_seis_grav_mod = 0;

   if(flag->index_mt != 0 && flag->index_tseis != 0)
		flag->link_seis_mt_mod = 1; /*The resistiviy model will be calculated by means of the seismic model in the joint inversion*/
   else if(flag->index_mt != 0 && flag->index_grav != 0)
		flag->link_seis_mt_mod = 2; /*The resistiviy model will be calculated by means of the gravity model in the joint inversion*/
   else
		flag->link_seis_mt_mod = 0;


	/*Specifing the seismic tomography related flags*/
   flag->kind_of_rays = par.ray_index;
   flag->balance_para_rays = par.para_rays_index; 
   flag->ray_density = par.ray_density_index;
   flag->ray_density2 = par.ray_density_index2;

   /*Specifying the MT related flags*/
		/*1-D MT modeling*/
   if(par.dimension_mt == 1)
		flag->dimension_mt = 1;
		/*2-D MT modeling*/
   else if(par.dimension_mt == 2)
		flag->dimension_mt = 2;
		/*2-D MT forward modeling and 1-D inversion*/
   else if(par.dimension_mt == 2)
	   flag->dimension_mt = 3;
		/*2-D MT forward modeling and 1-D inversion using RRI (see Smith and Booker; Journal of Geophysical Research; 1991)*/
   else
	   flag->dimension_mt = 4;


   flag->kind_of_data_mt = par.kind_of_data_mt;
   flag->index_mt_weighting = par.index_mt_weighting; 

   if(flag->dimension_mt != 1)
   {
   		flag->direc_2D_mt = par.direc_2D_mt;
		flag->nr_cells_2d_mt = par.nr_cells_2d_mt;
   }

   /*Specifying the inversion related flags*/
   if(par.damp == 0.0 || par.smoothfac == 0.0)
	   flag->do_damp = 0;
   else 
	   flag->do_damp = 1;
   if(par.damp == 0.0 || par.smoothfac == 1.0) 
	   flag->do_smooth = 0;
   else 
	   flag->do_smooth = 1;
   flag->flex_grid_size = par.flex_grid_size;
   
   /*Jin Chen changed new*/
   /*Specify the value of do_extra_smooth, yes=1, no=0.*/
   if(par.extra_smooth == 0)
	   flag->do_extra_smooth = 0;
   else
	   flag->do_extra_smooth = 1;

   flag->write_sens_out = par.write_sens_out;
   flag->sens_in_perc = par.sens_in_perc;

	return(1);
}

/*---------------------------------------------------------------------------*/
/* Set regularization from the input parameters*/
/* parameters:  par       := Structure of input parameters */                                     
/*				*regu	  := Regularization structure*/

int set_regu(PAR par,REGU_STRUCT  *regu)
{
	/*Jin Chen changed new*/
	regu->damp = par.damp;
	regu->smooth = par.smoothfac;
	regu->val_extra_smooth = par.val_extra_smooth; /*Set the extra smoothing value which the user can gives.*/
	/*Jin Chen changed 02.2007*/
	regu->dist_extra_smooth = par.dist_extra_smooth; /*Set the distance value.*/
	/*extra cell*/
	regu->w_extra_cell = par.extra_cell; /*Set the weighting value of the extra cell.*/
	/*BJOERN_MOD5*/
	regu->smoothratio_y = par.smoothratio_y;
	regu->smoothratio_z = par.smoothratio_z;

	return(1);
}


/*---------------------------------------------------------------------------*/
/* Set inversion parameters from the input parameters*/
/* parameters:  par       := Structure of input parameters */                                     
/*				*inv	  := Inversion structure*/

int set_inv(PAR par, INV *inv)
{
	inv->ninv = 0;
	inv->ninv_max = par.inv_max;
	
	/*Velocity*/
	if(par.vmina <= 0.0)
		inv->vmin = 0.0;
	else
		inv->vmin = par.vmina;
	
	if(par.vmaxa <= 0.0)
		inv->vmax = 1.0e+12;
	else
		inv->vmax = par.vmaxa;
	
	if(par.maxadj <= 0.0)
		inv->maxadj = 1.0e+12;
	else
		inv->maxadj = par.maxadj;

	inv->angle_para = par.angle_para; /*Angle between two rays (in degree), where they are considered as "parallel"*/

	/*Density*/
	if(par.gmina <= 0.0)
		inv->gmin = 0.0;
	else
		inv->gmin = par.gmina;

	if(par.gmaxa <= 0.0)
		inv->gmax = 1.0e+12;
	else
		inv->gmax = par.gmaxa;

	/*Resistivity*/
	if(par.rmina <= 0.0)
		inv->rmin = 0.0;
	else
		inv->rmin = par.rmina;

	if(par.rmaxa <= 0.0)
		inv->rmax = 1.0e+20;
	else
		inv->rmax = par.rmaxa;


	/*set the relative scaling of the data sets to default*/
	inv->rel_scal_seis = 1;
	inv->rel_scal_grav = 1;
	inv->rel_scal_mt = 1;
	
	/*extra cell for the gravity inversion (BJOERN_MOD)*/
	inv->dens_extra_cell = 0; /*Set the initial density value of the extra cell.*/

	return(1);
}
