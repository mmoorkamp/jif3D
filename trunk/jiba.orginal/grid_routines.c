#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "inv3d.h"
#include "fkt.h"

/*----------------------------------------------------------------*/
/*Assign the variables to the structure GRID_STRUCT*/
/* Parameter:  par			:= Input parameter in the parameter file*/
/*			  *grid_struct	:= Pointer on grid structure that will be modified here*/

int MakeGridStruct(PAR par,GRID_STRUCT *grid_struct)
{
	(grid_struct->nx) = par.nx;
	(grid_struct->ny) = par.ny;
	(grid_struct->nz) = par.nz;
	(grid_struct->h)  = par.h;
	(grid_struct->org[0]) = par.org[0];
	(grid_struct->org[1]) = par.org[1];
	(grid_struct->org[2]) = par.org[2];
	(grid_struct->topo_index) = par.topo_index;
	(grid_struct->grid_cell_factor[0]) = par.grid_cell_factor[0];
	(grid_struct->grid_cell_factor[1]) = par.grid_cell_factor[1];
	(grid_struct->grid_cell_factor[2]) = par.grid_cell_factor[2];
	grid_struct->border_index = (int *)memory(NULL,1,sizeof(int),"MakeGridStruct");

	/*seismic parameters*/
	(grid_struct->nborder) = par.nborder;
	(grid_struct->vborder) = par.vborder;
	(grid_struct->v_air_water) = par.v_air_water;
	grid_struct->slow = (double *)memory(NULL,1,sizeof(double),"MakeGridStruct");

	/*gravity parameters*/
	(grid_struct->dens_air_water) = par.g_air_water;
	grid_struct->dens = (double *)memory(NULL,1,sizeof(double),"MakeGridStruct");

	/*Electrical parameters*/
	(grid_struct->res_air_water) = par.r_air_water;
	grid_struct->res = (double *)memory(NULL,1,sizeof(double),"MakeGridStruct");


	return(1);
}

/*----------------------------------------------------------------*/
/* Assign the variables to the structure GRADIENT (can be used for the velocity, density or resistivity field)*/
/* Parameter: 	tmp_index_topo          			:= This index governs if the field is determined by the surface topography(1) or */
/*													   by a defined gridpoint "org"*/
/*				tmp_gx,tmp_gy,tmp_gz    			:= Gradient of velocity (m/s)/m, density (g/cm^3)/m or resistivity */
/*		        tmp_orgx,tmp_orgy,tmp_orgz          := Position where "value_0" is defined*/
/*				tmp_0								:= Starting value in m/s (velocity), g/cm^3 (density) or ohmm (resistivity)  (at the defined gridpoint "org_vo" or at the surface)*/						
/*              tmp_min,tmp_max						:= Minimum and maximum value (velocity[m/s], density[g/cm^3], resistivity[ohmm]) for the starting model*/

int MakeGradientStruct(short tmp_index_topo, float tmp_gx,float tmp_gy,float tmp_gz,float tmp_orgx, float tmp_orgy,float tmp_orgz,float tmp_0, float tmp_min,float tmp_max ,GRADIENT *grad)
{
	(grad->index_topo) = tmp_index_topo;
	(grad->g[0]) = tmp_gx;
	(grad->g[1]) = tmp_gy;
	(grad->g[2]) = tmp_gz;
	(grad->org[0]) = tmp_orgx;
	(grad->org[1]) = tmp_orgy;
	(grad->org[2]) = tmp_orgz;
	(grad->value_0) = tmp_0;
	(grad->min) = tmp_min;
	(grad->max) = tmp_max;

	return(1);
}