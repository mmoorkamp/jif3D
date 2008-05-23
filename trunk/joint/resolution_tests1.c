#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "nrutil.h"
#include "meschachd_no_complex/iter.h"
#include "meschachd_no_complex/sparse.h"
#include "inv3d.h"
#include "fkt.h"

/*------------------------------------------------------------*/
/*Calculate the hit matrix and the DWS (derivative weighted sum) for rays*/
/*Parameter:    *raypath		:= Raypath structure*/
/*				*inv_cell		:= Structure of the inversion cells*/
/*				 nray			:= Number of shot-receiver combinations*/
/*				 ncell			:= Number of inversion cells*/
/*				 ninv			:= Number of the iteration*/
/*              *weight_ray     := weightings of the rays*/
/*				 grid			:= Grid structure (of forward cells)*/

int RayDensity(GRID_STRUCT grid, RP_STRUCT *raypath, BIDX_STRUCT *inv_cell, long nray, long ncell, double *weight_ray, int ninv)
{
	long i,j;
	long *nrays_per_cell;	/*Number of rays in an inversion cell*/
	double *length;			/*sum of the normalized lengths of all rays in the inversion cell*/
	double *total_length;	/*total ray length*/
	double total_weight;	/*average weight of all rays*/

	length= (double *)memory(NULL,ncell,sizeof(double),"RayDensity");
	if(nray != 0)
		total_length= (double *)memory(NULL,nray,sizeof(double),"RayDensity");
	else
		total_length= (double *)memory(NULL,1,sizeof(double),"RayDensity");
	nrays_per_cell = (long *)memory(NULL,ncell,sizeof(long),"RayDensity");

	for(i=0;i<ncell;i++)
	{
		length[i]=0.0;
		nrays_per_cell[i]=0;
	}

	total_weight = 0.0;

	/*Determine the total length of the rays for normalization*/
	/*Loop over all rays*/
	for(i=0;i<nray;i++)
	{
		total_weight = total_weight + weight_ray[i];
		total_length[i] = 0.0;
		
		/*Loop over all ray segments*/
		for(j=0;j<raypath[i].nray;j++)
		{
			total_length[i] = total_length[i] + raypath[i].len[j];
		}

	}

	if(nray <= 0)
	{
		printf("WARNING!! The number of rays is zero or negative while determing the DWS-matrix\n");
	}
	else
	{

		total_weight = total_weight/nray;

		if(total_weight == 0.0)
		{
					printf("The cummulative weight of all rays is zero\n");
					printf("This makes no sense!!!\n");
					exit(0);
		}

		/*Loop over all rays*/
		for(i=0;i<nray;i++)
		{
			/*Loop over all ray segments*/
			for(j=0;j<raypath[i].nray;j++)
			{

				if(raypath[i].ele[j] >= ncell)
				{
					printf("NOT enough memory is allocated: used %d, allocated %d\n",raypath[i].ele[j]+1, ncell);
					exit(0);
				}

				if(total_length[i] > 0)

					length[raypath[i].ele[j]] = ((weight_ray[i] * raypath[i].len[j])/total_length[i]) + length[raypath[i].ele[j]];
				nrays_per_cell[raypath[i].ele[j]]++;

			}

		}

		/*Loop over all cells*/
		for(i=0;i<ncell;i++)
		{
			length[i] = length[i]/total_weight;
		}

	}

	/*Write out the ray densities*/
	WriteRayDenseOut(grid, ninv, ncell, nrays_per_cell, length, inv_cell);

	free(nrays_per_cell);
	free(length);
	free(total_length);

	return(1);
}



/*------------------------------------------------------------*/
/*Calculate the hit matrix and the DWS (derivative weighted sum) for fat-rays*/
/*Parameter:    *fatray			:= Fatray structure*/
/*				*inv_cell		:= Structure of the inversion cells*/
/*				 nray			:= Number of shot-receiver combinations*/
/*				 ncell			:= Number of inversion cells*/
/*              *weight_ray     := weightings of the rays*/
/*				 ninv			:= Number of the iteration*/
/*				 grid			:= Grid structure (of forward cells)*/

int FatRayDensity(GRID_STRUCT grid, F_RP_STRUCT *fatray, BIDX_STRUCT *inv_cell, long nray, long ncell, double *weight_ray, int ninv)
{
	long i,j;
	long *nrays_per_cell;	/*Number of fat-rays in an inversion cell*/
	double *over_all_size;	/*sum of all normalized fatray weights in a inversion cell*/
	double *total_size; /*The total sum of all weights*/
	double total_weight;	/*average weight of all fat-rays*/

	over_all_size= (double *)memory(NULL,ncell,sizeof(double),"FatRayDensity");
	if(nray != 0)
		total_size= (double *)memory(NULL,nray,sizeof(double),"FatRayDensity");
	else
		total_size= (double *)memory(NULL,1,sizeof(double),"FatRayDensity");
	nrays_per_cell = (long *)memory(NULL,ncell,sizeof(long),"FatRayDensity");

	for(i=0;i<ncell;i++)
	{
		over_all_size[i]=0.0;
		nrays_per_cell[i]=0;
	}

	total_weight = 0.0;
	/*Determine the total weight of the fatrays for the normalization*/
	/*Loop over all rays*/
	for(i=0;i<nray;i++)
	{
		total_weight = total_weight + weight_ray[i];
		total_size[i] = 0.0;
		
		/*Loop over all ray segments*/
		for(j=0;j<fatray[i].ncell;j++)
		{
			total_size[i] = total_size[i] + fatray[i].weight[j];
		}

	}


	if(nray <= 0)
	{
		printf("WARNING!! The number of fat-rays is zero or negative while determing the DWS-matrix\n");
	}
	else
	{
		total_weight = total_weight/nray;

		if(total_weight == 0.0)
		{
					printf("The cummulative weight of all rays is zero\n");
					printf("This makes no sense!!!\n");
					exit(0);
		}


		/*Loop over all fatrays*/
		for(i=0;i<nray;i++)
		{
			/*Loop over all fatray segments*/
			for(j=0;j<fatray[i].ncell;j++)
			{

				if(fatray[i].ele[j] >= ncell)
				{
					printf("NOT enough memory is allocated: used %d, allocated %d\n",fatray[i].ele[j]+1, ncell);
					exit(0);
				}

				if(total_size[i] >= 0)
					over_all_size[fatray[i].ele[j]] = ((weight_ray[i] *fatray[i].weight[j])/total_size[i]) + over_all_size[fatray[i].ele[j]];
				nrays_per_cell[fatray[i].ele[j]] = nrays_per_cell[fatray[i].ele[j]] + 1;

			}

		}

		/*Loop over all cells*/
		for(i=0;i<ncell;i++)
		{
			over_all_size[i] = over_all_size[i]/total_weight;
		}
	}

	/*Write out the ray densities*/
	WriteRayDenseOut( grid, ninv, ncell, nrays_per_cell, over_all_size, inv_cell);

	free(nrays_per_cell);
	free(over_all_size);
	free(total_size);

	return(1);
}



/*------------------------------------------------------------*/
/*Calculate the ray density tensor for rays*/
/*Parameter:    *raypath		:= Raypath structure*/
/*				*inv_cell		:= Structure of the inversion cells*/
/*				 nray			:= Number of shot-receiver combinations*/
/*				 ncell			:= Number of inversion cells*/
/*              *weight_ray     := weightings of the rays*/
/*				 ninv			:= Number of the iteration*/
/*				 grid			:= Grid structure (of forward cells)*/

int RayDensityTensor( GRID_STRUCT grid, RP_STRUCT *raypath, BIDX_STRUCT *inv_cell, long nray, long ncell, double *weight_ray, int ninv)
{
	long i,j,k;
	double ***ray_density_matrix; /*First index:  pos. of inversion cells*/
								  /*Second index: nr. of rows; the three components of the ray-segments*/
								  /*Third index:  nr. of columns; number of rays intersecting the cells*/
	double **A;					  /*Input matrix for SVD*/
	double **V;					  /*The second eigenvector matrix in the SVD*/
	double *total_length;		  /*Total length of the rays*/
	long *rays_per_cell;		  /*Number of rays in the inversion cells*/
	double eig[4];				  /*singular values resp. eigenvalues*/
	double *rel_minmax_eig;			  /*Index 1: Relation of the smallest/largest eigenvalue for all inversion cells*/
	double *rel_medmax_eig;			  /*Index 1: Relation of the smallest/largest eigenvalue for all inversion cells*/


	/*Allocate the memory for the matrices G used for the determination of the ray density tensor GG^t*/
	 ray_density_matrix = (double ***)memory(NULL,ncell,sizeof(double **),"RayDensityTensor");

	for(i=0;i<ncell;i++)
		ray_density_matrix[i] = (double **)memory(NULL,4,sizeof(double *),"RayDensityTensor"); 
		/*Remark: The svd-routine of num.rec. expects matrices with indices starting at 1*/

	for(i=0;i<ncell;i++)
		for(j=0;j<4;j++)
		{
			ray_density_matrix[i][j] = (double *)memory(NULL,1,sizeof(double),"RayDensityTensor");
		}

	for(i=0;i<ncell;i++)
		for(j=0;j<4;j++)
		{
				ray_density_matrix[i][j][0] = 0.0; 
		}

	/*Allocate memory for the number of rays per cell*/
	rays_per_cell = (long *)memory(NULL,ncell,sizeof(long),"RaydensityTensor");
	for(i=0;i<ncell;i++)
		rays_per_cell[i] = 0;

	/*Allocate memory the for total length of each ray*/
	if(nray != 0)
		total_length = (double *)memory(NULL,nray,sizeof(double),"RaydensityTensor");
	else
		total_length = (double *)memory(NULL,1,sizeof(double),"RaydensityTensor");

	/*Allocate memory for the relation of the different eigenvalues for each cell*/
	rel_minmax_eig = (double *)memory(NULL,ncell,sizeof(double),"RayDensityTensor");
	for(i=0;i<ncell;i++)
		rel_minmax_eig[i] = 0.0;

	rel_medmax_eig = (double *)memory(NULL,ncell,sizeof(double),"RayDensityTensor");
	for(i=0;i<ncell;i++)
		rel_medmax_eig[i] = 0.0;

	/*Determine the total length of the rays for normalization*/
	/*Loop over all rays*/
	for(i=0;i<nray;i++)
	{
		total_length[i] = 0.0;
		
		/*Loop over all ray segments*/
		for(j=0;j<raypath[i].nray;j++)
		{
			total_length[i] = total_length[i] + raypath[i].len[j];
		}

	}

	/*****************************************************************/
	/*Calculate the ray density matrices*/
	/*Loop over all rays*/

		for(i=0;i<nray;i++)
		{

			/*Loop over all ray segments*/
			for(j=0;j<raypath[i].nray;j++)
			{

			

				if(raypath[i].ele[j] >= ncell)
				{
					printf("NOT enough memory is allocated: used %d, allocated %d\n",raypath[i].ele[j]+1, ncell);
					exit(0);
				}

				/*Reallocate the memory for G*/
				rays_per_cell[raypath[i].ele[j]]++;

				ray_density_matrix[raypath[i].ele[j]][0] = (double *)memory((char *)ray_density_matrix[raypath[i].ele[j]][0],rays_per_cell[raypath[i].ele[j]]+1,sizeof(double),"RaydensityTensor");
				ray_density_matrix[raypath[i].ele[j]][1] = (double *)memory((char *)ray_density_matrix[raypath[i].ele[j]][1],rays_per_cell[raypath[i].ele[j]]+1,sizeof(double),"RaydensityTensor");
				ray_density_matrix[raypath[i].ele[j]][2] = (double *)memory((char *)ray_density_matrix[raypath[i].ele[j]][2],rays_per_cell[raypath[i].ele[j]]+1,sizeof(double),"RaydensityTensor");
				ray_density_matrix[raypath[i].ele[j]][3] = (double *)memory((char *)ray_density_matrix[raypath[i].ele[j]][3],rays_per_cell[raypath[i].ele[j]]+1,sizeof(double),"RaydensityTensor");
				/*Remark: The svd-routine of num.rec. expects matrices with indices starting at 1*/


				ray_density_matrix[raypath[i].ele[j]][0][rays_per_cell[raypath[i].ele[j]]] = 0.0;
				/*Remark: The svd-routine of num.rec. expects matrices with indices starting at 1*/


				if(total_length[i] >= 0)
				{	
					ray_density_matrix[raypath[i].ele[j]][1][rays_per_cell[raypath[i].ele[j]]] = (raypath[i].x[j]*weight_ray[i])/total_length[i];
					ray_density_matrix[raypath[i].ele[j]][2][rays_per_cell[raypath[i].ele[j]]] = (raypath[i].y[j]*weight_ray[i])/total_length[i];
					ray_density_matrix[raypath[i].ele[j]][3][rays_per_cell[raypath[i].ele[j]]] = (raypath[i].z[j]*weight_ray[i])/total_length[i];
				}
				else
				{	
					ray_density_matrix[raypath[i].ele[j]][1][rays_per_cell[raypath[i].ele[j]]] = 0.0;
					ray_density_matrix[raypath[i].ele[j]][2][rays_per_cell[raypath[i].ele[j]]] = 0.0;
					ray_density_matrix[raypath[i].ele[j]][3][rays_per_cell[raypath[i].ele[j]]] = 0.0;
				
				}

			}
		}

	/*****************************************************************/
	/*Perform the singular value decomposition*/
	for(i=0;i<ncell;i++)
	{
		/*Switch matrix indeces of the input matrix*/
		A = (double **)memory(NULL,rays_per_cell[i]+1,sizeof(double *),"RayDensityTensor");
		for(j=0;j<rays_per_cell[i]+1;j++)
		{
			A[j] = (double *)memory(NULL,4,sizeof(double),"RayDensityTensor");

			for(k=0;k<4;k++)
				A[j][k] = ray_density_matrix[i][k][j];
		}

		/*Allocate memory for the second eigenvector matrix:*/
		V = (double **)memory(NULL,4,sizeof(double *),"RayDensityTensor");
		for(j=0;j<4;j++)
		{
			V[j] = (double *)memory(NULL,4,sizeof(double),"RayDensityTensor");
		}

		/*calculate singular values*/
		dsvdcmp(A,rays_per_cell[i],3,eig,V);


		for(j=1;j<4;j++)
		{
			/*Calculate the eigenvalues from the singular values*/
			eig[j] = sqrt(eig[j]);
		}

		/*Sort the eigenvalues*/
		sort2(3,eig,eig);

		/*Calculate the relation of the smallest(medium)/largest eigenvalue*/
		if(eig[3]!=0.0)
		{
			rel_minmax_eig[i] = eig[1]/eig[3];
			rel_medmax_eig[i] = eig[2]/eig[3];

			if(rel_minmax_eig[i] > 1.0 || rel_medmax_eig[i] > 1.0)
			{
				printf("The relation of the eigenvales is not correct!\n");
				printf("First EV:   %f\n",eig[3]);
				printf("Second EV:  %f\n",eig[2]);
				printf("Third EV:	%f\n",eig[1]);
				exit(0);
			}
		}
		
		/*Free temporary matrices*/
		for(j=0;j<rays_per_cell[i]+1;j++)
			free(A[j]);
		free(A);

		for(j=0;j<4;j++)
			free(V[j]);
		free(V);

	}

	/*****************************************************************/
	/*Write out the ray density tensor*/
	WriteRayDensityTensorOut( grid, ninv, ncell, rel_minmax_eig, rel_medmax_eig, inv_cell);

	/*free the memory*/
   for(i=0;i<ncell;i++)
		for(j=0;j<4;j++)
		{
			free(ray_density_matrix[i][j]);
		}

	for(i=0;i<ncell;i++)
		free(ray_density_matrix[i]);

	free(ray_density_matrix);
	free(rays_per_cell);
	free(total_length);

	free(rel_minmax_eig);
	free(rel_medmax_eig);


	return(1);
}


/*------------------------------------------------------------*/
/*Sort and write out the total (scaled) sensitivities for the conventional rays/fat-rays*/
/*Parameter:  *inv = Inversion structure*/
/*			*raypath = Raypath structure for the conventional rays*/
/*			*fatray =  Structure for fat rays*/
/*		kind_of_rays = Kind of the used cells*/
/*          *inv_cell = Inversion cell structure (relating the forward and the inversion cells to each other)*/
/*			   *grid  = Grid structure */
/*           *data    = data structure*/
/*             geo	  = Geometry*/
/*		sens_in_perc  = The sensitivities are expressed in percentage (yes=1/no=0)*/

int DetSensSeisRays(INV *inv, RP_STRUCT *raypath, F_RP_STRUCT *fatray, int kind_of_rays ,BIDX_STRUCT *inv_cell, GRID_STRUCT *grid, DATA_STRUCT *data, GEOMETRY geo, int sens_in_perc)
{
	long i,j,k,nxyz;
	double *sens_seis; /*sensitivities from the seismics*/

	/*Allocate the memory*/
	nxyz = (2*grid->nborder + grid->nx) * (2*grid->nborder + grid->ny) * (2*grid->nborder + grid->nz); 
	sens_seis = (double *)memory(NULL,nxyz,sizeof(double),"DetSensSeisRays");

	for(i=0;i<nxyz;i++)
		sens_seis[i] =0.0;

	if(kind_of_rays == 1) /*conventional rays*/
	{
		/*Loop over all rays*/
		for(i=0;i<data->ndata_seis;i++)
		{
			/*Loop over all ray segments*/
			for(j=0;j<raypath[i].nray;j++)
			{
				/*Inversion cell has to active*/
				if(inv_cell[raypath[i].ele[j]].use == 1)
				{
					for(k=0;k<inv_cell[raypath[i].ele[j]].nele;k++)
					{
						if(sens_in_perc == 0) /*NOT normalized*/
							/*Sum all sensitvities*/
							sens_seis[inv_cell[raypath[i].ele[j]].ele[k]] = sens_seis[inv_cell[raypath[i].ele[j]].ele[k]] + (inv->rel_scal_seis * fabs(raypath[i].len[j]));
						else /*normalized (sensitivities are expressed in percentage)*/
							/*Sum all sensitvities*/
							sens_seis[inv_cell[raypath[i].ele[j]].ele[k]] = sens_seis[inv_cell[raypath[i].ele[j]].ele[k]] + (100 * fabs((raypath[i].len[j] * inv_cell[raypath[i].ele[j]].val_slow)/data->tcalc[i]));
					}
				}

			}
		}
	}
	else /*Fat rays*/
	{
		/*Loop over all rays*/
		for(i=0;i<data->ndata_seis;i++)
		{
			/*Loop over all ray segments*/
			for(j=0;j<fatray[i].ncell;j++)
			{
				/*Inversion cell has to be active*/
				if(inv_cell[fatray[i].ele[j]].use == 1)
				{
					for(k=0;k<inv_cell[fatray[i].ele[j]].nele;k++)
					{
						if(sens_in_perc == 0) /*NOT normalized*/
							/*Sum all sensitvities*/
							sens_seis[inv_cell[fatray[i].ele[j]].ele[k]] = sens_seis[inv_cell[fatray[i].ele[j]].ele[k]] + (inv->rel_scal_seis * fabs(fatray[i].weight[j]));
						else /*normalized (sensitivities are expressed in percentage)*/
							/*Sum all sensitvities*/
							sens_seis[inv_cell[fatray[i].ele[j]].ele[k]] = sens_seis[inv_cell[fatray[i].ele[j]].ele[k]] + (100 * fabs((fatray[i].weight[j] * inv_cell[fatray[i].ele[j]].val_slow)/data->tcalc[i]));

					}
					
				}
				
			}
		}
	}

	/*Write the sensitivities in a binary model file*/
	WriteSensSeisOut(inv->ninv, *data, geo, *grid, sens_seis);
	
	free(sens_seis);
	
	return(1);
}

/*------------------------------------------------------------*/
/*Sort and write out the total (scaled) sensitivities for the gravity measurements*/
/*Parameter:  *inv = Inversion structure*/
/*			  *grav = Gravity structure*/
/*          *inv_cell = Inversion cell structure (relating the forward and the inversion cells to each other)*/
/*			   *grid  = Grid structure */
/*           *data    = data structure*/
/*             geo	  = Geometry*/
/*		sens_in_perc  = The sensitivities are expressed in percentage (yes=1/no=0)*/
/*		kind_of_para  = Kind of model parameter (1= slowness; 2= density)*/

int DetSensGrav(INV *inv, GRAV_STRUCT *grav,BIDX_STRUCT *inv_cell, GRID_STRUCT *grid, DATA_STRUCT *data, GEOMETRY geo, int sens_in_perc, int kind_of_para)
{
	long i,j,k,nxyz;
	double *sens_grav;
	
	/*Allocate the memory*/
	nxyz = (2*grid->nborder + grid->nx) * (2*grid->nborder + grid->ny) * (2*grid->nborder + grid->nz); 
	sens_grav = (double *)memory(NULL,nxyz,sizeof(double),"DetSensSeisRays");

	for(i=0;i<nxyz;i++)
		sens_grav[i] =0.0;

	/*Loop over all gravity stations*/
	for(i=0;i<data->ndata_grav;i++)
	{
		/*Loop over all cells*/
		for(j=0;j<grav[i].ncell;j++)
		{
				/*Inversion cell has to be active*/
				if(inv_cell[grav[i].ele[j]].use == 1)
				{
					for(k=0;k<inv_cell[grav[i].ele[j]].nele;k++)
					{
						if(sens_in_perc == 0) /*NOT normalized*/
							/*Sum all sensitvities*/
							sens_grav[inv_cell[grav[i].ele[j]].ele[k]] = sens_grav[inv_cell[grav[i].ele[j]].ele[k]] + (inv->rel_scal_grav * fabs(grav[i].deriv[j]));
						else /*normalized (sensitivities are expressed in percentage)*/
						{
							if(kind_of_para == 1)
								sens_grav[inv_cell[grav[i].ele[j]].ele[k]] = sens_grav[inv_cell[grav[i].ele[j]].ele[k]] + (100 * fabs((grav[i].deriv[j] * inv_cell[grav[i].ele[j]].val_slow)/data->calc_grav[i]));
							else
								sens_grav[inv_cell[grav[i].ele[j]].ele[k]] = sens_grav[inv_cell[grav[i].ele[j]].ele[k]] + (100 * fabs((grav[i].deriv[j] * inv_cell[grav[i].ele[j]].val_dens)/data->calc_grav[i]));								
						}
					}
					
				}
		}
	}

	/*Write the sensitivities in a binary model file*/
	WriteSensGravOut(inv->ninv, *data, geo, *grid, sens_grav);

	free(sens_grav);

	return(1);
}

/*------------------------------------------------------------*/
/*Sort and write out the total (scaled) sensitivities for the MT measurements/
/*Parameter:  *inv = Inversion structure*/
/*			  *mt = Structure structure*/
/*          *inv_cell = Inversion cell structure (relating the forward and the inversion cells to each other)*/
/*			   *grid  = Grid structure */
/*           *data    = data structure*/
/*             geo	  = Geometry*/
/*		sens_in_perc  = The sensitivities are expressed in percentage (yes=1/no=0)*/
/*		kind_of_para  = Kind of model parameter (1= slowness; 2= resistivity)*/
/*              flag  = structure including general flags*/

int DetSensMT(INV *inv, MT_STRUCT *mt,BIDX_STRUCT *inv_cell, GRID_STRUCT *grid, DATA_STRUCT *data, GEOMETRY geo, int sens_in_perc, int kind_of_para, FLAG_STRUCT *flag)
{
	long i,j,k,m,nxyz;
	double *sens_mt;
	
	/*Allocate the memory*/
	nxyz = (2*grid->nborder + grid->nx) * (2*grid->nborder + grid->ny) * (2*grid->nborder + grid->nz); 
	sens_mt = (double *)memory(NULL,nxyz,sizeof(double),"DetSensSeisRays");

	for(i=0;i<nxyz;i++)
		sens_mt[i] =0.0;

	/*Loop over all mt stations*/
	for(i=0;i<data->ndata_mt;i++)
	{
		/*Loop over all inv. cells*/
		for(j=0;j<mt[i].ncell;j++)
		{
			/*Inversion cell has to be active*/
			if(inv_cell[mt[i].ele[j]].use == 1)
			{
				/*Loop over all forward cells in the inversion cell*/
				for(k=0;k<inv_cell[mt[i].ele[j]].nele;k++)
				{

					/*********************************************/
					/*Both modes are used modes*/
					if(flag->nr_of_modes == 2 &&  flag->dimension_mt == 2)
					{
						for(m=0;m<4*mt[i].nfreq;m++)
						{
							if(sens_in_perc == 0) /*NOT normalized*/
								/*Sum all sensitvities*/
								sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (inv->rel_scal_mt * fabs(mt[i].deriv[j][m]));
							else/*normalized (sensitivities are expressed in percentage)*/
							{
								if(kind_of_para == 1)
								{
									if(m%4==0)
									{
										sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_slow)/data->calc_real_mt_TE[i][m/4]));
										sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m+2] * inv_cell[mt[i].ele[j]].val_slow)/data->calc_real_mt_TM[i][m/4]));
									}
									else if((m-1)%4 == 0)
									{
										sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_slow)/data->calc_imag_mt_TE[i][(m-1)/4]));
										sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m+2] * inv_cell[mt[i].ele[j]].val_slow)/data->calc_imag_mt_TM[i][(m-1)/4]));
									}
								}
								else
								{
									if(m%4==0)
									{
										sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_res)/data->calc_real_mt_TE[i][m/4]));
										sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m+2] * inv_cell[mt[i].ele[j]].val_res)/data->calc_real_mt_TM[i][m/4]));
									}
									else if((m-1)%4 == 0)
									{
										sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_res)/data->calc_imag_mt_TE[i][(m-1)/4]));
										sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m+2] * inv_cell[mt[i].ele[j]].val_res)/data->calc_imag_mt_TM[i][(m-1)/4]));
									}
								}
							}
						}
					}
					/*********************************************/
					/*Only one mode is used*/
					else
					{
						for(m=0;m<2*mt[i].nfreq;m++)
						{
							if(sens_in_perc == 0) /*NOT normalized*/
								/*Sum all sensitvities*/
								sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (inv->rel_scal_mt * fabs(mt[i].deriv[j][m]));
							else/*normalized (sensitivities are expressed in percentage)*/
							{
								if(kind_of_para == 1)
								{
									if(m%2==0)
									{
										/*Only TE-mode*/
										if(flag->kind_of_data_mt == 1 || flag->dimension_mt == 1)
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_slow)/data->calc_real_mt_TE[i][m/2]));
										/*Only TM-mode*/
										else if(flag->kind_of_data_mt == 2)
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_slow)/data->calc_real_mt_TM[i][m/2]));
										/*Berdichewsky average*/
										else
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs(mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_slow)/(0.5*fabs(data->calc_real_mt_TE[i][m/2])+ 0.5*fabs(data->calc_real_mt_TM[i][m/2])));
									}
									else if((m-1)%2 == 0)
									{
										/*Only TE-mode*/
										if(flag->kind_of_data_mt == 1 || flag->dimension_mt == 1)
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_slow)/data->calc_imag_mt_TE[i][(m-1)/2]));
										/*Only TM-mode*/
										else if(flag->kind_of_data_mt == 2)
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_slow)/data->calc_imag_mt_TM[i][(m-1)/2]));
										/*Berdichewsky average*/
										else
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs(mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_slow)/(0.5*fabs(data->calc_imag_mt_TE[i][(m-1)/2])+ 0.5*fabs(data->calc_imag_mt_TM[i][(m-1)/2])));
											
									}
								}
								else
								{
									if(m%2==0)
									{
										/*Only TE-mode*/
										if(flag->kind_of_data_mt == 1 || flag->dimension_mt == 1)
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_res)/data->calc_real_mt_TE[i][m/2]));
										/*Only TM-mode*/
										else if(flag->kind_of_data_mt == 2)
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_res)/data->calc_real_mt_TM[i][m/2]));
										/*Berdichewsky average*/
										else
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs(mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_res)/(0.5*fabs(data->calc_real_mt_TE[i][m/2])+ 0.5*fabs(data->calc_real_mt_TM[i][m/2])));
									}
									else if((m-1)%2 == 0)
									{
										/*Only TE-mode*/
										if(flag->kind_of_data_mt == 1 || flag->dimension_mt == 1)
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_res)/data->calc_imag_mt_TE[i][(m-1)/2]));
										/*Only TM-mode*/
										else if(flag->kind_of_data_mt == 2)
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs((mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_res)/data->calc_imag_mt_TM[i][(m-1)/2]));
										/*Berdichewsky average*/
										else
											sens_mt[inv_cell[mt[i].ele[j]].ele[k]] = sens_mt[inv_cell[mt[i].ele[j]].ele[k]] + (100 * fabs(mt[i].deriv[j][m] * inv_cell[mt[i].ele[j]].val_res)/(0.5*fabs(data->calc_imag_mt_TE[i][(m-1)/2])+ 0.5*fabs(data->calc_imag_mt_TM[i][(m-1)/2])));
									}
								}
							}
						}
					}
					/*********************************************/
				}
			}
		}
	}

	/*Write the sensitivities in a binary model file*/
	WriteSensMTOut(inv->ninv, *data, geo, *grid, sens_mt);

	free(sens_mt);

	return(1);
}
