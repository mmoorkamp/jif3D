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
/*Prepare and perform the inversion*/
/*Parameter:	geo		:= Geometrie structure*/
/*				*grid	:= Grid structure*/
/*				raypath	:= Input:	Structure including the raypath of the forward cells;*/
/*						   Output:	Structure including the raypath of the inversion cells;*/
/*              fatray	:=  Input:	Structure including the fatray parameters of the forward cells*/
/*						   Output:	Structure including the fatray parameters of the inversion cells;*/
/*				  grav	:= Gravity structure for the derivatives*/
/*				   mt   := MT structure for the derivatives*/
/*				inv		:= Inversion structure */
/*			    *flag   := Flag structure*/
/*				*regu	:= Regularisation structure*/
/*         *calc_2D_mt  := 2D MT structure*/
/*  nr_of_slices_2D_mt  := Number of slices for the 2D MT*/
/*				*fname	:= Filename of the parameterfile including the spatially dependent sizes of the inversion grid cells*/
/*      *fname_vel_dens  := Filename of the parameterfile including the relationship of velocity and density*/
/*      *fname_vel_res  := Filename of the parameterfile including the relationship of velocity and resistivity*/

int InvRoutine(GEOMETRY geo, GRID_STRUCT *grid, RP_STRUCT *raypath, F_RP_STRUCT *fatray, GRAV_STRUCT *grav, MT_STRUCT *mt, INV *inv, DATA_STRUCT *data, FLAG_STRUCT *flag, REGU_STRUCT *regu, CALC_2D_MT *calc_2D_mt, long nr_of_slices_2D_mt, char *fname, char *fname_vel_dens, char *fname_vel_res)
{
	long i,j,nxyz,nxyz1,a;
	int steps,limit,kind_of_para, num_border;
	VEC *v_mod;
	BIDX_STRUCT *inv_cell;
	FLEX_STRUCT flex;
	MT_STRUCT *forward_mt;
	REL_VDR_STRUCT rel_vdr, rel_vdr2;// rel_vdr3, rel_vdr4; // rel_vdr5, rel_vdr6;

	double *Rho1, *Rho2; /*Perturbations of the densities: Rho2-Rho1 in mgal*/
	double *R1, *R2; /*Perturbations of the resistivities: R2-R1 in ohmm*/
	double tol = 1.0e-7;
	double deltaSlow, sum;
	double *v_border;

	FILE *out;
	
	limit =1000;
	/****************************************************/
	/*Determine the slowness derivatives for the gravity and MT*/
	nxyz = grid->nx * grid->ny *grid->nz;
	nxyz1 = (grid->nx +2*grid->nborder) * (grid->ny +2*grid->nborder) * (grid->nz +2*grid->nborder); 

	/******************/
	/*Gravity*/
	Rho1 = (double *)memory(NULL,1,sizeof(double),"InvRoutine");
	Rho2 = (double *)memory(NULL,1,sizeof(double),"InvRoutine");

	if(flag->index_grav !=0)
	{
		
		printf("\nFor gravity:\n");
		printf("Start determining the derivatives drho/ds\n");
		printf("----------------\n");

		Rho1 = (double *)memory((char *)Rho1,nxyz,sizeof(double),"InvRoutine");
		Rho2 = (double *)memory((char *)Rho2,nxyz,sizeof(double),"InvRoutine");

		for(i=0;i<nxyz;i++)
		{
			Rho1[i] = grid->dens[i];
			Rho2[i] = grid->dens[i];
		}

		/*For joint inversion (seismic and/or resistivity is activated):*/
		/*Derivatives of rho(s); s=slowness in s/m; are calculated*/
		if(flag->index_joint != 0)
		{

			/*Read in the file including the relationship of velocity and density*/
			ReadRelVelDensRes(fname_vel_dens, &rel_vdr);

			/*Calculate the derivatives drho/ds[g/cm^3/s/m]*/
			deltaSlow = CalcDerivDensSlow(grid, Rho2, Rho1, rel_vdr);

			free(rel_vdr.tab_v);
			free(rel_vdr.tab_dr);
		}

		/*For "simple" gravity inversion*/
		else
		{
			/*Fill the derivatives "tmp_deriv" with density distortions*/
			FillwithValues(grid, Rho2, 0.0005, 2);
			FillwithValues(grid, Rho1, -0.0005, 2);
			deltaSlow = 0.001;
		}

		printf("Start determining the derivatives dG(rho(s))/ds\n");
		printf("----------------\n\n");

		/*Calculate the slowness derivatives for the gravity measurements dG/ds (for the forward cells)*/
		/*Remark: temporary files containing the derivatives are written out*/
		DerivativesGrav(geo, *grid, data, Rho2, Rho1, deltaSlow);

	}

	free(Rho1);
	free(Rho2);


	/******************/
	/*MT*/
	R1 = (double *)memory(NULL,1,sizeof(double),"InvRoutine");
	R2 = (double *)memory(NULL,1,sizeof(double),"InvRoutine");

	if(flag->index_mt !=0)
	{
		printf("\nFor MT:\n");
		printf("Start determining the derivatives drhoR/ds\n");
		printf("----------------\n");

		R1 = (double *)memory((char *)R1,nxyz,sizeof(double),"InvRoutine");
		R2 = (double *)memory((char *)R2,nxyz,sizeof(double),"InvRoutine");

		for(i=0;i<nxyz;i++)
		{
			R1[i] = grid->res[i];
			R2[i] = grid->res[i];
		}

		/*For joint inversion (seismic and/or gravity is activated):*/
		/*Derivatives of rhoR(s); s=slowness in s/m; are calculated*/
		if(flag->index_joint != 0)
		{
			/*Read in the file including the relationship of velocity and resistivity*/
			ReadRelVelDensRes(fname_vel_res, &rel_vdr2);

			/*Calculate the derivatives drhoR/ds[ohmm/s/m]*/
			deltaSlow = CalcDerivDensSlow(grid, R2, R1, rel_vdr2);

			free(rel_vdr2.tab_v);
			free(rel_vdr2.tab_dr);
		}
		/*For "simple" MT inversion*/
		else
		{
			/*Fill the derivatives "tmp_deriv" with resistivity distortions*/
			FillwithValues(grid, R2, 0.0005, 3);
			FillwithValues(grid, R1, -0.0005, 3);
			deltaSlow = 0.001;
		}
	

		printf("Start determining the derivatives dZ(f,rhoR(s))/ds\nfor the real and the imaginary parts\n");
		printf("----------------\n\n");

		/*Calculate the slowness derivatives for the MT measurements dZ/ds (for the forward cells)*/
		if(data->ndata_mt != 0)
			forward_mt =(MT_STRUCT *)memory(NULL,data->ndata_mt,sizeof(MT_STRUCT),"InvRoutine");
		else
			forward_mt =(MT_STRUCT *)memory(NULL,1,sizeof(MT_STRUCT),"InvRoutine");

		if(flag->dimension_mt != 2 && flag->dimension_mt != 4)
			/*1D*/
			Derivatives1DMT(geo, *grid, data, R2, R1, deltaSlow, forward_mt); 
		else if(flag->dimension_mt == 2) 
			/*2D-analytical*/
			Derivavtives2DMTAnaly(calc_2D_mt, nr_of_slices_2D_mt, *grid, data, R2, R1, deltaSlow, forward_mt, flag);
		else
			/*2D-RRI*/
			; //HIER FORTSETZEN !!!

	}

	free(R1);
	free(R2);

	/*****************************************************/
	/*Set up grid for inversion*/

	printf("\n\nStart with building the inversion structures:\n");
	printf("----------------\n");

	/*Make a spatial invariant grid*/

	/*Determine the size of the inversion grid*/
	DetNumInvCells(*grid, inv);

	/*Build the spatial invariant grid structure:*/
	inv_cell = (BIDX_STRUCT *)memory(NULL,inv->nvel,sizeof(BIDX_STRUCT),"InvRoutine");
	inv->nvel_used = MakeRegularInvStruct( *grid, inv_cell, *inv, *flag);

	if(flag->index_tseis != 0)
	{

		/*Conventional rays*/
		if(flag->kind_of_rays ==1)
		{
			/*Make the ray structure for the inversion cells*/
			InvRayStruct(*grid,raypath,inv, *data);
		}
		/*Fat rays*/
		else
		{
			/*Make the fat-ray structure for the inversion cells*/
			InvFatRayStruct(*grid, fatray, inv, *data);
		}
	}

	/*Modify the inversion grid (spatial variable grid size)*/

	if(flag->flex_grid_size ==1)
	{
		/*Read in the parameter for the variable inversion grid*/
		flex = ReadFlexGridIn(fname);

		/*Rebuild the spatial variant grid structure*/
		MakeIrregularStructuresI(flex,inv_cell,*grid, inv, *flag);
		/*Rebuild the spatial variant raypath structures*/
		if(flag->index_tseis != 0)
			MakeIrregularRayStructures(raypath, fatray, inv_cell, inv->nvel , data->ndata_seis, flag->kind_of_rays);

		/*Free the memory of the flex structure*/
		free(flex.xpos);
		free(flex.ypos);
		free(flex.zpos);

		for(i=0;i<3;i++)
		{
			free(flex.xratio[i]);
			free(flex.yratio[i]);
			free(flex.zratio[i]);
		}

	}

	/*****************************************************/
	/*Determine the gravity structure for the inversion grid*/

	if(flag->index_grav !=0)
		/*Build the gravity structure*/
		InvGravStruct(geo, data, grav, inv_cell, inv->nvel_used, inv->nvel);

	/*****************************************************/
	/*Determine the MT structure for the inversion grid*/
	if(flag->index_mt != 0)
	{
		/*Build the MT inversion structure*/
		InvMTStruct(data, mt, forward_mt, inv_cell, inv->nvel, nxyz1, flag->kind_of_data_mt, flag->dimension_mt);

		for(i=0;i<data->ndata_mt;i++)
			{
				if(flag->dimension_mt != 2)
				{
					for(j=0;j<forward_mt[i].ncell;j++)
						free(forward_mt[i].deriv[j]);
				}

				free(forward_mt[i].ele);
				free(forward_mt[i].n);
				free(forward_mt[i].deriv);
			}

		free(forward_mt);
	}

	/*****************************************************/
	/*Down-scaling parallel rays in the inversion (following the idea of Hansruedi)*/
	if(flag->kind_of_rays == 1 && flag->index_tseis != 0 && (flag->balance_para_rays == 1) || (flag->balance_para_rays == 2))
		ScalParallelRays(raypath, inv_cell, inv->nvel, data->ndata_seis, inv->angle_para, data->weigth_seis, flag->balance_para_rays);

	/*****************************************************/
	/*Calculate resolution measures:*/
	/*Seismic*/
	if(flag->index_tseis != 0)
	{
		/*Calculate the ray density (hit matrix and Derivative weighted sum)*/
		if(flag->kind_of_rays ==1 && flag->ray_density !=0)
			RayDensity(*grid,raypath,inv_cell,data->ndata_seis,inv->nvel,data->weigth_seis,inv->ninv);
		else if(flag->ray_density != 0)
		/*Calculate the fatray density (hit matrix and Derivative weighted sum)*/
			FatRayDensity(*grid,fatray,inv_cell,data->ndata_seis,inv->nvel,data->weigth_seis,inv->ninv);
		/*Calculate the ray density tensor*/
		if(flag->kind_of_rays == 1 && flag->ray_density2 !=0)
		  RayDensityTensor(*grid,raypath,inv_cell,data->ndata_seis,inv->nvel,data->weigth_seis,inv->ninv);
	}


	/*Jin Chen changed new*/
	/*Calculate the border cell numbers of inversion cells.*/
	CalcNumBorder(inv,inv_cell);

	/*****************************************************/
	/*Determine the size of the matrix:*/
	SetupInvMatrix(data, grav, mt, inv, flag);

	/*Number of columns in the matrix*/
	inv->ncol = inv->nvel_used;

	/*Add an extra cell accounting for shifts in the gravity data*/
	/*Allocate memory for the sparse matrix*/
	if(flag->index_grav != 0 && regu->w_extra_cell > 0.0)
		inv->g = sp_get(inv->nrow+1,inv->ncol+1,0);
	else
		inv->g = sp_get(inv->nrow,inv->ncol,0);

	/*Assign values to the sparse matrix*/
	FillSparseMatrix(*flag,fatray,raypath,grav,mt,inv_cell,*data,inv->g);

	/*****************************************************/
	/*extra cell*/
	/*Determine the original model vector */
	if(flag->index_grav != 0 && regu->w_extra_cell > 0.0)
		inv->mx = (double *)memory(NULL,inv->ncol + 1, sizeof(double),"InvRoutine");
	else
		inv->mx = (double *)memory(NULL,inv->ncol, sizeof(double),"InvRoutine");
	
	for(i=0;i<inv->nvel;i++)
	{
		if(inv_cell[i].use == 1)
		{
			/*Case I: The model parameters are slowness values*/
			if(flag->index_tseis != 0 || (flag->index_mt != 0 && flag->index_grav != 0))
				inv->mx[inv_cell[i].used_nr] = inv_cell[i].val_slow;
			/*Case II: The model parameters are density values*/
			else if(flag->index_grav != 0)
				inv->mx[inv_cell[i].used_nr] = inv_cell[i].val_dens;
			/*Case III: The model parameters are resistivity values*/
			else if(flag->index_mt != 0)
				inv->mx[inv_cell[i].used_nr] = inv_cell[i].val_res;
			else
			{
				printf("The combination of inversion parameters is unknown\n");
				exit(0);
			}

			if(inv_cell[i].used_nr >= inv->ncol)
			{
				printf("The inversion cell number is larger than the number of columns\n");
				exit(0);
			}
		}
	}

	/*****************************************************/
	/*Calculate the residual data vector*/
	/*Add an extra cell accounting for shifts in the gravity data*/
	if(flag->index_grav != 0 && regu->w_extra_cell > 0.0)
	{
		inv->v_res = v_get(inv->nrow+1);
		memset(inv->v_res->ve,0,sizeof(double)*(inv->nrow+1));
	}
	else
	{
		inv->v_res = v_get(inv->nrow);
		memset(inv->v_res->ve,0,sizeof(double)*inv->nrow);
	}
		
	FillResdDataVector(*flag,inv->v_res,*data,raypath,fatray,grav,mt);
	
	/*****************************************************/
	/*Balance the MT data equally for different frequencies*/
// CODE IST NUR TEMPOR�R !!!
//	if(flag->index_mt != 0 && flag->index_mt_weighting != 0)
//		BalanceMTfrequencies(inv, inv->g, inv->v_res->ve);

	/*Improve the balance between the different data sets in the sparse matrix*/
	/*BJOERN_MOD5*/
	if(flag->index_joint != 0)
		BalanceSparseMatrix(inv, *flag, inv->g, inv->v_res->ve);

	/*Assign the external weights to the matrix*/
	ReweightSparseMatrix(*flag,fatray,raypath,grav,mt,inv_cell,*data,inv->g);

	/*Assign the external weights to the residual data vector*/
	ReweightResdDataVector(*flag,inv->v_res,*data,raypath,fatray,grav,mt);

	/*****************************************************/

	/*Calculate the sensitivities*/
	if(flag->write_sens_out == 1)
	{
		/*seismic*/
		if(flag->index_tseis != 0)
			DetSensSeisRays(inv, raypath, fatray, flag->kind_of_rays, inv_cell, grid, data, geo, flag->sens_in_perc);

		/*gravity*/
		if(flag->index_grav != 0)
		{
			if(flag->index_tseis != 0 || flag->index_mt != 0)
			/*Parameter is slowness*/
				kind_of_para = 1;
			else
			/*Parameter is density*/
				kind_of_para = 0;

			DetSensGrav(inv, grav, inv_cell, grid, data, geo, flag->sens_in_perc, kind_of_para);
		}

		/*mt*/
		if(flag->index_mt != 0)
		{
			if(flag->index_tseis != 0 || flag->index_grav != 0)
			/*Parameter is slowness*/
				kind_of_para = 1;
			else
			/*Parameter is resistivity*/
				kind_of_para = 0;

			DetSensMT(inv, mt, inv_cell, grid, data, geo, flag->sens_in_perc, kind_of_para, flag);
		}
	}

	/*****************************************************/
	/*Determine the regularisation (damping and smoothing)*/
	SetRegu(inv,inv_cell, inv->g, *flag, regu, inv->v_res->ve);
	/*****************************************************/

/*Jin Chen changed new*/
	/*Determine the extra smoothing*/
	if(flag->index_grav != 0 && flag->do_smooth != 0 && flag->do_extra_smooth != 0)
		SetRegu1(inv,inv_cell, inv->g, *flag, regu); /*Determine the extra smoothing equations*/
		                                                             /*and output the tmp_config1.bin*/
		                                                             	                                                             

	/*extra cell*/
	/*Allocate the memory for the model vector*/
	if(flag->index_grav != 0 && regu->w_extra_cell > 0.0)
		v_mod = v_get(inv->ncol+1);
	else
	    v_mod = v_get(inv->ncol);

	/*extra cell*/
	if(flag->index_grav != 0 && regu->w_extra_cell > 0.0)
		InvExtraCell(grid,inv,grav,inv->g,data,regu, inv->v_res->ve, *flag, inv->rel_scal_grav); /*BJOERN_MOD5*/
		
	/* compact g matrix */
	inv->g = sp_compact(inv->g,tol);

	/* Perform LSQR */
	mem_stat_mark(1);
	v_mod = iter_splsqr(inv->g,inv->v_res,tol,v_mod,limit,&steps);	//<-NOCH CHECKEN WIE DIE TOLERANZ SEIN MUSS !!!
	mem_stat_free(1);
		
	/*Copy results to inv->para_mod (parameters that is directly determined by the inversion process)*/
	inv->para_mod = (double *)memory(NULL,inv->nvel_used,sizeof(double),"InvRoutine");
	if(flag->index_tseis != 0 || (flag->index_grav != 0 && flag->index_mt != 0))
	{	/*Determine the second parameter from joint inversion (gravity OR resistivity)*/
		if(flag->index_grav != 0 || flag->index_mt != 0)
		{
			inv->para_mod2 = (double *)memory(NULL,inv->nvel_used,sizeof(double),"InvRoutine");
		}
			/*Determine the third parameter from joint inversions (resistivity)*/
		if(flag->index_grav != 0 && flag->index_mt != 0)
		{
			inv->para_mod3 = (double *)memory(NULL,inv->nvel_used,sizeof(double),"InvRoutine");
		}
	}

	/*Jin, 23.11.2006*/
	/*Copy results to inv->para_mod (parameters that is directly determined by the inversion process)*/
	CopyResultsToInvParaMod(inv, flag, v_mod, fname_vel_dens, fname_vel_res);

	/*Transfer the results from the extra cell to the inv structure*/
	if(flag->index_grav != 0 && regu->w_extra_cell > 0.0) /*BJOERN_MOD5*/
		TransExtraCellValue(inv, flag, v_mod, inv->ncol);

	sum = 0;
	num_border = 0;
	
	v_border = (double *)memory(NULL,1,sizeof(double),"InvRoutine");

	for(a=0;a<inv->nvel;a++)
	{
		if(inv_cell[a].use == 1)
		{
			if(inv_cell[a].border_inv_cell_x_right == 1)
			{
				num_border++;
				sum = sum + inv->para_mod[inv_cell[a].used_nr];
				v_border = (double *)memory((char *)v_border,num_border, sizeof(double),"InvRoutine");
				v_border[num_border-1] = inv->para_mod[inv_cell[a].used_nr];
			}
		}
	}
   
	out = fopen("v_border.txt","w");
	for(i=0;i<num_border;i++)
		fprintf(out,"%f\n",v_border[i]);
	fclose(out);

	/*****************************************************/
	/*Assign the modified slowness/density and resistivity values to the forward cells*/	
	AssignValueToForward(flag, grid, inv, inv_cell);
	/****************************************************/

	inv->ninv++;

	printf("\n----------------\n");
	printf("Iteration %d: The inversion is performed\n",inv->ninv);
	printf("----------------\n\n");

	/*Free the memory*/
	for(i=0;i<inv->nvel;i++)
	{
		free(inv_cell[i].ele);
		free(inv_cell[i].ele_inv);
	}
	free(inv_cell);

	free(inv->mx);
	SP_FREE(inv->g);
	V_FREE(inv->v_res);
	free(inv->para_mod);
	
	free(v_border);

	if(flag->index_joint != 0)
		free(inv->para_mod2);
	if(flag->index_grav != 0 && flag->index_mt != 0)
		free(inv->para_mod3);

	V_FREE(v_mod);


	return(1);
}

/*------------------------------------------------------------*/
/*Determine the size of the Jacobian matrix and determine the positions of the different entries/frechets*/
/*for the different (joint) inversion combinations*/
/*Parameter:	*data	:= Data structure*/
/*				*flag   := Structure that includes the general program flags*/
/*				grav	:= Gravity structure for the derivatives*/
/*				mt   := MT structure for the derivatives*/
/*				inv		:= Inversion structure */

int SetupInvMatrix(DATA_STRUCT *data, GRAV_STRUCT *grav, MT_STRUCT *mt, INV *inv, FLAG_STRUCT *flag)
{
	long nr_of_data_inv,i,j;

		/*Determine the size of the matrix:*/
	nr_of_data_inv = 0;

	if(flag->index_tseis != 0 && flag->index_mt == 0 && flag->index_grav == 0) /*only seismic inversion*/
		nr_of_data_inv =  data->ndata_seis;
	else if(flag->index_tseis == 0 && flag->index_mt == 0 && flag->index_grav != 0)  /*only gravity inversion*/
		nr_of_data_inv = data->ndata_grav;
	else if(flag->index_tseis == 0 && flag->index_mt != 0 && flag->index_grav == 0)  /*only MT inversion*/
	{
		/*Both modes will be used*/
		if(flag->nr_of_modes == 2)
			nr_of_data_inv = (4*data->ndata_freq_mt);
		/*Only one mode will be used*/
		else
			nr_of_data_inv = (2*data->ndata_freq_mt);
	}
	else if(flag->index_tseis_grav != 0)  /* joint inversion of seismic and gravity*/
		nr_of_data_inv = data->ndata_seis + data->ndata_grav;
	else if(flag->index_tseis_mt != 0)  /* joint inversion of seismic and MT*/
	{
		/*Both modes will be used*/
		if(flag->nr_of_modes == 2)
			nr_of_data_inv = data->ndata_seis + (4*data->ndata_freq_mt);
		/*Only one mode will be used*/
		else
			nr_of_data_inv = data->ndata_seis + (2*data->ndata_freq_mt);
	}
	else if(flag->index_grav_mt != 0)  /* joint inversion of gravity and MT*/
	{
		/*Both modes will be used*/
		if(flag->nr_of_modes == 2)
			nr_of_data_inv = data->ndata_grav + (4*data->ndata_freq_mt);
		/*Only one mode will be used*/
		else
			nr_of_data_inv = data->ndata_grav + (2*data->ndata_freq_mt);
	}
	else if(flag->index_tseis_grav_mt != 0)  /* joint inversion of seismic, gravity and MT*/
	{
		/*Both modes will be used*/
		if(flag->nr_of_modes == 2)
			nr_of_data_inv = data->ndata_seis + data->ndata_grav + (4*data->ndata_freq_mt);
		/*Only one mode will be used*/
		else
			nr_of_data_inv = data->ndata_seis + data->ndata_grav + (2*data->ndata_freq_mt);
	}
	else
	{
		printf("NO inversion is activated !!! This makes NO sense !!!\n");
		printf("Please change the parameter file !!!\n");
		exit(0);
	}

		/*Number of rows in the matrix*/

	/*Jin Chen changed new*/
	/*If the gravity data will be inverted.*/
     if(flag->index_grav != 0 && flag->do_smooth != 0 && flag->do_extra_smooth !=0)
	 { 
		if(flag->do_damp != 0) /*with damping.*/
				inv->nrow = nr_of_data_inv + 2*inv->nvel_used + inv->x_left_border_nr + inv->x_right_border_nr + inv->y_front_border_nr + inv->y_back_border_nr;
		else  /*without damping.*/
				inv->nrow = nr_of_data_inv + inv->nvel_used + inv->x_left_border_nr + inv->x_right_border_nr + inv->y_front_border_nr + inv->y_back_border_nr;
	}
	else if(flag->do_damp != 0 && flag->do_smooth != 0)	/*damping AND smoothing*/       /*           */
		inv->nrow = nr_of_data_inv + 2*inv->nvel_used;		                            /*  This is  */
	else if(flag->do_damp == 0 && flag->do_smooth != 0)	/*ONLY smoothing*/              /*    the    */ 
		inv->nrow = nr_of_data_inv + inv->nvel_used;                                    /* original  */
	else if(flag->do_damp != 0 && flag->do_smooth == 0)	/*ONLY damping*/                /*    part.  */
		inv->nrow = nr_of_data_inv + inv->nvel_used;                                    /*           */
	else if(flag->do_damp == 0 && flag->do_smooth == 0)	/*NO smoothing NOR damping*/    /*           */
		inv->nrow = nr_of_data_inv;                                                     /*           */
	
	/************************/
	inv->seis_first_row = 0;
	inv->seis_nr_row = 0;
	inv->grav_first_row = 0;
	inv->grav_nr_row = 0;
	inv->mt_first_row = 0;
	inv->mt_nr_row = 0;

	/*For simple inversions (determine the position of the different data in the matrix):*/
	/************************/
		/*seismic inversion:*/
	if(flag->index_tseis != 0 && flag->index_mt == 0 && flag->index_grav == 0)
	{
		inv->seis_first_row = 0;				/*First seismic row*/
		inv->seis_nr_row = data->ndata_seis;	/*Nr. of seismic rows*/
	}

	/************************/
		/*gravity inversion:*/
	if(flag->index_tseis == 0 && flag->index_mt == 0 && flag->index_grav != 0)
	{
		inv->grav_first_row = 0;				/*First gravity row*/
		inv->grav_nr_row = data->ndata_grav;	/*Nr. of gravity rows*/
	}

	/************************/
		/*MT inversion:*/
	if(flag->index_tseis == 0 && flag->index_mt != 0 && flag->index_grav == 0)
	{
		inv->mt_first_row = 0;					/*First MT row*/
		inv->mt_nr_row = 0;	

			/*Both modes will be used*/
		if(flag->nr_of_modes == 2)
		{
			for(i=0;i<data->ndata_mt;i++)
				inv->mt_nr_row = (4*data->nfreq_mt[i]) + inv->mt_nr_row;	/*Nr. of MT rows*/
		}
			/*Only one mode will be used*/
		else
		{
			for(i=0;i<data->ndata_mt;i++)
				inv->mt_nr_row = (2*data->nfreq_mt[i]) + inv->mt_nr_row;	/*Nr. of MT rows*/
		}

	}


	/************************/
	/*Adjust the row numbers of the gravity and MT entries for joint inversion in the inversion matrix*/
	/************************/
		/* joint inversion of seismic and gravity*/
	if(flag->index_tseis_grav != 0)
	{
		inv->seis_first_row = 0;				/*First seismic row*/
		inv->seis_nr_row = data->ndata_seis;	/*Nr. of seismic rows*/

		inv->grav_first_row = data->ndata_seis;	/*First gravity row*/
		inv->grav_nr_row = data->ndata_grav;	/*Nr. of gravity rows*/

		for(i=0;i<data->ndata_grav;i++)
			grav[i].n = grav[i].n + data->ndata_seis;
	}

	/************************/
		/* joint inversion of seismic and MT*/
	if(flag->index_tseis_mt != 0)
	{
		inv->seis_first_row = 0;				/*First seismic row*/
		inv->seis_nr_row = data->ndata_seis;	/*Nr. of seismic rows*/

		inv->mt_first_row = data->ndata_seis;	/*First MT row*/
		inv->mt_nr_row = 0;	
		
			/*Both modes will be used*/
		if(flag->nr_of_modes == 2)
		{
			for(i=0;i<data->ndata_mt;i++)
				inv->mt_nr_row = (4*data->nfreq_mt[i]) + inv->mt_nr_row;	/*Nr. of MT rows*/
		}
			/*Only one mode will be used*/
		else
		{
			for(i=0;i<data->ndata_mt;i++)
				inv->mt_nr_row = (2*data->nfreq_mt[i]) + inv->mt_nr_row;	/*Nr. of MT rows*/
		}

		for(i=0;i<data->ndata_mt;i++)
		{
			/*Both modes will be used*/
			if(flag->nr_of_modes == 2)
			{
				for(j=0;j<(4*mt[i].nfreq);j++)
					mt[i].n[j] = mt[i].n[j] + data->ndata_seis;
			}
			/*Only one mode will be used*/
			else
			{
				for(j=0;j<(2*mt[i].nfreq);j++)
					mt[i].n[j] = mt[i].n[j] + data->ndata_seis;
			}
		}
	}

	/************************/
		/* joint inversion of gravity and MT*/
	if(flag->index_grav_mt != 0)
	{

		inv->grav_first_row = 0;				/*First gravity row*/
		inv->grav_nr_row = data->ndata_grav;	/*Nr. of gravity rows*/

		inv->mt_first_row = data->ndata_grav;	/*First MT row*/
		inv->mt_nr_row = 0;	

		/*Both modes will be used*/
		if(flag->nr_of_modes == 2)
		{
			for(i=0;i<data->ndata_mt;i++)
				inv->mt_nr_row = (4*data->nfreq_mt[i]) + inv->mt_nr_row;	/*Nr. of MT rows*/
		}
		/*Only one mode will be used*/
		else
		{
			for(i=0;i<data->ndata_mt;i++)
				inv->mt_nr_row = (2*data->nfreq_mt[i]) + inv->mt_nr_row;	/*Nr. of MT rows*/
		}

		for(i=0;i<data->ndata_mt;i++)
		{
			/*Both modes will be used*/
			if(flag->nr_of_modes == 2)
			{
				for(j=0;j<(4*mt[i].nfreq);j++)
					mt[i].n[j] = mt[i].n[j] + data->ndata_grav;
			}
			/*Only one mode will be used*/
			else
			{
				for(j=0;j<(2*mt[i].nfreq);j++)
					mt[i].n[j] = mt[i].n[j] + data->ndata_grav;
			}
		}
	}


		/* joint inversion of seismic, gravity and MT*/
	if(flag->index_tseis_grav_mt != 0)
	{

		inv->seis_first_row = 0;				/*First seismic row*/
		inv->seis_nr_row = data->ndata_seis;	/*Nr. of seismic rows*/

		inv->grav_first_row = data->ndata_seis;	/*First gravity row*/
		inv->grav_nr_row = data->ndata_grav;	/*Nr. of gravity rows*/

		inv->mt_first_row = data->ndata_seis + data->ndata_grav;	/*First MT row*/
		inv->mt_nr_row = 0;

			/*Both modes will be used*/
		if(flag->nr_of_modes == 2)
		{
			for(i=0;i<data->ndata_mt;i++)
				inv->mt_nr_row = (4*data->nfreq_mt[i]) + inv->mt_nr_row;	/*Nr. of MT rows*/
		}
			/*Only one mode will be used*/
		else
		{
			for(i=0;i<data->ndata_mt;i++)
				inv->mt_nr_row = (2*data->nfreq_mt[i]) + inv->mt_nr_row;	/*Nr. of MT rows*/
		}

		for(i=0;i<data->ndata_grav;i++)
			grav[i].n = grav[i].n + data->ndata_seis;

		for(i=0;i<data->ndata_mt;i++)
		{
				/*Both modes will be used*/
			if(flag->nr_of_modes == 2)
			{
				for(j=0;j<(4*mt[i].nfreq);j++)
					mt[i].n[j] = mt[i].n[j] + data->ndata_seis + data->ndata_grav;
			}
				/*Only one mode will be used*/
			else
			{	
				for(j=0;j<(2*mt[i].nfreq);j++)
					mt[i].n[j] = mt[i].n[j] + data->ndata_seis + data->ndata_grav;
			}
		}
	}


	if(nr_of_data_inv == 0)
	{
		printf("NO data exists for the inversion !!!\n");
		printf("Change the data input files!!!\n");
		exit(0);
	}

	/*Specify the number of data rows*/
	inv->ndata = nr_of_data_inv;

	return(1);
}

/*------------------------------------------------------------*/
/*Fill the sparse matrix */
/*Parameter:    *fatray			:= Fat-ray structure*/
/*				*raypath		:= Ray structure*/
/*				*grav			:= Gravity structure*/
/*				*mt				:= MT structure*/
/*				*inv_cell		:= Structure of the inversion cells*/
/*				 nray			:= Number of shot-receiver combinations*/
/*				 *A				:= Sparse matrix*/
/*				flag			:= General flags*/

int FillSparseMatrix(FLAG_STRUCT flag, F_RP_STRUCT *fatray,RP_STRUCT *raypath, GRAV_STRUCT *grav, MT_STRUCT *mt, BIDX_STRUCT *inv_cell, DATA_STRUCT data, SPMAT *A)
{
	long i,k,j;

	/*Seismic*/
	if(flag.index_tseis != 0)
	{
		if(flag.kind_of_rays == 1)
		{
			/*for conventional rays*/
			/*Loop over all rays*/
			for(i=0;i<data.ndata_seis;i++)
				/*Loop over all ray segments*/
				for(j=0;j<raypath[i].nray;j++)
				{
					/*Fill the sparse matrix*/
					sp_set_val(A,raypath[i].n,inv_cell[raypath[i].ele[j]].used_nr,raypath[i].len[j]);

				}
		}
		else
		{
			/*for fat-rays*/
			/*Loop over all fatrays*/
			for(i=0;i<data.ndata_seis;i++)
				/*Loop over all fat-ray segments*/
				for(j=0;j<fatray[i].ncell;j++)
				{
					/*Fill the sparse matrix*/
					sp_set_val(A,fatray[i].n,inv_cell[fatray[i].ele[j]].used_nr,fatray[i].weight[j]);
				}	
		}
	}


	/*Gravity*/
	if(flag.index_grav != 0)
	{
		/*Loop over all stations*/
			for(i=0;i<data.ndata_grav;i++)
			{
				/*Loop over all inversion cells*/
				for(j=0;j<grav[i].ncell;j++)
				{
					/*Fill the sparse matrix*/
					sp_set_val(A,grav[i].n,inv_cell[grav[i].ele[j]].used_nr,grav[i].deriv[j]);
				}
			}
	}

	/*MT*/
	if(flag.index_mt != 0)
	{
		/*Loop over all stations*/
		for(i=0;i<data.ndata_mt;i++)
		{
			/*Loop over all cells that affects the MT measurements*/
			for(j=0;j<mt[i].ncell;j++)
			{
				/*BOTH modes are used (2-D modelling)*/
				if(flag.nr_of_modes == 2  && flag.dimension_mt == 2 )
				{
					/*Loop over all frequencies*/
					for(k=0;k<(4*mt[i].nfreq);k++)
					{
						/*Fill the sparse matrix*/
						sp_set_val(A,mt[i].n[k],inv_cell[mt[i].ele[j]].used_nr,mt[i].deriv[j][k]);
					}
				}
				/*BOTH modes are used (2-D forward modelling + 1-D inversion)*/
				else if(flag.nr_of_modes == 2 && flag.dimension_mt == 3)
				{
					/*Loop over all frequencies*/
					for(k=0;k<mt[i].nfreq;k++)
					{
						/*Fill the sparse matrix*/
							/*TE*/
						sp_set_val(A,mt[i].n[4*k],inv_cell[mt[i].ele[j]].used_nr,mt[i].deriv[j][2*k]); /*real*/
						sp_set_val(A,mt[i].n[4*k+1],inv_cell[mt[i].ele[j]].used_nr,mt[i].deriv[j][2*k+1]); /*imag*/
							/*TM*/
						sp_set_val(A,mt[i].n[4*k+2],inv_cell[mt[i].ele[j]].used_nr,(-1)*mt[i].deriv[j][2*k]); /*real*/
						sp_set_val(A,mt[i].n[4*k+3],inv_cell[mt[i].ele[j]].used_nr,(-1)*mt[i].deriv[j][2*k+1]); /*imag*/
					}
					
				}
				/*ONLY one mode is used*/
				/*(Special case of 2-D forward + 1-D inversion of the TM-mode; the sign has to be switched)*/
				else if(flag.dimension_mt == 3 && flag.kind_of_data_mt == 2)
				{
					/*Loop over all frequencies*/
					for(k=0;k<(2*mt[i].nfreq);k++)
					{
						/*Fill the sparse matrix*/
						sp_set_val(A,mt[i].n[k],inv_cell[mt[i].ele[j]].used_nr,(-1)*mt[i].deriv[j][k]);
					}
				}
				/*ONLY one mode is used*/
				else
				{
					/*Loop over all frequencies*/
					for(k=0;k<(2*mt[i].nfreq);k++)
					{
						/*Fill the sparse matrix*/
						sp_set_val(A,mt[i].n[k],inv_cell[mt[i].ele[j]].used_nr,mt[i].deriv[j][k]);
					}
				}
			}
		}
	}


	return(1);
}

/*------------------------------------------------------------*/
/*Fill the residual data vector d_obs -d_calc */
/*Remark: For the MT measurements d_obs can vary dependent on the data used (TE-Mode,TM-Mode, Berdichewsky average) */
/*Parameter:    *fatray			:= Fat-ray structure*/
/*				*raypath		:= Ray structure*/
/*				*grav			:= Gravity structure*/
/*				*mt				:= MT structure*/
/*				data			:= Data structure*/
/*				flag			:= General flags*/
/*				vec				:= Structure for the residual vector*/

int FillResdDataVector(FLAG_STRUCT flag,VEC *vec, DATA_STRUCT data, RP_STRUCT *raypath, F_RP_STRUCT *fatray, GRAV_STRUCT *grav, MT_STRUCT *mt)
{
	long i,j;
	double berd_real, berd_imag;

	/*Seismic: traveltimes in s*/
	if(flag.index_tseis != 0)
	{
		/*Conventional rays*/
		if(flag.kind_of_rays == 1)
		{
			for(i=0;i<data.ndata_seis;i++) 
				vec->ve[raypath[i].n] = (data.tobs[i] - data.tcalc[i])/1000;
		}
		else
		{
			for(i=0;i<data.ndata_seis;i++) 
				vec->ve[fatray[i].n] = (data.tobs[i] - data.tcalc[i])/1000;
		}
	}

	/*Gravity: acceleration in mgal*/
	if(flag.index_grav != 0)
	{
		for(i=0;i<data.ndata_grav;i++)
			vec->ve[grav[i].n] = (data.obs_grav[i] - data.calc_grav[i]);
	}

	/*MT:Real- and imaginary part of the impedance in ohm*/
	if(flag.index_mt != 0)
	{
		for(i=0;i<data.ndata_mt;i++)
		{
			/*Loop over the frequencies*/
			for(j=0;j<(data.nfreq_mt[i]);j++)
			{

				/*Using different kind of observed data*/

				/****************/
				/*		1D		*/
				/****************/
				if(flag.dimension_mt == 1 )
				{
						/*Using TE-mode data*/
					if(flag.kind_of_data_mt == 1)
					{
						/*Real part*/
						vec->ve[mt[i].n[2*j]] = (data.real_mt_TE[i][j] - data.calc_real_mt_TE[i][j]);
						/*Imaginary part*/
						vec->ve[mt[i].n[2*j+1]] = (data.imag_mt_TE[i][j] - data.calc_imag_mt_TE[i][j]);
					}
						/*Using TM-mode data*/
					else if(flag.kind_of_data_mt == 2)
					{
						/*Real part*/
						vec->ve[mt[i].n[2*j]] = ((-1)*data.real_mt_TM[i][j] - data.calc_real_mt_TE[i][j]);
						/*Imaginary part*/
						vec->ve[mt[i].n[2*j+1]] = ((-1)*data.imag_mt_TM[i][j] - data.calc_imag_mt_TE[i][j]);
					}
						/*Using Berdichewsky average*/
					else
					{
						berd_real =  (data.real_mt_TE[i][j] - data.real_mt_TM[i][j])/2;
						berd_imag =  (data.imag_mt_TE[i][j] - data.imag_mt_TM[i][j])/2;
	
						/*Real part*/
						vec->ve[mt[i].n[2*j]] =   (berd_real - data.calc_real_mt_TE[i][j]);
						/*Imaginary part*/
						vec->ve[mt[i].n[2*j+1]] = (berd_imag - data.calc_imag_mt_TE[i][j]);

					}
				}
				/****************/
				/*		2D (or 2D forward +1D inversion)*/
				/****************/
				else
				{
					/*Using TE-mode data*/
					if(flag.kind_of_data_mt == 1)
					{
						/*Real part*/
						vec->ve[mt[i].n[2*j]] = (data.real_mt_TE[i][j] - data.calc_real_mt_TE[i][j]);
						/*Imaginary part*/
						vec->ve[mt[i].n[2*j+1]] = (data.imag_mt_TE[i][j] - data.calc_imag_mt_TE[i][j]);
					}
					/*Using TM-mode data*/
					else if(flag.kind_of_data_mt == 2)
					{
						/*Real part*/
						vec->ve[mt[i].n[2*j]] = (data.real_mt_TM[i][j] - data.calc_real_mt_TM[i][j]);
						/*Imaginary part*/
						vec->ve[mt[i].n[2*j+1]] = (data.imag_mt_TM[i][j] - data.calc_imag_mt_TM[i][j]);
					}
					/*Using both TE and TM mode*/
					else
					{
						/*********/
						/*TE-mode*/
						/*********/
						/*Real part*/
						vec->ve[mt[i].n[4*j]] = (data.real_mt_TE[i][j] - data.calc_real_mt_TE[i][j]);
						/*Imaginary part*/
						vec->ve[mt[i].n[4*j+1]] = (data.imag_mt_TE[i][j] - data.calc_imag_mt_TE[i][j]);

						/*********/
						/*TM-mode*/
						/*********/
						/*Real part*/
						vec->ve[mt[i].n[4*j+2]] = (data.real_mt_TM[i][j] - data.calc_real_mt_TM[i][j]);
						/*Imaginary part*/
						vec->ve[mt[i].n[4*j+3]] = (data.imag_mt_TM[i][j] - data.calc_imag_mt_TM[i][j]);
					}
				}

			}
		}
	}

	return(1);
}

/*------------------------------------------------------------*/
/*Determine the regularisation terms (damping and smoothing) for the inversion*/
/*Parameter: *inv		:= inversion structure*/
/*			*inv_cell	:= Cell structures of the inversion cells*/
/*				*A		:= Sparse matrix*/
/*			  flags		:= Structure with flags*/
/*			 *regu		:= Regularization structure*/
/*			*res		:= Residual vector*/

int SetRegu(INV *inv, BIDX_STRUCT *inv_cell, SPMAT *A, FLAG_STRUCT flag, REGU_STRUCT *regu, double *res)
{
	long a,b,i,j,counter,num_neigh,*idx_neigh;
	long sti,pos,*pos_neigh,nr_inv_cells, c;
	double sum_damp, sum2, *idx; /*BJOERN_MOD5*/
	double *xa1,*xa2,*ya1,*ya2,*za1,*za2;
	double idxsum,*smooth,*smooth_neigh;

	char s1[32];
	FILE *ouf,*inf;

	printf("Start with building the regularization:\n");
	printf("----------------\n");

	/*Determine the slowness regularization factor*/
	if(flag.do_damp != 0 || flag.do_smooth != 0)
		regu->damp2 = RegFactor(A,0,inv->ndata,0,inv->nvel_used,1,&c);
	printf("Regularization factor is %lf\n",regu->damp2);

	/*******************************************************/
	/*NO smoothing, BUT damping*/
	/* Step length minimization ("creeping")*/
	if(flag.do_damp != 0 && flag.do_smooth == 0)
	{
		/*Jin Chen changed new*/
		sum_damp = regu->damp * regu->damp2; /*BJOERN_MOD5*/
		printf("Effective damping is %lf\n",sum_damp);

		/*Determine the damping constrain:*/
		/*Loop over all inversion cells*/
		for(a=0;a<inv->nvel;a++)
		{
			if(inv_cell[a].use == 1)
			{

				if(inv_cell[a].used_nr >= inv->ncol)
				{
					printf("The inversion cell number is larger than the number of columns\n");
					exit(0);
				}

				/*For the pure Marquardt inversion the creeping scheme will be used, therefore res = 0*/
				sp_set_val(A,inv_cell[a].used_nr+inv->ndata,inv_cell[a].used_nr,sum_damp*inv_cell[a].bdamp); /*BJOERN_MOD5*/
			}
		}
		printf("The damping term is introduced into the matrix!\n\n");
	}


	/******************************************************/
	/*WITH smoothing*/
	if(flag.do_smooth != 0)
	{
		/**************************************/
		/* Calculate modified residual vector */
			/*Loop over all rays*/
		for(a=0;a<inv->ndata;a++)
		{
				/*Loop over all inversion cells*/
			for(b=0;b<inv->ncol;b++)
			{
				res[a] = res[a] + (inv->mx[b] * sp_get_val(A,a,b));
			}
		}

		/***************************************/
		if(flag.do_damp != 0) /*and WITH damping*/
		{
			/*Jin Chen changed new*/
			sum_damp = regu->damp * regu->damp2;
			printf("Effective damping is %lf\n",sum_damp); /*BJOERN_MOD5*/

			/*Determine the damping constrain:*/
			/*Loop over all inversion cells*/
			for(a=0;a<inv->nvel;a++)
			{
				if(inv_cell[a].use == 1)
				{

					if(inv_cell[a].used_nr >= inv->ncol)
					{
						printf("The inversion cell number is larger than the number of columns\n");
						exit(0);
					}

					/*Adjust the damping constraints in the res. vector*/
					res[inv_cell[a].used_nr + inv->ndata] = inv->mx[inv_cell[a].used_nr]*sum_damp*inv_cell[a].bdamp; /*BJOERN_MOD5*/
					/*Fill the matrix with the damping constrains*/
					sp_set_val(A,inv_cell[a].used_nr + inv->ndata,inv_cell[a].used_nr,sum_damp*inv_cell[a].bdamp); /*BJOERN_MOD5*/
				}
			}
			printf("The damping term is introduced into the matrix!\n\n");
		}
		/***************************************/

		/*Setup the smoothing equation*/

		if(flag.do_damp != 0 && flag.do_smooth != 0)
			sti = inv->ndata + inv->nvel_used; /*First row of the smoothing constrain in the matrix*/
		else
			sti = inv->ndata;

		sum2 = regu->damp2 * regu->smooth;
		if(sum2 <= 0)
		{
			printf("The smoothing factor have to be larger than 0.0\n");
			exit(0);
		}

		/***************************************/
		sprintf(s1,"tmp_config.bin");

		make_config:; /*Because the former "tmp_config.bin" file does not fit, it has to be calculated once again*/

		if(access(s1,0) == -1) /*Check, if the temporary file "tmp_config.bin" exists*/ 
		{
		    ouf = fopen(s1,"wb"); /*Generate the temporary file "tmp_config.bin" with the spatial
			                        relationships of the cells*/

			fwrite(&(inv->nvel_used),sizeof(long),1,ouf); /*Write out the number of active inversion cells*/

			idx = (double *)memory(NULL,inv->nvel_used, sizeof(double),"SetRegu"); 

			counter = 0;

			/**************************************/
			/*Determine the boundaries of the inversion cells*/

			/*Memory allocation for the boundaries of the inversion cells*/
			xa1 = (double *)memory(NULL,inv->nvel, sizeof(double),"SetRegu");
			xa2 = (double *)memory(NULL,inv->nvel, sizeof(double),"SetRegu");
			ya1 = (double *)memory(NULL,inv->nvel, sizeof(double),"SetRegu");
			ya2 = (double *)memory(NULL,inv->nvel, sizeof(double),"SetRegu");
			za1 = (double *)memory(NULL,inv->nvel, sizeof(double),"SetRegu");
			za2 = (double *)memory(NULL,inv->nvel, sizeof(double),"SetRegu");


			/*Loop over all inversion cells*/
			for(a=0;a<inv->nvel;a++)
			{
				if(inv_cell[a].use == 1)	/*Only ACTIVE inversion cells will be considered*/
				{
					/*Determine the boundaries of the inversion cells*/
					xa1[a] = inv_cell[a].xo - 0.5*inv_cell[a].xdim;
					xa2[a] = inv_cell[a].xo + 0.5*inv_cell[a].xdim;
					ya1[a] = inv_cell[a].yo - 0.5*inv_cell[a].ydim;
					ya2[a] = inv_cell[a].yo + 0.5*inv_cell[a].ydim;
					za1[a] = inv_cell[a].zo - 0.5*inv_cell[a].zdim;
					za2[a] = inv_cell[a].zo + 0.5*inv_cell[a].zdim;
				}
			}


			/*Loop over all inversion cells*/
			for(a=0;a<inv->nvel;a++)
			{
				if(inv_cell[a].use == 1)	/*Only ACTIVE inversion cells will be considered*/
				{

					for(b=0;b<inv->nvel_used;b++)
						idx[b] = 0.0;

					counter++;
			

					for(b=0;b<inv->nvel;b++)
					{
						if(inv_cell[b].use == 1)
						{
							if( ((0.5*(inv_cell[a].xdim + inv_cell[b].xdim)) < fabs(inv_cell[a].xo - inv_cell[b].xo)) ||
								((0.5*(inv_cell[a].ydim + inv_cell[b].ydim)) < fabs(inv_cell[a].yo - inv_cell[b].yo)) ||
								((0.5*(inv_cell[a].zdim + inv_cell[b].zdim)) < fabs(inv_cell[a].zo - inv_cell[b].zo)) )
								;
							else
							{
							/*Find horizontal neighbors in x-direction*/
								if((ya2[b] > ya1[a]) && (ya1[b] < ya2[a]) && (za2[b] > za1[a]) && (za1[b] < za2[a]) && (a!=b))
									idx[inv_cell[b].used_nr] = sum2;
								
							/*Find horizontal neighbors in y-direction*/
								if((xa2[b] > xa1[a]) && (xa1[b] < xa2[a]) && (za2[b] > za1[a]) && (za1[b] < za2[a]) && (a!=b))
									idx[inv_cell[b].used_nr] = sum2 * regu->smoothratio_y;
							

							/*Find vertical neighbors (in z-direction)*/
								if((xa2[b] > xa1[a]) && (xa1[b] < xa2[a]) && (ya2[b] > ya1[a]) && (ya1[b] < ya2[a]) && (a!=b))
									idx[inv_cell[b].used_nr] = sum2 * regu->smoothratio_z;
						
							}

						}

					}

				
					
					idxsum = 0.0;
					num_neigh = 0;

					idx_neigh = (long *)memory(NULL,1, sizeof(long),"SetRegu"); 
					smooth = (double *)memory(NULL,1, sizeof(double),"SetRegu"); 


					/*Sum of the values of all neighboring cells*/
					for(b=0;b<inv->nvel_used;b++)
					{
						idxsum = idxsum + idx[b];

						if(idx[b] != 0)
						{
							num_neigh++;

							idx_neigh = (long *)memory((char *)idx_neigh,num_neigh, sizeof(long),"SetRegu");
							idx_neigh[num_neigh-1] = b;

							smooth = (double *)memory((char *)smooth,num_neigh, sizeof(double),"SetRegu");
							smooth[num_neigh-1] = idx[b]/sum2;

						}
					}

					/*Write neighboring relationships to the "tmp_config.bin" file*/

					fwrite(&(inv_cell[a].used_nr),sizeof(long),1,ouf);	/*Write out the index of the active inversion cell*/
					fwrite(&num_neigh,sizeof(long),1,ouf);				/*Write out the number of neighboring cells*/

					for(i=0;i<num_neigh;i++)
					{
						fwrite(&(idx_neigh[i]),sizeof(long),1,ouf);		/*Write out the indeces of the active neighboring inversion cells*/
						fwrite(&(smooth[i]),sizeof(double),1,ouf);		/*Write out the corresponding smooth ratios*/
					}


	
					free(idx_neigh);
					free(smooth);

					/*Find the cells without any neighboring cells*/
					if(idxsum == 0.0)
					{
						printf("Warning! Element %d (x= %lf, y= %lf, z= %lf) has NO neighbors\n", inv_cell[a].used_nr, inv_cell[a].xo, inv_cell[a].yo, inv_cell[a].zo);
						printf("Therefore this element will be NOT smoothed!!");
					}
					else
					{
						/*Fill in the matrix the smoothing constrain*/
						 sp_set_val(A, inv_cell[a].used_nr+sti, inv_cell[a].used_nr, sum2 + sp_get_val(A, inv_cell[a].used_nr+sti, inv_cell[a].used_nr));

						 for(b=0;b<inv->nvel;b++)
						 {
							 if(inv_cell[b].use == 1 && idx[inv_cell[b].used_nr] != 0.0)
								sp_set_val(A, inv_cell[a].used_nr+sti,inv_cell[b].used_nr, -idx[inv_cell[b].used_nr]/idxsum * sum2);
						 }
					}

					if(counter%1000 == 1)
					printf("For %d of %d inversion cells the smoothing is inserted in the matrix\n",counter,inv->nvel_used);
				}

			}


			printf("The smoothing term is introduced into the matrix!\n\n");

			free(idx);

			free(xa1);
			free(xa2);
			free(ya1);
			free(ya2);
			free(za1);
			free(za2);

			fclose(ouf);

		}
		/**************************************/
	    else /*The smoothing relations are already calculated in a former iteration, therefore the relations 
			  are taken from the temporary file "tmp_config.bin"*/
		{

			inf = fopen(s1,"rb");

	        if (inf == NULL)
			{
				fprintf(stderr,"Unable to open %s\n",s1);
				exit(0);
			}

			/*Read in the binary file*/
			fread(&nr_inv_cells,sizeof(long),1,inf);								/*Read in the number of active inversion cells*/

			if((inv->nvel_used) != nr_inv_cells)
			{
				printf("The tmp_config.bin file does not fit to the inversion grid!\n");
				printf("Number of active inversion cells in the grid: %d\n",inv->nvel_used);
				printf("Number of active inversion cells stored in the file: %d\n\n", nr_inv_cells);
				printf("=> the smoothing relations will be re-calculated\n\n");

				/*Close and remove the "tmp_config.bin" file*/
				fclose(inf);
				remove(s1);

				/*Make a new "tmp_config.bin" file*/
				goto make_config;
			}

			for(i=0;i<inv->nvel_used;i++)
			{
					fread(&pos,sizeof(long),1,inf);					/*Read in the index of the active inversion cell*/
					fread(&num_neigh,sizeof(long),1,inf);				/*Read in the number of neighboring cells*/

					idxsum = 0.0;
					
					pos_neigh = (long *)memory(NULL,1, sizeof(long),"SetRegu"); 
					smooth_neigh = (double *)memory(NULL,1, sizeof(double),"SetRegu"); 

					for(j=0;j<num_neigh;j++)
					{
						pos_neigh = (long *)memory((char *)pos_neigh,j+1, sizeof(long),"SetRegu"); 
						smooth_neigh = (double *)memory((char *)smooth_neigh,j+1, sizeof(double),"SetRegu"); 

						fread(&pos_neigh[j],sizeof(long),1,inf);					/*Read in the indeces of the active neighboring inversion cells*/
						fread(&smooth_neigh[j],sizeof(double),1,inf);				/*Read in the corresponding smooth ratios*/

						/*Calculate the smoothing ratio*/
						smooth_neigh[j] = sum2*smooth_neigh[j];

						idxsum = idxsum + smooth_neigh[j];
						
					}

					/*Fill in the smoothing constraints into the sparse matrix*/
					sp_set_val(A, pos + sti, pos, sum2 + sp_get_val(A, pos + sti, pos));

					for(j=0;j<num_neigh;j++)
					{
						sp_set_val(A, pos + sti, pos_neigh[j], - smooth_neigh[j]/idxsum * sum2);
					}

					free(pos_neigh);
					free(smooth_neigh);

			}

			printf("The smoothing is introduced into the matrix! \n\n");

			fclose(inf);

		}

	}
	return(1);
}


/*Jin Chen changed new*/
/*------------------------------------------------------------*/
/*Determine the extra smoothing for the inversion*/
/*REMARK: regu->damp2 MUST BE CALCULATED FROM SETREGU.!!!*/
/*Parameter: *inv		:= inversion structure*/
/*			*inv_cell	:= Cell structures of the inversion cells*/
/*				*A		:= Sparse matrix*/
/*			  flags		:= Structure with flags*/
/*			 *regu		:= Regularization structure*/
/*			*res		:= Residual vector*/

int SetRegu1(INV *inv, BIDX_STRUCT *inv_cell, SPMAT *A, FLAG_STRUCT flag, REGU_STRUCT *regu)
{
	long a,b,i,j,counter,*idx_x_left_border,*idx_x_right_border,*idx_y_front_border,*idx_y_back_border;
	long pos,sti,nr_new_row,nr_new_row1,num_x_left,num_x_right,num_y_front,num_y_back,counter_row;
	long *idx_x_left_layer_cell,*idx_x_right_layer_cell,*idx_y_front_layer_cell,*idx_y_back_layer_cell;
	/*Number of the cells in the same layer*/
	long num_x_left_border,num_x_right_border,num_y_front_border,num_y_back_border, nr_inv_cells;
	double *dist_x_left,*dist_x_right,*dist_y_front,*dist_y_back;
	int kind_of_border;                                 /*1--x_left_border, 2--x_right_border, 3--y_front_border, 4--y_back_border*/
    double sum_dist,sum3;                    /*The sum distance of the same layer cells and the extra smooth entry.*/


	char s1[32];
	FILE *ouf,*inf;

	printf("Start with building the extra smoothing:\n");
	printf("----------------\n");

	sum3 = regu->damp2 * regu->val_extra_smooth;

	/*Setup the smoothing equation*/

	if(flag.do_damp != 0)
		sti = inv->ndata + 2*inv->nvel_used; /*First row of the extra smoothing constrain in the matrix*/
	else
		sti = inv->ndata + inv->nvel_used;

	/***************************************/
	sprintf(s1,"tmp_config1.bin");
	make_config1:; /*Because the former "tmp_config.bin" file does not fit, it has to be calculated once again*/

	if(access(s1,0) == -1) /*Check, if the temporary file "tmp_config.bin" exists*/ 
	{
	    ouf = fopen(s1,"wb"); /*Generate the temporary file "tmp_config.bin" with the spatial
		                        relationships of the cells*/

		fwrite(&(inv->nvel_used),sizeof(long),1,ouf); /*Write out the number of active inversion cells*/

		nr_new_row = inv->x_left_border_nr + inv->x_right_border_nr + inv->y_front_border_nr + inv->y_back_border_nr;
		fwrite(&(nr_new_row),sizeof(long),1,ouf); /*Write out the number of active inversion cells*/

		idx_x_left_border = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
		idx_x_right_border = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
		idx_y_front_border = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
		idx_y_back_border = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 

		counter = 0;
		kind_of_border = 0;
		counter_row = 0;


		/*Find the x- and y- border active cells and the active cells in the same layer.*/
		num_x_left_border = 0;
		num_x_right_border = 0;
		num_y_front_border = 0;
		num_y_back_border = 0;

		/*Loop over all inversion cells*/
		for(a=0;a<inv->nvel;a++)
		{
			if(inv_cell[a].use == 1)	/*Only ACTIVE inversion cells will be considered*/
			{
				counter++;

				if(inv_cell[a].border_inv_cell_x_left == 1 && inv_cell[a].border_inv_cell_x_right != 1)
				{
					counter_row++;
					num_x_left_border++;
					idx_x_left_border = (long *)memory((char *)idx_x_left_border,num_x_left_border, sizeof(long),"SetRegu1");
					idx_x_left_border[num_x_left_border-1]=inv_cell[a].used_nr;
					kind_of_border = 1;

					idx_x_left_layer_cell = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
					dist_x_left = (double *)memory(NULL,1, sizeof(double),"SetRegu1");
					num_x_left = 0; /*The number of the inv_cells in the same layer with this left border cell*/

					for(b=0;b<inv->nvel;b++)
					{
						if(inv_cell[b].use == 1)
						{
							/*Find horizontal layer cells in x-direction*/
							if( ((0.5*(inv_cell[a].ydim + inv_cell[b].ydim)) > fabs(inv_cell[a].yo - inv_cell[b].yo)) &&
								((0.5*(inv_cell[a].zdim + inv_cell[b].zdim)) > fabs(inv_cell[a].zo - inv_cell[b].zo)) &&
								(a != b) && inv_cell[b].border_inv_cell_x_right != 1
								&& fabs(inv_cell[a].xo - inv_cell[b].xo) <= regu->dist_extra_smooth)
							{
                                num_x_left++;
								idx_x_left_layer_cell = (long *)memory((char *)idx_x_left_layer_cell,num_x_left, sizeof(long),"SetRegu1");
								dist_x_left = (double *)memory((char *)dist_x_left,num_x_left, sizeof(double),"SetRegu1");
									
								idx_x_left_layer_cell[num_x_left-1] = inv_cell[b].used_nr;
								dist_x_left[num_x_left-1] = fabs(inv_cell[a].xo - inv_cell[b].xo);
							}		
						}
					}

					/*Write layer cell relationships to the "tmp_config.bin" file*/

					fwrite(&(inv_cell[a].used_nr),sizeof(long),1,ouf);	/*Write out the index of the active left border inversion cell*/
					fwrite(&kind_of_border, sizeof(int),1,ouf);         /*Write out the kind of border cell.*/
					fwrite(&num_x_left,sizeof(long),1,ouf);				/*Write out the number of in same layer cells*/
	
					sum_dist = 0.0;
					for(i=0;i<num_x_left;i++)
					{	/*Calculate the sum value of the same layer cells*/
						sum_dist = sum_dist + dist_x_left[i];

						fwrite(&(idx_x_left_layer_cell[i]),sizeof(long),1,ouf);		/*Write out the indeces of the active inversion cells in same layer*/
						fwrite(&(dist_x_left[i]),sizeof(double),1,ouf);		/*Write out the corresponding distance from the border cell*/
					}

					sp_set_val(A, counter_row-1+sti, inv_cell[a].used_nr, sum3 + sp_get_val(A, counter_row-1+sti, inv_cell[a].used_nr));

					for(i=0;i<num_x_left;i++)
						sp_set_val(A, counter_row-1+sti,idx_x_left_layer_cell[i],-(2.0/num_x_left*sum_dist-dist_x_left[i])/sum_dist*sum3);

					free(idx_x_left_layer_cell);
					free(dist_x_left);

				}

				if(inv_cell[a].border_inv_cell_x_right == 1 && inv_cell[a].border_inv_cell_x_left != 1)
				{
					counter_row++;
					num_x_right_border++;
					idx_x_right_border = (long *)memory((char *)idx_x_right_border,num_x_right_border, sizeof(long),"SetRegu1");
					idx_x_right_border[num_x_right_border-1]=inv_cell[a].used_nr;
					kind_of_border = 2;

					idx_x_right_layer_cell = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
					dist_x_right = (double *)memory(NULL,1, sizeof(double),"SetRegu1");
					num_x_right = 0;

					for(b=0;b<inv->nvel;b++)
					{
						if(inv_cell[b].use == 1)
						{
							if( ((0.5*(inv_cell[a].ydim + inv_cell[b].ydim)) > fabs(inv_cell[a].yo - inv_cell[b].yo)) &&
								((0.5*(inv_cell[a].zdim + inv_cell[b].zdim)) > fabs(inv_cell[a].zo - inv_cell[b].zo)) &&
								(a != b) && inv_cell[b].border_inv_cell_x_left != 1
								&& fabs(inv_cell[a].xo - inv_cell[b].xo) <= regu->dist_extra_smooth)
							{
                                num_x_right++;
								idx_x_right_layer_cell = (long *)memory((char *)idx_x_right_layer_cell,num_x_right, sizeof(long),"SetRegu1");
								dist_x_right = (double *)memory((char *)dist_x_right,num_x_right, sizeof(double),"SetRegu1");
									
								idx_x_right_layer_cell[num_x_right-1] = inv_cell[b].used_nr;
								dist_x_right[num_x_right-1] = fabs(inv_cell[a].xo - inv_cell[b].xo);
							}		
						}
					}

					/*Write layer cell relationships to the "tmp_config.bin" file*/

					fwrite(&(inv_cell[a].used_nr),sizeof(long),1,ouf);	/*Write out the index of the active right border inversion cell*/
					fwrite(&kind_of_border, sizeof(int),1,ouf);         /*Write out the kind of border cell.*/
					fwrite(&num_x_right,sizeof(long),1,ouf);				/*Write out the number of in same layer cells*/

					sum_dist = 0.0;
					for(i=0;i<num_x_right;i++)
					{	/*Calculate the sum value of the same layer cells*/
						sum_dist = sum_dist + dist_x_right[i];

						fwrite(&(idx_x_right_layer_cell[i]),sizeof(long),1,ouf);		/*Write out the indeces of the active inversion cells in same layer*/
						fwrite(&(dist_x_right[i]),sizeof(double),1,ouf);		/*Write out the corresponding distance from the border cell*/
					}

					sp_set_val(A, counter_row-1+sti, inv_cell[a].used_nr, sum3 + sp_get_val(A, counter_row-1+sti, inv_cell[a].used_nr));

					for(i=0;i<num_x_right;i++)
						sp_set_val(A, counter_row-1+sti,idx_x_right_layer_cell[i], -(2.0/num_x_right*sum_dist-dist_x_right[i])/sum_dist*sum3);


					free(idx_x_right_layer_cell);
					free(dist_x_right);
				}

				if(inv_cell[a].border_inv_cell_y_front == 1 && inv_cell[a].border_inv_cell_y_back != 1)
				{
					counter_row++;
					num_y_front_border++;
					idx_y_front_border = (long *)memory((char *)idx_y_front_border,num_y_front_border, sizeof(long),"SetRegu1");
					idx_y_front_border[num_y_front_border-1]=inv_cell[a].used_nr;
					kind_of_border = 3;

					idx_y_front_layer_cell = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
					dist_y_front = (double *)memory(NULL,1, sizeof(double),"SetRegu1");
					num_y_front = 0;

					for(b=0;b<inv->nvel;b++)
					{
						if(inv_cell[b].use == 1)
						{
							if( ((0.5*(inv_cell[a].xdim + inv_cell[b].xdim)) > fabs(inv_cell[a].xo - inv_cell[b].xo)) &&
								((0.5*(inv_cell[a].zdim + inv_cell[b].zdim)) > fabs(inv_cell[a].zo - inv_cell[b].zo)) &&
								(a != b) && inv_cell[b].border_inv_cell_y_back != 1
								&& fabs(inv_cell[a].yo - inv_cell[b].yo) <= regu->dist_extra_smooth)
							{
                                num_y_front++;
								idx_y_front_layer_cell = (long *)memory((char *)idx_y_front_layer_cell,num_y_front, sizeof(long),"SetRegu1");
								dist_y_front = (double *)memory((char *)dist_y_front,num_y_front, sizeof(double),"SetRegu1");
									
								idx_y_front_layer_cell[num_y_front-1] = inv_cell[b].used_nr;
								dist_y_front[num_y_front-1] = fabs(inv_cell[a].yo - inv_cell[b].yo);
							}		
						}
					}

					/*Write layer cell relationships to the "tmp_config.bin" file*/

					fwrite(&(inv_cell[a].used_nr),sizeof(long),1,ouf);	/*Write out the index of the active front border inversion cell*/
					fwrite(&kind_of_border, sizeof(int),1,ouf);         /*Write out the kind of border cell.*/
					fwrite(&num_y_front,sizeof(long),1,ouf);				/*Write out the number of in same layer cells*/

					sum_dist = 0.0;
					for(i=0;i<num_y_front;i++)
					{	/*Calculate the sum value of the same layer cells*/
						sum_dist = sum_dist + dist_y_front[i];

						fwrite(&(idx_y_front_layer_cell[i]),sizeof(long),1,ouf);		/*Write out the indeces of the active inversion cells in same layer*/
						fwrite(&(dist_y_front[i]),sizeof(double),1,ouf);		/*Write out the corresponding distance from the border cell*/
					}


					sp_set_val(A, counter_row-1+sti, inv_cell[a].used_nr, sum3 + sp_get_val(A, counter_row-1+sti, inv_cell[a].used_nr));

					for(i=0;i<num_y_front;i++)
						sp_set_val(A, counter_row-1+sti,idx_y_front_layer_cell[i], -(2.0/num_y_front*sum_dist-dist_y_front[i])/sum_dist*sum3);


					free(idx_y_front_layer_cell);
					free(dist_y_front);
				}

				if(inv_cell[a].border_inv_cell_y_back == 1 && inv_cell[a].border_inv_cell_y_front != 1)
				{
					counter_row++;
					num_y_back_border++;
					idx_y_back_border = (long *)memory((char *)idx_y_back_border,num_y_back_border, sizeof(long),"SetRegu1");
					idx_y_back_border[num_y_back_border-1]=inv_cell[a].used_nr;
					kind_of_border = 4;

					idx_y_back_layer_cell = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
					dist_y_back = (double *)memory(NULL,1, sizeof(double),"SetRegu1");
					num_y_back = 0;

					for(b=0;b<inv->nvel;b++)
					{
						if(inv_cell[b].use == 1)
						{
							if( ((0.5*(inv_cell[a].xdim + inv_cell[b].xdim)) > fabs(inv_cell[a].xo - inv_cell[b].xo)) &&
								((0.5*(inv_cell[a].zdim + inv_cell[b].zdim)) > fabs(inv_cell[a].zo - inv_cell[b].zo)) &&
								(a != b) && inv_cell[b].border_inv_cell_y_front != 1
								&& fabs(inv_cell[a].yo - inv_cell[b].yo) <= regu->dist_extra_smooth)
							{
                                num_y_back++;
								idx_y_back_layer_cell = (long *)memory((char *)idx_y_back_layer_cell,num_y_back, sizeof(long),"SetRegu1");
								dist_y_back = (double *)memory((char *)dist_y_back,num_y_back, sizeof(double),"SetRegu1");
								
								idx_y_back_layer_cell[num_y_back-1] = inv_cell[b].used_nr;
								dist_y_back[num_y_back-1] = fabs(inv_cell[a].yo - inv_cell[b].yo);
							}		
						}
					}

					/*Write layer cell relationships to the "tmp_config.bin" file*/

					fwrite(&(inv_cell[a].used_nr),sizeof(long),1,ouf);	/*Write out the index of the active front border inversion cell*/
					fwrite(&kind_of_border, sizeof(int),1,ouf);         /*Write out the kind of border cell.*/
					fwrite(&num_y_back,sizeof(long),1,ouf);				/*Write out the number of in same layer cells*/

					sum_dist = 0.0;
					for(i=0;i<num_y_back;i++)
					{	/*Calculate the sum value of the same layer cells*/
						sum_dist = sum_dist + dist_y_back[i];

						fwrite(&(idx_y_back_layer_cell[i]),sizeof(long),1,ouf);		/*Write out the indeces of the active inversion cells in same layer*/
						fwrite(&(dist_y_back[i]),sizeof(double),1,ouf);		/*Write out the corresponding distance from the border cell*/
					}

					sp_set_val(A, counter_row-1+sti, inv_cell[a].used_nr, sum3 + sp_get_val(A, counter_row-1+sti, inv_cell[a].used_nr));

					for(i=0;i<num_y_back;i++)
						sp_set_val(A, counter_row-1+sti,idx_y_back_layer_cell[i], -(2.0/num_y_back*sum_dist-dist_y_back[i])/sum_dist*sum3);


					free(idx_y_back_layer_cell);
					free(dist_y_back);

				}

			}
	
			if(counter%1000 == 1)
				printf("For %d of %d inversion cells the extra smoothing is inserted in the matrix\n",counter,inv->nvel_used);

		}

		/*Jin Chen changed new*/
		if(counter_row != nr_new_row)
			printf("The number of the new rows is not equal to the countered number, please check it !!\n");

		printf("The smoothing term is introduced into the matrix!\n\n");

		free(idx_x_left_border);
		free(idx_x_right_border);
		free(idx_y_front_border);
		free(idx_y_back_border);

		fclose(ouf);

	}
		/**************************************/
    else /*The smoothing relations are already calculated in a former iteration, therefore the relations 
		  are taken from the temporary file "tmp_config.bin"*/
	{

		inf = fopen(s1,"rb");

        if (inf == NULL)
		{
			fprintf(stderr,"Unable to open %s\n",s1);
			exit(0);
		}

		/*Read in the binary file*/
		fread(&nr_inv_cells,sizeof(long),1,inf);								/*Read in the number of active inversion cells*/
		fread(&nr_new_row,sizeof(long),1,inf);								    /*Read in the number of extra smoothing inversion cells*/

		nr_new_row1 = inv->x_left_border_nr + inv->x_right_border_nr + inv->y_front_border_nr + inv->y_back_border_nr;

		/*Read in the binary file*/


		if((nr_new_row1 != nr_new_row) && (inv->nvel_used != nr_inv_cells))
		{
			printf("The tmp_config1.bin file does not fit to the inversion grid!\n");
			printf("Number of extra smoothing inversion cells in the grid: %d\n",nr_new_row1);
			printf("Number of extra smoothing inversion cells stored in the file: %d\n\n", nr_new_row);
			printf("Number of active inversion cells in the grid: %d\n",inv->nvel_used);
			printf("Number of active inversion cells stored in the file: %d\n\n", nr_inv_cells);
			printf("=> the extra smoothing relations will be re-calculated\n\n");

			/*Close and remove the "tmp_config.bin" file*/
			fclose(inf);
			remove(s1);

			/*Make a new "tmp_config.bin" file*/
			goto make_config1;
		}

        counter_row = 0;
		for(i=0;i<inv->nvel;i++)
		{
			if(inv_cell[i].use == 1)	/*Only ACTIVE inversion cells will be considered*/
			{
				if(inv_cell[i].border_inv_cell_x_left == 1 && inv_cell[i].border_inv_cell_x_right != 1)
				{
					counter_row++;
					fread(&pos,sizeof(long),1,inf);					/*Read in the index of the active inversion cell*/
					fread(&kind_of_border,sizeof(long),1,inf);		/*Read in the index of the kind of border*/		
					fread(&num_x_left,sizeof(long),1,inf);          /*Read in the number of neighboring cells*/

					sum_dist = 0.0;
					
					idx_x_left_layer_cell = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
					dist_x_left = (double *)memory(NULL,1, sizeof(double),"SetRegu1"); 

					for(j=0;j<num_x_left;j++)
					{
						idx_x_left_layer_cell = (long *)memory((char *)idx_x_left_layer_cell,j+1, sizeof(long),"SetRegu1"); 
						dist_x_left = (double *)memory((char *)dist_x_left,j+1, sizeof(double),"SetRegu1"); 

						fread(&idx_x_left_layer_cell[j],sizeof(long),1,inf);					/*Read in the indeces of the active neighboring inversion cells*/
						fread(&dist_x_left[j],sizeof(double),1,inf);				/*Read in the distance of the same layer*/

						sum_dist = sum_dist + dist_x_left[j];
						
					}

					/*Fill in the smoothing constraints into the sparse matrix*/
					sp_set_val(A, counter_row-1 + sti, pos, sum3 + sp_get_val(A, counter_row-1 + sti, pos));

					for(j=0;j<num_x_left;j++)
					{
						sp_set_val(A, counter_row-1 + sti, idx_x_left_layer_cell[j], -(2.0/num_x_left*sum_dist-dist_x_left[j])/sum_dist*sum3);
					}

					free(idx_x_left_layer_cell);
					free(dist_x_left);
				}

				if(inv_cell[i].border_inv_cell_x_right == 1 && inv_cell[i].border_inv_cell_x_left != 1)
				{
					counter_row++;
					fread(&pos,sizeof(long),1,inf);					/*Read in the index of the active inversion cell*/
					fread(&kind_of_border,sizeof(long),1,inf);		/*Read in the index of the kind of border*/		
					fread(&num_x_right,sizeof(long),1,inf);          /*Read in the number of neighboring cells*/

					sum_dist = 0.0;
					
					idx_x_right_layer_cell = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
					dist_x_right = (double *)memory(NULL,1, sizeof(double),"SetRegu1"); 

					for(j=0;j<num_x_right;j++)
					{
						idx_x_right_layer_cell = (long *)memory((char *)idx_x_right_layer_cell,j+1, sizeof(long),"SetRegu1"); 
						dist_x_right = (double *)memory((char *)dist_x_right,j+1, sizeof(double),"SetRegu1"); 

						fread(&idx_x_right_layer_cell[j],sizeof(long),1,inf);					/*Read in the indeces of the active neighboring inversion cells*/
						fread(&dist_x_right[j],sizeof(double),1,inf);				/*Read in the distance of the same layer*/

						sum_dist = sum_dist + dist_x_right[j];
						
					}

					/*Fill in the smoothing constraints into the sparse matrix*/
					sp_set_val(A, counter_row-1 + sti, pos, sum3 + sp_get_val(A, counter_row-1 + sti, pos));

					for(j=0;j<num_x_right;j++)
					{
						sp_set_val(A, counter_row-1 + sti, idx_x_right_layer_cell[j], -(2.0/num_x_right*sum_dist-dist_x_right[j])/sum_dist*sum3);
					}

					free(idx_x_right_layer_cell);
					free(dist_x_right);
				}

				if(inv_cell[i].border_inv_cell_y_front == 1 && inv_cell[i].border_inv_cell_y_back != 1)
				{
					counter_row++;
					fread(&pos,sizeof(long),1,inf);					/*Read in the index of the active inversion cell*/
					fread(&kind_of_border,sizeof(long),1,inf);		/*Read in the index of the kind of border*/		
					fread(&num_y_front,sizeof(long),1,inf);          /*Read in the number of neighboring cells*/

					sum_dist = 0.0;
					
					idx_y_front_layer_cell = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
					dist_y_front = (double *)memory(NULL,1, sizeof(double),"SetRegu1"); 

					for(j=0;j<num_y_front;j++)
					{
						idx_y_front_layer_cell = (long *)memory((char *)idx_y_front_layer_cell,j+1, sizeof(long),"SetRegu1"); 
						dist_y_front = (double *)memory((char *)dist_y_front,j+1, sizeof(double),"SetRegu1"); 

						fread(&idx_y_front_layer_cell[j],sizeof(long),1,inf);					/*Read in the indeces of the active neighboring inversion cells*/
						fread(&dist_y_front[j],sizeof(double),1,inf);				/*Read in the distance of the same layer*/

						sum_dist = sum_dist + dist_y_front[j];
						
					}

					/*Fill in the smoothing constraints into the sparse matrix*/
					sp_set_val(A, counter_row-1 + sti, pos, sum3 + sp_get_val(A, counter_row-1 + sti, pos));

					for(j=0;j<num_y_front;j++)
					{
						sp_set_val(A, counter_row-1 + sti, idx_y_front_layer_cell[j], -(2.0/num_y_front*sum_dist-dist_y_front[j])/sum_dist*sum3);
					}

					free(idx_y_front_layer_cell);
					free(dist_y_front);
				}

				if(inv_cell[i].border_inv_cell_y_back == 1 && inv_cell[i].border_inv_cell_y_front != 1)
				{
					counter_row++;
					fread(&pos,sizeof(long),1,inf);					/*Read in the index of the active inversion cell*/
					fread(&kind_of_border,sizeof(long),1,inf);		/*Read in the index of the kind of border*/		
					fread(&num_y_back,sizeof(long),1,inf);          /*Read in the number of neighboring cells*/

					sum_dist = 0.0;
					
					idx_y_back_layer_cell = (long *)memory(NULL,1, sizeof(long),"SetRegu1"); 
					dist_y_back = (double *)memory(NULL,1, sizeof(double),"SetRegu1"); 

					for(j=0;j<num_y_back;j++)
					{
						idx_y_back_layer_cell = (long *)memory((char *)idx_y_back_layer_cell,j+1, sizeof(long),"SetRegu1"); 
						dist_y_back = (double *)memory((char *)dist_y_back,j+1, sizeof(double),"SetRegu1"); 

						fread(&idx_y_back_layer_cell[j],sizeof(long),1,inf);					/*Read in the indeces of the active neighboring inversion cells*/
						fread(&dist_y_back[j],sizeof(double),1,inf);				/*Read in the distance of the same layer*/

						sum_dist = sum_dist + dist_y_back[j];
						
					}

					/*Fill in the smoothing constraints into the sparse matrix*/
					sp_set_val(A, counter_row-1 + sti, pos, sum3 + sp_get_val(A, counter_row-1 + sti, pos));

					for(j=0;j<num_y_back;j++)
					{
						sp_set_val(A, counter_row-1 + sti, idx_y_back_layer_cell[j], -(2.0/num_y_back*sum_dist-dist_y_back[j])/sum_dist*sum3);
					}

					free(idx_y_back_layer_cell);
					free(dist_y_back);
				}
			}
		}
	
		/*Jin Chen changed new*/
		if(counter_row != nr_new_row1)
			printf("The number of the new row is not equal to the countered number, please check it !!\n");

		printf("The extra smoothing is introduced into the matrix! \n\n");

		fclose(inf);

	}

	return(1);
}
/*------------------------------------------------------------*/
/*Calculate the normalisation (sum of all sensitivities/sum of all inversion cells unequal 0)*/
/*Parameter:	*A		:= Sparse matrix*/
/*			first_row	:= First row of A*/
/*			last_row	:= Last row of A*/
/*			first_col	:= First column of A*/
/*			last_col	:= Last column of A*/
/*          kind_of_reg := Specify the regulraziation (1 = Divide by the number of active columns, 2 = Divide by (the number of active columns * number of data), 3 = Divide by the number of rows)*/
/*Output: Value for the normalisation*/
/*        c: the number of the columns of the nonzero entries. */

double RegFactor(SPMAT *A, int first_row, int last_row, int first_col, int last_col, int kind_of_reg, long *c)
{
	long a,b,d;
	long double sum = 0.0;	/*Sum of all sensitivities*/
	long double sum2;		/*Sum of all sensitivities in an inversion cell*/

	d = 0;


	/*Loop over all columns*/
	for(a=first_col;a<last_col;a++)
	{
		sum2 = 0.0;

		/*Loop over all rows*/
		for(b=first_row;b<last_row;b++)
			sum2 = fabs(sp_get_val(A,b,a)) + sum2;

		if(sum2 != 0.0) /*Inversion cell has non-zero entries*/
		{
			sum = sum2 + sum;
			d++;
		}
	}

	(*c)=d;

	/*Normalize by the number of active columns*/
	if(d != 0 && kind_of_reg == 1)
	{
//		printf("Sum of all matrix entries: %lf\n",sum);
//		printf("Number of all non-zero matrix columns: %d\n",c);

		sum = sum/d;
	}
	/*Normalize by the number of active columns multiplied with the number of rows*/
	else if((d * (last_row - first_row)) != 0 && kind_of_reg == 2)
	{
//		printf("Sum of all matrix entries: %lf\n",sum);
		sum = sum/(d * (last_row - first_row));
	}
	/*Normalize by the number of rows*/
	else if ((last_row - first_row) != 0 && kind_of_reg != 2 && kind_of_reg != 1)
	{
//		printf("Sum of all matrix entries: %lf\n",sum);
		sum = sum/(last_row - first_row);
	}
	else
	{
		printf("WARNING !!! All entries are empty in the range: rows %d to %d; columns %d to %d\n\n", first_row, last_row, first_col, last_col);
		sum = 0.0;
	}

	 return((double)sum);
}



/*------------------------------------------------------------*/
/*All velocities below/above the min/max allowed velocity are assigned to*/
/*the min/max allowed velocity*/
/*Parameter: inv := inversion structure*/

int	MinMaxConstrainsSeis(INV *inv)
{
	long i;

	for(i=0;i<inv->nvel_used;i++)
	{
		/*Readjust the velocity that has changed too much during a iteration*/
		if(fabs((1/inv->para_mod[i]) - (1/inv->mx[i])) > inv->maxadj && inv->para_mod[i] > inv->mx[i])
			inv->para_mod[i] = 1/((1/inv->mx[i]) - inv->maxadj);
		else if(fabs((1/inv->para_mod[i]) - (1/inv->mx[i])) > inv->maxadj && inv->para_mod[i] < inv->mx[i])
			inv->para_mod[i] = 1/((1/inv->mx[i]) + inv->maxadj);
		/*Reset the velocities that are smaller than the min. allowed velocity*/
		if((1/inv->para_mod[i]) < inv->vmin)
			inv->para_mod[i] = (1/inv->vmin);
		/*Reset the velocities that are higher than the max. allowed velocity*/
		if((1/inv->para_mod[i]) > inv->vmax)
			inv->para_mod[i] = (1/inv->vmax);
	}

	return(1);
}

/*------------------------------------------------------------*/
/*All densities below/above the min/max allowed density are assigned to*/
/*the min/max allowed density*/
/*Parameter: inv := inversion structure*/

int MinMaxConstrainsGrav(INV *inv)
{
	long i;

	for(i=0;i<inv->nvel_used;i++)
	{
		/*Reset the densities that are smaller than the min.allowed density*/
		if(inv->para_mod[i] < inv->gmin)
			inv->para_mod[i] = inv->gmin;
		/*Reset the densities that are higher than the max. allowed density*/
		if(inv->para_mod[i] > inv->gmax)
			inv->para_mod[i] = inv->gmax;
	}

	return(1);
}

/*------------------------------------------------------------*/
/*All resistivities below/above the min/max allowed resistivity are assigned to*/
/*the min/max allowed resistivity*/
/*Parameter: inv := inversion structure*/

int MinMaxConstrainsRes(INV *inv)
{
	long i;

	for(i=0;i<inv->nvel_used;i++)
	{
		/*Reset the resistivities that are smaller than the min.allowed resistivity*/
		if(inv->para_mod[i] < inv->rmin)
			inv->para_mod[i] = inv->rmin;
		/*Reset the resistivities that are higher than the max. allowed resistivity*/
		if(inv->para_mod[i] > inv->rmax)
			inv->para_mod[i] = inv->rmax;
	}

	return(1);
}

/*------------------------------------------------------------*/
/*Calculate the density or resistivity from the seismic model vector and afterwards */
/*adjust the seismic model vector to the constraints of the density or resistivity values*/
/*Parameter: *seis_mod := Seismic model vector*/
/*			      mod2 := corresponding density/resistivity values*/
/*				nvel   := Number of active inversion cells*/
/*			vmin,vmax  := Lower and upper threshold for the density or resistivity values*/

int CalcResDensfromVel(double *seis_mod, double *mod2, long nvel, double vmin, double vmax, REL_VDR_STRUCT rel_vdr)
{
	double eps;
	int index1,index2;
	long i;

	index1 = 0;
	index2 = 0;
	eps = 1E-11; /*ensure that the values are always slightly positive*/

	/*Loop over all active inversion cells*/
	for(i=0;i<nvel;i++)
		/*Determine the density/resistivity from the velocity*/
		mod2[i] = VeltoDensRes((1/seis_mod[i]), rel_vdr);		


	/*Sort the densities in ascending order for tabular links*/
	if(rel_vdr.index == 9)
		sort2(rel_vdr.tab_nr_pairs, rel_vdr.tab_dr, rel_vdr.tab_v);

	/*Loop over all active inversion cells*/
	for(i=0;i<nvel;i++)
	{
		/*Find density/resistivity values that are below the lower boundary value*/
		if(mod2[i] < vmin)
		{
			index1 = 1;
			/*Set the densities/resistivity to the value of the lower boundary*/
			mod2[i] = vmin + eps;
			/*... and re-adjust the velocities*/
			seis_mod[i] = 1/(DensRestoVel(mod2[i],rel_vdr));
		}

		/*Find density/resistivity values that are above the upper boundary value*/
		if(mod2[i] > vmax)
		{
			index2 = 1;
			/*Set the densities/resistivity to the value of the upper boundary*/
			mod2[i] = vmax;
			/*... and re-adjust the velocities*/
			seis_mod[i] = 1/(DensRestoVel(mod2[i],rel_vdr));
		}
	}

	if(index1 == 1)
	{
		printf("REMARK: Density or resistivity values <%f g/cm^3 or ohmm were found\nafter the inversion\n",vmin);
		printf("These values were set to: %f\n", vmin+eps);
	}

	if(index2 == 1)
	{
		printf("REMARK: Density or resistivity values >%f g/cm^3 or ohmm were found\nafter the inversion\n",vmax);
		printf("These values were set to: %f\n", vmax);
	}

	return(1);
}

/*------------------------------------------------------------*/
/*Assign the new slowness values to the data structure*/
/*Parameter: *grid  := Grid structure*/
/*			 nvel	:= Number of inversion cells*/
/*			 *mod   := model vector (velocity)*/
/*			 *inv_cell := Inversion cell structure*/

int	MakeNewSlowMod(GRID_STRUCT *grid,long nvel, double *mod, BIDX_STRUCT *inv_cell)
{
	long i,j,a,b,c;
	long nx2,ny2,nz2,nyz2,nxyz2;
	long ny3,nz3,nyz3;
	double *tmp_slow;

	#define tslow(x,y,z) tmp_slow[nyz2*(x) + nz2*(y) + (z)]
	#define slow(x,y,z) grid->slow[nyz3*(x) + nz3*(y) + (z)] 

	nx2 = grid->nx + 2*grid->nborder;
	ny2 = grid->ny + 2*grid->nborder;
	nz2 = grid->nz + 2*grid->nborder;
	nyz2 = ny2*nz2;
	nxyz2 = nx2*ny2*nz2;

	ny3 = grid->ny + 2*grid->nborder+1;
	nz3 = grid->nz + 2*grid->nborder+1;
	nyz3 = ny3*nz3;
	
	tmp_slow = (double *)memory(NULL,nxyz2,sizeof(double),"MakeNewSlowMod");
	for(i=0;i<nxyz2;i++)
		tmp_slow[i] = -99999.9;

	/*Loop over all inversion cells*/
	for(i=0;i<nvel;i++)
	{
		if(inv_cell[i].use == 1)
		{
			/*Loop over all active foward cells in the inversion cell*/
			for(j=0;j<inv_cell[i].nele;j++)
			{
				/*Write the modified slownesses in a temporary slowness file*/
				tmp_slow[inv_cell[i].ele[j]] = mod[inv_cell[i].used_nr] * grid->h;
			}
		}
	}

	/*Modify the grid structure*/
	for(a=0;a<nx2;a++)
		for(b=0;b<ny2;b++)
			for(c=0;c<nz2;c++)
			{
				if(tslow(a,b,c) != -99999.9)
				{
					slow(a,b,c) = tslow(a,b,c);
				}
			}




	#undef slow
	#undef tslow

	free(tmp_slow);

	return(1);
}


/*------------------------------------------------------------*/
/*Assign the new density/resistivity values to the data structure*/
/*Parameter: *grid  := Grid structure*/
/*			 nvel	:= Number of inversion cells*/
/*			 *mod   := model vector (density or resistivity)*/
/*			 *inv_cell := Inversion cell structure*/
/*		 kind_of_model := (density = 1; resistivity =2)*/

int	MakeNewDensResMod(GRID_STRUCT *grid,long nvel, double *mod, BIDX_STRUCT *inv_cell, int kind_of_model)
{
	long i,j,a,b,c,a1,b1,c1;
	long nx,ny,nz,nyz,nxyz;
	long nx2,ny2,nz2,nyz2,nxyz2;
	double *tmp_res_dens;

	#define tmod(x,y,z) tmp_res_dens[nyz2*(x) + nz2*(y) + (z)]
	#define dens(x,y,z) grid->dens[nyz*(x) + nz*(y) + (z)] 
	#define res(x,y,z) grid->res[nyz*(x) + nz*(y) + (z)]

	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nyz = ny*nz;
	nxyz = nx*ny*nz;

	nx2 = nx + 2*grid->nborder;
	ny2 = ny + 2*grid->nborder;
	nz2 = nz + 2*grid->nborder;
	nyz2 = ny2*nz2;
	nxyz2 = nx2*ny2*nz2;
	
	tmp_res_dens = (double *)memory(NULL,nxyz2,sizeof(double),"MakeNewDensResMod");
	for(i=0;i<nxyz2;i++)
		tmp_res_dens[i] = -99999.9;

	/*Loop over all inversion cells*/
	for(i=0;i<nvel;i++)
	{
		if(inv_cell[i].use == 1)
		{
			/*Loop over all active foward cells in the inversion cell*/
			for(j=0;j<inv_cell[i].nele;j++)
			{
				/*Write the modified density/resistivity in a temporary density/resistivity file*/
				if(kind_of_model == 1)
					tmp_res_dens[inv_cell[i].ele[j]] = mod[inv_cell[i].used_nr];	/*density*/					
				else
					tmp_res_dens[inv_cell[i].ele[j]] = mod[inv_cell[i].used_nr];	/*resistivity*/
			}
		}
	}

	/*Modify the grid structure*/
	for(a=0;a<nx;a++)
		for(b=0;b<ny;b++)
			for(c=0;c<nz;c++)
			{
				a1 = a + grid->nborder;
				b1 = b + grid->nborder;
				c1 = c + grid->nborder; 

				if(tmod(a1,b1,c1) != -99999.9)
				{
					if(kind_of_model == 1)
						dens(a,b,c) = tmod(a1,b1,c1);
					else
						res(a,b,c) = tmod(a1,b1,c1);
				}
			}


	#undef tmod
	#undef res 
	#undef dens

	free(tmp_res_dens);

	return(1);
}


/*------------------------------------------------------------*/
	/*Copy results to inv->para_mod (parameters that is directly determined by the inversion process)*/
/*Parameter:  inv := inversion structure*/
/*           flag := flag structure*/
/*          v_mod := VEC */
/*      *fname_vel_dens  := Filename of the parameterfile including the relationship of velocity and density*/
/*      *fname_vel_res  := Filename of the parameterfile including the relationship of velocity and resistivity*/
/*		 rel_vdr3 := Structure that organizes the relationship of velocity and density and velocity and resistivity*/
/*  	 rel_vdr4 := Structure that organizes the relationship of velocity and density and velocity and resistivity*/

int	CopyResultsToInvParaMod(INV *inv, FLAG_STRUCT *flag, VEC *v_mod, char *fname_vel_dens, char *fname_vel_res)
{
	int i;
	REL_VDR_STRUCT  rel_vdr3, rel_vdr4;

	for(i=0;i<inv->nvel_used;i++)
	{
		/*Consider the "creeping" case (without smoothing)*/
		if(flag->do_smooth == 0)
			inv->para_mod[i] = v_mod->ve[i] + inv->mx[i];

		/*Consider the "jumping" case (with smoothing)*/
		else
			inv->para_mod[i] = v_mod->ve[i];
	}

	/*****************************************************/
	if(flag->index_tseis != 0 || (flag->index_grav != 0 && flag->index_mt != 0))
	{
		/*Remove unexspected values (for velocities)*/
		MinMaxConstrainsSeis(inv);
		
		/***************/
		/*Determine the second parameter from joint inversion (gravity OR resistivity)*/
		if(flag->index_grav != 0 || flag->index_mt != 0)
		{
			if(flag->index_grav != 0)
			{
				/*Read in the file including the relationship of velocity and density*/
				ReadRelVelDensRes(fname_vel_dens, &rel_vdr3);

				/*Determine the density from the velocity AND modify the velocity by means of the density constraints*/
				CalcResDensfromVel(inv->para_mod, inv->para_mod2, inv->nvel_used, inv->gmin, inv->gmax, rel_vdr3);
			}
			else
			{
				/*Read in the file including the relationship of velocity and resistivity*/
				ReadRelVelDensRes(fname_vel_res, &rel_vdr3);

				/*Determine the resistivity from the velocity AND modify the velocity by means of the resistivity constraints*/
				CalcResDensfromVel(inv->para_mod, inv->para_mod2, inv->nvel_used, inv->rmin, inv->rmax, rel_vdr3);
			}


			free(rel_vdr3.tab_v);
			free(rel_vdr3.tab_dr);
		}
	
		/***************/
		/*Determine the third parameter from joint inversions (resistivity)*/
		if(flag->index_grav != 0 && flag->index_mt != 0)
		{
			/*Read in the file including the relationship of velocity and resistivity*/
			ReadRelVelDensRes(fname_vel_res, &rel_vdr4);

			/*Determine the resistivity from the velocity AND modify the velocity by means of the resistivity constraints*/
			CalcResDensfromVel(inv->para_mod, inv->para_mod3, inv->nvel_used, inv->rmin, inv->rmax, rel_vdr4);

			free(rel_vdr4.tab_v);
			free(rel_vdr4.tab_dr);
		}
		
	}
	else if(flag->index_grav != 0)
	{
		/*Remove unexspected values (for density)*/
		MinMaxConstrainsGrav(inv);
	}
	else if(flag->index_mt != 0)
	{
		/*Remove unexspected values (for resistivity)*/
		MinMaxConstrainsRes(inv);
	}

	return(1);

}

/*------------------------------------------------------------*/
	/*Assign the modified slowness/density and resistivity values to the forward cells*/
/*Parameter:  flag := flag structure*/
/*            grid := grid structure*/
/*             inv := inversion structure*/
/*        inv_cell := BIDX_STRUCTURE*/

int AssignValueToForward(FLAG_STRUCT *flag, GRID_STRUCT *grid, INV *inv, BIDX_STRUCT *inv_cell)
{
	if(flag->index_tseis != 0 || (flag->index_grav != 0 && flag->index_mt != 0))
	{
		MakeNewSlowMod(grid,inv->nvel,inv->para_mod,inv_cell); /*velocity in simple or joint inversion*/

		/*joint inversions; second parameter*/
		if(flag->index_grav != 0)
			MakeNewDensResMod(grid,inv->nvel,inv->para_mod2,inv_cell,1); /*density in joint inversion*/
		else if(flag->index_mt != 0)
			MakeNewDensResMod(grid,inv->nvel,inv->para_mod2,inv_cell,2); /*resistivity in joint inversion*/

		/*joint inversions; third parameter*/
		if(flag->index_grav != 0 && flag->index_mt != 0)
			MakeNewDensResMod(grid,inv->nvel,inv->para_mod3,inv_cell,2); /*resistivity in joint inversion*/
	}
	else if(flag->index_grav != 0) /*density in simple inversion*/
		MakeNewDensResMod(grid,inv->nvel,inv->para_mod,inv_cell,1);
	else if(flag->index_mt != 0) /*resistivity in simple MT inversion*/
		MakeNewDensResMod(grid,inv->nvel,inv->para_mod,inv_cell,2);

	return(1);
}


/*------------------------------------------------------------*/
	/*Jin Chen changed new*/
	/*Calculate the border cell numbers of inversion cells.*/
/*Parameter: 
/*         inv := inversion structure*/
/*         inv_cell := BIDX_STRUCTURE*/

int CalcNumBorder(INV *inv, BIDX_STRUCT *inv_cell)
{
	int i;
	inv->x_left_border_nr = 0;
	inv->x_right_border_nr = 0;
	inv->y_front_border_nr = 0;
	inv->y_back_border_nr = 0;

	for(i=0;i<inv->nvel;i++)
	{
		if(inv_cell[i].use == 1)
		{
			if(inv_cell[i].border_inv_cell_x_left == 1 && inv_cell[i].border_inv_cell_x_right !=1)
				inv->x_left_border_nr++;
			if(inv_cell[i].border_inv_cell_x_right == 1 && inv_cell[i].border_inv_cell_x_left != 1)
				inv->x_right_border_nr++;
			if(inv_cell[i].border_inv_cell_y_front == 1 && inv_cell[i].border_inv_cell_y_back != 1)
				inv->y_front_border_nr++;
			if(inv_cell[i].border_inv_cell_y_back == 1 && inv_cell[i].border_inv_cell_y_front != 1)
				inv->y_back_border_nr++;
		}
	}
	return(1);
}


/*------------------------------------------------------------*/
/*Perform the extra cell inversion; this means that an additional cell will be introduced for the gravity*/
/*Parameter:	*grid	:= Grid structure*/
/*				inv		:= Inversion structure */
/*				*grav			:= Gravity structure*/
/*				 *A				:= Sparse matrix*/
/*				data			:= Data structure*/
/*				*regu	:= Regularisation structure*/
/*				*res	:= Residual vector*/
/*				flag    := Flag structure*/
/*				weight_grav := ext.weighting of the gravity (due to balancing of the different data sets)*/

int InvExtraCell(GRID_STRUCT *grid, INV *inv, GRAV_STRUCT *grav, SPMAT *A, DATA_STRUCT *data, REGU_STRUCT *regu, double *res, FLAG_STRUCT  flag, double weight_grav) /*BJOERN_MOD5*/
{
	long i;
	double rho1, rho2; /*Perturbations of the densities: Rho2-Rho1 in mgal*/
	double deltaSlow, grav1,grav2, deriv_grav_extra;

	#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/

	/*Calculate the gravity derivative of the extra cell*/
	deltaSlow = 0.0001;

	rho2 = inv->dens_extra_cell + (deltaSlow/2.0);
	rho1 = inv->dens_extra_cell - (deltaSlow/2.0);

	grav2 = 2*PI*G_KONST*rho2*(grid->h * grid->nz); /*infinite sheet*/
	grav1 = 2*PI*G_KONST*rho1*(grid->h * grid->nz); /*infinite sheet*/

	deriv_grav_extra = 1.0E8 * (grav2 -grav1)/deltaSlow; /*Calculate the total derivatives*/

	/*Assign values to the last column of the sparse matrix*/
	/*Loop over all stations*/
	for(i=0;i<data->ndata_grav;i++)
	{
		/*Fill the last column of the sparse matrix*/
		sp_set_val(A,grav[i].n,inv->nvel_used,deriv_grav_extra*weight_grav); /*BJOERN_MOD5*/
	}

	/*Determine the last parameter of the original model vector*/
	inv->mx[inv->nvel_used] = inv->dens_extra_cell;

	/*Fill the residual vector*/
	/*Loop over all rays*/
	if(flag.do_smooth != 0)
	{
		for(i=0;i<data->ndata_grav;i++)
			res[grav[i].n] = res[grav[i].n] + (inv->mx[inv->nvel_used] * sp_get_val(A,grav[i].n,inv->nvel_used));

		res[inv->nrow] = (regu->w_extra_cell*deriv_grav_extra*data->ndata_grav*weight_grav)*inv->mx[inv->nvel_used]; 
	}
	

	/*Set damping factor to the extra cell*/	
	sp_set_val(A,inv->nrow,inv->ncol,regu->w_extra_cell*deriv_grav_extra*data->ndata_grav*weight_grav); /*BJOERN_MOD5*/

	printf("\nThe rubbish bin for the gravity is added\n\n");
	printf("----------------\n");

	#undef G_KONST

	return(1);
}


/*------------------------------------------------------------*/
/*Transfer the value of the cell controlling the gravity shift in the data to the inversion structure*/
/*Parameter:  inv := inversion structure*/
/*           flag := flag structure*/
/*          v_mod := (modified) model vector */
/*          index := index that specify the position of the column that includes the extra cell value*/

int TransExtraCellValue(INV *inv, FLAG_STRUCT *flag, VEC *v_mod, int index)
{
	/*Consider the "creeping" case (without smoothing)*/
	if(flag->do_smooth == 0)
		inv->dens_extra_cell = v_mod->ve[index] + inv->mx[index];
	
	/*Consider the "jumping" case (with smoothing)*/
	else
		inv->dens_extra_cell = v_mod->ve[index];
	
	return(1);
}
