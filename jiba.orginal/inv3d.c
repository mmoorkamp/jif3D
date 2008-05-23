/* Program is developed by Bjoern Heincke */
#define _CRTDBG_MAP_ALLOC

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "meschachd_no_complex/iter.h"
#include "meschachd_no_complex/matrix2.h"
#include "meschachd_no_complex/sparse.h"
#include "inv3d.h"
#include "fkt.h"

//#include <crtdbg.h>

main(int argc, char *argv[])
{ 
  int i,j,ngeo_file;	/*number of geometry files*/ 
  int nmod_file;    /*number of model files*/
  long nr_of_slices_2D_mt, tmp_nr_of_slices_2D_mt; /*Number of 2D MT slices*/
  int mt_data_format; /*Found MT data format*/
  char fname[128];

  long *sample_topo;        /*Indices of the cells that specify the topography*/

  PAR par;					/*Input parameter from parameter file*/
  GRID_STRUCT grid;			/*Structure including all information about the grid (inversion)*/ 
  GEOMETRY *geo,			/*unmerged geometries*/ 
			geo_all;		/*all geometries merged*/
  TOPO_STRUCT topo;			/*Topography structure*/
  DATA_STRUCT *data,		/*unmerged data information*/ 
			   data_all;	/*all data structures merged */
  INV inv;					/*Structure for inversion*/
  RP_STRUCT	*raypath;		/*Raypath structures including ray-positions and lengths*/
  F_RP_STRUCT *fat_rays;	/*Fat_ray structure including the fat-ray locations and their weight*/
  GRAV_STRUCT *grav;		/*Structure including the gravity derivatives*/
  MT_STRUCT *mt;			/*Structure including the MT derivatives*/
  CALC_2D_MT *calc_2D_mt;    /*Structure for the 2-D MT forward modelling*/
  REGU_STRUCT regu;			/*Parameters of the regularistion*/
  FLAG_STRUCT flag;			/*Flags*/

  FILE *out_rms;

  time_t time_start; /*Computing time*/
  
	/*Start memory check*/
   //  _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );


   /*--------------------------------------------------------------*/
   /* Start the timer:*/		
	time(&time_start);

   /* Read parameter-file*/
   if (argc < 2) strcpy(fname,"inv3d.inp"); else strcpy(fname,argv[1]);
   ReadInput(fname,&par,&ngeo_file,&nmod_file); 

   /*Set the flags*/
   set_flags(par,&flag);
 
   /*Set some regularisation parameters*/
   set_regu(par, &regu);

	/*Set some of the inversion parameters*/
   set_inv(par, &inv);

   /*--------------------------------------------------------------*/

	/*Make the geometry and data structure*/
	geo = (GEOMETRY *)memory(NULL,ngeo_file,sizeof(GEOMETRY),"Main");
	data = (DATA_STRUCT *)memory(NULL,ngeo_file,sizeof(DATA_STRUCT),"Main");

	mt_data_format = MakeGeoandDataStructure(geo,par.files_geo,&geo_all,data,par.files_dat,&data_all,ngeo_file,par.eps, &flag.kind_of_data_mt);

	/*specify the number of mt modes that will be used for inversion*/
	if(flag.dimension_mt != 1 && flag.kind_of_data_mt != 1 && flag.kind_of_data_mt != 2)
		flag.nr_of_modes = 2;
	else
		flag.nr_of_modes = 1;

	/*Write out the merged data file*/
	WriteMergedDatOut(data_all, mt_data_format);

	/*-----------------------------------------------------------*/
	if(par.binary_index != 1) /*The starting model will be built by means of the parameters*/
	{
		free(par.binary_mod_dat);

		/* Initialize the grid structure and the gradient structures for the starting models*/
		MakeGridStruct(par,&grid);
		if(flag.index_tseis != 0) /*for the velocity model*/
			MakeGradientStruct(par.index_velocity_field,par.vg[0],par.vg[1],par.vg[2],par.org_vo[0],par.org_vo[1],par.org_vo[2],par.vo,par.vmin,par.vmax,&(grid.grad_vel));
		if(flag.index_grav != 0) /*for the gravity model*/
			MakeGradientStruct(par.index_density_field,par.gg[0],par.gg[1],par.gg[2],par.org_go[0],par.org_go[1],par.org_go[2],par.go,par.gmin,par.gmax,&(grid.grad_dens));
		if(flag.index_mt != 0) /*for the elec.resistivity model*/
			MakeGradientStruct(par.index_resistivity_field,par.rg[0],par.rg[1],par.rg[2],par.org_ro[0],par.org_ro[1],par.org_ro[2],par.ro,par.rmin,par.rmax,&(grid.grad_res));

		/*-----------------------------------------------------------*/
		/* Initialize velocity/density/resistivity model*/
		MakeInitModels(&grid, &geo_all, flag, par.file_seis_grav_rel->name, par.file_seis_res_rel->name, &topo);
	}
	else /*use a binary file to determine the model*/
	{

		ReadFromBinary(nmod_file, &grid, par, &flag, geo_all, par.file_seis_grav_rel->name);

		for(i=0;i<3;i++)
			grid.grid_cell_factor[i] = par.grid_cell_factor[i];


		/*Read out topography information*/
		sample_topo = TopoCalc(&geo_all, &grid, &topo);
		free(sample_topo);

		free(par.binary_mod_dat);
	}

	/*-----------------------------------------------------------*/
	/*Check if all shots and receivers are inside the cube*/
	 CheckCoordinates(geo_all, data_all, grid, flag);

	/*Check if no shots and receivers are in the "air"*/
	 CheckCoordinatesTopo(geo_all, data_all, grid);
   /*-----------------------------------------------------------*/

   	/*Write out the binary file of the velocity model*/
	if(flag.index_tseis != 0 || (flag.index_grav != 0 && flag.index_mt != 0))
		WriteModSeisOut(inv.ninv, data_all, geo_all, grid);
	if(flag.index_grav != 0)
		WriteModGravOut(inv.ninv, data_all, geo_all, grid);
	if(flag.index_mt != 0)
		WriteModResOut(inv.ninv, data_all, geo_all, grid);

	/*Write out the binary file of the data*/
	if(flag.index_tseis != 0)
		WriteDatSeisOut(data_all,inv.ninv,geo_all);
	if(flag.index_grav != 0)
		WriteDatGravOut(data_all,inv.ninv,geo_all);
	if(flag.index_mt != 0)
		WriteDatMTOut(data_all,inv.ninv,geo_all, flag.dimension_mt, flag.kind_of_data_mt);
   /*-----------------------------------------------------------*/

	out_rms = fopen("rms.txt","wt");

	/*HERE STARTS THE ITERATION-LOOP:*/
	while(inv.ninv < inv.ninv_max)
	{
		fprintf(out_rms,"Number of inversion: %d\n", inv.ninv);
		fflush(out_rms);

		/*--------------Forward modeling----------------------------*/
		
		printf("\n\n\n\n\n");
		printf("Start forward modeling: After iteration: %d\n\n", inv.ninv);

		/*------------------ Seismic --------------------------------*/
		raypath =(RP_STRUCT *)memory(NULL,1,sizeof(RP_STRUCT),"Main");
		fat_rays = (F_RP_STRUCT *)memory(NULL,1,sizeof(F_RP_STRUCT),"Main");

		if(flag.index_tseis != 0)
		{
			printf("Start seismic modeling\n---------------------\n");

			fprintf(out_rms,"Seismic traveltimes:\n");
			fprintf(out_rms,"RMS-value[ms] Act.rays All rays\n");
			fflush(out_rms);

			/*Forward algorithm (Podvin&Lecomte):*/
			if(flag.kind_of_rays == 1)	/*Using Conventional rays*/
			{
				if(data_all.ndata_seis != 0)
					raypath =(RP_STRUCT *)memory((char *)raypath,data_all.ndata_seis,sizeof(RP_STRUCT),"Main");
				else
					raypath =(RP_STRUCT *)memory((char *)raypath,1,sizeof(RP_STRUCT),"Main");
				fat_rays = (F_RP_STRUCT *)memory((char *)fat_rays,1,sizeof(F_RP_STRUCT),"Main");

				ForwardModRay(geo_all, grid, &data_all, raypath, time_start);
			}
			else /*Using Fat-rays*/
			{
				if(data_all.ndata_seis != 0)
					fat_rays = (F_RP_STRUCT *)memory((char *)fat_rays,data_all.ndata_seis,sizeof(F_RP_STRUCT),"Main");
				else
					fat_rays = (F_RP_STRUCT *)memory((char *)fat_rays,1,sizeof(F_RP_STRUCT),"Main");
				raypath =(RP_STRUCT *)memory((char *)raypath,1,sizeof(RP_STRUCT),"Main");

				for(i=0;i<data_all.ndata_seis;i++) 
				{
					fat_rays[i].fatb = par.fatb;
					fat_rays[i].fatthres = par.fatthres;
				}
				ForwardModFatRay(geo_all, grid, &data_all, fat_rays, time_start);
			}

			/*Write out the binary file of the seismic data*/
			WriteDatSeisOut(data_all,inv.ninv,geo_all);
			/*Write out the rms-values*/

			fprintf(out_rms," %f %d %d\n", data_all.rms_seis, data_all.ndata_seis_act, data_all.ndata_seis);
			fflush(out_rms);
		}

		/*------------------ Gravity --------------------------------*/
		grav =(GRAV_STRUCT *)memory(NULL,1,sizeof(GRAV_STRUCT),"Main");

		if(flag.index_grav != 0)
		{
			printf("Start gravity modeling\n---------------------\n");

				if(data_all.ndata_grav != 0)
					grav =(GRAV_STRUCT *)memory((char *)grav,data_all.ndata_grav,sizeof(GRAV_STRUCT),"Main");
				else
					grav =(GRAV_STRUCT *)memory((char *)grav,1,sizeof(GRAV_STRUCT),"Main");

			fprintf(out_rms,"Gravity:\n");
			fprintf(out_rms,"RMS-value[mgal] All measurements\n");
			fflush(out_rms);

			ForwardModGrav(geo_all, grid, &data_all, topo, inv);/*Jin,12_03*/
		
			/*Write out the binary file of the gravity data*/
			WriteDatGravOut(data_all,inv.ninv,geo_all);

			fprintf(out_rms,"%f %d\n",data_all.rms_grav, data_all.ndata_grav);
			fflush(out_rms);
		}

		/*----------------------- MT --------------------------------*/
		mt =(MT_STRUCT *)memory(NULL,1,sizeof(MT_STRUCT),"Main");
		calc_2D_mt = (CALC_2D_MT *)memory(NULL,1,sizeof(CALC_2D_MT),"Main");
		nr_of_slices_2D_mt = 0;

		if(flag.index_mt != 0)
		{
			printf("Start MT modeling\n---------------------\n");

				if(data_all.ndata_mt != 0)
					mt =(MT_STRUCT *)memory((char *)mt,data_all.ndata_mt,sizeof(MT_STRUCT),"Main");
				else
					mt =(MT_STRUCT *)memory((char *)mt,1,sizeof(MT_STRUCT),"Main");


			fprintf(out_rms,"MT:\n");
			fprintf(out_rms,"RMS-value of Re(Z)    Im(Z)  Nr.of data\n");
			fflush(out_rms);

			if(flag.dimension_mt == 1)
			/*1D modelling*/
				ForwardModMT_1D(geo_all, grid, &data_all, flag.kind_of_data_mt);
			else
			{
				/*Number of 2-D mt-slices*/
				if(flag.direc_2D_mt != 2) 
					nr_of_slices_2D_mt = (long) ceil((double)grid.ny/(double)(flag.nr_cells_2d_mt)); 
				else
					nr_of_slices_2D_mt = (long) ceil((double)grid.nx/(double)(flag.nr_cells_2d_mt));

				calc_2D_mt = (CALC_2D_MT *)memory((char *)calc_2D_mt,nr_of_slices_2D_mt,sizeof(CALC_2D_MT),"Main");

			/*2D modelling*/
				tmp_nr_of_slices_2D_mt = ForwardModMT_2D(geo_all, grid, &data_all, calc_2D_mt, flag, topo);

				if(tmp_nr_of_slices_2D_mt != nr_of_slices_2D_mt)
				{
					printf("WARNING!!! The number slices used for the 2D MT do not fit!!");
					nr_of_slices_2D_mt = tmp_nr_of_slices_2D_mt;
				}

			}

			/*Write out the binary file of the MT data*/
			WriteDatMTOut(data_all,inv.ninv,geo_all, flag.dimension_mt, flag.kind_of_data_mt);
			
			fprintf(out_rms,"           %10.7f %10.7f   %d\n",data_all.rms_mt_real, data_all.rms_mt_imag, data_all.ndata_freq_mt);
			fflush(out_rms);

		}

		/*-----------------------------------------------------------*/
/******************/
//HIER FORTSETZEN !!
/******************/

		/*Inversion*/
 		InvRoutine(geo_all, &grid, raypath, fat_rays, grav, mt, &inv, &data_all, &flag, &regu, calc_2D_mt, nr_of_slices_2D_mt, par.file_grid_size[0].name, par.file_seis_grav_rel[0].name, par.file_seis_res_rel[0].name);
		printf("\n\n\n\n\n");

		/*-----------------------------------------------------------*/
		/*Write out the binary file of the velocity model*/
		if(flag.index_tseis != 0 || (flag.index_grav != 0 && flag.index_mt != 0))
			WriteModSeisOut(inv.ninv, data_all, geo_all, grid);
		if(flag.index_grav != 0)
			WriteModGravOut(inv.ninv, data_all, geo_all, grid);
		if(flag.index_mt != 0)
			WriteModResOut(inv.ninv, data_all, geo_all, grid);
		/*-----------------------------------------------------------*/

		 /* Free memory*/
		if(flag.index_tseis != 0)
		{
			if(flag.kind_of_rays == 1)
			{	
				for(i=0;i<data_all.ndata_seis;i++)
				{
					free(raypath[i].len);
					free(raypath[i].ele);
					free(raypath[i].x);
					free(raypath[i].y);
					free(raypath[i].z);
				}
			}
			else
			{
				for(i=0;i<data_all.ndata_seis;i++)
				{
					free(fat_rays[i].weight);
					free(fat_rays[i].ele);
				}
			}
		}

		free(raypath);
		free(fat_rays);

		if(flag.index_grav != 0)
		{
			for(i=0;i<data_all.ndata_grav;i++)
			{
				free(grav[i].ele);
				free(grav[i].deriv);
			}
		}

		free(grav);

		if(flag.index_mt!= 0)
		{
			for(i=0;i<data_all.ndata_mt;i++)
			{
				for(j=0;j<mt[i].ncell;j++)
					free(mt[i].deriv[j]);

				free(mt[i].ele);
				free(mt[i].n);
				free(mt[i].deriv);
			}
		}

		free(mt);

		if(flag.index_mt != 0 && flag.dimension_mt != 1)
		{
			for(i=0;i<nr_of_slices_2D_mt;i++)
			{
				for(j=0;j<calc_2D_mt[i].nfreq;j++)
				{
					free(calc_2D_mt[i].freq_data[j]);
					free(calc_2D_mt[i].freq_stat[j]);
					free(calc_2D_mt[i].hx[j]);
					free(calc_2D_mt[i].hz[j]);
					free(calc_2D_mt[i].res_slice[j]);
					free(calc_2D_mt[i].index_slice[j]);

					free(calc_2D_mt[i].Ex_r[j]);
					free(calc_2D_mt[i].Ex_i[j]);
					free(calc_2D_mt[i].Ey_r[j]);
					free(calc_2D_mt[i].Ey_i[j]);
					free(calc_2D_mt[i].Hx_r[j]);
					free(calc_2D_mt[i].Hx_i[j]);
					free(calc_2D_mt[i].Hy_r[j]);
					free(calc_2D_mt[i].Hy_i[j]);

				}
		
				free(calc_2D_mt[i].freq);
				free(calc_2D_mt[i].freq_ndata);
				free(calc_2D_mt[i].freq_data);
				free(calc_2D_mt[i].freq_stat);
				free(calc_2D_mt[i].res_slice);
				free(calc_2D_mt[i].index_slice);
				free(calc_2D_mt[i].mts);
				free(calc_2D_mt[i].x);
				free(calc_2D_mt[i].z);
				free(calc_2D_mt[i].nx);
				free(calc_2D_mt[i].nz);
				free(calc_2D_mt[i].nz_atm);
				free(calc_2D_mt[i].nz_ionos);
				free(calc_2D_mt[i].nz_basement);
				free(calc_2D_mt[i].nx_left_border);
				free(calc_2D_mt[i].nx_right_border);
				free(calc_2D_mt[i].hx);
				free(calc_2D_mt[i].hz);
				free(calc_2D_mt[i].org[0]);
				free(calc_2D_mt[i].org[1]);
				free(calc_2D_mt[i].cell_scal_factor);

				free(calc_2D_mt[i].Ex_r);
				free(calc_2D_mt[i].Ex_i);
				free(calc_2D_mt[i].Ey_r);
				free(calc_2D_mt[i].Ey_i);
				free(calc_2D_mt[i].Hx_r);
				free(calc_2D_mt[i].Hx_i);
				free(calc_2D_mt[i].Hy_r);
				free(calc_2D_mt[i].Hy_i);

			}
		}

		free(calc_2D_mt);
	 
	/*HERE ENDS THE ITERATION-LOOP:*/
	}
							
	fclose(out_rms);
   /*-----------------------------------------------------------*/
   /* Free memory*/
		free(par.file_grid_size);
		free(par.file_seis_grav_rel);
		free(par.file_seis_res_rel);

		/*  for Grid structure*/			
		free(grid.slow);
		free(grid.dens);
		free(grid.res);
		free(grid.border_index);

		/* for structures GEOMETRY*/
		free(geo_all.x);
		free(geo_all.y);
		free(geo_all.z);
		free(geo_all.coor_info);

		/*for structure TOPO*/
		free(topo.x);
		free(topo.y);
		free(topo.z);
		free(topo.index);
   
		/* for structures DATA_STRUCT */
		/*Seismic*/
		free(data_all.sno);
		free(data_all.rno);
		free(data_all.tobs);
		free(data_all.tcalc);
		free(data_all.weigth_seis);
		free(data_all.xdist);
		free(data_all.ydist);
		free(data_all.zdist);
		free(data_all.shots);
		free(data_all.recs);
		free(data_all.lshots);
		free(data_all.lrecs);

		/*Gravimetry*/
		free(data_all.obs_grav);
		free(data_all.calc_grav);
		free(data_all.gno);
		free(data_all.weigth_grav);
		free(data_all.gravs);

		/*MT*/
		free(data_all.mno);
		free(data_all.nfreq_mt);
		free(data_all.mts);

		for(i=0;i<data_all.ndata_mt;i++)
		{
			free(data_all.freq_mt[i]);
			free(data_all.real_mt_TE[i]);
			free(data_all.imag_mt_TE[i]);
			free(data_all.calc_real_mt_TE[i]);
			free(data_all.calc_imag_mt_TE[i]);
			free(data_all.real_mt_TM[i]);
			free(data_all.imag_mt_TM[i]);
			free(data_all.calc_real_mt_TM[i]);
			free(data_all.calc_imag_mt_TM[i]);
			free(data_all.weigth_mt[i]);
		}
		free(data_all.freq_mt);
		free(data_all.real_mt_TE);
		free(data_all.imag_mt_TE);
		free(data_all.calc_real_mt_TE);
		free(data_all.calc_imag_mt_TE);
		free(data_all.real_mt_TM);
		free(data_all.imag_mt_TM);
		free(data_all.calc_real_mt_TM);
		free(data_all.calc_imag_mt_TM);
		free(data_all.weigth_mt);

		
		
   /*-----------------------------------------------------------*/
		/*finish memory check*/
//		_CrtDumpMemoryLeaks();

printf("----------------\n");
printf("program INV3D terminated without an error!\n");
printf("----------------\n\n");

return(1);
}