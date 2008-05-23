#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "inv3d.h"
#include "fkt.h"


/* Useful working tools: */

/*-------------------------------------------------------------------------*/
/* Read in comment lines from the parameter file; consider a line as a comment line, when it starts with a # */
/* Parameters: *line      := Line which has to be checked*/
/*             *inf       := Pointer on input file location*/

int GetNext(char *line,FILE *inf)
{
   line[0] = 0;
   do fgets(line,128,inf); while(line[0] == '#');
   return(0);
}

/*-------------------------------------------------------------------------*/
/* re-/allocate memory in a flexible way; */
/* Parameters: *prev_addr := Pointer on the field for which the memory is reallocated */
/*                           (if no memory is allocated before use the NULL pointer) */
/*             n          := Number of objects, for which memory should be allocated */
/*             size       := Size of objects                                         */
/*             *programme := Part of the program, where the allocation is performed  */
/*                                                                                   */
/* Output: new pointer on the field                                                  */    

char *memory (char *prev_addr,int n,int size,char *progname)
{
   char *tmp;

//   if (n == 0) return((char *) NULL); /* Take care of n = 0 */
	if (n == 0)
	{
		printf ("\n\n!!!WARNING!!!: In %s, no memory was allocated, n = %d\n\n", progname, n);
		return((char *) NULL); /* Take care of n = 0 */
	}
   
   if (prev_addr) 
   {
      if ((tmp = realloc ((char *) prev_addr, (unsigned) (n * size))) == NULL)
      {
	 printf ("Fatal Error: %s could not reallocate more memory, n = %d\n", progname, n);
	 exit (-1);
      }
   }
   else 
   {
      if ((tmp = calloc ((unsigned) n, (unsigned) size)) == NULL) 
      {
		printf ("Fatal Error: %s could not allocate memory, n = %d\n", progname, n);
		exit (-1);
      }
   }
   return (tmp);
}

/*---------------------------------------------------------------------------*/
/* read in the parameter file*/
/* parameters: *fname  := Pointer on the Inputfilename */
/*             *par    := Pointer on the Input-parameters*/
/*             *num_geo:= Pointer on the number of GEO-files */
/*			   *num_mod:= Pointer on the number of MOD-files */

int ReadInput(char *fname,PAR *par,int *num_geo, int *num_mod)
{
	int i,j,num_dat,zero_line,pos_string_geo,pos_string_mod, pos_string_dat;
	FILE *inf;
    char line[128];

   inf = fopen(fname,"rt");
   if (inf == NULL)
   {
      fprintf(stderr,"Unable to open %s\n",fname);
      exit(0);
   }

   /* Read in parameterfile */
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->binary_index);
   	/*binary MOD files*/
   pos_string_mod = ftell(inf);
   GetNext(line,inf);
   i=0;
   while(line[0] != '#' ) 
   {
	   fgets(line,128,inf);
	   i++;
   }
   *num_mod = i;
   par->binary_mod_dat = (FILENAME *)memory(NULL,*num_mod,sizeof(FILENAME),"ReadInput");
   fseek(inf,pos_string_mod,SEEK_SET);

   GetNext(line,inf);
   i=0;
   while(line[0] != '#')
   {	 
	   zero_line = (int)strlen(line);
	   for(j=0;!isgraph(line[j]);j++);
	   if(j-1 == zero_line)
	   {
		   printf("No filename exists; Format-error!");
		   exit(0);
	   }
	   memcpy((*par).binary_mod_dat[i].name, line+j, strlen(line+j)-1);
	   fgets(line,128,inf);
	   i++;
   }
   GetNext(line,inf);
   sscanf(line,"%d %d %d\n", &par->index_tseis, &par->index_grav, &par->index_mt);
   GetNext(line,inf);
   par->file_seis_grav_rel = (FILENAME *)memory(NULL,1,sizeof(FILENAME),"ReadInput");
   for(j=0;!isgraph(line[j]);j++);
     memcpy((*par).file_seis_grav_rel[0].name, line+j, strlen(line+j)-1);
   GetNext(line,inf);
   par->file_seis_res_rel = (FILENAME *)memory(NULL,1,sizeof(FILENAME),"ReadInput");
   for(j=0;!isgraph(line[j]);j++);
     memcpy((*par).file_seis_res_rel[0].name, line+j, strlen(line+j)-1);
   GetNext(line,inf);
   sscanf(line,"%d %d %d %f\n", &par->nx, &par->ny, &par->nz, &par->h);
   GetNext(line,inf);
   sscanf(line,"%f %f %f\n", &par->org[0], &par->org[1], &par->org[2]);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->topo_index);
   GetNext(line,inf);
   sscanf(line,"%d %d %d\n", &par->grid_cell_factor[0], &par->grid_cell_factor[1], &par->grid_cell_factor[2]);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->flex_grid_size);
   GetNext(line,inf);
   par->file_grid_size = (FILENAME *)memory(NULL,1,sizeof(FILENAME),"ReadInput");
   for(j=0;!isgraph(line[j]);j++);
     memcpy((*par).file_grid_size[0].name, line+j, strlen(line+j)-1);

   GetNext(line,inf);
   sscanf(line,"%d\n", &par->index_velocity_field);
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->vo);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->vmin, &par->vmax);
   GetNext(line,inf);
   sscanf(line,"%f %f %f\n", &par->vg[0],  &par->vg[1],   &par->vg[2]);
   GetNext(line,inf);
   sscanf(line,"%f %f %f\n", &par->org_vo[0],  &par->org_vo[1],   &par->org_vo[2]);
   GetNext(line,inf);
   sscanf(line,"%d %f\n", &par->nborder, &par->vborder);
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->v_air_water);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->ray_index);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->fatthres, &par->fatb);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->ray_density_index);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->ray_density_index2);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->vmina, &par->vmaxa);
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->maxadj);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->para_rays_index);
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->angle_para);
   
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->index_density_field);
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->go);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->gmin, &par->gmax);
   GetNext(line,inf);
   sscanf(line,"%f %f %f\n", &par->gg[0], &par->gg[1], &par->gg[2]);
   GetNext(line,inf);
   sscanf(line,"%f %f %f\n", &par->org_go[0], &par->org_go[1], &par->org_go[2]);
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->g_air_water);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->gmina, &par->gmaxa);

   /*extra cell*/
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->extra_cell);

   GetNext(line,inf);
   sscanf(line,"%d\n", &par->dimension_mt);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->direc_2D_mt);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->nr_cells_2d_mt);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->kind_of_data_mt);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->index_resistivity_field);
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->ro);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->rmin, &par->rmax);
   GetNext(line,inf);
   sscanf(line,"%f %f %f\n", &par->rg[0], &par->rg[1], &par->rg[2]);
   GetNext(line,inf);
   sscanf(line,"%f %f %f\n", &par->org_ro[0], &par->org_ro[1], &par->org_ro[2]);
   GetNext(line,inf);
   sscanf(line,"%f\n", &par->r_air_water);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->rmina, &par->rmaxa);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->index_mt_weighting);

   GetNext(line,inf);
   sscanf(line,"%f\n", &par->eps);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->damp, &par->smoothfac);
   GetNext(line,inf);
   sscanf(line,"%f %f\n", &par->smoothratio_y, &par->smoothratio_z);

   /*Jin Chen changed new*/
   GetNext(line,inf);             /*Read in the parameter which governs if the extra_smooth will be done from "XXX.inp" file.*/
   sscanf(line,"%d\n", &par->extra_smooth);
   GetNext(line,inf);             /**Read in the extra smooth parameter from "XXX.inp" file.*/
   sscanf(line,"%f\n", &par->val_extra_smooth);
   /*Jin Chen changed 02.2007*/
   GetNext(line,inf);             /**Read in the extra smooth distant parameter from "XXX.inp" file.*/
   sscanf(line,"%f\n", &par->dist_extra_smooth);

   GetNext(line,inf);
   sscanf(line,"%d\n", &par->write_sens_out);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->sens_in_perc);
   GetNext(line,inf);
   sscanf(line,"%d\n", &par->inv_max);

		/*Geo-files:*/
   pos_string_geo = ftell(inf);
   GetNext(line,inf);
   i=0;
   while(line[0] != '#' ) 
   {
	   fgets(line,128,inf);
	   i++;
   }
   *num_geo = i;
   par->files_geo = (FILENAME *)memory(NULL,*num_geo,sizeof(FILENAME),"ReadInput");
   fseek(inf,pos_string_geo,SEEK_SET);

   GetNext(line,inf);
   i=0;
   while(line[0] != '#')
   {	 
	   zero_line = (int)strlen(line);
	   for(j=0;!isgraph(line[j]);j++);
	   if(j-1 == zero_line)
	   {
		   printf("No filename exists; Format-error!");
		   exit(0);
	   }
	   memcpy((*par).files_geo[i].name, line+j, strlen(line+j)-1);
	   fgets(line,128,inf);
	   i++;
   }

   		/*Dat-files:*/
   pos_string_dat = ftell(inf);
   GetNext(line,inf);
   i=0;
   while(line[0] != '#' && !feof(inf))
   {
	   fgets(line,128,inf);
	   i++;
   }
   num_dat = i;
   par->files_dat = (FILENAME *)memory(NULL,num_dat,sizeof(FILENAME),"ReadInput");
   fseek(inf,pos_string_dat,SEEK_SET);

   GetNext(line,inf);
   i=0;
   while(line[0] != '#'  && !feof(inf))
   {	
	   zero_line = (int)strlen(line);
	   for(j=0;!isgraph(line[j]);j++);
	   if(j-1 == zero_line)
	   {
		   printf("No filename exists; Format-error!");
		   exit(0);
	   }
	   memcpy((*par).files_dat[i].name, line+j, strlen(line+j)-1);
	   fgets(line,128,inf);
	   i++;
   }

   fclose(inf);

   /*********************************************************************************************/
   /*Error messages*/

   if(*num_geo != num_dat)
   {   
	   printf("The number of GEO and DAT-Files is NOT equal\n");
	   exit(0);
   }

   if((par->binary_index ==1) && (*num_mod > 3))
   {
	   printf("The maximum number of binary MOD Files is 3\n");
	   exit(0);
   }

 
   if(par->index_tseis == 0 && par->index_grav == 0 && par->index_mt == 0)
   {
	   printf("NO inversion is choosen;\n check your index choice in the parameter file!! \n");
	   exit(0);
   }

   if(par->vmin < 0.0)
   {
	   printf("The minimum velocity HAVE TO be positive!!\n");
	   exit(0);
   }

   if(par->vmin > par->vmax)
   {
	   printf("The minimum velocity HAVE TO be smaller than the maximum velocity!!\n");
	   exit(0);
   }


   if(par->vborder < 0.0)
   {
	   printf("The minimum velocity of the border HAVE TO be positive!!\n");
	   exit(0);
   }

   if(par->v_air_water <0.0)
   {
	   printf("The minimum velocity of the air/water HAVE TO be positive!!\n");
	   exit(0);
   }


   if(par->gmin < 0.0)
   {
	   printf("The minimum density HAVE TO be positive!!\n");
	   exit(0);
   }

   if(par->gmin > par->gmax)
   {
	   printf("The minimum density HAVE TO be smaller than the maximum density!!\n");
	   exit(0);
   }

   if(par->g_air_water <0.0)
   {
	   printf("The minimum density of the air/water HAVE TO be positive or zero!!\n");
	   exit(0);
   }


   if(par->nr_cells_2d_mt < 1)
   {
	   printf("The thickness of the 2-D MT layers HAS TO larger equal 1 grid cell!!\n");
	   exit(0);
   }

   if(par->rmin < 0.0)
   {
	   printf("The minimum elec.resistivity HAVE TO be positive!!\n");
	   exit(0);
   }

   if(par->rmin > par->rmax)
   {
	   printf("The minimum elec.resistivity HAVE TO be smaller than the maximum elec.resistivity!!\n");
	   exit(0);
   }

   if(par->r_air_water <0.0)
   {
	   printf("The minimum elec.resistivity of the air/water HAVE TO be positive or zero!!\n");
	   exit(0);
   }


   for(i=0;i<3;i++)
   {
		if(par->grid_cell_factor[i] <= 0)
		{
			printf("The grid cell dimensions for the forward modeling HAVE TO integer value\n of the grid cell dimensions for the inversion!\n");
			exit(0);
		}
   }

   if(par->fatb < 0)
   {
	   printf("The exponent for the fat-rays have to have values that are >=0\n");
	   exit(0);
   }

   if(par->smoothfac < 0.0 || par->smoothfac >1.0)
   {
	   printf("The factor governing the relationship between damping and smoothing \nhave to be between 0 and 1 !!\n");
	   printf("Actual value:%f\n", par->smoothfac);
	   exit(0);
   }

   if(par->damp < 0.0)
   {
	   printf("The damping factor have to be >=0 !!\n");
	   printf("Actual value:%f\n", par->damp);
	   exit(0);
   }

   if(par->inv_max < 0)
   {
	   printf("The max. number of iteration %d makes NO sense\n", par->inv_max);
	   exit(0);
   }

   if(par->vmina > par->vmaxa)
   {
	   printf("The lower and upper velocity boundaries are not set correctly!!\n");
	   printf("v_min: %f m/s\n",par->vmina);
	   printf("v_max: %f m/s\n",par->vmaxa);
	   exit(0);
   }

   if(par->gmina > par->gmaxa)
   {
	   printf("The lower and upper density boundaries are not set correctly!!\n");
	   printf("g_min: %f g/cm^3\n",par->gmina);
	   printf("g_max: %f g/cm^3\n",par->gmaxa);
	   exit(0);
   }

   if(par->rmina > par->rmaxa)
   {
	   printf("The lower and upper resistivity boundaries are not set correctly!!\n");
	   printf("r_min: %f ohmm\n",par->rmina);
	   printf("r_max: %f ohmm\n",par->rmaxa);
	   exit(0);
   }

   
	if(par->angle_para < 0 || par->angle_para > 180)
	{
		printf("The angle to identify parallel rays in the seismic tomography\nHAVE to be between 0 and 180 degree\n");
		printf("Used angle: %f degree\n", par->angle_para);
		exit(0);
	}

/*********************************************************************************************/

   if(par->binary_index ==1)
   {
	   par->index_tseis = 1;
	   par->index_grav = 1;
	   par->index_mt = 1;
   }

   /*Print parameters:*/
   printf("Parameter-file:\n\n");

   if(par->binary_index !=1)
   {
		printf("Kind of inversion:\n");
		if(par->index_tseis != 0 && par->index_grav == 0 && par->index_mt == 0)
		   printf("ONLY inversion of the seismic traveltime picks\n\n");
		else if(par->index_tseis == 0 && par->index_grav != 0 && par->index_mt == 0)
		   printf("ONLY inversion of gravity data\n\n");
		else if(par->index_tseis == 0 && par->index_grav == 0 && par->index_mt != 0)
		   printf("ONLY inversion of the MT data\n\n");
		else if(par->index_tseis != 0 && par->index_grav != 0 && par->index_mt == 0)
		   printf("Seismic traveltimes and gravity\n\n");
		else if(par->index_tseis != 0 && par->index_grav == 0 && par->index_mt != 0)
		   printf("Seismic traveltimes and MT\n\n");
		else if(par->index_tseis == 0 && par->index_grav != 0 && par->index_mt != 0)
		   printf("Gravity and MT data\n\n");
		else
		   printf("Seismic traveltimes, gravity and MT\n\n");
   }

   if(par->binary_index ==1)
   {
	   printf("Binary model input:\n");

	   for(i=0;i<*num_mod;i++)
	   {
			printf("Filename of the %d binary model file: %s\n",i+1, par->binary_mod_dat[i].name);
	   }

	   printf("\nREMARK: The kind of considered methods will be specified\nlateron by the binary modell\n");
   }

   printf("\n------------------- Grid parameters: -------------------\n\n");
   
   if(par->binary_index !=1)
   {
		printf("nx is %6d\nny is %6d\nnz is %6d\nh is %8.3fm\n\n", par->nx, par->ny, par->nz, par->h);
		printf("Origin of the grid:\n");
		printf("x0 is %9.3fm\ny0 is %9.3fm\nz0 is %9.3fm\n\n", par->org[0], par->org[1], par->org[2]);

		if(par->topo_index == 0)
			printf("The topography information will be ignored\n\n");
		else
			printf("The topography will be incorporated by means of triangulation\n\n");
   }

   printf("\nThe grid cell dimensions for the forward modeling\ncompared to the grid cell dimensions for the inversion: x=%d, y=%d, z=%d\n\n", par->grid_cell_factor[0],par->grid_cell_factor[1],par->grid_cell_factor[2]);
   if(par->flex_grid_size != 1)   
	   printf("The grid cell size of the inversion cells is spatially invariant\n\n");
   else
   {
	   printf("The grid cell size of the inversion cells is dependent on their location.\n");
	   printf(" Parameter-file, where spatially varying sizes of the inversion cells\n are specified: %s\n\n",par->file_grid_size[0].name);
   }


   if(par->index_tseis != 0)
   {
	   printf("------------------- Seismic parameters: -------------------\n\n");

		if(par->binary_index  != 1)
		{
			

			if(par->index_velocity_field != 1) 
				printf("Velocity at the origin of the gradient velocity model:\n");
			else 
				printf("Velocity at the surface:\n");
			printf("velocity is %8.2fm/s:\n\n",par->vo);
			printf("The lowest velocity is %8.2fm/s\nThe highest velocity is %8.2fm/s\n\n",par->vmin,par->vmax);
			if(par->index_velocity_field != 1)
			{
				printf("Velocity-gradient:\n");
				printf("dv/dr is (%8.2f)1/s\n         (%8.2f)1/s\n         (%8.2f)1/s:\n\n",par->vg[0],par->vg[1],par->vg[2]);
			}
			else
			{
				printf("Abs.-value of the velocity-gradient:\n");
				printf("|dv/dr| is (%8.2f)1/s\n\n",(float)sqrt((par->vg[0])*(par->vg[0])+(par->vg[1])*(par->vg[1])+(par->vg[2])*(par->vg[2])));
			}

			if(par->index_velocity_field != 1) 
			{
				printf("The origin of the gradient velocity model:\n");
				printf("x is %9.3fm\ny is %9.3fm\nz is %9.3fm\n\n", par->org_vo[0], par->org_vo[1], par->org_vo[2]);
			}

			printf("Properties of the boundaries:\n");
			printf("Number of cells is %3d\nvelocity is %8.2fm/s\n\n",par->nborder, par->vborder);

			printf("Velocity of the air/water: %8.2fm/s\n\n", par->v_air_water);
		}
   
		if(par->ray_index == 1)
			printf("Normal rays are used!\n\n");
		else
		{
			printf("Fat rays are used!\n");
			printf("Threshold:%5.2f ms\nDecay:    %5.2f 1/ms\n\n",par->fatthres,par->fatb);
		}

		if(par->vmina <= 0.0)
			printf("The allowed minimum velocity is: 0.0 m/s\n");
		else
			printf("The allowed minimum velocity is: %8.2f m/s\n",par->vmina);

		if(par->vmaxa <= 0.0)
			printf("The allowed maximum velocity is infinity\n\n");
		else
			printf("The allowed maximum velocity is %8.2f m/s\n\n",par->vmaxa);

		if(par->maxadj <= 0.0)
			printf("The velocity changes per iteration are NOT constrained\n\n");
		else
			printf("The max. velocity change per iteration is %8.2f m/s\n\n",par->maxadj);

		if(par->ray_density_index == 1)
			printf("The hit matrix and the DWS (derivative weighted sum) will be calculated\nand written out in a binary file\n\n");
		if(par->ray_density_index2 == 1)
			printf("The ray density tensor will be calculated and written out in a binary file\n\n");

		if(par->ray_index == 1 && par->para_rays_index != 0)
		{
			printf("Parallel rays will be downscaled !! ");
			printf("Rays that have an angle of <%5.2f \nbetween each other will be considered as 'parallel'\n\n", par->angle_para);
		}

   }

   if(par->index_grav != 0)
   {

	   printf("------------------- Gravity parameters: -------------------\n\n");

	   if(par->binary_index !=1)
	   {
			if(par->index_tseis != 0 || par->index_mt != 0)
				printf("The starting model will be determined by a link with the velocity \n model by using the parameter in the file %s\n\n", par->file_seis_grav_rel);
			else
			{
		   		if(par->index_density_field != 1) 
					printf("Density at the origin of the gradient density model:\n");
				else 
					printf("Density at the surface:\n");
				printf("density is %8.2fg/cm^3:\n\n",par->go);
				printf("The lowest density is  %8.2fg/cm^3\nThe highest density is %8.2fg/cm^3\n\n",par->gmin,par->gmax);
				if(par->index_density_field != 1)
				{
					printf("Density-gradient:\n");
					printf("dg/dr is (%8.2f)g/cm^3/m\n         (%8.2f)g/cm^3/m\n         (%8.2f)g/cm^3/m:\n\n",par->gg[0],par->gg[1],par->gg[2]);
				}
				else
				{
					printf("Abs.-value of the density-gradient:\n");
					printf("|dg/dr| is (%8.2f)g/cm^3/m\n\n",(float)sqrt((par->gg[0])*(par->gg[0])+(par->gg[1])*(par->gg[1])+(par->gg[2])*(par->gg[2])));
				}

				if(par->index_density_field != 1) 
				{
					printf("The origin of the gradient density model:\n");
					printf("x is %9.3fm\ny is %9.3fm\nz is %9.3fm\n\n", par->org_go[0], par->org_go[1], par->org_go[2]);
				}
			}
	   }

	   printf("Density of the air/water: %8.2fg/cm^3\n\n", par->g_air_water);

		if(par->gmina <= 0.0)
			printf("The allowed minimum density is: 0.0 g/cm^3\n");
		else
			printf("The allowed minimum density is: %8.2f g/cm^3\n",par->gmina);

		if(par->gmaxa <= 0.0)
			printf("The allowed maximum density is: infinity\n\n");
		else
			printf("The allowed maximum density is: %8.2f g/cm^3\n\n",par->gmaxa);
			
		/*extra cell (BJOERN_MOD)*/
   		if(par->extra_cell > 0.0)
   		{
	   		printf("The extra cell under the model will be involved, and the weighting is: %f\n",par->extra_cell);
   		}

   }


   if(par->index_mt != 0)
   {
	   printf("---------------------- MT parameters: -------------------\n\n");

	   if(par->dimension_mt == 1)
		   printf("The MT modelling is 1-D\n\n");
	   else
	   {
		   if(par->dimension_mt == 2)
				printf("The MT modelling is 2-D\n\n");
		   else
		   {
			   if(par->dimension_mt == 3)
				   printf("MT forward modelling is 2-D and inversion 1-D\n\n");
			   else
			   {
				   printf("MT forward modelling is 2-D and inversion 1-D\n");
				   printf("The RRI-method is used for the inversion\n\n");
			   }
				
		   }

           if(par->direc_2D_mt !=2)
			   printf("The 2-D modelling will be performed along the x-direction\n");
		   else
			   printf("The 2-D modelling will be performed along the y-direction\n");
		   printf("Width of one 2-D layer:    %d cells\n\n",par->nr_cells_2d_mt);
	   }

	   if(par->kind_of_data_mt == 1)
		   printf("Data from the TE-Mode will be used\n\n");
	   else if(par->kind_of_data_mt == 2)
		   printf("Data from the TM-Mode will be used\n\n");
	   else if(par->dimension_mt == 1)
		   printf("Berdichewsky average will be used for the 1-D calculation\n\n");
	   else
		   printf("Data both from the TE- and TM-Mode will be used\n\n");




	   if(par->binary_index !=1)
	   {
		   	if(par->index_tseis != 0 || par->index_grav != 0)
				printf("The starting model will be determined by a link with the velocity \n model by using the parameter in the file %s\n\n", par->file_seis_grav_rel);
			else
			{
		   		if(par->index_resistivity_field != 1) 
					printf("Elec.resistivity at the origin of the gradient model:\n");
				else 
					printf("Elec.resistivity at the surface:\n");
				printf("Elec.resistivity is %8.2fohmm:\n\n",par->ro);
				printf("The lowest elec.resistivity is  %8.2fohmm\nThe highest elec.resistivity is %8.2fohmm\n\n",par->rmin,par->rmax);
				if(par->index_resistivity_field != 1)
				{
					printf("Elec.resistivity-gradient:\n");
					printf("dR/dr is (%8.2f)ohmm/m\n         (%8.2f)ohmm/m\n         (%8.2f)ohmm/m:\n\n",par->rg[0],par->rg[1],par->rg[2]);
				}
				else
				{
					printf("Abs.-value of the elec.resistivity-gradient:\n");
					printf("|dR/dr| is (%8.2f)ohmm/m\n\n",(float)sqrt((par->rg[0])*(par->rg[0])+(par->rg[1])*(par->rg[1])+(par->rg[2])*(par->rg[2])));
				}
	
				if(par->index_resistivity_field != 1) 
				{
					printf("The origin of the gradient elec.resistivity model:\n");
					printf("x is %9.3fm\ny is %9.3fm\nz is %9.3fm\n\n", par->org_ro[0], par->org_ro[1], par->org_ro[2]);
				}
			}
	   }

	    printf("Elec.resistivity of the air/water: %8.2fohmm\n\n", par->r_air_water);

		if(par->rmina <= 0.0)
			printf("The allowed minimum resistivity is: 0.0 ohmm\n");
		else
			printf("The allowed minimum resistivity is: %8.2f ohmm\n",par->rmina);

		if(par->rmaxa <= 0.0)
			printf("The allowed maximum resistivity is infinity\n\n");
		else
			printf("The allowed maximum resistivity is: %8.2f ohmm\n\n",par->rmaxa);

		if(par->index_mt_weighting != 0)
			printf("The MT impedances will weighted individually for different frequencies\nin the inversion\n\n");

   }


   printf("------------------- Other parameters: -------------------\n\n");

//(BJOERN_MOD:start changes)
   
   printf("Maximum coordinate distance Eps for the same shot/receiver location: %6.3fm\n", fabs(par->eps));
   if(fabs(par->eps) > 5.0) 
	printf("!!Attention!! Eps is very large\n");
   printf("\n");

   printf("Regularisation:\n");
   printf("Damping parameter:     %4.2f\n",par->damp);
   printf("Rel. damping/smooting: %4.2f\n",par->smoothfac);
   if(par->smoothfac == 0.0)
   	printf("NO damping !!\n");
   else if(par->smoothfac == 1.0)
	printf("NO smoothing !!\n\n");

   if(par->smoothfac != 1.0)
   {
	printf("Smoothing ratios: y/x = %4.2f\n",par->smoothratio_y);
	printf("                  z/x = %4.2f\n\n",par->smoothratio_z);
   }
		
   /*Jin Chen changed new*/
   if(par->extra_smooth == 0)
	printf("NO extra smoothing !!\n");
   else
   {
	printf("Extra smoothing will be done in inversion !!\n");
	printf("Extra smoothing: %4.2f\n",par->val_extra_smooth);
	/*Jin Chen changed 02.2007*/
	printf("Inner cells, which is not far away than %8.2f\n m from the border cell, will be considered in extra smooth.",par->dist_extra_smooth);
   }

//(BJOERN_MOD:end changes)
    if(par->write_sens_out == 1)
		printf("The sensitivities for the different methods are written out in binary files!\n\n");
	if(par->sens_in_perc == 1)
		printf("The sensitivities are expressed in percentages (100*|(dD/dm)/(D/m)|)\n\n");

	printf("GEO-inputfiles:\n");
	for(i=0;i<*num_geo;i++)
		printf("%s\n",par->files_geo[i].name);
	printf("DAT-inputfiles:\n");
	for(i=0;i<num_dat;i++)
		printf("%s\n",par->files_dat[i].name);

   printf("\nMax. number of iteration: %d\n\n",par->inv_max);

   
   printf("----------------\n\n");


   return(1);
}


/*---------------------------------------------------------------------------*/
/* Write out the binary file of seismic velocity model*/
/* parameters:   data       := Data structure*/ 
/*				 geo        := Geometry structure */
/*               grid       := Grid structure of the model used for INVERSION*/
/*               ninv       := Number of inversions  */  


int WriteModSeisOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid)
{
	int index_seis;
	long nx2,ny2,nz2,nyz2,nborder;
	long ny3,nz3,nyz3;
	long a,b,c;
	float inv_help;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"modseis%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	nborder = grid.nborder;

	index_seis = 1; /*Specifying that this is a seismic file*/

	nx2 = grid.nx + 2*nborder;
	ny2 = grid.ny + 2*nborder;
	nz2 = grid.nz + 2*nborder;
	nyz2=ny2*nz2;

	ny3 = ny2+1;
	nz3 = nz2+1;
	nyz3 = ny3*nz3;

	/*Write out the binary header*/
	fwrite(&index_seis,sizeof(int),1,ouf);							/*Specifying that this is a seismic file*/

	fwrite(&grid.nx,sizeof(int),1,ouf);								/*Number of grid cells in x-direction*/
	fwrite(&grid.ny,sizeof(int),1,ouf);								/*Number of grid cells in y-direction*/
	fwrite(&grid.nz,sizeof(int),1,ouf);								/*Number of grid cells in z-direction*/
	fwrite(&grid.h,sizeof(float),1,ouf);							/*Size of grid cell*/
	fwrite(&grid.nborder,sizeof(int),1,ouf);						/*Number of grid cells of the boundaries*/

	for(a=0;a<3;a++) 
		fwrite(&grid.org[a],sizeof(float),1,ouf);					/*Origin(x,y,z) of grid*/

	fwrite(&geo.nstat,sizeof(int),1,ouf);							/*Number of active stations*/

	b=0;
	for(a=0;a<(geo.nshot);a++) if(data.lshots[a] >= 0) b++;
	fwrite(&b,sizeof(int),1,ouf);									/*Number of active shots*/

	b=0;
	for(a=0;a<(geo.nrec);a++) if(data.lrecs[a] >= 0) b++;
	fwrite(&b,sizeof(int),1,ouf);									/*Number of active receivers*/

	for(a=0;a<geo.nshot;a++)
	{
		if(data.lshots[a] >= 0)
		{
			fwrite(&(data.shots[a]),sizeof(int),1,ouf);				/*Location number of the shots*/
		}

	}

	for(a=0;a<geo.nrec;a++)
	{
		if(data.lrecs[a] >= 0)
		{
			fwrite(&(data.recs[a]),sizeof(int),1,ouf);				/*Location number of the receivers*/
		}

	}

	for(a=0;a<geo.nshot;a++)
	{
		/*Memory check*/
		if(data.shots[a]-1 >= geo.nstat)
		{
			printf("NOT enough memory is allocated: used %d, allocated %d\n",data.shots[a]-1, geo.nstat);
			exit(0);
		}

		if(data.lshots[a] >= 0)
		{
			fwrite(&(geo.x[data.shots[a]-1]),sizeof(float),1,ouf);		/*x-coordinates of shot positions*/
		}

	}

	for(a=0;a<geo.nshot;a++)
	{
		if(data.lshots[a] >= 0)
		{
			fwrite(&(geo.y[data.shots[a]-1]),sizeof(float),1,ouf);		/*y-coordinates of shot positions*/	
		}

	}


	for(a=0;a<geo.nshot;a++)
	{
		if(data.lshots[a] >= 0)
		{
			fwrite(&(geo.z[data.shots[a]-1]),sizeof(float),1,ouf);		/*z-coordinates of shot positions*/
		}

	}


	for(a=0;a<geo.nrec;a++)
	{
		/*Memory check*/
		if(data.recs[a]-1 >= geo.nstat)
		{
			printf("NOT enough memory is allocated: used %d, allocated %d\n",data.recs[a]-1, geo.nstat);
			exit(0);
		}

		if(data.lrecs[a] >= 0)
		{
			fwrite(&(geo.x[data.recs[a]-1]),sizeof(float),1,ouf);		/*x-coordinates of receiver positions*/
		}

	}

	for(a=0;a<geo.nrec;a++)
	{
		if(data.lrecs[a] >= 0)
		{
			fwrite(&(geo.y[data.recs[a]-1]),sizeof(float),1,ouf);		/*y-coordinates of receiver positions*/	
		}

	}

	
	for(a=0;a<geo.nrec;a++)
	{
		if(data.lrecs[a] >= 0)
		{
			fwrite(&(geo.z[data.recs[a]-1]),sizeof(float),1,ouf);		/*z-coordinates of receiver positions*/
		}

	}


	for(a=0;a<nx2;a++)
		for(b=0;b<ny2;b++)
			for(c=0;c<nz2;c++)
			{
				inv_help = grid.h/(float)grid.slow[nyz3*a + nz3*b +c];
				fwrite(&inv_help,sizeof(float),1,ouf);					/*3-D velocity model*/
			}

	for(a=0;a<nx2;a++)
		for(b=0;b<ny2;b++)
			for(c=0;c<nz2;c++)
			{
				fwrite(&grid.border_index[nyz2*a + nz2*b +c],sizeof(int),1,ouf);	/*3-D "boundary model*/
			}


	fclose(ouf);

	printf("----------------\n");
	printf("The binary velocity model file modseis%03d.dat is written out\n",ninv);
	printf("----------------\n\n");

	return(1);
}


/*---------------------------------------------------------------------------*/
/* Write out the binary file of density model*/
/* parameters:   data       := Data structure*/ 
/*				 geo        := Geometry structure */
/*               grid       := Grid structure of the model used for INVERSION*/
/*               ninv       := Number of inversions   */  


int WriteModGravOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid)
{
	int index_grav;
	long ny2,nz2,nyz2,nborder;
	long a,b,c;
	float inv_help;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"modgrav%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	nborder = 0; /*!!Attention!! Is set to zero for density data*/

	index_grav = 2; /*Specifying that this is a gravity file*/

	ny2 = grid.ny + 2*grid.nborder;
	nz2 = grid.nz + 2*grid.nborder;
	nyz2=ny2*nz2;

	/*Write out the binary header*/
	fwrite(&index_grav,sizeof(int),1,ouf);							/*Specifying that this is a gravity file*/

	fwrite(&grid.nx,sizeof(int),1,ouf);								/*Number of grid cells in x-direction*/
	fwrite(&grid.ny,sizeof(int),1,ouf);								/*Number of grid cells in y-direction*/
	fwrite(&grid.nz,sizeof(int),1,ouf);								/*Number of grid cells in z-direction*/
	fwrite(&grid.h,sizeof(float),1,ouf);							/*Size of grid cell*/
	fwrite(&nborder,sizeof(int),1,ouf);								/*Number of grid cells of the boundaries is set to 0*/

	for(a=0;a<3;a++) 
		fwrite(&grid.org[a],sizeof(float),1,ouf);					/*Origin(x,y,z) of grid*/

	fwrite(&geo.nstat,sizeof(int),1,ouf);							/*Number of active stations*/

	fwrite(&geo.nstat_grav,sizeof(int),1,ouf);						/*Number of gravity stations*/

	b=0;
	fwrite(&b,sizeof(int),1,ouf);									/*NOT USED*/

	for(a=0;a<geo.nstat_grav;a++)
			fwrite(&(data.gravs[a]),sizeof(int),1,ouf);				/*Location number of the gravity stations*/


	for(a=0;a<geo.nstat_grav;a++)
	{
		/*Memory check*/
		if(data.gravs[a]-1 >= geo.nstat)
		{
			printf("NOT enough memory is allocated: used %d, allocated %d\n",data.gravs[a]-1, geo.nstat);
			exit(0);
		}

		fwrite(&(geo.x[data.gravs[a]-1]),sizeof(float),1,ouf);		/*x-coordinates of gravity station*/
	}

	for(a=0;a<geo.nstat_grav;a++)
		fwrite(&(geo.y[data.gravs[a]-1]),sizeof(float),1,ouf);		/*y-coordinates of gravity station*/	


	for(a=0;a<geo.nstat_grav;a++)
		fwrite(&(geo.z[data.gravs[a]-1]),sizeof(float),1,ouf);		/*z-coordinates of gravity station*/

	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				inv_help = (float)grid.dens[grid.ny*grid.nz*a + grid.nz*b + c];
				fwrite(&inv_help,sizeof(float),1,ouf);					/*3-D density model*/
			}

	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				fwrite(&grid.border_index[nyz2*(a + grid.nborder) + nz2*(b+grid.nborder) + (c+grid.nborder)],sizeof(int),1,ouf);	/*3-D "boundary model*/
			}																		/*Attention!! Size of the grid differs from the seismic model*/						


	fclose(ouf);

	printf("----------------\n");
	printf("The binary density model file modgrav%03d.dat is written out\n",ninv);
	printf("----------------\n\n");

	return(1);
}


/*---------------------------------------------------------------------------*/
/* Write out the binary file of resistivity model*/
/* parameters:   data       := Data structure*/ 
/*				 geo        := Geometry structure */
/*               grid       := Grid structure of the model used for INVERSION*/
/*               ninv       := Number of inversions   */  


int WriteModResOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid)
{
	int index_res;
	long ny2,nz2,nyz2,nborder;
	long a,b,c;
	float inv_help;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"modres%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	nborder = 0; /*!!Attention!! Is set to zero for resistivity data*/

	index_res = 3; /*Specifying that this is a MT file*/

	ny2 = grid.ny + 2*grid.nborder;
	nz2 = grid.nz + 2*grid.nborder;
	nyz2=ny2*nz2;

	/*Write out the binary header*/
	fwrite(&index_res,sizeof(int),1,ouf);							/*Specifying that this is a MT file*/

	fwrite(&grid.nx,sizeof(int),1,ouf);								/*Number of grid cells in x-direction*/
	fwrite(&grid.ny,sizeof(int),1,ouf);								/*Number of grid cells in y-direction*/
	fwrite(&grid.nz,sizeof(int),1,ouf);								/*Number of grid cells in z-direction*/
	fwrite(&grid.h,sizeof(float),1,ouf);							/*Size of grid cell*/
	fwrite(&nborder,sizeof(int),1,ouf);								/*Number of grid cells of the boundaries is set to 0*/

	for(a=0;a<3;a++) 
		fwrite(&grid.org[a],sizeof(float),1,ouf);					/*Origin(x,y,z) of grid*/

	fwrite(&geo.nstat,sizeof(int),1,ouf);							/*Number of active stations*/

	fwrite(&geo.nstat_mt,sizeof(int),1,ouf);						/*Number of MT stations*/

	b=0;
	fwrite(&b,sizeof(int),1,ouf);									/*NOT USED*/

	for(a=0;a<geo.nstat_mt;a++)
			fwrite(&(data.mts[a]),sizeof(int),1,ouf);				/*Location number of the MT stations*/


	for(a=0;a<geo.nstat_mt;a++)
	{
		/*Memory check*/
		if(data.mts[a]-1 >= geo.nstat)
		{
			printf("NOT enough memory is allocated: used %d, allocated %d\n",data.mts[a]-1, geo.nstat);
			exit(0);
		}

		fwrite(&(geo.x[data.mts[a]-1]),sizeof(float),1,ouf);		/*x-coordinates of MT station*/
	}

	for(a=0;a<geo.nstat_mt;a++)
		fwrite(&(geo.y[data.mts[a]-1]),sizeof(float),1,ouf);		/*y-coordinates of MT station*/	


	for(a=0;a<geo.nstat_mt;a++)
		fwrite(&(geo.z[data.mts[a]-1]),sizeof(float),1,ouf);		/*z-coordinates of MT station*/

	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				inv_help = (float)grid.res[grid.ny*grid.nz*a + grid.nz*b + c];
				fwrite(&inv_help,sizeof(float),1,ouf);					/*3-D ele. resistivity model*/
			}

	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				fwrite(&grid.border_index[nyz2*(a + grid.nborder) + nz2*(b+grid.nborder) + (c+grid.nborder)],sizeof(int),1,ouf);	/*3-D "boundary model*/
			}																		/*Attention!! Size of the grid differs from the resistivity model*/						


	fclose(ouf);

	printf("----------------\n");
	printf("The binary resistivity model file modres%03d.dat is written out\n",ninv);
	printf("----------------\n\n");

	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out the binary file of the total seismic sensitivity*/
/* parameters:   data       := Data structure*/ 
/*				 geo        := Geometry structure */
/*               grid       := Grid structure of the model used for INVERSION*/
/*				 sens		:= Total sensitivities for the cells*/
/*               ninv       := Number of inversions   */
/*Remark: The format is the same like for the model files*/

int WriteSensSeisOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid, double *sens)
{
	int index_seis;
	long nx2,ny2,nz2,nyz2,nborder;
	long a,b,c;
	float inv_help;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"senseis%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	nborder = grid.nborder;

	index_seis = 1; /*Specifying that this is a seismic file*/

	nx2 = grid.nx + 2*nborder;
	ny2 = grid.ny + 2*nborder;
	nz2 = grid.nz + 2*nborder;
	nyz2=ny2*nz2;

	/*Write out the binary header*/
	fwrite(&index_seis,sizeof(int),1,ouf);							/*Specifying that this is a seismic file*/

	fwrite(&grid.nx,sizeof(int),1,ouf);								/*Number of grid cells in x-direction*/
	fwrite(&grid.ny,sizeof(int),1,ouf);								/*Number of grid cells in y-direction*/
	fwrite(&grid.nz,sizeof(int),1,ouf);								/*Number of grid cells in z-direction*/
	fwrite(&grid.h,sizeof(float),1,ouf);							/*Size of grid cell*/
	fwrite(&grid.nborder,sizeof(int),1,ouf);						/*Number of grid cells of the boundaries*/

	for(a=0;a<3;a++) 
		fwrite(&grid.org[a],sizeof(float),1,ouf);					/*Origin(x,y,z) of grid*/

	fwrite(&geo.nstat,sizeof(int),1,ouf);							/*Number of active stations*/

	b=0;
	for(a=0;a<(geo.nshot);a++) if(data.lshots[a] >= 0) b++;
	fwrite(&b,sizeof(int),1,ouf);									/*Number of active shots*/

	b=0;
	for(a=0;a<(geo.nrec);a++) if(data.lrecs[a] >= 0) b++;
	fwrite(&b,sizeof(int),1,ouf);									/*Number of active receivers*/

	for(a=0;a<geo.nshot;a++)
	{
		if(data.lshots[a] >= 0)
		{
			fwrite(&(data.shots[a]),sizeof(int),1,ouf);				/*Location number of the shots*/
		}

	}

	for(a=0;a<geo.nrec;a++)
	{
		if(data.lrecs[a] >= 0)
		{
			fwrite(&(data.recs[a]),sizeof(int),1,ouf);				/*Location number of the receivers*/
		}

	}

	for(a=0;a<geo.nshot;a++)
	{
		/*Memory check*/
		if(data.shots[a]-1 >= geo.nstat)
		{
			printf("NOT enough memory is allocated: used %d, allocated %d\n",data.shots[a]-1, geo.nstat);
			exit(0);
		}

		if(data.lshots[a] >= 0)
		{
			fwrite(&(geo.x[data.shots[a]-1]),sizeof(float),1,ouf);		/*x-coordinates of shot positions*/
		}

	}

	for(a=0;a<geo.nshot;a++)
	{
		if(data.lshots[a] >= 0)
		{
			fwrite(&(geo.y[data.shots[a]-1]),sizeof(float),1,ouf);		/*y-coordinates of shot positions*/	
		}

	}


	for(a=0;a<geo.nshot;a++)
	{
		if(data.lshots[a] >= 0)
		{
			fwrite(&(geo.z[data.shots[a]-1]),sizeof(float),1,ouf);		/*z-coordinates of shot positions*/
		}

	}


	for(a=0;a<geo.nrec;a++)
	{
		/*Memory check*/
		if(data.recs[a]-1 >= geo.nstat)
		{
			printf("NOT enough memory is allocated: used %d, allocated %d\n",data.recs[a]-1, geo.nstat);
			exit(0);
		}

		if(data.lrecs[a] >= 0)
		{
			fwrite(&(geo.x[data.recs[a]-1]),sizeof(float),1,ouf);		/*x-coordinates of receiver positions*/
		}

	}

	for(a=0;a<geo.nrec;a++)
	{
		if(data.lrecs[a] >= 0)
		{
			fwrite(&(geo.y[data.recs[a]-1]),sizeof(float),1,ouf);		/*y-coordinates of receiver positions*/	
		}

	}

	
	for(a=0;a<geo.nrec;a++)
	{
		if(data.lrecs[a] >= 0)
		{
			fwrite(&(geo.z[data.recs[a]-1]),sizeof(float),1,ouf);		/*z-coordinates of receiver positions*/
		}

	}


	for(a=0;a<nx2;a++)
		for(b=0;b<ny2;b++)
			for(c=0;c<nz2;c++)
			{
				inv_help = (float)sens[nyz2*a + nz2*b +c];
				fwrite(&inv_help,sizeof(float),1,ouf);					/*seismic sensitivities in 3-D*/
			}

	for(a=0;a<nx2;a++)
		for(b=0;b<ny2;b++)
			for(c=0;c<nz2;c++)
			{
				fwrite(&grid.border_index[nyz2*a + nz2*b +c],sizeof(int),1,ouf);	/*3-D "boundary model*/
			}


	fclose(ouf);

	printf("----------------\n");
	printf("The seismic sensitivities are written out in the file senseis%03d.dat\n",ninv);
	printf("----------------\n\n");


	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out the binary file of the total gravity sensitivity*/
/* parameters:   data       := Data structure*/ 
/*				 geo        := Geometry structure */
/*               grid       := Grid structure of the model used for INVERSION*/
/*				 sens		:= Total sensitivities for the cells*/
/*               ninv       := Number of inversions   */
/*Remark: The format is the same like for the model files*/

int WriteSensGravOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid, double *sens)
{

	int index_grav;
	long ny2,nz2,nyz2,nborder;
	long a,b,c;
	float inv_help;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"sengrav%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	nborder = 0; /*!!Attention!! Is set to zero for density data*/

	index_grav = 2; /*Specifying that this is a gravity file*/

	ny2 = grid.ny + 2*grid.nborder;
	nz2 = grid.nz + 2*grid.nborder;
	nyz2=ny2*nz2;

	/*Write out the binary header*/
	fwrite(&index_grav,sizeof(int),1,ouf);							/*Specifying that this is a gravity file*/

	fwrite(&grid.nx,sizeof(int),1,ouf);								/*Number of grid cells in x-direction*/
	fwrite(&grid.ny,sizeof(int),1,ouf);								/*Number of grid cells in y-direction*/
	fwrite(&grid.nz,sizeof(int),1,ouf);								/*Number of grid cells in z-direction*/
	fwrite(&grid.h,sizeof(float),1,ouf);							/*Size of grid cell*/
	fwrite(&nborder,sizeof(int),1,ouf);								/*Number of grid cells of the boundaries is set to 0*/

	for(a=0;a<3;a++) 
		fwrite(&grid.org[a],sizeof(float),1,ouf);					/*Origin(x,y,z) of grid*/

	fwrite(&geo.nstat,sizeof(int),1,ouf);							/*Number of active stations*/

	fwrite(&geo.nstat_grav,sizeof(int),1,ouf);						/*Number of gravity stations*/

	b=0;
	fwrite(&b,sizeof(int),1,ouf);									/*NOT USED*/

	for(a=0;a<geo.nstat_grav;a++)
			fwrite(&(data.gravs[a]),sizeof(int),1,ouf);				/*Location number of the gravity stations*/


	for(a=0;a<geo.nstat_grav;a++)
	{
		/*Memory check*/
		if(data.gravs[a]-1 >= geo.nstat)
		{
			printf("NOT enough memory is allocated: used %d, allocated %d\n",data.gravs[a]-1, geo.nstat);
			exit(0);
		}

		fwrite(&(geo.x[data.gravs[a]-1]),sizeof(float),1,ouf);		/*x-coordinates of gravity station*/
	}

	for(a=0;a<geo.nstat_grav;a++)
		fwrite(&(geo.y[data.gravs[a]-1]),sizeof(float),1,ouf);		/*y-coordinates of gravity station*/	


	for(a=0;a<geo.nstat_grav;a++)
		fwrite(&(geo.z[data.gravs[a]-1]),sizeof(float),1,ouf);		/*z-coordinates of gravity station*/

	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				inv_help = (float)sens[nyz2*(a + grid.nborder) + nz2*(b+grid.nborder) + (c+grid.nborder)];
				fwrite(&inv_help,sizeof(float),1,ouf);					/*Sensitivity model for the gravity*/
			}

	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				fwrite(&grid.border_index[nyz2*(a + grid.nborder) + nz2*(b+grid.nborder) + (c+grid.nborder)],sizeof(int),1,ouf);	/*3-D "boundary model*/
			}																		/*Attention!! Size of the grid differs from the seismic model*/						


	fclose(ouf);

	printf("----------------\n");
	printf("The gravity sensitivities are written out in the file sengrav%03d.dat\n",ninv);
	printf("----------------\n\n");

	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out the binary file of the total MT sensitivity*/
/* parameters:   data       := Data structure*/ 
/*				 geo        := Geometry structure */
/*               grid       := Grid structure of the model used for INVERSION*/
/*				 sens		:= Total sensitivities for the cells*/
/*               ninv       := Number of inversions   */
/*Remark: The format is the same like for the model files*/

int WriteSensMTOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid, double *sens)
{

	int index_res;
	long ny2,nz2,nyz2,nborder;
	long a,b,c;
	float inv_help;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"senres%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	nborder = 0; /*!!Attention!! Is set to zero for resistivity data*/

	index_res = 3; /*Specifying that this is a MT file*/

	ny2 = grid.ny + 2*grid.nborder;
	nz2 = grid.nz + 2*grid.nborder;
	nyz2=ny2*nz2;

	/*Write out the binary header*/
	fwrite(&index_res,sizeof(int),1,ouf);							/*Specifying that this is a MT file*/

	fwrite(&grid.nx,sizeof(int),1,ouf);								/*Number of grid cells in x-direction*/
	fwrite(&grid.ny,sizeof(int),1,ouf);								/*Number of grid cells in y-direction*/
	fwrite(&grid.nz,sizeof(int),1,ouf);								/*Number of grid cells in z-direction*/
	fwrite(&grid.h,sizeof(float),1,ouf);							/*Size of grid cell*/
	fwrite(&nborder,sizeof(int),1,ouf);								/*Number of grid cells of the boundaries is set to 0*/

	for(a=0;a<3;a++) 
		fwrite(&grid.org[a],sizeof(float),1,ouf);					/*Origin(x,y,z) of grid*/

	fwrite(&geo.nstat,sizeof(int),1,ouf);							/*Number of active stations*/

	fwrite(&geo.nstat_mt,sizeof(int),1,ouf);						/*Number of MT stations*/

	b=0;
	fwrite(&b,sizeof(int),1,ouf);									/*NOT USED*/

	for(a=0;a<geo.nstat_mt;a++)
			fwrite(&(data.mts[a]),sizeof(int),1,ouf);				/*Location number of the MT stations*/


	for(a=0;a<geo.nstat_mt;a++)
	{
		/*Memory check*/
		if(data.mts[a]-1 >= geo.nstat)
		{
			printf("NOT enough memory is allocated: used %d, allocated %d\n",data.mts[a]-1, geo.nstat);
			exit(0);
		}

		fwrite(&(geo.x[data.mts[a]-1]),sizeof(float),1,ouf);		/*x-coordinates of MT station*/
	}

	for(a=0;a<geo.nstat_mt;a++)
		fwrite(&(geo.y[data.mts[a]-1]),sizeof(float),1,ouf);		/*y-coordinates of MT station*/	


	for(a=0;a<geo.nstat_mt;a++)
		fwrite(&(geo.z[data.mts[a]-1]),sizeof(float),1,ouf);		/*z-coordinates of MT station*/

	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				inv_help = (float)sens[nyz2*(a + grid.nborder) + nz2*(b+grid.nborder) + (c+grid.nborder)];
				fwrite(&inv_help,sizeof(float),1,ouf);					/*Sensitivity model for the MT*/
			}

	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				fwrite(&grid.border_index[nyz2*(a + grid.nborder) + nz2*(b+grid.nborder) + (c+grid.nborder)],sizeof(int),1,ouf);	/*3-D "boundary model*/
			}																		/*Attention!! Size of the grid differs from the resistivity model*/						


	fclose(ouf);

	printf("----------------\n");
	printf("The MT sensitivities are written out in the file sengrav%03d.dat\n",ninv);
	printf("----------------\n\n");

	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out the binary file of the seismic data structure*/
/* parameters:   data       := Data structure*/
/*               geo		:= Geometry structure*/ 
/*               ninv        := Number of the inversion*/

int WriteDatSeisOut(DATA_STRUCT data, int ninv, GEOMETRY geo)
{
	int a, kind_of_data;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"resdseis%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	kind_of_data = 1;

	/*Write out the binary header*/
	fwrite(&kind_of_data,sizeof(int),1,ouf);							/*Specified that this are seismic data*/

	fwrite(&data.ndata_seis,sizeof(long),1,ouf);						/*Number of shot-receiver combinations*/
	fwrite(&geo.nshot,sizeof(int),1,ouf);								/*Number of shots*/
	fwrite(&geo.nrec,sizeof(int),1,ouf);								/*Number of receivers*/

	for(a=0;a<geo.nshot;a++)
		fwrite(&data.shots[a],sizeof(int),1,ouf);						/*Shot location indeces*/
	for(a=0;a<geo.nrec;a++)
		fwrite(&data.recs[a],sizeof(int),1,ouf);						/*Receiver location indeces*/

	for(a=0;a<geo.nshot;a++)
		fwrite(&data.lshots[a],sizeof(int),1,ouf);							/*Nr of receivers per shot*/
	for(a=0;a<geo.nrec;a++)
		fwrite(&data.lrecs[a],sizeof(int),1,ouf);							/*Nr of shots per receiver*/

	
	for(a=0;a<data.ndata_seis;a++)
		fwrite(&data.sno[a],sizeof(int),1,ouf);							/*Shot location index of the shot-rec. combination a*/
	for(a=0;a<data.ndata_seis;a++)
		fwrite(&data.rno[a],sizeof(int),1,ouf);							/*Receiver location index of the shot-rec. combination a*/
	for(a=0;a<data.ndata_seis;a++)
		fwrite(&data.tobs[a],sizeof(double),1,ouf);						/*Obs. traveltimes of the shot-rec. combination a in ms*/
	for(a=0;a<data.ndata_seis;a++)
		fwrite(&data.tcalc[a],sizeof(double),1,ouf);					/*Calc. traveltimes of the shot-rec. combination a in ms*/
	for(a=0;a<data.ndata_seis;a++)
		fwrite(&data.xdist[a],sizeof(float),1,ouf);						/*Shot-receiver distance in the x-direction in m*/
	for(a=0;a<data.ndata_seis;a++)
		fwrite(&data.ydist[a],sizeof(float),1,ouf);						/*Shot-receiver distance in the y-direction in m*/
	for(a=0;a<data.ndata_seis;a++)
		fwrite(&data.zdist[a],sizeof(float),1,ouf);						/*Shot-receiver distance in the z-direction in m*/
	for(a=0;a<data.ndata_seis;a++)
		fwrite(&data.weigth_seis[a],sizeof(double),1,ouf);					/*Weighting of the rays in the inversion*/

	fclose(ouf);

	printf("----------------\n");
	printf("The binary data file resdseis%03d.dat is written out\n",ninv);
	printf("----------------\n\n");

	return(1);
}


/*---------------------------------------------------------------------------*/
/* Write out the binary file of the gravity data structure*/
/* parameters:   data       := Data structure*/ 
/*               ninv        := Number of the inversion*/
/*               geo		:= Geometry of the stations*/

int WriteDatGravOut(DATA_STRUCT data, int ninv, GEOMETRY geo)
{
	int kind_of_data,a;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"resdgrav%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	kind_of_data = 2;

	/*Write out the binary header*/
	fwrite(&kind_of_data,sizeof(int),1,ouf);							/*Specified that this are gravity data*/

	fwrite(&data.ndata_grav,sizeof(long),1,ouf);						/*Number of gravity data*/
	fwrite(&geo.nstat_grav,sizeof(int),1,ouf);							/*Number of gravity stations*/

	for(a=0;a<geo.nstat_grav;a++)
		fwrite(&data.gravs[a],sizeof(int),1,ouf);						/*Gravity station indeces*/

	for(a=0;a<geo.nstat_grav;a++)
		fwrite(&geo.x[data.gravs[a] - 1],sizeof(float),1,ouf);		    /*x coordinate of the stations*/

	for(a=0;a<geo.nstat_grav;a++)
		fwrite(&geo.y[data.gravs[a] - 1],sizeof(float),1,ouf);		    /*y coordinate of the stations*/

	for(a=0;a<data.ndata_grav;a++)
		fwrite(&data.gno[a],sizeof(int),1,ouf);							/*Gravity location index of the data measurement a*/
	for(a=0;a<data.ndata_grav;a++)
		fwrite(&data.obs_grav[a],sizeof(double),1,ouf);					/*Obs. densities of the data measurement a in mgal*/
	for(a=0;a<data.ndata_grav;a++)
		fwrite(&data.calc_grav[a],sizeof(double),1,ouf);				/*Calc. densities of the data measurement a in mgal*/
	for(a=0;a<data.ndata_grav;a++)
		fwrite(&data.weigth_grav[a],sizeof(double),1,ouf);				/*Weighting of the gravity measurements in the inversion*/


	fclose(ouf);

	printf("----------------\n");
	printf("The binary data file resdgrav%03d.dat is written out\n",ninv);
	printf("----------------\n\n");


	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out the binary file of the MT data structure*/
/* parameters:   data       := Data structure*/ 
/*               ninv        := Number of the inversion*/
/*               geo		:= Geometry of the stations*/
/*			dimension_mt	:= dimesion for the MT-modeling*/
/*          kind_of_data_mt    : Specify the kind of data used for the inversion*/
/*								1= TE-Mode, 2= TM-Mode, 3= Both modes for 2-D, Berdichewsky average for 1-D*/

int WriteDatMTOut(DATA_STRUCT data, int ninv, GEOMETRY geo, int dimension_mt, int kind_of_data_mt)
{
	int kind_of_data,a,b;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"resdmt%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	kind_of_data = 3;

	/*Write out the binary header*/
	fwrite(&kind_of_data,sizeof(int),1,ouf);							/*Specified that this are MT data*/

	fwrite(&dimension_mt,sizeof(int),1,ouf);							/*Dimension of the forward modeling (1=1D or 2=2D or 3=2D forward calculation and 1D inversion, 4= rapid relaxation inversion (RRI))*/

	fwrite(&kind_of_data_mt, sizeof(int),1,ouf);						/*Specify the kind of data used (TE-,TM-mode,both modes, Berdichewsky average*/

	fwrite(&data.ndata_mt,sizeof(long),1,ouf);							/*Number of MT data points*/
	fwrite(&geo.nstat_mt,sizeof(int),1,ouf);							/*Number of MT stations*/

	for(a=0;a<geo.nstat_mt;a++)
		fwrite(&geo.x[data.mts[a] - 1],sizeof(float),1,ouf);		    /*x coordinate of the stations*/

	for(a=0;a<geo.nstat_mt;a++)
		fwrite(&geo.y[data.mts[a] - 1],sizeof(float),1,ouf);		    /*y coordinate of the stations*/

	for(a=0;a<geo.nstat_mt;a++)
		fwrite(&data.mts[a],sizeof(int),1,ouf);							/*MT station indeces*/

	for(a=0;a<data.ndata_mt;a++)
		fwrite(&data.mno[a],sizeof(int),1,ouf);							/*MT location index of the data measurement point a*/

	for(a=0;a<data.ndata_mt;a++)
		fwrite(&data.nfreq_mt[a],sizeof(int),1,ouf);					/*Nr of frequencies for every measurement*/

	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.freq_mt[a][b],sizeof(double),1,ouf);			/*List of frequencies (in MHz) for each data point*/
		}
	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.real_mt_TE[a][b],sizeof(double),1,ouf);			/*Obs.real part of TE-Mode impedance E_y/H_x (Ohmm) for each data point*/
		}
	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.imag_mt_TE[a][b],sizeof(double),1,ouf);			/*Obs.imag.part of TE-Mode impedance E_y/H_x (Ohmm) for each data point*/
		}
	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.real_mt_TM[a][b],sizeof(double),1,ouf);			/*Obs.real part of TM-Mode impedance E_x/H_y (Ohmm) for each data point*/
		}
	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.imag_mt_TM[a][b],sizeof(double),1,ouf);			/*Obs.imag.part of TM-Mode impedance E_x/H_y (Ohmm) for each data point*/
		}
	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.calc_real_mt_TE[a][b],sizeof(double),1,ouf);			/*Calc.real part of TE-Mode impedance E_y/H_x (Ohmm) for each data point*/
		}
	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.calc_imag_mt_TE[a][b],sizeof(double),1,ouf);			/*Calc.imag.part of TE-Mode impedance E_y/H_x (Ohmm) for each data point*/
		}
	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.calc_real_mt_TM[a][b],sizeof(double),1,ouf);			/*Calc.real part of TM-Mode impedance E_x/H_y (Ohmm) for each data point*/
		}
	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.calc_imag_mt_TM[a][b],sizeof(double),1,ouf);			/*Calc.imag.part of TM-Mode impedance E_x/H_y (Ohmm) for each data point*/
		}

	for(a=0;a<data.ndata_mt;a++)
		for(b=0;b<data.nfreq_mt[a];b++)
		{
			fwrite(&data.weigth_mt[a][b],sizeof(double),1,ouf);			/*Weighting of the different frequency dependent measurements for each data point*/
		}

	fclose(ouf);

	printf("----------------\n");
	printf("The binary data file resdmt%03d.dat is written out\n",ninv);
	printf("----------------\n\n");

	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write the merged dat-file out (ascii-File) including infomation about the used shot receiver-combination,*/
/* gravity and MT stations */
/* Parameter:   data :=  Data structure*/
/*		mt_data_format:= 0=no MT data, 1= only one impedance per frequency, 2= two impedances per frequency*/

int WriteMergedDatOut(DATA_STRUCT data, int  mt_data_format)
{
	long i,j,k;
	FILE *out;

			/* Write out opt. Data File: */
			 out = fopen("merged.dat","wt");
			 /*seismic data*/
			 
			 for(i=0;i<data.ndata_seis;i++)
			   	fprintf(out,"%6d %6d %6d %1d   %7f   %7f\n", i+1, data.sno[i], data.rno[i], 1, data.tobs[i], data.weigth_seis[i]);
			 /*gravity data*/
			 j=i;
			 for(i=0;i<data.ndata_grav;i++)
				fprintf(out,"%6d %6d %6d %1d   %7f   %7f\n", i+j+1, data.gno[i], -99999, 2, data.obs_grav[i], data.weigth_grav[i]); 
			 /*MT data*/
			 j=i+j;
			 for(i=0;i<data.ndata_mt;i++)
			 {
				 fprintf(out,"%6d %6d %6d %1d %6d\n", i+j+1, data.mno[i], -88888, 3, data.nfreq_mt[i]);
				 for(k=0;k<data.nfreq_mt[i];k++)
				 {
					if(mt_data_format == 1)	/*MT data format consists only of one impedance per frequency*/
						fprintf(out,"      %7f   %7f   %7f   %7f\n", data.freq_mt[i][k], data.real_mt_TE[i][k], data.imag_mt_TE[i][k], data.weigth_mt[i][k]);
					else
						fprintf(out,"      %7f   %7f   %7f   %7f   %7f   %7f\n", data.freq_mt[i][k], data.real_mt_TE[i][k], data.imag_mt_TE[i][k], data.real_mt_TM[i][k], data.imag_mt_TM[i][k], data.weigth_mt[i][k]);
				 }
			 }

			 fclose(out);

	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out the calculated raypaths of model*/
/* parameters:   rp         := Vector of the ray path structures*/ 
/*				 data		:= Data structure*/	
/*               grid       := Grid structure (of the FORWARD MODEL)*/ 

int WriteRayOut(RP_STRUCT *rp, DATA_STRUCT data, GRID_STRUCT grid)
{
	int b,c;
	int nx2,ny2,nz2,nborder;
	float x,y,z;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"rayout.dat");
	ouf = fopen(fname,"wb");

	nborder = grid.nborder;


	nx2 = grid.nx + 2*nborder +1;
	ny2 = grid.ny + 2*nborder +1;
	nz2 = grid.nz + 2*nborder +1;

	/*Write out the binary file*/
	fwrite(&data.ndata_seis,sizeof(long),1,ouf);						 /*Number of rays*/

	for(b=0;b<data.ndata_seis;b++)
	{
		fwrite(&rp[b].nray,sizeof(long),1,ouf);					/*Number of segments per ray*/

		for(c=0;c<rp[b].nray+1;c++)
		{
			x=(float)(rp[b].x[c]*grid.h - grid.nborder*grid.h + grid.org[0] - 0.5*grid.h);
			fwrite(&x,sizeof(float),1,ouf);				/*Write out the x-positions of the ray segment corners*/
		}
		for(c=0;c<rp[b].nray+1;c++)
		{
			y=(float)(rp[b].y[c]*grid.h - grid.nborder*grid.h + grid.org[1] - 0.5*grid.h);
			fwrite(&y,sizeof(float),1,ouf);				/*Write out the y-positions of the ray segment corners*/
		}
		for(c=0;c<rp[b].nray+1;c++)
		{
			z=(float)(rp[b].z[c]*grid.h - grid.nborder*grid.h + grid.org[2] - 0.5*grid.h);
			fwrite(&z,sizeof(float),1,ouf);				/*Write out the z-positions of the ray segment corners*/

		}
	}
	fclose(ouf);

	printf("----------------\n");
	printf("The  determined raypaths are written out in rayout.dat\n\n");


	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out the calculated fat-rays of model*/
/* parameters:   fr         := Vector of the fat-ray structures*/ 
/*				 data		:= Data structure*/	
/*               grid       := Grid structure (of the FORWARD MODEL)*/ 

int WriteFatRayOut(F_RP_STRUCT *fr,DATA_STRUCT data,GRID_STRUCT grid)
{
	long nx2,ny2,nz2,nyz2,a,b,c;
	float org[3],dummy;
	char fname[32];
	FILE *ouf;

	sprintf(fname,"fat_rayout.dat");
	ouf = fopen(fname,"wb");

	nx2 = grid.nx;
	ny2 = grid.ny;
	nz2 = grid.nz;
	nyz2 = ny2*nz2;

	org[0] =(float) grid.org[0];
	org[1] =(float) grid.org[1];
	org[2] =(float) grid.org[2];

	/*Write out the binary file*/
	fwrite(&nx2,sizeof(int),1,ouf);	    						 /*Number of grid cells in x-direction*/
	fwrite(&ny2,sizeof(int),1,ouf);	 	    					 /*Number of grid cells in y-direction*/
	fwrite(&nz2,sizeof(int),1,ouf); 							 /*Number of grid cells in z-direction*/
	fwrite(&grid.h,sizeof(float),1,ouf);						 /*Size of grid cell*/
	fwrite(&grid.nborder,sizeof(int),1,ouf);					 /*Number of grid cells of the boundaries*/

	for(a=0;a<3;a++) 
		fwrite(&org[a],sizeof(float),1,ouf);					 /*Origin(x,y,z) of grid*/

	fwrite(&data.ndata_seis,sizeof(long),1,ouf);	     				 /*Number of fat-rays*/

	for(b=0;b<data.ndata_seis;b++)
	{

		fwrite(&fr[b].ncell,sizeof(long),1,ouf);				/*Number of segments per fat-ray*/

		for(c=0;c<fr[b].ncell;c++)
		{
			fwrite(&fr[b].ele[c],sizeof(long),1,ouf);			/*Write out the position of the cell*/
		}

		for(c=0;c<fr[b].ncell;c++)
		{
			dummy = (float)fr[b].weight[c];
			fwrite(&dummy,sizeof(float),1,ouf);					/*Write out the weight of the cell*/
		}


	}

	fclose(ouf);

	printf("----------------\n");
	printf("The determined fat-rays are written out in fat_rayout.dat\n\n");


	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out (temporary) binary traveltime models*/
/* parameters:   tt         := Traveltimes in s*/ 
/*               loc        := Location number of the shot/receiver position(where the traveltimes starts;t=0)*/ 
/*                nx,ny,nz  := Number of the grid cell edges of the model*/

/*Remark: nx,ny,nz are the number of the grid cell NODES*/

int WriteTravelTimeModOut(float *tt, int loc,int nx,int ny,int nz)
{
	char fname[32];
	FILE *ouf;
	long a,nxyz;

	nxyz= nx*ny*nz;

	sprintf(fname,"tmp_traveltimes%d.dat",loc);	/*The number in the file-name correspond to the shot/receiver location-nr.*/
	ouf = fopen(fname,"wb");

	/*Write out the binary header*/
	fwrite(&nx,sizeof(int),1,ouf);								/*Number of grid cell edges in x-direction*/
	fwrite(&ny,sizeof(int),1,ouf);								/*Number of grid cell edges in y-direction*/
	fwrite(&nz,sizeof(int),1,ouf);								/*Number of grid cell edges in z-direction*/
	
	for(a=0;a<nxyz;a++)
		fwrite(&(tt[a]),sizeof(float),1,ouf);					/*3-D traveltime model*/

	fclose(ouf);

	printf("The binary file tmp_traveltimes%d.dat containing the traveltimes\n of shot/receiver loc. %d is written out\n\n",loc,loc);

	return(1);
}
	

/*---------------------------------------------------------------------------*/
/* Read in (temporary) the traveltime model from a binary file*/
/* parameters:   fname         := Filename of the binary-file*/ 
/*
/* Output = traveltimes taken from the binary file */

float *ReadTravelTimeModIn(char *fname)
{
	int nx,ny,nz;
	long a,nxyz;
	float *tt;

	FILE *inf;
	
	inf = fopen(fname,"rb");
	if (inf == NULL)
	{
      fprintf(stderr,"Unable to open %s\n",fname);
      exit(0);
	}

	/*Read in the binary file*/
	fread(&nx,sizeof(int),1,inf);								/*Number of grid cell edges in x-direction*/
	fread(&ny,sizeof(int),1,inf);								/*Number of grid cell edges in y-direction*/
	fread(&nz,sizeof(int),1,inf);								/*Number of grid cell edges in z-direction*/
	nxyz = nx*ny*nz;

	tt = (float *)memory(NULL,nxyz,sizeof(float),"ReadTravelTimeModIn");

	for(a=0;a<nxyz;a++)
		fread(&(tt[a]),sizeof(float),1,inf);					/*Traveltimes of the model*/

	fclose(inf);


	printf("The binary file %s is read in\n\n",fname);


	return(tt);
}

/*---------------------------------------------------------------------------*/
/* Write out the ray-densities for the inversion cells*/
/* parameters:   ninv		:= Iteration of inversion */
/*				 ninv_cell  := Number of inversion cells*/ 
/*               nrays		:= Number of rays per inversion cell*/ 
/*				 length		:= The over-all length of all rays in the inversion cells*/
/*               inv_cell	:= Structures of the inversion cells*/

int WriteRayDenseOut(GRID_STRUCT grid, int ninv, long ninv_cell, long *nrays, double *length, BIDX_STRUCT *inv_cell)
{
	long a,b;
	float org[3];
	char fname[32];
	FILE *ouf;

	sprintf(fname,"ray_density%03d.dat",ninv);
	ouf = fopen(fname,"wb");


	/*Write out the binary header*/
	fwrite(&grid.nx,sizeof(int),1,ouf);							/*Number of cells in the x-direction of the forward grid*/
	fwrite(&grid.ny,sizeof(int),1,ouf);							/*Number of cells in the y-direction of the forward grid*/
	fwrite(&grid.nz,sizeof(int),1,ouf);							/*Number of cells in the z-direction of the forward grid*/
	fwrite(&grid.nborder,sizeof(int),1,ouf);					/*Width of border in number of cells*/
	fwrite(&grid.h,sizeof(float),1,ouf);						/*Size of the forward cells*/

	for(a=0;a<3;a++) 
		fwrite(&grid.org[a],sizeof(float),1,ouf);				/*Origin(x,y,z) of grid*/

	fwrite(&ninv_cell,sizeof(long),1,ouf);						/*Number of inversion cells*/

	for(a=0;a<ninv_cell;a++)
	{
		org[0] = (float)inv_cell[a].xo;							/*x-position of the inversion cell center*/
		fwrite(&org[0],sizeof(float),1,ouf);
	}											

	for(a=0;a<ninv_cell;a++)
	{
		org[1] = (float)inv_cell[a].yo;							/*y-position of the inversion cell center*/
		fwrite(&org[1],sizeof(float),1,ouf);
	}											

	for(a=0;a<ninv_cell;a++)
	{
		org[2] = (float)inv_cell[a].zo;							/*z-position of the inversion cell center*/
		fwrite(&org[2],sizeof(float),1,ouf);
	}											

	for(a=0;a<ninv_cell;a++)
		fwrite(&nrays[a],sizeof(long),1,ouf);					/*Number of rays/fat-rays in a inversion cell*/

	for(a=0;a<ninv_cell;a++)
		fwrite(&length[a],sizeof(double),1,ouf);				/*normalized over-all length of rays(weight of fatrays) in a inversion cell*/


	for(a=0;a<ninv_cell;a++)
		fwrite(&inv_cell[a].nele,sizeof(long),1,ouf);			/*Number of forward cells in a inversion cell*/

	for(a=0;a<ninv_cell;a++)
		for(b=0;b<inv_cell[a].nele;b++)
		{
			fwrite(&inv_cell[a].ele[b],sizeof(long),1,ouf);		/*Location numbers of the forward cells in the inversion cells*/
		}



	fclose(ouf);

	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out the ray densities tensor for the inversion cells*/
/* parameters:   ninv		:= Iteration of inversion */
/*				 ninv_cell  := Number of inversion cells*/ 
/*               rel_min_max:= The quotient of the smallest and the largest eigenvalue of the ray density tensor*/
/*               rel_min_max:= The quotient of the medium and the largest eigenvalue of the ray density tensor*/
/*				 length		:= The over-all length of all rays in the inversion cells*/
/*               inv_cell	:= Structures of the inversion cells*/

int	WriteRayDensityTensorOut(GRID_STRUCT grid, int ninv, long ninv_cell, double *rel_minmax_eig, double *rel_medmax_eig, BIDX_STRUCT *inv_cell)
{
	long a,b;
	float org[3];
	char fname[32];
	FILE *ouf;

	sprintf(fname,"ray_density_tensor%03d.dat",ninv);
	ouf = fopen(fname,"wb");

	
	/*Write out the binary header*/
	fwrite(&grid.nx,sizeof(int),1,ouf);							/*Number of cells in the x-direction of the forward grid*/
	fwrite(&grid.ny,sizeof(int),1,ouf);							/*Number of cells in the y-direction of the forward grid*/
	fwrite(&grid.nz,sizeof(int),1,ouf);							/*Number of cells in the z-direction of the forward grid*/
	fwrite(&grid.nborder,sizeof(int),1,ouf);					/*Width of border in number of cells*/
	fwrite(&grid.h,sizeof(float),1,ouf);						/*Size of the forward cells*/

	for(a=0;a<3;a++) 
		fwrite(&grid.org[a],sizeof(float),1,ouf);				/*Origin(x,y,z) of grid*/

	fwrite(&ninv_cell,sizeof(long),1,ouf);						/*Number of inversion cells*/

	for(a=0;a<ninv_cell;a++)
	{
		org[0] = (float)inv_cell[a].xo;							/*x-position of the inversion cell center*/
		fwrite(&org[0],sizeof(float),1,ouf);
	}											

	for(a=0;a<ninv_cell;a++)
	{
		org[1] = (float)inv_cell[a].yo;							/*y-position of the inversion cell center*/
		fwrite(&org[1],sizeof(float),1,ouf);
	}											

	for(a=0;a<ninv_cell;a++)
	{
		org[2] = (float)inv_cell[a].zo;							/*z-position of the inversion cell center*/
		fwrite(&org[2],sizeof(float),1,ouf);
	}				

	for(a=0;a<ninv_cell;a++)
		fwrite(&rel_minmax_eig[a],sizeof(double),1,ouf);		/*Relation of the smallest and the largest eigenvalue in the ray density tensor*/

	for(a=0;a<ninv_cell;a++)
		fwrite(&rel_medmax_eig[a],sizeof(double),1,ouf);		/*Relation of the medium and the largest eigenvalue in the ray density tensor*/

	for(a=0;a<ninv_cell;a++)
		fwrite(&inv_cell[a].nele,sizeof(long),1,ouf);			/*Number of forward cells in a inversion cell*/

	for(a=0;a<ninv_cell;a++)
		for(b=0;b<inv_cell[a].nele;b++)
		{
			fwrite(&inv_cell[a].ele[b],sizeof(long),1,ouf);		/*Location numbers of the forward cells in the inversion cells*/
		}


	fclose(ouf);

	return(1);
}

/*---------------------------------------------------------------------------*/
/* Write out (temporary) binary resistivity/Ex/Ey/Hy/Hx model from the 2D mt calculation*/
/* Parameter:  layer_i = nr of the considered 2-D slice*/
/*                freq = used frequency*/
/*          nsamples   = number of cells (nx*nz)*/
/*         *parameter  = vector of the physical parameter which will be written out*/
/*   kind_of_parameter = Kind of parameter:  0:= resistivity; 1:= Ex(Real); 2:= Ex(Imag) ; 3:= Ey(Real); 4:= Ey(Imag)*/
/*                                                            5:= Hx(Real); 6:= Hx(Imag) ; 7:= Hy(Real); 8:= Hy(Imag)*/
/*                                                            9:= Ez(Real); 10:= Ez(Imag)*/


int WriteMTResistivityOut(long layer_i, double freq, long nsamples, double *parameter, int kind_of_parameter)
{
	long i;
	double tmp_para;
	char fname[40];
	FILE *ouf;

	/********************************************/
	/*Write the resistivities of field components in a temporary file*/
	/*The numbers in the file-name correspond to the layer nr and the frequency*/
	
	if(kind_of_parameter == 1) /*Ex(Real) in V/m*/
		sprintf(fname,"tmp_ExR_layer%d_freq%f.dat",layer_i,freq);	
	else if(kind_of_parameter == 2) /*Ex(Imag) in V/m*/
		sprintf(fname,"tmp_ExI_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 3) /*Ey(Real) in V/m*/
		sprintf(fname,"tmp_EyR_layer%d_freq%f.dat",layer_i,freq);	
	else if(kind_of_parameter == 4) /*Ey(Imag) in V/m*/
		sprintf(fname,"tmp_EyI_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 5) /*Hx(Real) in T*/
		sprintf(fname,"tmp_HxR_layer%d_freq%f.dat",layer_i,freq);	
	else if(kind_of_parameter == 6) /*Hx(Imag) in T*/
		sprintf(fname,"tmp_HxI_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 7) /*Hy(Real) in T*/
		sprintf(fname,"tmp_HyR_layer%d_freq%f.dat",layer_i,freq);	
	else if(kind_of_parameter == 8) /*Hy(Imag) in T*/
		sprintf(fname,"tmp_HyI_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 9) /*Ez(Real) in T*/
		sprintf(fname,"tmp_EzR_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 10) /*Ez(Imag) in T*/
		sprintf(fname,"tmp_EzI_layer%d_freq%f.dat",layer_i,freq);
	else /*resistivity in ohmm*/
		sprintf(fname,"tmp_res_layer%d_freq%f.dat",layer_i,freq);

	ouf = fopen(fname,"wb");

	for(i=0;i<nsamples;i++)
	{
		tmp_para = (double)parameter[i];
		fwrite(&(tmp_para),sizeof(double),1,ouf);					/*parameters of the 2-D slice*/
	}

	fclose(ouf);

	return(0);
}

/*---------------------------------------------------------------------------*/
/* Write out (temporary) binary index model from the 2D mt calculation*/
/* Parameter:  layer_i = nr of the considered 2-D slice*/
/*                freq = used frequency*/
/*          nsamples   = number of cells (nx*nz)*/
/*         *parameter  = vector of the physical parameter which will be written out*/

int WriteMTResistivityIndexOut(long layer_i, double freq, long nsamples, int *parameter)
{
	long i;
	int tmp_para;
	char fname[40];
	FILE *ouf;


	/********************************************/
	/*Write the indeces of field components in a temporary file*/
	/*The numbers in the file-name correspond to the layer nr and the frequency*/
	sprintf(fname,"tmp_indeces_layer%d_freq%f.dat",layer_i,freq);

	ouf = fopen(fname,"wb");

	for(i=0;i<nsamples;i++)
	{
		tmp_para = (int)parameter[i];
		fwrite(&(tmp_para),sizeof(int),1,ouf);					/*parameters of the 2-D slice*/
	}

	fclose(ouf);

	return(0);
}

/*---------------------------------------------------------------------------*/
/* ReadIn (temporary) binary resistivity/Ex/Ey/Hy/Hx model from the 2D mt calculation*/
/* Parameter:  layer_i = nr of the considered 2-D slice*/
/*                freq = used frequency*/
/*          nsamples   = number of cells (nx*nz)*/
/*         *parameter  = vector of the physical parameter which will be written out*/
/*   kind_of_parameter = Kind of parameter:  0:= resistivity; 1:= Ex(Real); 2:= Ex(Imag) ; 3:= Ey(Real); 4:= Ey(Imag)*/
/*                                                            5:= Hx(Real); 6:= Hx(Imag) ; 7:= Hy(Real); 8:= Hy(Imag)*/
/*                                                            9:= Ez(Real); 10:= Ez(Imag)*/

int ReadTResistivityIn(long layer_i, double freq, long nsamples, double *parameter, int kind_of_parameter)
{
	long a;
	double tmp_para;
	char fname[40];
	FILE *inf;

	/********************************************/
	/*Read the resistivities or field components from a temporary file*/
	/*The numbers in the file-name correspond to the layer nr and the frequency*/
	
	if(kind_of_parameter == 1) /*Ex(Real) in V/m*/
		sprintf(fname,"tmp_ExR_layer%d_freq%f.dat",layer_i,freq);	
	else if(kind_of_parameter == 2) /*Ex(Imag) in V/m*/
		sprintf(fname,"tmp_ExI_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 3) /*Ey(Real) in V/m*/
		sprintf(fname,"tmp_EyR_layer%d_freq%f.dat",layer_i,freq);	
	else if(kind_of_parameter == 4) /*Ey(Imag) in V/m*/
		sprintf(fname,"tmp_EyI_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 5) /*Hx(Real) in T*/
		sprintf(fname,"tmp_HxR_layer%d_freq%f.dat",layer_i,freq);	
	else if(kind_of_parameter == 6) /*Hx(Imag) in T*/
		sprintf(fname,"tmp_HxI_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 7) /*Hy(Real) in T*/
		sprintf(fname,"tmp_HyR_layer%d_freq%f.dat",layer_i,freq);	
	else if(kind_of_parameter == 8) /*Hy(Imag) in T*/
		sprintf(fname,"tmp_HyI_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 9) /*Ez(Real) in T*/
		sprintf(fname,"tmp_EzR_layer%d_freq%f.dat",layer_i,freq);
	else if(kind_of_parameter == 10) /*Ez(Imag) in T*/
		sprintf(fname,"tmp_EzI_layer%d_freq%f.dat",layer_i,freq);
	else /*resistivity in ohmm*/
		sprintf(fname,"tmp_res_layer%d_freq%f.dat",layer_i,freq);

	inf = fopen(fname,"rb");

	if (inf == NULL)
	{
      fprintf(stderr,"Unable to open %s\n",fname);
      exit(0);
	}

	for(a=0;a<nsamples;a++)
	{
		fread(&(tmp_para),sizeof(double),1,inf);					/*Read in the resistivities or the field components*/
		parameter[a] = tmp_para;
	}

	fclose(inf);

	return(0);
}