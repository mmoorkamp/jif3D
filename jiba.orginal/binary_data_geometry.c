#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "inv3d.h"
#include "fkt.h"

/*-------------------------------------------------------------*/
/* Build the grid structure using the binary model */
/* terminates) */
/* Parameter:	*/
/*				nmmod_file := Number of the binary model files*/
/*              *grid := Pointer on the grid structure*/
/*				par   := Input parameter structure*/
/*				*fname:= File name including the velocity-resistivity relationship*/
/*                geo := Geometry structure*/
/* Remark: The kind of inversions are determined by the binary model files that are read in*/

int ReadFromBinary(int nmod_file, GRID_STRUCT *grid, PAR par, FLAG_STRUCT *flag, GEOMETRY geo, char *fname)
{
	int i, index, count; /*Number of the inversion*/
	int kind_of_mod[4]; /*Indeces to check, what kind of models are read in (seismic, gravity, and/or MT)*/
	int tmp_nborder, *tmp_border_index;
	long a,b,c,nx2,ny2,nz2,nyz2, nxyz2;
	long nx2_tmp, ny2_tmp, nz2_tmp, nyz2_tmp, nxyz2_tmp;
	long  nyz3,nxyz3,nz3;
	double vel_value;

	REL_VDR_STRUCT rel_vdr;
	GRID_STRUCT tmp_grid[3];

	if(nmod_file == 0) /*Check if MOD-files exist*/
	{
		printf("NO MOD-File is written in !!!\n");
		exit(0);
	}

	for(i=0;i<4;i++)
		kind_of_mod[i] = 0;

	nxyz2_tmp = 0;
	nx2_tmp = 0; ny2_tmp = 0; nz2_tmp = 0;
	nyz2_tmp = 0;
	tmp_nborder = 0;
	count = 0;

	tmp_border_index = (int *)memory(NULL,1,sizeof(int),"ReadFromBinary");

	grid->border_index = (int *)memory(NULL,1,sizeof(int),"ReadfromBinary");
	grid->slow = (double *)memory(NULL,1,sizeof(double),"ReadFromBinary");
	grid->dens = (double *)memory(NULL,1,sizeof(double),"ReadFromBinary");
	grid->res = (double *)memory(NULL,1,sizeof(double),"ReadFromBinary");


	for(i=0;i<nmod_file;i++)
	{
	    /*Read in binary files for the model*/
		index = ReadModIn(grid, par.binary_mod_dat[i].name);

		/****************************************************/
		/*Check, if the model files have the right specification*/
		if(index != 1 && index !=2 && index != 3)
		{
			printf("The specification of the model type within the MOD-file %s is not\ncorrect !", par.binary_mod_dat[i].name);
			exit(0);

		}
		/****************************************************/
		/*Check,if the same models are not read in two times:*/
		if(kind_of_mod[index] == 1 && index == 1)
		{
			printf("TWO MOD files with a seismic model are read in!\nRemove one of the two files!\n");
			exit(0);
		}

		if(kind_of_mod[index] == 1 && index == 2)
		{
			printf("TWO MOD files with a gravity model are read in!\nRemove one of the two files!\n");
			exit(0);
		}

		if(kind_of_mod[index] == 1 && index == 3)
		{
			printf("TWO MOD files with a MT model are read in!\nRemove one of the two files!\n");
			exit(0);
		}

		kind_of_mod[index] = 1;
		/****************************************************/

		tmp_grid[i] = *grid;

		/*Take the number of border cells from the seismic file:*/
		if(index ==1)
			tmp_nborder = grid->nborder;

		nx2 = grid->nx + 2*grid->nborder;
		ny2 = grid->ny + 2*grid->nborder;
		nz2 = grid->nz + 2*grid->nborder;
		nyz2 = ny2*nz2;
		nxyz2 = nx2*ny2*nz2;

			
		/*Take the border_index from the largest grid size*/
		if(nxyz2_tmp < nxyz2)
		{
			tmp_border_index = (int *)memory((char *)tmp_border_index,nxyz2,sizeof(int),"ReadfromBinary");

			for(a=0;a<nx2;a++)
				for(b=0;b<ny2;b++)
					for(c=0;c<nz2;c++)
					{
					   tmp_border_index[nyz2*a + nz2*b + c] = grid->border_index[nyz2*a + nz2*b + c];
					}

			nx2_tmp = nx2; ny2_tmp = ny2; nz2_tmp = nz2;
			nyz2_tmp = nyz2;
			nxyz2_tmp = nxyz2;
		}


		count++;
	}

	/*Check, if the GRID structure is the same in all models*/
	for(i=0;i<count;i++)
	{
		if(grid->nx != tmp_grid[i].nx ||
		   grid->ny != tmp_grid[i].ny ||
		   grid->nz != tmp_grid[i].nz ||
		   grid->h != tmp_grid[i].h ||
		   grid->org[0] != tmp_grid[i].org[0] ||
		   grid->org[1] != tmp_grid[i].org[1] ||
		   grid->org[2] != tmp_grid[i].org[2])
		{
			printf("The grid parameter are different in the different MOD-Files!!\n");
			exit(0);
		}
	}

	/*Reset the number of border cells*/
	grid->nborder = tmp_nborder;

	/*Set the topography index*/
	grid->topo_index = 0;

	/*Remark: The topography information will be used for the refinement of the gravity grid and/or MT grid*/
	for(i=0;i<(geo.nstat);i++)
	{
		if(geo.coor_info[i] == 1)
		{
			grid->topo_index = 1;
		}
	}

	if(grid->topo_index == 1)
	{
		printf("\nRemark: Because topography points are found in the GEO-file\n");
		printf("The topography is activated\n\n");
	}


	/*Reset the border_index*/
	grid->border_index = (int *)memory((char *)grid->border_index,nxyz2_tmp,sizeof(int),"ReadfromBinary");

	for(a=0;a<nx2_tmp;a++)
		for(b=0;b<ny2_tmp;b++)
			for(c=0;c<nz2_tmp;c++)
			{
				grid->border_index[nyz2_tmp*a + nz2_tmp*b + c] = tmp_border_index[nyz2_tmp*a + nz2_tmp*b + c];
			}


	/*Reset flags what kind of model is read in and what kind of inversions will be performed*/
	if(kind_of_mod[1] == 1)
		flag->index_tseis = 1;
	else
		flag->index_tseis = 0;

	if(kind_of_mod[2] == 1)
	{
		grid->dens_air_water = par.g_air_water; /*Set the density of water from the input file; required for the refinement of the topography*/
		flag->index_grav = 1;
	}
	else
		flag->index_grav = 0;
	
	if(kind_of_mod[3] == 1)
	{
		grid->res_air_water = par.r_air_water; /*Set the resistivity of water from the input file; required for the refinement of the topography*/
		flag->index_mt = 1;
	}
	else
		flag->index_mt = 0;

	if(kind_of_mod[1] == 1 && kind_of_mod[2] == 1 && kind_of_mod[3] != 1)
		flag->index_tseis_grav = 1;
	else
		flag->index_tseis_grav = 0;

	if(kind_of_mod[1] == 1 && kind_of_mod[3] == 1 && kind_of_mod[2] != 1)
		flag->index_tseis_mt = 1;
	else
		flag->index_tseis_mt = 0;

	if(kind_of_mod[2] == 1 && kind_of_mod[3] == 1 && kind_of_mod[1] != 1)
		flag->index_grav_mt = 1;
	else
		flag->index_grav_mt = 0;

	if(kind_of_mod[1] == 1 && kind_of_mod[2] == 1 && kind_of_mod[3] == 1)
		flag->index_tseis_grav_mt = 1;
	else
		flag->index_tseis_grav_mt = 0;

	if(flag->index_tseis_grav != 0 || flag->index_tseis_mt != 0 || flag->index_grav_mt != 0 || flag->index_tseis_grav_mt != 0)
	   flag->index_joint = 1;
    else	
	   flag->index_joint = 0;


	/*Make the slowness model for the case that only MT and gravity were calculated*/
	nz3 = nz2_tmp+1;
	nyz3 = (nz2_tmp+1)*(ny2_tmp+1);
	nxyz3 = (nx2_tmp+1)*(ny2_tmp+1)*(nz2_tmp+1);

	if(flag->index_mt != 0 && flag->index_grav != 0 && flag->index_tseis == 0)
	{
		grid->slow = (double *)memory((char *) grid->slow,nxyz3,sizeof(double),"ReadFromBinary");

		

		rel_vdr.tab_v = (double *)memory(NULL,1,sizeof(double),"ReadFromBinary");
		rel_vdr.tab_dr = (double *)memory(NULL,1,sizeof(double),"ReadFromBinary");


		/*Read in the file including the relationship of velocity and density*/
		ReadRelVelDensRes(fname, &rel_vdr);

		printf("The velocity model is determined by linking it to the density model\n\n");

		for(a=0;a<nx2_tmp;a++)
			for(b=0;b<ny2_tmp;b++)
				for(c=0;c<nz2_tmp;c++)
				{
					if(grid->border_index[nyz2_tmp*a + nz2_tmp*b + c] != 0)
					{
						/*Determine the corresponding velocity value from the density model*/
						vel_value = DensRestoVel(grid->dens[nyz2_tmp*a + nz2_tmp*b + c],rel_vdr);
					}
					else
					{
						/*ATTENTION!!! This is the only case where the velocities of the air have to be filled in the vel. model*/
						vel_value = par.v_air_water;
					}

					/*Make the velocity grid*/
					grid->slow[nyz3*a + nz3*b + c] = grid->h/vel_value;

				}

		free(rel_vdr.tab_v);
		free(rel_vdr.tab_dr);
	}


	free(tmp_border_index);

	return(1);
}



/*---------------------------------------------------------------------------*/
/* Read in a model from a binary file*/
/* parameters:   *grid       := Grid structure of the model*/
/*				 *fname		 := Filename of the binary input file*/
/* Output:		index that specifies the kind of model(seismic,gravity,MT) */

int ReadModIn(GRID_STRUCT *grid, char *fname)
{
	int index;
	int nstat1, nstat2, nstat_all;
	int recs,shots;
	float xs,ys,zs,xr,yr,zr;
	int inv_help2;
	long a,b,c,nx2,ny2,nz2,nyz2,nxyz2;
	long nx3,ny3,nz3,nyz3,nxyz3;
	float inv_help;
	FILE *inf;

   inf = fopen(fname,"rb");
   if (inf == NULL)
   {
      fprintf(stderr,"Unable to open %s\n",fname);
      exit(0);
   }

   /*Read in the binary header*/
	fread(&index,sizeof(int),1,inf);									/*Index specifying the kind of model*/

	fread(&(grid->nx),sizeof(int),1,inf);								/*Number of grid cells in x-direction*/
	fread(&(grid->ny),sizeof(int),1,inf);								/*Number of grid cells in y-direction*/
	fread(&(grid->nz),sizeof(int),1,inf);								/*Number of grid cells in z-direction*/
	fread(&(grid->h),sizeof(float),1,inf);								/*Size of grid cell*/
	fread(&(grid->nborder),sizeof(int),1,inf);							/*Number of grid cells of the boundaries*/

	for(a=0;a<3;a++) 
		fread(&(grid->org[a]),sizeof(float),1,inf);						/*Origin(x,y,z) of grid*/

	fread(&nstat_all,sizeof(int),1,inf);									/*Number of all stations*/

/*************************************************************************/

		fread(&nstat1,sizeof(int),1,inf);									/*Number of active shots, gravity or MT stations*/
		fread(&nstat2,sizeof(int),1,inf);									/*Number of active receivers*/

		for(a=0;a<nstat1;a++)
				fread(&shots,sizeof(int),1,inf);						/*Indices of the shot/gravity or MT positions*/
		for(a=0;a<nstat2;a++)
				fread(&recs,sizeof(int),1,inf);						/*Indices of the receiver positions*/

		
		for(a=0;a<nstat1;a++)
				fread(&xs,sizeof(float),1,inf);							/*x-coordinates of shot, gravity station or MT station positions*/
		for(a=0;a<nstat1;a++)
				fread(&ys,sizeof(float),1,inf);							/*x-coordinates of shot, gravity station or MT station positions*/
		for(a=0;a<nstat1;a++)
				fread(&zs,sizeof(float),1,inf);							/*x-coordinates of shot, gravity station or MT station positions*/

		for(a=0;a<nstat2;a++)
			fread(&xr,sizeof(float),1,inf);								/*x-coordinates of receiver, gravity station or MT station positions*/
		for(a=0;a<nstat2;a++)
			fread(&yr,sizeof(float),1,inf);								/*x-coordinates of receiver, gravity station or MT station positions*/
		for(a=0;a<nstat2;a++)
			fread(&zr,sizeof(float),1,inf);								/*x-coordinates of receiver, gravity station or MT station positions*/

		nx2 = grid->nx + 2*grid->nborder;
		ny2 = grid->ny + 2*grid->nborder;
		nz2 = grid->nz + 2*grid->nborder;
		nyz2 = ny2*nz2;
		nxyz2 = nx2*ny2*nz2;

		nx3 = nx2 +1;
		ny3 = ny2 +1;
		nz3 = nz2 +1;
		nyz3 = ny3*nz3;
		nxyz3 = nx3*ny3*nz3;

		if(index == 1) /*seismic velocity model*/
		{
			grid->slow = (double *)memory((char *)grid->slow,nxyz3,sizeof(double),"ReadModin");


			for(a=0;a<nx2;a++)
				for(b=0;b<ny2;b++)
					for(c=0;c<nz2;c++)
					{
						fread(&inv_help,sizeof(float),1,inf);					/*3-D velocity model*/
						grid->slow[nyz3*a + nz3*b + c] =  grid->h/inv_help;
					}

		}
		else if(index == 2) /*density model*/
		{
			grid->dens = (double *)memory((char *)grid->dens,nxyz2,sizeof(double),"ReadModin");

			for(a=0;a<nx2;a++)
				for(b=0;b<ny2;b++)
					for(c=0;c<nz2;c++)
					{
						fread(&inv_help,sizeof(float),1,inf);					/*3-D density model*/
						grid->dens[nyz2*a + nz2*b + c] =  inv_help;
					}

		}
		else if(index ==3) /*resistivity model*/
		{

			grid->res = (double *)memory((char *)grid->res,nxyz2,sizeof(double),"ReadModin");

			for(a=0;a<nx2;a++)
				for(b=0;b<ny2;b++)
					for(c=0;c<nz2;c++)
					{
						fread(&inv_help,sizeof(float),1,inf);					/*3-D resistivity model*/
						grid->res[nyz2*a + nz2*b + c] =  inv_help;
					}
		}


		grid->border_index = (int *)memory((char *)grid->border_index,nxyz2,sizeof(int *),"ReadModin");

		for(a=0;a<nx2;a++)
			for(b=0;b<ny2;b++)
				for(c=0;c<nz2;c++)
				{
					fread(&inv_help2,sizeof(int),1,inf);					/*3-D "grid index model*/
					grid->border_index[nyz2*a + nz2*b + c] = inv_help2;
				}


    fclose(inf);

	printf("----------------\n");
	printf("The binary model file %s is read in\n",fname);
	printf("----------------\n");

	return(index); /*<-write out the index that specify the kind of model(seismic,gravity,MT)*/
}