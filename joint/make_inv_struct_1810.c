/*Contain some changes about to determine the border of inversion cells by Jin on 5th, July, 2006.*/

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
/*Determine the number of inversion cells*/
/*Parameter:	grid	:= Grid structure*/
/*				inv		:= Inversion structure */
int DetNumInvCells(GRID_STRUCT grid, INV *inv)
{
long nx2,ny2,nz2,nyz2;

	nx2 = grid.nx + 2*grid.nborder;
	ny2 = grid.ny + 2*grid.nborder;
	nz2 = grid.nz + 2*grid.nborder;
	nyz2 =ny2*nz2;

	/*Build the inversion structure (by means of a regular forward grid)*/
	inv->nx = (long)(ceil((double)nx2/grid.grid_cell_factor[0]));
	inv->ny = (long)(ceil((double)ny2/grid.grid_cell_factor[1]));
	inv->nz = (long)(ceil((double)nz2/grid.grid_cell_factor[2]));
	inv->nvel = inv->nx*inv->ny*inv->nz; /*Number of cells*/

	return(0);
}


/*------------------------------------------------------------*/
/*Build the inversion structures by means of the forward grid*/
/*Parameter:     grid			:= Grid structure*/
/*				*inv_cell		:= Structure of the inversion cells*/
/*				inv				:= Inversion structure */
/*				flag			:= Structure including flags*/

/*Output: number of active inversion cells*/


long MakeRegularInvStruct(GRID_STRUCT grid,BIDX_STRUCT *inv_cell,INV inv, FLAG_STRUCT flag)
{
	long i,a,b,c;
	long ia,ib,ic;
	long nx,ny,nz,nyz;
	long nx1,ny1,nz1,nyz1;
	long ny2,nz2,nyz2;
	long nyz_inv, nxyz_inv;
	long nxyz_inv_used;
	long ix_inv,iy_inv,iz_inv;
	long nr_cell,nr_cell_active;

	double xpos,ypos,zpos;

	#define Coor(x,y,z) nyz_inv*(x)+(inv.nz)*(y)+(z)

	/*-----------------------------------------------------------*/
	/*Regular grid:*/

	/*Dimensions of the forward grid:*/
	nx = grid.nx;
	ny = grid.ny;
	nz = grid.nz;
	nyz = ny*nz;

	nx1 = grid.nx + 2*grid.nborder;
	ny1 = grid.ny + 2*grid.nborder;
	nz1 = grid.nz + 2*grid.nborder;
	nyz1 = ny1*nz1;

	nz2 = nz1+1;
	ny2 = ny1+1;
	nyz2 = ny2*nz2;

	nyz_inv  = inv.ny*inv.nz;
	nxyz_inv = inv.nvel;

	/*Initialize the inversion grid*/
	for(i=0;i<nxyz_inv;i++)
	{
	    inv_cell[i].bdamp = 0;
		inv_cell[i].val_slow = 0;
		inv_cell[i].val_dens = 0;
		inv_cell[i].val_res = 0;
		inv_cell[i].nele = 0;
		inv_cell[i].use = 0;
		inv_cell[i].used_nr = -1;

		inv_cell[i].ele = (long *)memory(NULL,1,sizeof(long),"MakeInvStruct");
		inv_cell[i].ele_inv = (long *)memory(NULL,1,sizeof(long),"MakeInvStruct");

		/*Jin, 5th, July, 2006.*/
		inv_cell[i].border_inv_cell_x_left = 0; 
		inv_cell[i].border_inv_cell_x_right = 0; 

		inv_cell[i].border_inv_cell_y_front = 0; 
		inv_cell[i].border_inv_cell_y_back = 0; 

		inv_cell[i].border_inv_cell_z = 0; 

	}


	nxyz_inv_used = 0; /*Active inversion cells*/

	/*Loop over the inversion cells*/
	for(ix_inv=0;ix_inv<inv.nx;ix_inv++)
		for(iy_inv=0;iy_inv<inv.ny;iy_inv++)
			for(iz_inv=0;iz_inv<inv.nz;iz_inv++)
			{
				
				if(Coor(ix_inv,iy_inv,iz_inv) >= nxyz_inv)
				{
					printf("NOT enough memory is allocated: used %d, allocated %d\n",Coor(ix_inv,iy_inv,iz_inv)+1, nxyz_inv);
					exit(0);
				}

				
				xpos = 0.0;
			    ypos = 0.0;
			    zpos = 0.0;

				nr_cell_active = 0;
				nr_cell = 0;

				ia = 0;
				/*Loop over the forward cells*/
				for(a= ix_inv*grid.grid_cell_factor[0];a<nx1 && a<(ix_inv+1)*grid.grid_cell_factor[0]; a++)
				{
					ib = 0;
					for(b=iy_inv*grid.grid_cell_factor[1];b<ny1 && b<(iy_inv+1)*grid.grid_cell_factor[1]; b++)
					{
						ic = 0;
						for(c=iz_inv*grid.grid_cell_factor[2];c<nz1 && c<(iz_inv+1)*grid.grid_cell_factor[2]; c++)		
						{
							nr_cell++;

							/*Assign forward cells to the inversion cells; cells in the "air" are NOT assigned*/
							if(0 != grid.border_index[nyz1*a + nz1*b +c])
							{

								inv_cell[Coor(ix_inv,iy_inv,iz_inv)].nele++;
								nr_cell_active = inv_cell[Coor(ix_inv,iy_inv,iz_inv)].nele;/*number of forward cells in the inversion cell*/
								
								inv_cell[Coor(ix_inv,iy_inv,iz_inv)].ele = (long *)memory((char *)inv_cell[Coor(ix_inv,iy_inv,iz_inv)].ele,nr_cell_active,sizeof(long),"MakeInvStruct");
								inv_cell[Coor(ix_inv,iy_inv,iz_inv)].ele[nr_cell_active-1] = a*nyz1 + b*nz1 + c; /*location numbers of forward cells*/

								if(flag.index_tseis != 0 || (flag.index_grav != 0 && flag.index_mt != 0))
									inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_slow = inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_slow + (grid.slow[nyz2*a + nz2*b + c]/grid.h);

								if(flag.index_grav != 0)
								{
									if(a >= grid.nborder && b >= grid.nborder && c >= grid.nborder && a < (grid.nborder + nx) && b < (grid.nborder + ny) && c < (grid.nborder + nz))
										inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_dens = inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_dens + grid.dens[nyz*(a - grid.nborder) + nz*(b - grid.nborder) + (c - grid.nborder)];
								}

								if(flag.index_mt != 0)
								{
									if(a >= grid.nborder && b >= grid.nborder && c >= grid.nborder && a < (grid.nborder + nx) && b < (grid.nborder + ny) && c < (grid.nborder + nz))
										inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_res = inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_res + grid.res[nyz*(a - grid.nborder) + nz*(b - grid.nborder) + (c - grid.nborder)];
								}

								inv_cell[Coor(ix_inv,iy_inv,iz_inv)].use = 1; /*inversion cell is "active"*/ 

								/*Jin, 5th, July, 2006.*/
								if(c == grid.nborder || c == (nz1 - 1 - grid.nborder))
									inv_cell[Coor(ix_inv,iy_inv,iz_inv)].border_inv_cell_z = 1; /*Cells is close to the border of the grid in z direction*/
							}
							ic++;
							/*Calculate the center of the inversion cells:*/
							xpos = grid.org[0]+ grid.h*(a - grid.nborder) + xpos;
							ypos = grid.org[1]+ grid.h*(b - grid.nborder) + ypos;
							zpos = grid.org[2]+ grid.h*(c - grid.nborder) + zpos;
						}
						ib++;
						/*Calculate the dimension of the inv. cell in z-direction*/
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].zdim = grid.h*ic;

						/*Jin, 5th, July, 2006.*/
						/*If only one cell exists in y-direction, this cells will NOT be considered as cells at the border*/
						if(b == grid.nborder && grid.ny > 1)
							inv_cell[Coor(ix_inv,iy_inv,iz_inv)].border_inv_cell_y_front = 1; /*Cells is close to the border of the grid in y direction*/
						if(b == (ny1 - 1 - grid.nborder) && grid.ny > 1)
							inv_cell[Coor(ix_inv,iy_inv,iz_inv)].border_inv_cell_y_back = 1; /*Cells is close to the border of the grid in y direction*/

					}
					ia++;
					/*Calculate the dimension of the inv. cell in y-direction*/
					inv_cell[Coor(ix_inv,iy_inv,iz_inv)].ydim = grid.h*ib;

					/*Jin, 5th, July, 2006.*/
					/*If only one cell exists in x-direction, this cells will NOT be considered as cells at the border*/
					if(a == grid.nborder && grid.nx > 1)
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].border_inv_cell_x_left = 1; /*Cells is close to the border of the grid in x direction*/
					if(a == (nx1 - 1 - grid.nborder) && grid.nx > 1)
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].border_inv_cell_x_right = 1; /*Cells is close to the border of the grid in x direction*/

				}

				/*Calculate the dimension of the inv. cell in x-direction*/
				inv_cell[Coor(ix_inv,iy_inv,iz_inv)].xdim = grid.h*ia;

				if(nr_cell <= 0)
				{
					printf("ERROR: Number of forward cells in the inversion cell %d is 0\n", Coor(ix_inv,iy_inv,iz_inv));
					exit(0);
				}

				inv_cell[Coor(ix_inv,iy_inv,iz_inv)].xo = xpos/nr_cell; /*coordinates of the cell centers*/
				inv_cell[Coor(ix_inv,iy_inv,iz_inv)].yo = ypos/nr_cell;
				inv_cell[Coor(ix_inv,iy_inv,iz_inv)].zo = zpos/nr_cell;

				

				if(inv_cell[Coor(ix_inv,iy_inv,iz_inv)].use == 1) /*cell is used for inversion*/
				{

					if(nr_cell_active <= 0)
					{
						printf("ERROR: Number of active forward cells in the inversion cell %d is 0\n", Coor(ix_inv,iy_inv,iz_inv));
						exit(0);
					}

					if(flag.index_tseis != 0 || (flag.index_grav != 0 && flag.index_mt != 0))
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_slow = inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_slow/nr_cell_active; /*slowness(The slowness in the air is not considered)*/
					else
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_slow = 0.0;

					if(flag.index_grav != 0)
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_dens = inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_dens/nr_cell_active; /*density(The density in the air is not considered)*/
					else
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_dens = 0.0;
					
					if(flag.index_mt != 0)
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_res = inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_res/nr_cell_active; /*resistivity(The resistivity in the air is not considered)*/
					else
						inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_res = 0.0;

					inv_cell[Coor(ix_inv,iy_inv,iz_inv)].bdamp = (double)nr_cell_active/(double)nr_cell; /*damping*/

					inv_cell[Coor(ix_inv,iy_inv,iz_inv)].used_nr = nxyz_inv_used;
					nxyz_inv_used++;
				}
				else /*cell is not used for inversion*/
				{
					inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_slow = 0.0;
					inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_dens = 0.0;
					inv_cell[Coor(ix_inv,iy_inv,iz_inv)].val_res = 0.0;
					inv_cell[Coor(ix_inv,iy_inv,iz_inv)].bdamp = 0.0;
				}
			

			}

	#undef Coor

	printf("The inversion cell structure is determined\n");
	printf("%d of %d cells are used for inversion\n", nxyz_inv_used, inv.nvel);
	printf("----------------\n");

	return(nxyz_inv_used);
}


/*------------------------------------------------------------*/
/*Determine the ray-structs of the inversion cells*/
/*Parameter:	grid	:= Grid structure*/
/*				raypath	:= Input:	Structure including the raypath of the forward cells;*/
/*						   Output:	Structure including the raypath of the inversion cells;*/
/*							(Remark: -The ray-length in the output structures are NOT normalized any more)
/*									 -The x,y,z components are not used any more*/
/*									 -The shot-receiver index remains the same*/
/*									 -Rays in the air or at the borders are NOT considered*/
/*				inv		:= Inversion structure */

int InvRayStruct(GRID_STRUCT grid,RP_STRUCT *raypath, INV *inv, DATA_STRUCT data)
{
	int index;
	long i,j,k,a,nr_inv_cell;
	long ix,iy,iz,nx2,ny2,nz2,nyz2;
	long inv_ix,inv_iy,inv_iz,inv_nz,inv_nyz;
	long inv_cell;
	RP_STRUCT tmp_raypath;

	nx2 = grid.nx + 2*grid.nborder;
	ny2 = grid.ny + 2*grid.nborder;
	nz2 = grid.nz + 2*grid.nborder;
	nyz2 =ny2*nz2;

	inv_nz = inv->nz; 
	inv_nyz = inv->ny *inv->nz;

	/*************************************************************/
	/*Make a raypath structure for the inversion cells*/
	/*Loop over all rays*/
	for(i=0;i<data.ndata_seis;i++)
	{
		tmp_raypath.ele = (long *)memory(NULL,1,sizeof(long),"InvRayStruct");
		tmp_raypath.len = (double *)memory(NULL,1,sizeof(double),"InvRayStruct");
		tmp_raypath.x = (double *)memory(NULL,1,sizeof(double),"InvRayStruct");
		tmp_raypath.y = (double *)memory(NULL,1,sizeof(double),"InvRayStruct");
		tmp_raypath.z = (double *)memory(NULL,1,sizeof(double),"InvRayStruct");

		tmp_raypath.nray =0;
		nr_inv_cell=0;

		/*Loop over all ray-segments*/
		for(j=0;j<raypath[i].nray;j++)
		{
			 if(grid.border_index[raypath[i].ele[j]] != 0) /*Consider just cell that are not in the air*/
			 {
					/*position in the forward grid*/
				 	ix = (long)floor((double)(raypath[i].ele[j]/nyz2));
					iy = (long)floor((double)((raypath[i].ele[j] - ix*nyz2)/nz2));
					iz = raypath[i].ele[j] - ix*nyz2 - iy*nz2;

					/*position in the inversion grid*/
					inv_ix = (long)(floor((double)ix/grid.grid_cell_factor[0]));
					inv_iy = (long)(floor((double)iy/grid.grid_cell_factor[1]));
					inv_iz = (long)(floor((double)iz/grid.grid_cell_factor[2]));

					inv_cell = inv_nyz*inv_ix + inv_nz*inv_iy + inv_iz;

					if(inv_cell > inv->nvel)
					{
						printf("The inversion cell number %d does not exist;\nThe highest inversion cell number is %d\n",inv_cell,inv->nvel);
						exit(0);
					}

					index = 0;

					/*Check, if ray has run already through the inversion cell*/
					for(a=0;a<nr_inv_cell;a++)
					{
						if(tmp_raypath.ele[a] == inv_cell)
						{
							index = 1;
							break;
						}
					}

					/*Build the ray structure for the inversion cells*/
					if(index == 0)
					{
						tmp_raypath.ele = (long *)memory((char *)tmp_raypath.ele,nr_inv_cell+1,sizeof(long),"InvRayStruct");
						tmp_raypath.ele[nr_inv_cell] = inv_cell;

						tmp_raypath.len = (double *)memory((char *)tmp_raypath.len,nr_inv_cell+1,sizeof(double),"InvRayStruct");
						tmp_raypath.len[nr_inv_cell] = raypath[i].len[j]*grid.h; /*remove the normalization*/


						/********************************************************************************************************/
						/*ATTENTION: In the forward grid the x.y and.z component describe the position, where the ray enter the cell. 
						However, in the inversion grid it describe the components of the average x,y and z component of the ray */
						/********************************************************************************************************/

						tmp_raypath.x = (double *)memory((char *)tmp_raypath.x,nr_inv_cell+1,sizeof(double),"InvRayStruct");
					    tmp_raypath.x[nr_inv_cell] = (raypath[i].x[j+1] - raypath[i].x[j])*grid.h; /*remove the normalization*/

						tmp_raypath.y = (double *)memory((char *)tmp_raypath.y,nr_inv_cell+1,sizeof(double),"InvRayStruct");
					    tmp_raypath.y[nr_inv_cell] = (raypath[i].y[j+1] - raypath[i].y[j])*grid.h; /*remove the normalization*/

						tmp_raypath.z = (double *)memory((char *)tmp_raypath.z,nr_inv_cell+1,sizeof(double),"InvRayStruct");
					    tmp_raypath.z[nr_inv_cell] = (raypath[i].z[j+1] - raypath[i].z[j])*grid.h; /*remove the normalization*/

						nr_inv_cell++;
						tmp_raypath.nray = nr_inv_cell;
					}
					else
					{
						tmp_raypath.len[a] =  tmp_raypath.len[a] + raypath[i].len[j]*grid.h;

						tmp_raypath.x[a] = tmp_raypath.x[a] + ((raypath[i].x[j+1] - raypath[i].x[j])*grid.h);
						tmp_raypath.y[a] = tmp_raypath.y[a] + ((raypath[i].y[j+1] - raypath[i].y[j])*grid.h);
						tmp_raypath.z[a] = tmp_raypath.z[a] + ((raypath[i].z[j+1] - raypath[i].z[j])*grid.h);

					}
					 
			 }

			 
		}

		/*transfer the parameters to the ray-structure*/
		raypath[i].nray = tmp_raypath.nray;
		if(tmp_raypath.nray != 0)
		{
			raypath[i].ele = (long *)memory((char *)raypath[i].ele,tmp_raypath.nray,sizeof(long),"InvRayStruct");
			raypath[i].len = (double *)memory((char *)raypath[i].len,tmp_raypath.nray,sizeof(double),"InvRayStruct");
			raypath[i].x = (double *)memory((char *)raypath[i].x,tmp_raypath.nray,sizeof(double),"InvRayStruct");
			raypath[i].y = (double *)memory((char *)raypath[i].y,tmp_raypath.nray,sizeof(double),"InvRayStruct");
			raypath[i].z = (double *)memory((char *)raypath[i].z,tmp_raypath.nray,sizeof(double),"InvRayStruct");
		}
		else
		{
			raypath[i].ele = (long *)memory((char *)raypath[i].ele,1,sizeof(long),"InvRayStruct");
			raypath[i].len = (double *)memory((char *)raypath[i].len,1,sizeof(double),"InvRayStruct");
			raypath[i].x = (double *)memory((char *)raypath[i].x,1,sizeof(double),"InvRayStruct");
			raypath[i].y = (double *)memory((char *)raypath[i].y,1,sizeof(double),"InvRayStruct");
			raypath[i].z = (double *)memory((char *)raypath[i].z,1,sizeof(double),"InvRayStruct");
		}

		for(k=0;k<raypath[i].nray;k++)
		{
			raypath[i].ele[k] = tmp_raypath.ele[k];
			raypath[i].len[k] = tmp_raypath.len[k];
			raypath[i].x[k] = tmp_raypath.x[k];
			raypath[i].y[k] = tmp_raypath.y[k];
			raypath[i].z[k] = tmp_raypath.z[k];
		}

		free(tmp_raypath.ele);
		free(tmp_raypath.len);
		free(tmp_raypath.x);
		free(tmp_raypath.y);
		free(tmp_raypath.z);
	}

	printf("The ray structures for the inversion grid are determined\n");
	printf("----------------\n\n");

	return(1);
}

/*------------------------------------------------------------*/
/*Determine the fat-ray-structs of the inversion cells*/
/*Parameter:	grid	:= Grid structure*/
/*				fatray	:= Input:	Structure including the fatrays in the forward grid;*/
/*						   Output:	Structure including the fatrays in the inversion grid;*/
/*							(Remark: -The shot-receiver index remains the same*/
/*									 -Fat-ray contributions  in the air or at the borders are NOT considered*/
/*				inv		:= Inversion structure */

int InvFatRayStruct(GRID_STRUCT grid,F_RP_STRUCT *fatray, INV *inv, DATA_STRUCT data)
{
	int index;
	long i,j,k,a,nr_inv_cell;
	long ix,iy,iz,nx2,ny2,nz2,nyz2;
	long inv_ix,inv_iy,inv_iz,inv_nz,inv_nyz;
	long inv_cell;
	F_RP_STRUCT tmp_fatray;

	nx2 = grid.nx + 2*grid.nborder;
	ny2 = grid.ny + 2*grid.nborder;
	nz2 = grid.nz + 2*grid.nborder;
	nyz2 =ny2*nz2;

	/*Build the inversion structure (by means of a regular forward grid)*/
	inv_nz = inv->nz; 
	inv_nyz = inv->ny *inv->nz;

	/*************************************************************/
	/*Make a fat-ray structure for the inversion cells*/
	/*Loop over all fat-rays*/
	for(i=0;i<data.ndata_seis;i++)
	{
		tmp_fatray.ele = (long *)memory(NULL,1,sizeof(long),"InvFatRayStruct");
		tmp_fatray.weight = (double *)memory(NULL,1,sizeof(double),"InvFatRayStruct");

		tmp_fatray.ncell =0;
		nr_inv_cell=0;

		/*Loop over all fat-ray contributions*/
		for(j=0;j<fatray[i].ncell;j++)
		{
			 if(grid.border_index[fatray[i].ele[j]] != 0) /*Consider just cell that are not in the air*/
			 {
					/*position in the forward grid*/
				 	ix = (long)floor((double)(fatray[i].ele[j]/nyz2));
					iy = (long)floor((double)((fatray[i].ele[j] - ix*nyz2)/nz2));
					iz = fatray[i].ele[j] - ix*nyz2 - iy*nz2;

					/*position in the inversion grid*/
					inv_ix = (long)(floor((double)(ix/grid.grid_cell_factor[0])));
					inv_iy = (long)(floor((double)(iy/grid.grid_cell_factor[1])));
					inv_iz = (long)(floor((double)(iz/grid.grid_cell_factor[2])));

					inv_cell = inv_nyz*inv_ix + inv_nz*inv_iy + inv_iz;


					if(inv_cell > inv->nvel)
					{
						printf("The inversion cell number %d does not exist;\nThe highest inversion cell number is %d\n",inv_cell,inv->nvel);
						exit(0);
					}

					index = 0;

					/*Check, if fat-ray runs already through the inversion cell*/
					for(a=0;a<nr_inv_cell;a++)
					{
						if(tmp_fatray.ele[a] == inv_cell)
						{
							index = 1;
							break;
						}
					}

					/*Build the fat-ray structure for the inversion cells*/
					if(index == 0)
					{
						tmp_fatray.ele = (long *)memory((char *)tmp_fatray.ele,nr_inv_cell+1,sizeof(long),"InvFatRayStruct");
						tmp_fatray.ele[nr_inv_cell] = inv_cell;

						tmp_fatray.weight = (double *)memory((char *)tmp_fatray.weight,nr_inv_cell+1,sizeof(double),"InvFatRayStruct");
						tmp_fatray.weight[nr_inv_cell] = fatray[i].weight[j];

						nr_inv_cell++;
						tmp_fatray.ncell = nr_inv_cell;
					}
					else
					{
						tmp_fatray.weight[a] =  tmp_fatray.weight[a] + fatray[i].weight[j];
					}
					 
			 }
			 
		}

		/*transfer the parameters to the fat-ray-structure*/

		fatray[i].ncell = tmp_fatray.ncell;
		if(tmp_fatray.ncell != 0)
		{
			fatray[i].ele = (long *)memory((char *)fatray[i].ele,tmp_fatray.ncell,sizeof(long),"InvFatRayStruct");
			fatray[i].weight = (double *)memory((char *)fatray[i].weight,tmp_fatray.ncell,sizeof(double),"InvFatRayStruct");
		}
		else
		{
			fatray[i].ele = (long *)memory((char *)fatray[i].ele,1,sizeof(long),"InvFatRayStruct");
			fatray[i].weight = (double *)memory((char *)fatray[i].weight,1,sizeof(double),"InvFatRayStruct");
		}

		for(k=0;k<fatray[i].ncell;k++)
		{
			fatray[i].ele[k] = tmp_fatray.ele[k];
			fatray[i].weight[k] = tmp_fatray.weight[k];
		}

		free(tmp_fatray.ele);
		free(tmp_fatray.weight);

		if(i%100 == 1)
		printf("%d of %d fat-ray structures are determined\n",i,data.ndata_seis);
	}

	printf("The fat-ray structures for the inversion grid are determined\n");
	printf("----------------\n\n");
	return(1);
}


/*------------------------------------------------------------*/
/*Read in the spatially dependent inversion grid structure*/
/*Parameter: *fname := Filename */
/*Output			:= The structure governing the spatially dependency of the inversion grid*/

FLEX_STRUCT ReadFlexGridIn(char *fname)
{
	long ix,iy,iz,j,i;
	char direc[1];
	char line[128];
	FLEX_STRUCT flex;
	FILE *inf;

   inf = fopen(fname,"rt");
   if (inf == NULL)
   {
      fprintf(stderr,"Unable to open %s\n",fname);
      exit(0);
   }

   	printf("Readjust the inversion structures (spatially invariant)\n");
	printf("----------------\n");


   /*Read in parameter file*/
	flex.nx = 0;
	flex.ny = 0;
	flex.nz = 0;
	flex.xpos = (double *)memory(NULL,1,sizeof(double),"ReadFlexGridIn");
	flex.ypos = (double *)memory(NULL,1,sizeof(double),"ReadFlexGridIn");
	flex.zpos = (double *)memory(NULL,1,sizeof(double),"ReadFlexGridIn");
	for(i=0;i<3;i++)
	{
		flex.xratio[i] = (int *)memory(NULL,1,sizeof(int),"ReadFlexGridIn");
		flex.yratio[i] = (int *)memory(NULL,1,sizeof(int),"ReadFlexGridIn");
		flex.zratio[i] = (int *)memory(NULL,1,sizeof(int),"ReadFlexGridIn");
	}

	ix=0,iy=0,iz=0;
	GetNext(line,inf);
	/*Loop over all bodies*/
	while('#' != line[0])
	{
		for(j=0;!isgraph(line[j]);j++);
		memcpy(direc, line+j, 1);

		/*Case: x-direction:*/
		if(direc[0] == 'x')
		{
			flex.nx++;
			flex.xpos = (double *)memory((char *)flex.xpos,ix+1,sizeof(double),"ReadFlexGridIn");
			for(i=0;i<3;i++)
				flex.xratio[i] = (int *)memory((char *)flex.xratio[i],ix+1,sizeof(int),"ReadFlexGridIn");
			 
			fgets(line,128,inf);
			sscanf(line,"%lf\n", &(flex.xpos[ix]));
			fgets(line,128,inf);
			sscanf(line,"%d %d %d\n", &(flex.xratio[0][ix]),&(flex.xratio[1][ix]),&(flex.xratio[2][ix]));
			fgets(line,128,inf);

			for(i=0;i<ix;i++)
			{
				if(flex.xpos[ix] < flex.xpos[i])
				{
					printf("The sorting in x-direction in the parameter file %s is not correct\n",fname);
					exit(0);
				}
			}

			ix++;
		}
		/*Case: y-direction:*/
		else if(direc[0] == 'y')
		{
			flex.ny++;
			flex.ypos = (double *)memory((char *)flex.ypos,iy+1,sizeof(double),"ReadFlexGridIn");
			for(i=0;i<3;i++)
				flex.yratio[i] = (int *)memory((char *)flex.yratio[i],iy+1,sizeof(int),"ReadFlexGridIn");

			fgets(line,128,inf);
			sscanf(line,"%lf\n", &(flex.ypos[iy]));
			fgets(line,128,inf);
			sscanf(line,"%d %d %d\n", &(flex.yratio[0][iy]),&(flex.yratio[1][iy]),&(flex.yratio[2][iy]));
			fgets(line,128,inf);

			for(i=0;i<iy;i++)
			{
				if(flex.ypos[iy] < flex.ypos[i])
				{
					printf("The sorting in y-direction in the parameter file %s is not correct\n",fname);
					exit(0);
				}
			}

			iy++;
		}
		/*Case: z-direction:*/
		else if(direc[0] == 'z')
		{
			flex.nz++;
			flex.zpos = (double *)memory((char *)flex.zpos,iz+1,sizeof(double),"ReadFlexGridIn");
			for(i=0;i<3;i++)
				flex.zratio[i] = (int *)memory((char *)flex.zratio[i],iz+1,sizeof(int),"ReadFlexGridIn");

			fgets(line,128,inf);
			sscanf(line,"%lf\n", &(flex.zpos[iz]));
			fgets(line,128,inf);
			sscanf(line,"%d %d %d\n", &(flex.zratio[0][iz]),&(flex.zratio[1][iz]),&(flex.zratio[2][iz]));
			fgets(line,128,inf);

			for(i=0;i<iz;i++)
			{
				if(flex.zpos[iz] < flex.zpos[i])
				{
					printf("The sorting in z-direction in the parameter file %s is not correct\n",fname);
					exit(0);
				}
			}

			iz++;
		}
		else
		{
			printf("Specification %s is unknown\n",direc[0]);
			exit(0);
		}

	}

	fclose(inf);

	printf("The parameter file for the flexible inversion grid %s is read in\n",fname);
	printf("----------------\n");

   return(flex);
}


/*------------------------------------------------------------*/
/*Readjust the inversion grid structure to a spatially variant grid*/
/*Parameter: flex		:= Structure comprising the parameters of the spatially variant grid*/
/*			*inv_cell	:= Structure of the inversion cells (will be modified by this procedure)*/
/*			grid		:= Parameters of the forward grid*/
/*			inv			:= Inversion grid*/
/*			flag		:= Structures including the flags*/

int MakeIrregularStructuresI(FLEX_STRUCT flex,BIDX_STRUCT *inv_cell,GRID_STRUCT grid, INV *inv, FLAG_STRUCT flag)
{
	int i,j,k,last_i;
	int	xscal,yscal,zscal,flx_int;
	int nr_x,nr_y,nr_z,ix,iy,iz, max_nr_forw;
	int *tmp_nele_inv;
	long **tmp_ele_inv;
	long nr_inv_cell,nr_inv_cell_use,npos;
	float xdim,ydim,zdim,resx,resy,resz;
	double grid_min[3],grid_max[3];

	BIDX_STRUCT *tmp_inv_cell;

	tmp_inv_cell = (BIDX_STRUCT *)memory(NULL,1,sizeof(BIDX_STRUCT),"MakeIrregularStructures");
	/********************************************/

	/*Changes in z-direction*/
	/*REMARK: Variability in the x- and y-direction are not implemented until now*/

	/*Make a new grid structure for the inversion grid:*/
	nr_inv_cell = 0;
	last_i = 0;

	grid_min[2] = grid.org[2] - grid.nborder*grid.h;

	if(flex.nz == 0 || flex.zpos[0] > grid.org[2] + grid.h* (grid.nborder + grid.nz -1))
	{
		grid_max[2] = grid.org[2] + grid.h* (grid.nborder + grid.nz -1);
		last_i = 1;
	}
	else
	{
		flx_int = (int)((flex.zpos[0] - grid.org[2])/grid.h);
		grid_max[2] = flx_int*grid.h + grid.org[2];
	}

	xscal = 1;
	yscal = 1;
	zscal = 1;

	for(i=0;i<flex.nz+1;i++)
	{

		/*Cell size of the new inversion cells*/
	    xdim = xscal * grid.h * grid.grid_cell_factor[0];
		ydim = yscal * grid.h * grid.grid_cell_factor[1];
		zdim = zscal * grid.h * grid.grid_cell_factor[2];

		/*Number of cells in x,y,z direction*/
		nr_x = (int)ceil((double)(((2*grid.nborder + grid.nx)* grid.h)/xdim));
		nr_y = (int)ceil((double)(((2*grid.nborder + grid.ny)* grid.h)/ydim));
		nr_z = (int)ceil((double)((grid_max[2]-grid_min[2]+grid.h)/zdim));
		if(nr_z < 0)
			nr_z=0;

		/*Grid size at the upper boundaries*/
		resx = nr_x - (((2*grid.nborder + grid.nx)* grid.h)/xdim);
		resy = nr_y - (((2*grid.nborder + grid.ny)* grid.h)/ydim);
		resz = nr_z - ((float)(grid_max[2]-grid_min[2]+grid.h)/zdim);

			resx=(1-resx)*xdim;
			resy=(1-resy)*ydim;
			resz=(1-resz)*zdim;

		/*Calculate the positions of the new inversion cells:*/
	    for(ix=0;ix<nr_x;ix++)
		{
			for(iy=0;iy<nr_y;iy++)
			{
				for(iz=0;iz<nr_z;iz++)
				{
					tmp_inv_cell = (BIDX_STRUCT *)memory((char *)tmp_inv_cell,nr_inv_cell+1,sizeof(BIDX_STRUCT),"MakeIrregularStructures");			
					tmp_inv_cell[nr_inv_cell].ele = (long *)memory(NULL,1,sizeof(long),"MakeIrregularStructures");

					tmp_inv_cell[nr_inv_cell].nele = 0;
					tmp_inv_cell[nr_inv_cell].use = 0;
					tmp_inv_cell[nr_inv_cell].val_slow = 0.0;
					tmp_inv_cell[nr_inv_cell].val_dens = 0.0;
					tmp_inv_cell[nr_inv_cell].val_res = 0.0;
					tmp_inv_cell[nr_inv_cell].bdamp = 0.0;

					/*Jin, 06, July, 2006.*/
					tmp_inv_cell[nr_inv_cell].border_inv_cell_x_left = 0;
					tmp_inv_cell[nr_inv_cell].border_inv_cell_x_right = 0;
					tmp_inv_cell[nr_inv_cell].border_inv_cell_y_front = 0;
					tmp_inv_cell[nr_inv_cell].border_inv_cell_y_back = 0;
					tmp_inv_cell[nr_inv_cell].border_inv_cell_z = 0;

					if(ix +1 != nr_x)
					{
						tmp_inv_cell[nr_inv_cell].xo = (grid.org[0] - (0.5+grid.nborder )*grid.h) + (xdim/2) + ix*xdim; /*coordinates of the new inversion cells*/
						tmp_inv_cell[nr_inv_cell].xdim = xdim;
					}
					else /*Cells at the upper x-border*/
					{
						tmp_inv_cell[nr_inv_cell].xo = (grid.org[0] - (0.5+grid.nborder)*grid.h) + (resx/2) + ix*xdim;
						tmp_inv_cell[nr_inv_cell].xdim = resx;
					}

					if(iy +1 != nr_y)
					{
						tmp_inv_cell[nr_inv_cell].yo = (grid.org[1] - (0.5+grid.nborder)*grid.h) + (ydim/2) + iy*ydim;
						tmp_inv_cell[nr_inv_cell].ydim = ydim;
					}
					else /*Cells at the upper y-border*/
					{
						tmp_inv_cell[nr_inv_cell].yo = (grid.org[1] - (0.5+grid.nborder)*grid.h) + (resy/2) + iy*ydim;
						tmp_inv_cell[nr_inv_cell].ydim = resy;
					}

					if(iz +1 != nr_z || last_i != 1)
					{
						tmp_inv_cell[nr_inv_cell].zo = grid_min[2] - (0.5*grid.h) + (zdim/2) + iz*zdim;
						tmp_inv_cell[nr_inv_cell].zdim = zdim;
					}
					else /*Cells at the upper z-border*/
					{
						tmp_inv_cell[nr_inv_cell].zo = grid_min[2] - (0.5*grid.h) + (resz/2) + iz*zdim;
						tmp_inv_cell[nr_inv_cell].zdim = resz;
					}

					nr_inv_cell++; /*Number of new inversion cells*/
				}				

			}

		}

		/*Determine the lower and upper boundary for readjustments of the inversion cells*/

		if(flex.nz != i)
		{
			grid_min[2] = grid_min[2] + zdim*nr_z;

			if(i+1 == flex.nz || flex.zpos[i] > grid.org[2] + grid.h* (grid.nborder + grid.nz -1))
			{
				grid_max[2] = grid.org[2] + grid.h* (grid.nborder + grid.nz -1); /*Lower z-boundary of the forward grid*/
				last_i = 1;
			}
			else
			{
				flx_int = (int)((flex.zpos[i+1] - grid.org[2])/grid.h);
				grid_max[2] = flx_int*grid.h + grid.org[2];
			}

			xscal = flex.zratio[0][i];
			yscal = flex.zratio[1][i];
			zscal = flex.zratio[2][i];
		}


	}
	/********************************************/
	/*Assign the former inversion cells to the modified inversion cells*/

	max_nr_forw = 0; /*max. number of forw. cells in a modified inversion cell*/

	tmp_nele_inv = (int *)memory(NULL,nr_inv_cell,sizeof(int),"MakeIrregularStructures");
	tmp_ele_inv = (long **)memory(NULL,nr_inv_cell,sizeof(long *),"MakeIrregularStructures");

	for(i=0;i<nr_inv_cell;i++)
	{
		tmp_nele_inv[i] = 0;	/*Number of former inversion cells in the modified inversion cell*/
		tmp_ele_inv[i] = (long *)memory(NULL,1,sizeof(long),"MakeIrregularStructures"); /*and the corresponding indeces*/
	}



	/*Loop over the former inversion cells*/
	for(i=0;i<inv->nvel;i++)
	{
		/*Loop over the modified inversion cells*/
		for(j=0;j<nr_inv_cell;j++)
		{
		
			/*The former inversion cell belongs to the modified inversion cell*/
			if(inv_cell[i].xo >= tmp_inv_cell[j].xo -0.5*tmp_inv_cell[j].xdim && inv_cell[i].xo <= tmp_inv_cell[j].xo + 0.5*tmp_inv_cell[j].xdim &&
			   inv_cell[i].yo >= tmp_inv_cell[j].yo -0.5*tmp_inv_cell[j].ydim && inv_cell[i].yo <= tmp_inv_cell[j].yo + 0.5*tmp_inv_cell[j].ydim &&	
			   inv_cell[i].zo >= tmp_inv_cell[j].zo -0.5*tmp_inv_cell[j].zdim && inv_cell[i].zo <= tmp_inv_cell[j].zo + 0.5*tmp_inv_cell[j].zdim)
			{
				if(inv_cell[i].use == 1) /*Activate the modified inversion cell*/
				{
					tmp_inv_cell[j].use = 1;

					if(flag.index_tseis != 0 || (flag.index_grav != 0 && flag.index_mt != 0))
						tmp_inv_cell[j].val_slow = (inv_cell[i].nele*inv_cell[i].val_slow) + tmp_inv_cell[j].val_slow;			/*Slowness*/

					if(flag.index_grav != 0)
						tmp_inv_cell[j].val_dens = (inv_cell[i].nele*inv_cell[i].val_dens) + tmp_inv_cell[j].val_dens;			/*Density*/

					if(flag.index_mt != 0)
						tmp_inv_cell[j].val_res = (inv_cell[i].nele*inv_cell[i].val_res) + tmp_inv_cell[j].val_res;			/*Density*/

					tmp_inv_cell[j].bdamp = (inv_cell[i].nele*inv_cell[i].bdamp) + tmp_inv_cell[j].bdamp;	/*Scaling of the damping*/

					/*Changed by Jin 07.07.06*/
					if(inv_cell[i].border_inv_cell_x_left != 0)
						tmp_inv_cell[j].border_inv_cell_x_left = inv_cell[i].border_inv_cell_x_left;
					if(inv_cell[i].border_inv_cell_x_right != 0)
						tmp_inv_cell[j].border_inv_cell_x_right = inv_cell[i].border_inv_cell_x_right;

					if(inv_cell[i].border_inv_cell_y_front != 0)
						tmp_inv_cell[j].border_inv_cell_y_front = inv_cell[i].border_inv_cell_y_front;
					if(inv_cell[i].border_inv_cell_y_back != 0)
						tmp_inv_cell[j].border_inv_cell_y_back = inv_cell[i].border_inv_cell_y_back;

					if(inv_cell[i].border_inv_cell_z != 0)
						tmp_inv_cell[j].border_inv_cell_z = inv_cell[i].border_inv_cell_z;

				}

				tmp_nele_inv[j]++;	/*Number of former inversion cells in the modified inversion cell*/
				tmp_ele_inv[j] = (long *)memory((char *)tmp_ele_inv[j],tmp_nele_inv[j],sizeof(long),"MakeIrregularStructures");
				tmp_ele_inv[j][tmp_nele_inv[j]-1] = i;

				npos = tmp_inv_cell[j].nele;
				tmp_inv_cell[j].nele = tmp_inv_cell[j].nele + inv_cell[i].nele; /*Number of forward cells in the modified inversion cells*/

				if(tmp_inv_cell[j].nele > max_nr_forw)
					max_nr_forw = tmp_inv_cell[j].nele; /*Max. number of forward cells in a inversion cell*/ 

				if(tmp_inv_cell[j].nele != 0)
				{
					tmp_inv_cell[j].ele = (long *)memory((char *)tmp_inv_cell[j].ele,tmp_inv_cell[j].nele,sizeof(long),"MakeIrregularStructures");
					memcpy(&(tmp_inv_cell[j].ele[npos]),&(inv_cell[i].ele[0]),inv_cell[i].nele*sizeof(long));
				}

				break;

			}

		}
	}

	/*Modify the number of active inversion cells*/
	nr_inv_cell_use = 0;

	for(i=0;i<nr_inv_cell;i++)
	{
		if(tmp_inv_cell[i].use ==1)
		{
			if(tmp_inv_cell[i].nele <= 0)
			{
				printf("ERROR! The number of forward cells in the inversion cell %d is unsuited\n", i);
				exit(0);
			}

			if(flag.index_tseis != 0 || (flag.index_grav != 0 && flag.index_mt != 0))
				tmp_inv_cell[i].val_slow = tmp_inv_cell[i].val_slow/tmp_inv_cell[i].nele;
			else
				tmp_inv_cell[i].val_slow = 0.0;

			if(flag.index_grav != 0)
				tmp_inv_cell[i].val_dens = tmp_inv_cell[i].val_dens/tmp_inv_cell[i].nele;
			else
				tmp_inv_cell[i].val_dens = 0.0;

			if(flag.index_mt != 0)
				tmp_inv_cell[i].val_res = tmp_inv_cell[i].val_res/tmp_inv_cell[i].nele;
			else
				tmp_inv_cell[i].val_res = 0.0;

			tmp_inv_cell[i].bdamp = tmp_inv_cell[i].bdamp/tmp_inv_cell[i].nele; /* I: Weighting of the inversion cell is independent of its size (number of act forward cells in considered inv. cell/number of all forward cells in considered inv. cell)*/
		//	tmp_inv_cell[i].bdamp = tmp_inv_cell[i].bdamp/max_nr_forw;		    /*II: Weighting of the inversion cell is dependent of its size: (number of act. forward cells in considered inv. cell/ max. number of act. forward cells in a inv. cell)*/

			tmp_inv_cell[i].used_nr = nr_inv_cell_use;
			nr_inv_cell_use++;
		}
	}

	printf("The inversion cell structure is readjusted:\n");
	printf("Number of inversion cells: %d\nNumber of active inversion cells: %d\n", nr_inv_cell, nr_inv_cell_use);
	printf("----------------\n");


	/********************************************/
	/*Copy the modified into the former structures*/

	/*Grid structure*/
	for(i=nr_inv_cell;i<inv->nvel;i++)
	{
		free(inv_cell[i].ele);
		free(inv_cell[i].ele_inv);
	}

//	inv_cell = (BIDX_STRUCT *)memory((char *)inv_cell,nr_inv_cell,sizeof(BIDX_STRUCT),"MakeIrregularStructures");	

	for(i=0;i<nr_inv_cell;i++)
	{
		inv_cell[i].bdamp = tmp_inv_cell[i].bdamp;
		inv_cell[i].val_slow =	tmp_inv_cell[i].val_slow;
		inv_cell[i].val_dens = tmp_inv_cell[i].val_dens;
		inv_cell[i].val_res = tmp_inv_cell[i].val_res;
		inv_cell[i].xo =	tmp_inv_cell[i].xo;
		inv_cell[i].yo =	tmp_inv_cell[i].yo;
		inv_cell[i].zo =	tmp_inv_cell[i].zo;
		inv_cell[i].xdim =	tmp_inv_cell[i].xdim;
		inv_cell[i].ydim =	tmp_inv_cell[i].ydim;
		inv_cell[i].zdim =	tmp_inv_cell[i].zdim;
		inv_cell[i].use =	tmp_inv_cell[i].use;
		inv_cell[i].used_nr=tmp_inv_cell[i].used_nr;
		inv_cell[i].nele =  tmp_inv_cell[i].nele;
		inv_cell[i].nele_inv = tmp_nele_inv[i];
		
		if(inv_cell[i].nele >0)
		{
			inv_cell[i].ele = (long *)memory((char *)inv_cell[i].ele, inv_cell[i].nele, sizeof(long),"MakeIrregularStructures");
			for(k=0;k<inv_cell[i].nele;k++)
			{
				inv_cell[i].ele[k] = tmp_inv_cell[i].ele[k];
			}
		}
		else
			inv_cell[i].ele = (long *)memory((char *)inv_cell[i].ele, 1, sizeof(long),"MakeIrregularStructures");

		if(inv_cell[i].nele_inv >0)
		{
			inv_cell[i].ele_inv = (long *)memory((char *)inv_cell[i].ele_inv, inv_cell[i].nele_inv, sizeof(long),"MakeIrregularStructures");
	
			for(k=0;k<inv_cell[i].nele_inv;k++)
			{
				inv_cell[i].ele_inv[k] = tmp_ele_inv[i][k];
			}
		}
		else
			inv_cell[i].ele_inv = (long *)memory((char *)inv_cell[i].ele_inv, 1, sizeof(long),"MakeIrregularStructures");

	}

	inv->nvel = nr_inv_cell;
	inv->nvel_used  = nr_inv_cell_use;

	/********************************************/
	/*free the memory */
	for(i=0;i<nr_inv_cell;i++)
		free(tmp_inv_cell[i].ele);

	free(tmp_inv_cell);

	for(i=0;i<nr_inv_cell;i++)
		free(tmp_ele_inv[i]);

	free(tmp_nele_inv);
	free(tmp_ele_inv);

	return(0);
}


/*------------------------------------------------------------*/
/*Readjust the raypath/fatray structure to a spatially variant grid*/
/*Parameter: *raypath	:= Raypath structure (will be modified by this procedure)*/
/*			*fatray		:= Fatray structure (will be modified by this routine)*/
/*			nr_inv_cell	:= Number of inversion cells*/
/*			*inv_cell	:= Structure of the inversion cells (will be modified by this procedure)*/
/*			grid		:= Parameters of the forward grid*/
/*			inv			:= Inversion grid*/
/*			nrays		:= Number of rays*/
/*			kind_of_rays:= 1== conventional rays; 2== fat rays*/

int MakeIrregularRayStructures(RP_STRUCT *raypath,F_RP_STRUCT *fatray, BIDX_STRUCT *inv_cell, long nr_inv_cell, long nrays, int kind_of_rays)
{
	int i,j,k,l,m;
	int index,index2;
	double laenge;

	RP_STRUCT *tmp_raypath;
	F_RP_STRUCT *tmp_fatray;

	/********************************************/
/*Modify the raypath structure*/
	if(kind_of_rays == 1)
	{
		/*conventional rays*/
		if(nrays != 0)
			tmp_raypath = (RP_STRUCT *)memory(NULL,nrays,sizeof(RP_STRUCT),"MakeIrregularStructures");
		else
			tmp_raypath = (RP_STRUCT *)memory(NULL,1,sizeof(RP_STRUCT),"MakeIrregularStructures");

		/*Loop over all rays*/
		for(i=0;i<nrays;i++)
		{
			tmp_raypath[i].nray = 0;
			tmp_raypath[i].ele = (long *)memory(NULL,1, sizeof(long), "MakeIrregularStructures");
			tmp_raypath[i].len = (double *)memory(NULL,1, sizeof(double), "MakeIrregularStructures");
			tmp_raypath[i].x = (double *)memory(NULL,1, sizeof(double), "MakeIrregularStructures");
			tmp_raypath[i].y = (double *)memory(NULL,1, sizeof(double), "MakeIrregularStructures");
			tmp_raypath[i].z = (double *)memory(NULL,1, sizeof(double), "MakeIrregularStructures");

			for(j=0;j<raypath[i].nray;j++)
			{
			
				/*Loop over all modified inversion cells*/
				for(k=0;k<nr_inv_cell;k++)
				{
					laenge = 0.0;
					index=0;

					for(l=0;l<inv_cell[k].nele_inv;l++)
					{
						if(inv_cell[k].ele_inv[l]==raypath[i].ele[j])
						{
							laenge = raypath[i].len[j];
							index=1;
							break;
						}
					}

					if(index==1)
					{
						index2 = 0; /*Check if the modified inversion cell was already touched by the ray*/
						for(m=0;m<tmp_raypath[i].nray;m++)
						{
							if(k == tmp_raypath[i].ele[m])
							{
								index2 = 1;
								break;
							}
						}

						if(index2 == 0)
						{
							tmp_raypath[i].nray++;
							tmp_raypath[i].ele = (long *)memory((char *)tmp_raypath[i].ele,tmp_raypath[i].nray, sizeof(long), "MakeIrregularStructures");
							tmp_raypath[i].ele[tmp_raypath[i].nray - 1] = k;
							tmp_raypath[i].len = (double *)memory((char *)tmp_raypath[i].len,tmp_raypath[i].nray, sizeof(double), "MakeIrregularStructures");
							tmp_raypath[i].len[tmp_raypath[i].nray - 1] = laenge; /*Length of ray segment in the modified inversion cell*/

							tmp_raypath[i].x = (double *)memory((char *)tmp_raypath[i].x,tmp_raypath[i].nray, sizeof(double), "MakeIrregularStructures");
							tmp_raypath[i].x[tmp_raypath[i].nray - 1] = raypath[i].x[j];
							tmp_raypath[i].y = (double *)memory((char *)tmp_raypath[i].y,tmp_raypath[i].nray, sizeof(double), "MakeIrregularStructures");
							tmp_raypath[i].y[tmp_raypath[i].nray - 1] = raypath[i].y[j];
							tmp_raypath[i].z = (double *)memory((char *)tmp_raypath[i].z,tmp_raypath[i].nray, sizeof(double), "MakeIrregularStructures");
							tmp_raypath[i].z[tmp_raypath[i].nray - 1] = raypath[i].z[j];
						}
						else
						{
							tmp_raypath[i].len[m] = tmp_raypath[i].len[m] + laenge;

							tmp_raypath[i].x[m] = tmp_raypath[i].x[m] + raypath[i].x[j];
							tmp_raypath[i].y[m] = tmp_raypath[i].y[m] + raypath[i].y[j];
							tmp_raypath[i].z[m] = tmp_raypath[i].z[m] + raypath[i].z[j];
						}

						break;
					}

				}
			}

		}

	}
	else
	{
		/*Fat-rays*/
		tmp_fatray = (F_RP_STRUCT *)memory(NULL,nrays,sizeof(F_RP_STRUCT),"MakeIrregularStructures");

		/*Loop over all rays*/
		for(i=0;i<nrays;i++)
		{
			tmp_fatray[i].ncell = 0;
			tmp_fatray[i].ele = (long *)memory(NULL,1, sizeof(long), "MakeIrregularStructures");
			tmp_fatray[i].weight = (double *)memory(NULL,1, sizeof(double), "MakeIrregularStructures");

			for(j=0;j<fatray[i].ncell;j++)
			{
			
				/*Loop over all modified inversion cells*/
				for(k=0;k<nr_inv_cell;k++)
				{
					laenge = 0.0;
					index=0;

					for(l=0;l<inv_cell[k].nele_inv;l++)
					{
						if(inv_cell[k].ele_inv[l]==fatray[i].ele[j])
						{
							laenge = fatray[i].weight[j];
							index=1;
							break;
						}
					}
					if(index==1)
					{
						index2 = 0; /*Check if the modified inversion cell was already touched by the ray*/
						for(m=0;m<tmp_fatray[i].ncell;m++)
						{
							if(k == tmp_fatray[i].ele[m])
							{
								index2 = 1;
								break;
							}
						}

						if(index2 == 0)
						{
							tmp_fatray[i].ncell++;
							tmp_fatray[i].ele = (long *)memory((char *)tmp_fatray[i].ele,tmp_fatray[i].ncell, sizeof(long), "MakeIrregularStructures");
							tmp_fatray[i].ele[tmp_fatray[i].ncell - 1] = k;
							tmp_fatray[i].weight = (double *)memory((char *)tmp_fatray[i].weight,tmp_fatray[i].ncell, sizeof(double), "MakeIrregularStructures");
							tmp_fatray[i].weight[tmp_fatray[i].ncell - 1] = laenge; /*Length of ray segment in the modified inversion cell*/
						}
						else
						{
							tmp_fatray[i].weight[m] = tmp_fatray[i].weight[m] + laenge;
						}

						break;
					}

				}
			}

		}
	
	}


	/********************************************/
	/*Copy the modified into the former structures*/

	if(kind_of_rays == 1)
	{
		/*Ray structure*/
		/*Conventional rays*/
		for(i=0;i<nrays;i++)
		{
			raypath[i].nray = tmp_raypath[i].nray;
			if(raypath[i].nray > 0)
			{
				raypath[i].ele  = (long *)memory((char *)raypath[i].ele,tmp_raypath[i].nray,sizeof(long), "MakeIrregularStructures");
				raypath[i].len = (double *)memory((char *)raypath[i].len,tmp_raypath[i].nray,sizeof(double), "MakeIrregularStructures");
				raypath[i].x = (double *)memory((char *)raypath[i].x,tmp_raypath[i].nray,sizeof(double), "MakeIrregularStructures");
				raypath[i].y = (double *)memory((char *)raypath[i].y,tmp_raypath[i].nray,sizeof(double), "MakeIrregularStructures");
				raypath[i].z = (double *)memory((char *)raypath[i].z,tmp_raypath[i].nray,sizeof(double), "MakeIrregularStructures");

			}
			else
			{
				raypath[i].ele  = (long *)memory((char *)raypath[i].ele,1,sizeof(long), "MakeIrregularStructures");
				raypath[i].len = (double *)memory((char *)raypath[i].len,1,sizeof(double), "MakeIrregularStructures");
				raypath[i].x = (double *)memory((char *)raypath[i].x,1,sizeof(double), "MakeIrregularStructures");
				raypath[i].y = (double *)memory((char *)raypath[i].y,1,sizeof(double), "MakeIrregularStructures");
				raypath[i].z = (double *)memory((char *)raypath[i].z,1,sizeof(double), "MakeIrregularStructures");
			}

			for(k=0;k<raypath[i].nray;k++)
			{
				raypath[i].ele[k]=tmp_raypath[i].ele[k];
				raypath[i].len[k]=tmp_raypath[i].len[k];
				raypath[i].x[k]=tmp_raypath[i].x[k];
				raypath[i].y[k]=tmp_raypath[i].y[k];
				raypath[i].z[k]=tmp_raypath[i].z[k];
			}
		}
	}
	else
	{
		/*Ray structure*/
		/*Fatrays*/
		for(i=0;i<nrays;i++)
		{
			fatray[i].ncell = tmp_fatray[i].ncell;
			if(fatray[i].ncell > 0)
			{
				fatray[i].ele  = (long *)memory((char *)fatray[i].ele,tmp_fatray[i].ncell,sizeof(long), "MakeIrregularStructures");
				fatray[i].weight = (double *)memory((char *)fatray[i].weight,tmp_fatray[i].ncell,sizeof(double), "MakeIrregularStructures");
			}
			else
			{
				fatray[i].ele  = (long *)memory((char *)fatray[i].ele,1,sizeof(long), "MakeIrregularStructures");
				fatray[i].weight = (double *)memory((char *)fatray[i].weight,1,sizeof(double), "MakeIrregularStructures");
			}

			for(k=0;k<fatray[i].ncell;k++)
			{
				fatray[i].ele[k]=tmp_fatray[i].ele[k];
				fatray[i].weight[k]=tmp_fatray[i].weight[k];
			}
		}
	}

	printf("The raypath structure is readjusted\n");
	printf("----------------\n\n");

	/********************************************/
	/*free the memory */

	if(kind_of_rays == 1)
	{
		for(i=0;i<nrays;i++)
		{
			free(tmp_raypath[i].ele);
			free(tmp_raypath[i].len);
			free(tmp_raypath[i].x);
			free(tmp_raypath[i].y);
			free(tmp_raypath[i].z);
		}
		free(tmp_raypath);
	}
	else
	{
		for(i=0;i<nrays;i++)
		{
			free(tmp_fatray[i].ele);
			free(tmp_fatray[i].weight);
		}
		free(tmp_fatray);
	}

	return(1);
}



/*********************************************************************************/
/*Assign the calculated derivatives to the corresponding gravity structure (corresponds to determing the frechets)*/
/*Parameter:	geo  := Geometry structure  */
/*              *data  := Pointer on data structure*/
/*				grav := Gravity structure including the frechets*/
/*				inv_cell := Inversion cell structure*/
/*				nr_inv_cell := Number of all inversion cells*/
/*				nr_inv_cell_used := Number of the used inversion cells*/

int InvGravStruct(GEOMETRY geo, DATA_STRUCT *data, GRAV_STRUCT *grav, BIDX_STRUCT *inv_cell, long nr_inv_cell_used, long nr_inv_cell)
{
	long   a,b,c,i,j,k,n; 
	long   nx,ny,nz,nborder,nz2,nyz2,nxyz2;	/*Nr of samples*/
	long   counter;
	double tmp_para;
	double *f_deriv;

	char fname[40];
	FILE *inf;

	/*Assign the calculated derivatives to the corresponding gravity structure*/

	printf("Start assigning the calculated derivatives to the corresponding \ngravity structure\n");
	printf("----------------\n");

	/*Loop over all stations*/
	for(i=0;i<geo.nstat_grav;i++)
	{

		/*************************************************************/
		/*Read in temporary files containing the gravity derivavtives*/
		sprintf(fname,"tmp_grav_deriv%d.dat",i);

		inf = fopen(fname,"rb");

		if (inf == NULL)
		{
			fprintf(stderr,"Unable to open %s\n",fname);
			exit(0);
		}

		fread(&(nx),sizeof(long),1,inf);
		fread(&(ny),sizeof(long),1,inf);
		fread(&(nz),sizeof(long),1,inf);
		fread(&(nborder),sizeof(long),1,inf);

		nz2 = (nz+2*nborder);
		nyz2 = (ny+2*nborder)*(nz+2*nborder);
		nxyz2 = (nx+2*nborder)*(ny+2*nborder)*(nz+2*nborder);

		f_deriv = (double *)memory(NULL,(nx+2*nborder)*(ny+2*nborder)*(nz+2*nborder),sizeof(double),"InvGravStruc");

		for(a=0;a<nxyz2;a++)
			f_deriv[a] = 0.0;

		for(a=0;a<nx;a++)
			for(b=0;b<ny;b++)
				for(c=0;c<nz;c++)
				{
					fread(&(tmp_para),sizeof(double),1,inf);					/*Read in the derivatives*/
					f_deriv[(a+nborder)*nyz2 + (b+nborder)*nz2 + (c+nborder)] = (double)tmp_para;
				}
	
		fclose(inf);
		/*Delete the temporary files*/
		remove(fname);
		/*************************************************************/

		for(j=0;j<data->ndata_grav;j++)
		{
			
			if(data->gravs[i] == data->gno[j])
			{
				/*Index starts at 0*/
				grav[j].n = j; 
				grav[j].ncell = nr_inv_cell_used;

				grav[j].ele = (long *)memory(NULL, nr_inv_cell_used, sizeof(long),"DerivativesGrav");
				grav[j].deriv = (double *)memory(NULL, nr_inv_cell_used, sizeof(double),"DerivativesGrav");

				for(k=0;k<nr_inv_cell_used;k++)
				{
					grav[j].ele[k] = -1;
					grav[j].deriv[k] = 0.0; 
				}

				counter = 0;

				/*Loop over all inversion cells*/
				for(k=0;k<nr_inv_cell;k++)
				{
					if(inv_cell[k].use == 1)
					{
						grav[j].ele[counter] = k;

						/*Loop over all forward cells in the inversion cells*/
						for(n=0;n<inv_cell[k].nele;n++)
						{
							grav[j].deriv[counter] = grav[j].deriv[counter] + f_deriv[inv_cell[k].ele[n]];	
						}
						counter++;
					}
				}
			}
		}


		free(f_deriv);
	}

	printf("The gravity derivatives are determined\n");
	printf("----------------\n\n");

	return(0);
}

/*********************************************************************************/
/*Assign the calculated derivatives to the corresponding MT structure (corresponds to determing the frechets)*/
/*Parameter:    *data  := Pointer on data structure*/
/*				*mt := MT structure (of the inversion grid) including the frechets*/
/*				*forward_mt := MT structure of the forward grid*/
/*				nr_inv_cell := Number of all inversion cells*/
/*              nr_of_forward_cells := Number of forward cells*/
/*			    kind_of_data_mt := Specify the kind of input data used: 1= TE-Mode, 2= TM-Mode, 3= 1D: Berdichewsky average/2D: Both Modes*/
/*				dimension_mt := 1 == 1D and 2 == 2D*/

int InvMTStruct(DATA_STRUCT *data, MT_STRUCT *mt, MT_STRUCT *forward_mt, BIDX_STRUCT *inv_cell, long nr_inv_cell, long nr_of_forward_cells, int kind_of_data_mt, int dimension_mt)
{
	/*Assign the calculated derivatives to the corresponding MT structure*/

	int index;
	long i,j,k,n,m,p, *inv_cell_index;
	long pos_inv,*count,nr_row;
	double **derivTE, **derivTM;

	printf("Start assigning the calculated derivatives to the corresponding \nMT structure\n");
	printf("----------------\n");

	inv_cell_index = (long *)memory(NULL,nr_of_forward_cells,sizeof(long),"InvMTStruct");

	for(n=0;n<nr_of_forward_cells;n++)
		inv_cell_index[n] = -1;

	/************************/
	/*Determine the inversion cell index for each forward cell*/
	for(k=0;k<nr_inv_cell;k++)
	{
		/*Loop over all forward cells in the inversion cells*/
		for(n=0;n<inv_cell[k].nele;n++)
		{
			if(nr_of_forward_cells <= inv_cell[k].ele[n])
			{
				printf("The index for the forward cells (%d) is larger than the specified number\nof forward cells (%d)\n", inv_cell[k].ele[n], nr_of_forward_cells);
				exit(0);
			}

			if(inv_cell[k].use == 1)
			{
				inv_cell_index[inv_cell[k].ele[n]] = k;
			}
			else
			{
				/*If the cells are in the air, they will be NOT considered*/
				inv_cell_index[inv_cell[k].ele[n]] = -1;
			}
		}
	}

	/************************/
	/*Allocate memory for the MT (inversion) structure*/
	for(i=0;i<data->ndata_mt;i++)
	{
		mt[i].n = (long *)memory(NULL,1,sizeof(long),"InvMTStruct");
		mt[i].ele = (long *)memory(NULL,1,sizeof(long),"InvMTStruct");
		mt[i].deriv = (double **)memory(NULL,1,sizeof(double *),"InvMTStruct");

		mt[i].ncell = 0;
		mt[i].nfreq = 0;
	}

	/************************/
	if(data->ndata_mt != 0)
		count = (long *)memory(NULL,data->ndata_mt,sizeof(long),"InvMTStruct");
	else
		count = (long *)memory(NULL,1,sizeof(long),"InvMTStruct");

	for(i=0;i<data->ndata_mt;i++)
		count[i] = 0;
	/************************/

	/*Loop over all data*/
	for(j=0;j<data->ndata_mt;j++)
	{
		mt[j].nfreq = forward_mt[j].nfreq;

		if(dimension_mt == 2)
		{
			
			/*Allocate memory for the derivatives*/
			if(forward_mt[j].nfreq != 0)
			{
				derivTE = (double **)memory(NULL,forward_mt[j].nfreq,sizeof(double *),"InvMTStruct");
				derivTM = (double **)memory(NULL,forward_mt[j].nfreq,sizeof(double *),"InvMTStruct");
			}
			else
			{
				derivTE = (double **)memory(NULL,1,sizeof(double *),"InvMTStruct");
				derivTM = (double **)memory(NULL,1,sizeof(double *),"InvMTStruct");
			}

			for(k=0;k<forward_mt[j].nfreq;k++)
			{
				derivTE[k] = (double *)memory(NULL,2*forward_mt[j].ncell,sizeof(double),"InvMTStruct");
				derivTM[k] = (double *)memory(NULL,2*forward_mt[j].ncell,sizeof(double),"InvMTStruct");

				for(m=0;m<2*forward_mt[j].ncell;m++)
				{
					derivTE[k][m] = 0.0;
					derivTM[k][m] = 0.0;
				}

				/*Read in the derivatives from tmp-files*/
				ReadinTmpMTfile(derivTE[k],derivTM[k],j,k, forward_mt[j].ncell, kind_of_data_mt);
			}
		}

		/*Loop over all cells affecting the MT measurement*/
		for(m=0;m<forward_mt[j].ncell;m++)
		{
			if(inv_cell_index[forward_mt[j].ele[m]] >= 0)
			{

				/*Check if the inversion cell is found already before*/
				index = 0;
				pos_inv = -99999;

				for(i=0;i<mt[j].ncell;i++)
				{
					if(mt[j].ele[i] == inv_cell_index[forward_mt[j].ele[m]])
					{
						pos_inv = i;
						index = 1;
						break;
					}
				}

				if(index == 0)
				{
					/*Fill the MT inversion structure*/
					mt[j].ele = (long *)memory((char *)mt[j].ele, count[j]+1, sizeof(long),"InvMTStruct");
					mt[j].deriv = (double **)memory((char *)mt[j].deriv, count[j]+1, sizeof(double *),"InvMTStruct");
								
					/*Two modes are used (in the 2-D modelling)*/
					if(kind_of_data_mt != 1 && kind_of_data_mt != 2 && dimension_mt == 2 && mt[j].nfreq != 0)
						mt[j].deriv[count[j]] = (double *)memory(NULL, 4*mt[j].nfreq, sizeof(double),"InvMTStruct");
					/*Only one mode is used*/
					else if(mt[j].nfreq != 0)
						mt[j].deriv[count[j]] = (double *)memory(NULL, 2*mt[j].nfreq, sizeof(double),"InvMTStruct");
					else
						mt[j].deriv[count[j]] = (double *)memory(NULL, 1, sizeof(double),"InvMTStruct");

					mt[j].ele[count[j]] = inv_cell_index[forward_mt[j].ele[m]];
					mt[j].ncell = count[j] + 1;

					/*Make the frechets for the 1-D modelling*/
					/*and for the 2-D modelling + 1-D inversion*/
					if(dimension_mt !=2)
					{
						for(p=0;p<(2*mt[j].nfreq);p++)
							mt[j].deriv[count[j]][p] = forward_mt[j].deriv[m][p];
					}
					/*... for the 2-D modelling*/
					else
					{
						/*Only TE*/
						if(kind_of_data_mt == 1)
						{
							for(p=0;p<mt[j].nfreq;p++)
							{
								/*Real part*/
								mt[j].deriv[count[j]][2*p] = derivTE[p][2*m];
								/*Imaginary part*/
								mt[j].deriv[count[j]][2*p+1] = derivTE[p][2*m+1];
							}
						}
						/*Only TM*/
						else if(kind_of_data_mt == 2)
						{
							for(p=0;p<mt[j].nfreq;p++)
							{
								/*Real part*/
								mt[j].deriv[count[j]][2*p] = derivTM[p][2*m];
								/*Imaginary part*/
								mt[j].deriv[count[j]][2*p+1] = derivTM[p][2*m+1];
							}
						}
						/*BOTH TE and TM-mode*/
						else
						{
							for(p=0;p<mt[j].nfreq;p++)
							{
								/*Real part*/
								mt[j].deriv[count[j]][4*p] = derivTE[p][2*m];
								mt[j].deriv[count[j]][4*p+2] = derivTM[p][2*m];
								/*Imaginary part*/
								mt[j].deriv[count[j]][4*p+1] = derivTE[p][2*m+1];
								mt[j].deriv[count[j]][4*p+3] = derivTM[p][2*m+1];
							}

						}
					}

					count[j]++;
				}
				else
				{
					if(pos_inv == -99999)
					{
						printf("Something goes wrong during generating the MT inversion srtucture\n");
						exit(0);
					}
	
					/*Make the frechets for the 1-D modelling*/
					/*and for the 2-D modelling + 1-D inversion*/
					if(dimension_mt != 2)
					{
						for(p=0;p<(2*mt[j].nfreq);p++)
							mt[j].deriv[pos_inv][p] = mt[j].deriv[pos_inv][p] + forward_mt[j].deriv[m][p];
					}
					/*... for the 2-D modelling*/
					else
					{
						/*Only TE*/
						if(kind_of_data_mt == 1)
						{
							for(p=0;p<mt[j].nfreq;p++)
							{
								/*Real parts*/
								mt[j].deriv[pos_inv][2*p] = derivTE[p][2*m] + mt[j].deriv[pos_inv][2*p];
								/*Imaginary part*/
								mt[j].deriv[pos_inv][2*p+1] = derivTE[p][2*m+1] + mt[j].deriv[pos_inv][2*p+1];
							}
						}
						/*Only TM*/
						else if(kind_of_data_mt == 2)
						{
							for(p=0;p<mt[j].nfreq;p++)
							{
								/*Real parts*/
								mt[j].deriv[pos_inv][2*p] = derivTM[p][2*m] + mt[j].deriv[pos_inv][2*p];
								/*Imaginary part*/
								mt[j].deriv[pos_inv][2*p+1] = derivTM[p][2*m+1] + mt[j].deriv[pos_inv][2*p+1];
							}
						}
						/*BOTH TE and TM-mode*/
						else
						{
							for(p=0;p<mt[j].nfreq;p++)
							{
								/*Real parts*/
								mt[j].deriv[pos_inv][4*p] = derivTE[p][2*m] + mt[j].deriv[pos_inv][4*p];
								mt[j].deriv[pos_inv][4*p+2] = derivTM[p][2*m] + mt[j].deriv[pos_inv][4*p+2];
								/*Imaginary part*/
								mt[j].deriv[pos_inv][4*p+1] = derivTE[p][2*m+1] + mt[j].deriv[pos_inv][4*p+1];
								mt[j].deriv[pos_inv][4*p+3] = derivTM[p][2*m+1] + mt[j].deriv[pos_inv][4*p+3];
							}

						}
						
					}
				}
			}
							
		}

		if(dimension_mt == 2)
		{
			for(k=0;k<forward_mt[j].nfreq;k++)
			{
				free(derivTE[k]);
				free(derivTM[k]);
			}

			free(derivTE);
			free(derivTM);
		}

	}

	/*Specify the index that correspond to the row in the inversion matrix*/
	nr_row = 0;

	for(j=0;j<data->ndata_mt;j++)
	{
		/*Two mode used*/
		if(kind_of_data_mt != 1 && kind_of_data_mt != 2	&& ( dimension_mt == 2 || dimension_mt == 3))
		{
			mt[j].n = (long *)memory((char *)mt[j].n,(4*mt[j].nfreq),sizeof(long),"InvMTStruct");

			/*Index starts at 0*/
			for(k=0;k<(4*mt[j].nfreq);k++)
				mt[j].n[k] = k + nr_row;

			nr_row = nr_row + (4*mt[j].nfreq);
		}
		/*Only one mode used*/
		else
		{
			mt[j].n = (long *)memory((char *)mt[j].n,(2*mt[j].nfreq),sizeof(long),"InvMTStruct");

			/*Index starts at 0*/
			for(k=0;k<(2*mt[j].nfreq);k++)
				mt[j].n[k] = k + nr_row;

			nr_row = nr_row + (2*mt[j].nfreq);
		}
	}

	free(count);
	free(inv_cell_index);

	printf("The MT derivatives are determined\n");
	printf("----------------\n\n");

	return(0);
}


/*********************************************************************************/
/*Read in the MT-derivatives of the 2D modeling from the temporary files*/
/*Parameter:    deriv = MT derivatives of the TE and TM mode*/
/*              i_stat = station index*/
/*              j_freq  = frequency index*/
/*              nr_of_cells := Number of forward cells in the model*/
/*			    kind_of_data_mt := Specify the kind of input data used: 1= TE-Mode, 2= TM-Mode, 3= 1D: Berdichewsky average/2D: Both Modes*/

int ReadinTmpMTfile(double *derivTE, double *derivTM ,long i_stat, long j_freq, long nr_of_cells, int kind_of_data_mt)
{
	long i;

	double tmp_deriv;
	char fname[40];
	FILE *inf;

	/*Read in the derivatives of the TE-mode*/
	if(kind_of_data_mt != 2)
	{
		sprintf(fname,"tmp_mt_TE_deriv%d_freq%d.dat",i_stat,j_freq);

		inf = fopen(fname,"rb");

		if (inf == NULL)
		{
			fprintf(stderr,"Unable to open %s\n",fname);
			exit(0);
		}

		for(i=0;i<2*nr_of_cells;i++)
		{
			fread(&(tmp_deriv),sizeof(double),1,inf);		/*Read in the derivatives*/
			derivTE[i] = tmp_deriv;
		}

		fclose(inf);

		/*Remove the temporary file*/
		remove(fname);
	}

	/*Read in the derivatives of the TM-mode*/
	if(kind_of_data_mt != 1) 
	{
		sprintf(fname,"tmp_mt_TM_deriv%d_freq%d.dat",i_stat,j_freq);
		
		inf = fopen(fname,"rb");

		if (inf == NULL)
		{
			fprintf(stderr,"Unable to open %s\n",fname);
			exit(0);
		}

		for(i=0;i<2*nr_of_cells;i++)
		{
			fread(&(tmp_deriv),sizeof(double),1,inf);		/*Read in the derivatives*/
			derivTM[i] = tmp_deriv;
		}

		fclose(inf);

		/*Remove the temporary file*/
		remove(fname);
	}

	return(0);
}
