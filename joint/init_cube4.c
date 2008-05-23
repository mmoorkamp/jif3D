#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "inv3d.h"
#include "fkt.h"


/*------------------------------------------------------------*/
/*Make a topographic model and construct the initial velocity, density and/ elec.resistivity models */
/*Parameter:    *grid   := Pointer on the grid structure*/
/*              *geo    := Pointer on the structure of the receiver/shot locations */
/*               flag   := Global settings/flags*/
/*		   fname_dens   := File name for the density-velocity parameters*/
/*		   fname_res	:= File name for the resistivity-velocity parameters*/
/*             *topo    := Topography structure*/

/* Output: Indices of the cells that describe the topography*/

int MakeInitModels(GRID_STRUCT *grid, GEOMETRY *geo, FLAG_STRUCT flag, char *fname_dens, char *fname_res, TOPO_STRUCT *topo)
{ 
	 int i;
	 long nxyz,nxyz2;
	 long *sample_topo;

	 if(flag.index_tseis == 0)
			grid->nborder = 0; /*If the velocity model is NOT calculated, no borders will be used*/

	 nxyz  = (grid->nx)*(grid->ny)*(grid->nz);
	 nxyz2 = ((2*grid->nborder) + (grid->nx))*((2*grid->nborder) + (grid->ny))*((2*grid->nborder) + (grid->nz));

    /*************************************************************************************/ 
	/*					Calculation of the topography model:							 */
	/*************************************************************************************/
	
	 sample_topo = TopoCalc(geo,grid,topo);

	/***********************************************************************************/
	/***********************************************************************************/
	/* Start loop to determine the different models                                    */
	/*     first run:   velocity                                                       */
	/*     second run:  density                                                        */
	/*      third run:  resistivity													   */
	/***********************************************************************************/
	
	 /*Calculate the models and the "border-index"*/
		grid->border_index = (int *)memory((char *)grid->border_index,nxyz2,sizeof(int),"MakeIndividualModel");

	 for(i=0;i<3;i++)
	 {
		 /*Select the parameters for the corresponding starting model*/
		 if(i==0)
		 {
			 if(flag.index_tseis != 0)
			 {
				printf("SEISMIC VELOCITY:\n------------------\n");
				/*seismic model will be calculated by means of the input parameters*/
				MakeIndividualModels(grid,i+1, sample_topo, 0, 0, 0);
			 }
		 }

		 if(i==1)
		 {
			if(flag.index_grav != 0)
			{
				printf("DENSITY:\n------------------\n");
				/*density model will be calculated by means of the input parameters*/
				MakeIndividualModels(grid,i+1, sample_topo, flag.link_seis_grav_mod, fname_dens, 0);
			}
		 }

		 if(i==2)
		 {
			if(flag.index_mt != 0)
			{
				printf("RESISTIVITY MODEL:\n------------------\n");
				/*elec.resistivity model will be calculated*/
				MakeIndividualModels(grid,i+1, sample_topo, flag.link_seis_mt_mod, fname_res, fname_dens);
			}
		 }

	 }

 free(sample_topo);

  return(0);
 }

/*------------------------------------------------------------*/
/* Construct the initial velocity, density or elec.resistivity models */
/*Parameter:    *grid   := Pointer on the grid structure*/
/*        kind_of_model	:= Specify the kind of model (1= velocity, 2 = density, 3= resistivity)
/*       *sample_topo   := pointer on samples of the grid cell centers (in z-direction), where the topography starts; (The samples in z-direction can lie */ 
/*		 vel_dens_link	:= specify if the velocity and density/resistivity are linked to each other (0= NO/ 1= YES(ANALYTICAL)/ 2= YES(By a table))*/
/*				*fname	:= Filename including the relationsship of velocity and density resp. resistivity*/
/*				*fname2	:= Filename including the relationsship of velocity and density for the case that the joint inversion of gravity and MT data will be performed*/

int MakeIndividualModels(GRID_STRUCT *grid, int kind_of_model, long *sample_topo, int vel_dens_link, char *fname, char *fname2)
{
	long a,b,c;
	long nx,ny,nz,nyz,nxyz,nx2,ny2,nz2,nyz2,nxyz2,nx3,ny3,nz3,nyz3,nxyz3;
	long nborder;
	int *borderIndex, *borderIndex2;
	double *mod, *mod2;
	float h;
	REL_VDR_STRUCT rel_vdr, rel_vdr2;

	#define b_index(x,y,z) borderIndex[nyz*(x) + nz*(y) + (z)]

	h = grid->h;
	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nborder = grid->nborder;
    nyz = nz*ny;
	nxyz = nx*ny*nz;


	nx2 = 2*nborder + nx;
	ny2 = 2*nborder + ny;
	nz2 = 2*nborder + nz;
	nyz2 = ny2*nz2;
	nxyz2 = nx2*ny2*nz2;

	nx3= nx2+1;
	ny3= ny2+1;
	nz3= nz2+1;
	nyz3= ny3*nz3;
	nxyz3 = nx3*ny3*nz3;

	/* Allocate memory for the border_index: */
	borderIndex = (int *)memory(NULL,nxyz,sizeof(int),"MakeIndividualModels");
	/*Allocate memory for the models*/
	mod = (double *)memory(NULL,nxyz,sizeof(double),"MakeIndividualModel");

	/*Set "boundary-index" to 1 for all grid cells:*/
	for(a=0; a<nx; a++)
		for(b=0; b<ny;b++)
			for(c=0; c<nz;c++)
			{
				b_index(a,b,c) = 1; 
			}


	/*Allocate the memory for the models*/
	if(kind_of_model == 1) /*s.velocity (ATTENTION!! The seismic model is larger than the corresponding density and res.models)*/
		grid->slow = (double *)memory((char *)grid->slow,nxyz3,sizeof(double),"MakeIndividualModels");

	if(kind_of_model == 2) /*density*/
		grid->dens = (double *)memory((char *)grid->dens,nxyz,sizeof(double),"MakeIndividualModels");

	if(kind_of_model == 3) /*resistivity*/
		grid->res = (double *)memory((char *)grid->res,nxyz,sizeof(double),"MakeIndividualModels");

	/*************************************************************************************/
	/* Determine the model by means of the link with the velocity model					 */
	/*************************************************************************************/
	if((vel_dens_link != 0 && kind_of_model  == 2) ||  /*gravity*/
	   (vel_dens_link != 0 && kind_of_model  == 3 ))    /*resistivity*/
	{
		if(kind_of_model == 2)
			printf("The density starting model is determined by linking it to the velocity model\n\n");
		else if(kind_of_model == 3 && vel_dens_link == 1)
			printf("The resistivity starting model is determined by linking it to the velocity model\n\n");
		else if(kind_of_model == 3 && vel_dens_link == 2)
			printf("The resistivity starting model is determined by linking it to the gravity model\n\n");

		/*The starting density OR resistivity model is linked to the seismic model*/
		if(vel_dens_link == 1)
		{
			/*Read in the file including the relationship of velocity and density resp. resistivity*/
			ReadRelVelDensRes(fname, &rel_vdr);

			/*Calculate the density/resistivity model from the link to the velocity model*/
			CalcVelDepMod(grid, mod, rel_vdr);
		}

		/*The starting resistivity model is linked to the gravity model*/
		else if(vel_dens_link == 2 && kind_of_model == 3)
		{
			/*Read in the file including the relationship of velocity and resistivity*/
			ReadRelVelDensRes(fname, &rel_vdr);
			
			/*Read in the file including the relationship of velocity and density*/
			ReadRelVelDensRes(fname2, &rel_vdr2);

			/*Calculate the resistivity model from the link to the gravity model*/
			/*IMPORTANT !!! In this routine also the slowness model (grid->slow) will be filled */
			CalcGravDepMod(grid, mod, rel_vdr, rel_vdr2);

			free(rel_vdr2.tab_v);
			free(rel_vdr2.tab_dr);
		}

		else
		{
			printf("The starting models could not linked to each other !!!\n");
			exit(0);
		}

		/*Assign the determined topography to the model*/
		ModCubeCalc(mod,borderIndex, grid, sample_topo, kind_of_model);

		free(rel_vdr.tab_v);
		free(rel_vdr.tab_dr);
	}

	/*************************************************************************************/
	/* Determine the model by means of the topography and the mean of the gradient       */
	/*************************************************************************************/
	else if((grid->grad_vel.index_topo == 1 && kind_of_model == 1) ||	/*seismic*/
	   (grid->grad_dens.index_topo == 1 && kind_of_model == 2) ||	/*gravity*/
	   (grid->grad_res.index_topo == 1 && kind_of_model == 3))		/*MT*/
	{
		/*Calculate the gradient model (velocity, density or resistivity) (by means of the topography)*/
		CalcSurfDepMod(sample_topo,grid,mod,kind_of_model);

		/*Assign the determined topography to the model*/
		 ModCubeCalc(mod,borderIndex, grid, sample_topo, kind_of_model);

	}
	/*********************************************************************************/
	/* Determine the model by means of a defined origin and the gradient components  */
	/*********************************************************************************/
	else
	{
		/*Calculate the gradient model (velocity, density or resistivity) (by means of the defined origin)*/
		CalcOrgDepMod(grid, mod, kind_of_model);

		/*Assign the determined topography to the model*/
		 ModCubeCalc(mod,borderIndex, grid, sample_topo, kind_of_model);

	}

	/********************************************************************************/

	/* Add boundaries to data cube*/
  	mod2 = (double *)memory(NULL,(nx2*ny2*nz2),sizeof(double),"MakeIndividualModel");
	borderIndex2 = (int *)memory(NULL,(nx2*ny2*nz2),sizeof(int),"MakeIndividualModel");

	AddBoundaries(grid,mod,borderIndex,mod2,borderIndex2);


	for(a=0; a<nx2; a++)
		for(b=0; b<ny2;b++)
			for(c=0; c<nz2;c++)
			{
				grid->border_index[nyz2*(a) + nz2*(b) + (c)] = borderIndex2[nyz2*(a) + nz2*(b) + (c)];
			}

	
	if(kind_of_model == 1) /*make the final seismic starting model*/
	{
		for(a=0; a<nx2; a++)
			for(b=0; b<ny2;b++)
				for(c=0; c<nz2;c++)
				{
					grid->slow[nyz3*(a) + nz3*(b) + (c)] = (h/mod2[nyz2*(a) + nz2*(b) + (c)]); /*ATTENTION: The size of the seismic model is specified by the number of nodes*/
				}
	}

	if(kind_of_model == 2) /*make the final density starting model*/
	{
		for(a=0; a<nx; a++)
			for(b=0; b<ny;b++)
				for(c=0; c<nz;c++)
				{
					grid->dens[nyz*(a) + nz*(b) + (c)] = (mod[nyz*(a) + nz*(b) + (c)]);
				}
	}

	if(kind_of_model == 3) /*make the final density starting model*/
	{
		for(a=0; a<nx; a++)
			for(b=0; b<ny;b++)
				for(c=0; c<nz;c++)
				{
					grid->res[nyz*(a) + nz*(b) + (c)] = (mod[nyz*(a) + nz*(b) + (c)]);
				}
	}


	printf("----------------\n\n");

	free(mod); 
	free(mod2);
	free(borderIndex);
	free(borderIndex2);

	#undef b_index

	return(1);
}



/*-------------------------------------------------------------*/
/* Calculate the topography model out of the GEOMETRY structure (receiver/shot and other positional informations) */
/* by using the Delaunay triangulation */
/* Parameter:    *geo := pointer on the geometry structure including the positional information*/
/*               *grid:= pointer on the grid structure */
/*               *topo:= Topography structure*/
/* Output: Topography model for the complete area (the samples in z-direction, at which the surface topography starts are written out) */


long *TopoCalc(GEOMETRY *geo, GRID_STRUCT *grid, TOPO_STRUCT *topo)
 {
	  int i,j, num_topo_data			/*Number of the topographic points used for the triangulation*/;
	  int ntriangles;					/*Number of triangles determined by the triangulation*/
	  int *corner_triangles;			/*corner of the triangles; the values i of the vector correspond to the index of x(i), y(i) and z(i) coordinates; 
											all triangles are listed consecutively leading to (3*nr. of triangles) values */
	  long  *sample_topo;				/*samples in z-direction, at which the surface topography starts*/
	  int topo_index;
      long  nx,ny,nxy;

	  float orgx,orgy,orgz,h;

	  double *x,*y,*z;					/*coordinates used to determine the topography model*/
	  double *xgrid_edges, *ygrid_edges, *zgrid_edges; /*x and y coordinates of the edges of the grid cells; "zgrid_edges" are the calculated height 
														of the topography at the grid edge positions using triangulation*/

	  double dist;

	  FILE *out;
	  FILE *out1;

	  #define Xgrid(x,y) xgrid_edges[(ny+1)*(x)+(y)]
	  #define Ygrid(x,y) ygrid_edges[(ny+1)*(x)+(y)]
	  #define Zgrid(x,y) zgrid_edges[(ny+1)*(x)+(y)]
	  #define sample(x,y) sample_topo[ny*(x)+(y)]
	  #define EPS 0.01
	  
	  nx = grid->nx;
	  ny = grid->ny;
	  orgx = grid->org[0];
	  orgy = grid->org[1];
	  orgz = grid->org[2];
	  h = grid->h;
	  topo_index = grid->topo_index;
	  nxy = nx*ny; 


	  /***************************************************************************/
	  /*		Ignore the topography											 */
	  /***************************************************************************/
	  if(topo_index == 0)
	  {
		sample_topo = (long *)memory(NULL,nxy,sizeof(long),"TopoCalc");

		for(i=0;i<nx;i++)
			for(j=0;j<ny;j++)
			{
				sample(i,j)=0;
			}
         
		x = (double *)memory(NULL,1,sizeof(double),"TopoCalc");
		y = (double *)memory(NULL,1,sizeof(double),"TopoCalc");
		z = (double *)memory(NULL,1,sizeof(double),"TopoCalc");
		corner_triangles = (int *)memory(NULL,1,sizeof(int),"TopoCalc");

		 num_topo_data = 0;
		 ntriangles = 0;

		 printf("The grid cells of the 'topography' running along the upper border of the model\n");
		 printf("----------------\n\n");
	  }
	  /***************************************************************************/
	  /*		Determine the topography by means of triangulation				 */
	  /***************************************************************************/
	  else
	  {	 
		/* Determine the number of coordinates belonging to the topography*/
		num_topo_data = 0;
			for(i=0;i<(geo->nstat);i++)
			{
				if(geo->coor_info[i] == 1)
					num_topo_data++;
			}


		if(num_topo_data > 0)
		{
			x = (double *)memory(NULL,num_topo_data,sizeof(double),"TopoCalc");
			y = (double *)memory(NULL,num_topo_data,sizeof(double),"TopoCalc");
			z = (double *)memory(NULL,num_topo_data,sizeof(double),"TopoCalc");
		}
		else
		{
			x = (double *)memory(NULL,1,sizeof(double),"TopoCalc");
			y = (double *)memory(NULL,1,sizeof(double),"TopoCalc");
			z = (double *)memory(NULL,1,sizeof(double),"TopoCalc");
		}

	
		 num_topo_data = 0;
			for(i=0;i<(geo->nstat);i++)
			{
				if(geo->coor_info[i] == 1)
				{
					x[num_topo_data]=(double)(geo->x[i]);
					y[num_topo_data]=(double)(geo->y[i]);
					z[num_topo_data]=(double)(geo->z[i]);
					num_topo_data++;
				}
			}

		/*Remove similar coordinates (necessary, because the triangulation terminates if the coordinates are located too close to each other)*/
			for(i=0;i<num_topo_data;i++)
				for(j=i+1;j<num_topo_data;j++)
				{
					dist = sqrt((x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]));
					if(dist <= EPS)
					{
			
						num_topo_data = num_topo_data - 1;
						memmove(&(x[j]),&(x[j+1]),(num_topo_data -j)*sizeof(double));
						memmove(&(y[j]),&(y[j+1]),(num_topo_data -j)*sizeof(double));
						memmove(&(z[j]),&(z[j+1]),(num_topo_data -j)*sizeof(double));
					}
		  
				}

		/*Performing Delaunay triangulation:*/


		printf("Delaunay triangulation will be started:\n");
		printf("Number of used coordinates for triangulation: %d\n\n",num_topo_data);

		ntriangles = delaunay(x,y,num_topo_data,&(corner_triangles));
	
		printf("All triangles are calculated:\n");
		printf("Number of triangles: %d\n",ntriangles);
		printf("----------------\n\n");

		if(ntriangles == 0)
		{
			printf("!!! TOO less coordinate information for triangulation !!!\n");
			exit(0);
		}
				/***************************************/
				/*Write out the triangle positions in "triangle.txt"*/
				out = fopen("triangle.txt","wt");
				fprintf(out,"#Position of all triangles determined by the triangulation:\n");
				fprintf(out, "#     x1          x2          x3          y1         y2          y3          z1           z2          z3\n");
				for(i=0;i<3*ntriangles;i=i+3)
				{
					fprintf(out,"%f %f %f %f %f %f %f %f %f\n",x[corner_triangles[i]],x[corner_triangles[i+1]],x[corner_triangles[i+2]],y[corner_triangles[i]],y[corner_triangles[i+1]],y[corner_triangles[i+2]],z[corner_triangles[i]],z[corner_triangles[i+1]],z[corner_triangles[i+2]]);
				}
				fclose(out);
				/***************************************/
				/*Write the topography points out*/
				out1 = fopen("topo.txt","wt");
				for(i=0;i<num_topo_data;i++)
				{
					fprintf(out,"%d %f %f %f\n",i+1,x[i],y[i],z[i]);
				}
				fclose(out1);

		/*Determine the x-y coordinates of the edges of the grid cells*/
		xgrid_edges = (double *)memory(NULL,(nx+1)*(ny+1),sizeof(double),"TopoCalc");
		ygrid_edges = (double *)memory(NULL,(nx+1)*(ny+1),sizeof(double),"TopoCalc");
		zgrid_edges = (double *)memory(NULL,(nx+1)*(ny+1),sizeof(double),"TopoCalc");

		for(i=0;i<(nx+1);i++)
			for(j=0;j<(ny+1);j++)
			{
				Xgrid(i,j) = orgx - (h/2.0) + h*i;
				Ygrid(i,j) = orgy - (h/2.0) + h*j;
			}


	   /*Interpolation by means of triangulation/ using barycentric coordinates (see in: Contouring (1992) from D.F. Watson; Pergamon Press; pages 76-78)*/
	      TriangleInterpol(xgrid_edges,ygrid_edges,zgrid_edges,(nx+1)*(ny+1),x,y,z,num_topo_data,corner_triangles,ntriangles);
		  
			   /***************************************/
			   /*Write out the gridded coordinates in "topography.txt"*/
				printf("\n\nThe topography heights for all coordinates are calculated\n");
				printf("----------------\n\n");

				out = fopen("topography.txt","wt");
				fprintf(out,"# Topography modell:\n");
				fprintf(out,"#     x          y         z\n");
				for(i=0;i<(nx+1)*(ny+1);i++)
				{
						fprintf(out,"%f %f %f\n",xgrid_edges[i],ygrid_edges[i],zgrid_edges[i]);
				}
				fclose(out);
				/***************************************/


	    /*Fit the grid to the topography*/
			sample_topo = FitGridtoTopo(zgrid_edges,grid,*geo);

			printf("The grid cells along the topography are determined\n");
			printf("----------------\n\n");

			/***************************************/
			/*out = fopen("topogridd.txt","wt");
			for(i=0;i<nxy;i++)
				fprintf(out,"%d\n", sample_topo[i]);
			fclose(out);*/
			/***************************************/


		free(xgrid_edges);
		free(ygrid_edges);
		free(zgrid_edges);
	  }


	  /***************************************/
	  /***************************************/
		/*Make the topography structure*/
	  /***************************************/
	  /***************************************/
		topo->nr_of_topo_points = num_topo_data;
		topo->nr_of_triangles = ntriangles;

		if(num_topo_data != 0)
		{
			topo->x = (double *)memory(NULL, num_topo_data, sizeof(double),"TopoCalc");
			topo->y = (double *)memory(NULL, num_topo_data, sizeof(double),"TopoCalc");
			topo->z = (double *)memory(NULL, num_topo_data, sizeof(double),"TopoCalc");
		}
		else
		{
			topo->x = (double *)memory(NULL, 1, sizeof(double),"TopoCalc");
			topo->y = (double *)memory(NULL, 1, sizeof(double),"TopoCalc");
			topo->z = (double *)memory(NULL, 1, sizeof(double),"TopoCalc");
		}

		for(i=0;i<num_topo_data;i++)
		{
			topo->x[i] = x[i];
			topo->y[i] = y[i];
			topo->z[i] = z[i];
		}

		if(ntriangles != 0)
		    topo->index = (int *)memory(NULL, 3*ntriangles, sizeof(int),"TopoCalc");
		else
			topo->index = (int *)memory(NULL, 1, sizeof(int),"TopoCalc");

		for(i=0;i<(3*ntriangles);i++)
			topo->index[i] = corner_triangles[i];
		/***************************************/
		/***************************************/

		free(x);
		free(y);
		free(z);
		free(corner_triangles);


	  #undef Xgrid
	  #undef Ygrid
 	  #undef Zgrid
	  #undef sample
	  #undef EPS

	  return(sample_topo);
   }



/*-------------------------------------------------------------*/
/* Calculate the samples of the (nx*ny) grid centers (in z-direction), where the topography starts*/
/* Parameter:    *Z := pointer on the z-coordinates of the topography at the (nx+1)*(ny+1) grid cell edges*/
/*               *grid:= pointer on the grid structure*/
/*				  geo:= Geometry structure */
/* Output:       samples of the (nx*ny) grid cell centers (in z-direction) where the topography starts; (The samples in z-direction can lie */ 
/*               outside of the subsequent grid)*/



long *FitGridtoTopo(double *Z, GRID_STRUCT *grid, GEOMETRY geo)
{
    long i,j,smin,*sample_topo; /*sample of the topography at the grid cell center*/
	int ix,iy,iz,index_x,index_y,index_z;
	long *sample_topo_edge; /*sample of the topography at the grid cell edges*/
	long nx,ny,nz;
	float h, orgx, orgy, orgz; 
	double iX,iY,iZ,deltaX,deltaY,deltaZ;

	#define Z(x,y) Z[(ny+1)*(x)+(y)]
	#define sampleE(x,y) sample_topo_edge[(ny+1)*(x)+(y)]
	#define sample(x,y) sample_topo[ny*(x)+(y)]

	nx = (grid->nx);
	ny = (grid->ny);
	nz = (grid->nz);
	h  = (grid->h);
	orgx = (grid->org[0]);
	orgy = (grid->org[1]);
	orgz = (grid->org[2]);
	
	sample_topo_edge = (long *)memory(NULL,(nx+1)*(ny+1),sizeof(long),"FitGridtoTopo");
	sample_topo = (long *)memory(NULL,nx*ny,sizeof(long),"FitGridtoTopo");


	/* Determine the sample of the topography for the grid cell edges*/
	for(i=0;i<(nx+1);i++)
		for(j=0;j<(ny+1);j++)
		{
			/*For cells without topography informations*/
			  if(Z(i,j)== -99999.9)
			  {
				  sampleE(i,j) = nz;
			  }
			 /*For cells with topography informations*/
			  else
			  {
				  sampleE(i,j) =  (int) floor((Z(i,j)-(orgz-(h/2.0)))/h);
			  }

		}

	/* Determine the sample of the topography for the grid cell centers*/
	for(i=0;i<nx;i++)
		for(j=0;j<ny;j++)
		{
			/*Looking for the edge with the highest topography value at the edges*/
			smin = sampleE(i,j);
			if(sampleE(i+1,j) < smin)
				smin = sampleE(i+1,j);
			if(sampleE(i,j+1) < smin)
				smin = sampleE(i,j+1);
			if(sampleE(i+1,j+1) < smin)
				smin = sampleE(i+1,j+1);

			sample(i,j) = smin;
		}

	/*Introduce the exception that the height of the GEOMETRY points is higher than the*/
	/*heights at the surrounding grid cell edges*/

		for(i=0;i<geo.nstat;i++)
		{
			if(geo.coor_info[i] == 1)
			{
				ix = (int)((geo.x[i] - orgx)/h);
				iy = (int)((geo.y[i] - orgy)/h);
				iz = (int)((geo.z[i] - orgz)/h);

				iX = (double)((geo.x[i]- orgx)/h); /*normalized distances between Grid-origins and coordinates*/
				iY = (double)((geo.y[i]- orgy)/h);
				iZ = (double)((geo.z[i]- orgz)/h);

		deltaX = iX-ix; /*distance from the shot coordintes to next lower grid cell center*/
		deltaY = iY-iy;
		deltaZ = iZ-iz;

		if(deltaX <= 0.5) 
			index_x = ix;	/*index of next grid cell center in x-direction*/
		else
			index_x = ix+1;

		if(deltaY <= 0.5) 
			index_y = iy;
		else
			index_y = iy+1;

		if(deltaZ <= 0.5) 
			index_z = iz;
		else
			index_z = iz+1;
				
		if(index_x >=0 && index_x <nx && index_y >= 0 && index_y <ny)
					if(sample(index_x,index_y) > index_z) /*The case where the heights of the Geometry point is larger than those from the surrounding grid cell edges*/
						sample(index_x,index_y) = index_z;

			}
		}

	free(sample_topo_edge);

	#undef Z
	#undef sampleE
	#undef sample

	return(sample_topo);
}


/*-------------------------------------------------------------*/
/* Introduce the topography into the model*/
/* Parameter:    *mod := pointer on the model*/
/*				 *border_index := marking cells that belong to the boundaries or to the air*/
/*               *grid:= pointer on the grid structure */
/*               *sample_topo:= pointer on samples of the grid cell centers (in z-direction), where the topography starts; (The samples in z-direction can lie */ 
/*								outside of the subsequent grid)*/
/*				 kind_of_model:= Identify the kind of the model(velocity =1; density =2; resistivity =3)*/
 
int ModCubeCalc(double *mod,int *border_index, GRID_STRUCT *grid ,long *sample_topo, int kind_of_model)
{
	int a,b,c;
	long nx,ny,nyz,nz;
	int Upper_boundary,Lower_boundary;
	double mborder, m_air_water;

	#define m(x,y,z) mod[nyz*(x) + nz*(y) + (z)]
	#define b_index(x,y,z) border_index[nyz*(x) + nz*(y) + (z)]
	#define sample(x,y) sample_topo[ny*(x)+(y)]
	#define NOT_REACHED 0
	#define REACHED 1

	# define EPS 0.000000000000000000001

	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nyz=ny*nz;

	
	/*Choose the model*/
	if(kind_of_model == 1)		/*s.velocity*/
	{
		mborder = grid->vborder;
		m_air_water = grid->v_air_water;
	}
	else if(kind_of_model == 2) /*density*/
	{
		mborder = 0.0;		/*<- If the topography is not specified in a region, the cells in this region will be not considered (g == 0)*/
		m_air_water = grid->dens_air_water;
	}
	else if(kind_of_model == 3)	/*resistivity*/
	{
		mborder = 1/EPS;
		m_air_water = grid->res_air_water; /*<- If the topography is not specified in a region, the cells in this region will be not considered (res == inf)*/
	}
	else
	{
		printf("Choosen model is unknown\n");
		exit(0);
	}



		Upper_boundary = NOT_REACHED;
		Lower_boundary = NOT_REACHED;

		for(a=0; a<nx; a++)
			for(b=0; b<ny; b++)
		{
			if(sample(a,b) < 0)
				Upper_boundary = REACHED;

			if(sample(a,b) >= nz)
			{
				Lower_boundary = REACHED;

				for(c=0;c<nz;c++)
				{
					/*Assign boundary velocity in regions where no topography exists*/
					m(a,b,c) = mborder;
					b_index(a,b,c) = 0;
				}
			}
			else
			{
				for(c=0;c<sample(a,b);c++)
				{
					/*Assign air/water velocity above the surface topography*/
					m(a,b,c) = m_air_water;
					b_index(a,b,c) = 0;
				}
			}
				
		}

		if(kind_of_model == 1)
			printf("The topography is introduced into the velocity model!\n");
		if(kind_of_model == 2)
			printf("The topography is introduced into the density model!\n");
		if(kind_of_model == 3)
			printf("The topography is introduced into the resistivity model!\n");

		if(Upper_boundary == REACHED)
			 printf("!!ATTENTION!! The topography touch the upper boundary of the grid!\n----------------\n\n");
		if(Lower_boundary == REACHED && kind_of_model == 1)
			 printf("!!ATTENTION!! The topography touch the lower boundary of the grid!\n----------------\n\n");
		if(Lower_boundary == REACHED && kind_of_model == 2)
		{
			 printf("!!ATTENTION!! The topography touch the lower boundary of the grid!\n");
			 printf("The density is set to g = 0.0 g/cm^3, where the topography is not specified !!\n----------------\n\n");
		}
		if(Lower_boundary == REACHED && kind_of_model == 3)
		{
			 printf("!!ATTENTION!! The topography touch the lower boundary of the grid!\n");
			 printf("The electrical res. is set to r = inf ohmm, where the topography is not specified !!\n----------------\n\n");
		}

		#undef m
		#undef b_index
		#undef sample
		#undef NOT_REACHED
		#undef REACHED
		#undef EPS

	return(1);
}

/*---------------------------------------------------------------------*/
/*Calculate the resistivity model by linking it to the starting density model (via the velocity relations)*/
/*Parameters    *grid   := Pointer on the grid structure*/
/*               *mod   := resistivity model */
/*			   rel_vdr	:= Structure that specifying the relation of velocity and gravity*/
/*			   rel_vdr2	:= Structure that specifying the relation of velocity and resistivity*/
/*Remark: In this routine also the slowness model will be calculated, since it is required lateron in the inversion*/

int CalcGravDepMod(GRID_STRUCT *grid, double *mod, REL_VDR_STRUCT rel_vdr, REL_VDR_STRUCT rel_vdr2)
{
	long a,b,c,nx,ny,nz,nyz,nx1,ny1,nz1,nyz1;
	double vel_value;

	#define m(x,y,z) mod[nyz*(x) + nz*(y) + (z)]
	#define dens(x,y,z) grid->dens[nyz*(x) + nz*(y) + (z)]

	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nyz = ny*nz;

	nx1 = nx + 1;
	ny1 = ny + 1;
	nz1 = nz + 1;
	nyz1 = ny1*nz1;

	/*Sort the densities in ascending order for tabular links*/
	if(rel_vdr2.index == 9)
		sort2(rel_vdr2.tab_nr_pairs,rel_vdr2.tab_dr,rel_vdr2.tab_v);

	grid->slow = (double *)memory((char *) grid->slow,nx1*ny1*nz1,sizeof(double),"CalcGravDepMod");

	/*Loop over all cells in the model*/
	for(a=0; a<nx; a++)
		for(b=0; b<ny;b++)
			for(c=0;c<nz;c++)
			{
				/*Determine the corresponding velocity value from the density model*/
				vel_value = DensRestoVel(dens(a,b,c),rel_vdr2);
				/*Calculate the corresponding resistivity value from the velocity value*/
				m(a,b,c) = VeltoDensRes(vel_value,rel_vdr);

				if(m(a,b,c) < 0)
				{
					printf("After determing the resistivities from the velocities,\n");
					printf("some restivity values occur that are smaller than 0 !!\n");
					printf("Please adjust the formulars in the corresponding parameter files !!\n");
					exit(0);
				}

				/*Make a slowness model (since the data will be inverted for slownesses for the combination of gravity and MT data)*/
				grid->slow[nyz1*a + nz1*b +c] = (grid->h)/vel_value;
			}

	printf("model is calculated\n");
	printf("----------------\n");

	#undef m
	#undef dens 

	return(0);
}

/*---------------------------------------------------------------------*/
/*Calculate the density/resistivity model by linking it to the starting velocity model*/
/*Parameters    *grid   := Pointer on the grid structure*/
/*               *mod   := model (dependent on the parameter the density or resistivity model will be determined)*/
/*			   rel_vd	:= Structure that specifying the relation of velocity and density/resistivity*/

int CalcVelDepMod(GRID_STRUCT *grid, double *mod, REL_VDR_STRUCT rel_vdr)
{
	long a,b,c,a1,b1,c1;
	long nx,ny,nz,nyz;
	long ny3,nz3,nyz3;
	double vel_value;

	#define m(x,y,z) mod[nyz*(x) + nz*(y) + (z)]
    #define slow(x,y,z) grid->slow[nyz3*(x) + nz3*(y) + (z)]

	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nyz = ny*nz;

	ny3 = ny + (2*grid->nborder) + 1;
	nz3 = nz + (2*grid->nborder) + 1;
	nyz3 = ny3*nz3;

	/*Loop over all cells in the model*/
	for(a=0; a<nx; a++)
		for(b=0; b<ny;b++)
			for(c=0;c<nz;c++)
			{
				a1 = a + grid->nborder;
				b1 = b + grid->nborder;
				c1 = c + grid->nborder;

				/*Determine the corresponding velocity value of he model*/
				vel_value = grid->h/slow(a1,b1,c1); /*ATTENTION !!! The slowness model has a different number of cells than the density and the resistivity model*/

				/*Calculate the corresponding density/resistivity value from the velocity value*/
				m(a,b,c) = VeltoDensRes(vel_value,rel_vdr);

				/*Terminate the program if density/resistivity values are < 0*/
				if( m(a,b,c) < 0)
				{
					printf("After determing the resistivities/densities from the velocities,\n");
					printf("some restivity OR density values occur that are smaller than 0 !!\n");
					printf("Please adjust the formulars in the corresponding parameter files !!\n");
					exit(0);
				}
			}
	
	printf("model is calculated\n");
	printf("----------------\n");

	#undef m
	#undef slow

	return(0);
}

/*---------------------------------------------------------------------*/
/*Calculate the surface-topography dependent gradient model (velocity, density or resistivity) v by determine the shortest distance from the surface (d,e,sample(d,e)) to every point */
/*in the volume (a,b,c) (below the surface) */
/*Parameter: *sample_topo:= pointer on samples of the grid cell centers (in z-direction), where the topography starts; (The samples in z-direction can lie */ 
/*							outside of the subsequent grid)*/
/*              *grid   := Pointer on the grid structure*/
/*               *mod   := model (dependent on the parameter the velocity, density or resistivity model will be determined)*/
/*       kind_of_model  := Identify the kind of the model(velocity =1; density =2; resistivity =3)*/


int CalcSurfDepMod(long *sample_topo,GRID_STRUCT *grid,double *mod, int kind_of_model)
{
	long a,b,c,d,e;
	long nx,ny,nz,nyz;
	int max_int,dmin,emin;
	float *distance, abstand,h;
	float abs_gradient;
	GRADIENT grad;

	#define m(x,y,z) mod[nyz*(x) + nz*(y) + (z)]
	#define sample(x,y) sample_topo[ny*(x)+(y)]
	#define dist(x,y,z) distance[nyz*(x) + nz*(y) + (z)]

	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nyz = ny*nz;

	/*Choose the model*/
	if(kind_of_model == 1)		/*s.velocity*/
		grad = grid->grad_vel;
	else if(kind_of_model == 2) /*density*/
		grad = grid->grad_dens;
	else if(kind_of_model == 3)	/*resistivity*/
		grad = grid->grad_res;
	else
	{
		printf("Choosen model is unknown\n");
		exit(0);
	}

	h =grid->h;

	/*Calculate the abs. value of the gradient:*/
	abs_gradient = (float)sqrt((grad.g[0])*(grad.g[0]) + (grad.g[1])*(grad.g[1]) + (grad.g[2])*(grad.g[2]));

	distance = (float *)memory(NULL,nx*ny*nz,sizeof(float),"CalcSurfDepMod");
	
	for(a=0; a<nx; a++)
		for(b=0; b<ny;b++)
			for(c=0;c<nz;c++)
			{
				dist(a,b,c) = (float)9999999.9;
			}

	/*Determine the shortest distance from the surface to every point in the subsurface:*/
	printf("Start calculating the gradient model:\n");
	printf("----------------\n");

		for(a=0; a<nx; a++)
		{
			if(a%10 == 0) printf("For %d of %d slices the gradient model is determined\n",a,nx);

			for(b=0; b<ny;b++)
				for(c=0;c<nz;c++)
				{
					/*Consider only samples below the topography*/
					if(c - sample(a,b) >= 0)
					{
						max_int = c - sample(a,b);
						if(a-max_int>0) dmin= a-max_int;
						else dmin = 0;
						if(b-max_int>0) emin= b-max_int;
						else emin =0;

					/*Consider only useful (x,y) surface points combinations*/
					for(d=dmin; d<(a+max_int+1) && d<nx; d++)
						for(e=emin; e<(b+max_int+1) && e<ny; e++) 
						{
							/*Use only surface points WITH topography information*/
						   if(sample(d,e) < nz)
								abstand = (float)((d-a)*(d-a) + (e-b)*(e-b) + (sample(d,e)-c)*(sample(d,e)-c));
						   else
							    abstand = (float)9999999.9;

						   if(dist(a,b,c) > abstand) dist(a,b,c)=abstand;
						}
					}
				}
		}


	/*Calculate the model:*/
	for(a=0; a<nx; a++)
		for(b=0; b<ny;b++)
			for(c=0;c<nz;c++)
			{
				m(a,b,c) = grad.value_0 + abs_gradient*h*(float)sqrt(dist(a,b,c));
					if(m(a,b,c) < grad.min) m(a,b,c) = grad.min;
					if(m(a,b,c) > grad.max) m(a,b,c) = grad.max;
			}
	
	printf("model is calculated\n");
	printf("----------------\n");

	free(distance);

	#undef m
	#undef sample
	#undef dist

	return(1);
}

/*---------------------------------------------------------------------*/
/*Calculate the starting model by means of a origin and gradient parameters*/
/*Parameter:     *grid   := Pointer on the grid structure*/
/*				 *mod    := 
/*           kind_of_mod :=  Identify the kind of the model(velocity =1; density =2; resistivity =3)*/

int	CalcOrgDepMod(GRID_STRUCT *grid, double *mod, int kind_of_model)
{
	long a,b,c,nz,nyz;
	float h;
	GRADIENT grad;

	#define m(x,y,z) mod[nyz*(x) + nz*(y) + (z)]

	/*Choose the model*/
	if(kind_of_model == 1)		/*s.velocity*/
		grad = grid->grad_vel;
	else if(kind_of_model == 2) /*density*/
		grad = grid->grad_dens;
	else if(kind_of_model == 3)	/*resistivity*/
		grad = grid->grad_res;
	else
	{
		printf("Choosen model is unknown\n");
		exit(0);
	}

	h = grid->h;
	nz = grid->nz;
	nyz = (grid->ny) * (grid->nz);

		/*Calculate the model:*/
		for(a=0; a<(grid->nx); a++)
			for(b=0; b<(grid->ny);b++)
				for(c=0; c<(grid->nz);c++)
				{
					m(a,b,c) = (grad.value_0) + ((grid->org[0])+a*h-(grad.org[0]))*(grad.g[0]) + ((grid->org[1])+b*h-(grad.org[1]))*(grad.g[1]) + ((grid->org[2])+c*h-(grad.org[2]))*(grad.g[2]);
					if(m(a,b,c) < (grad.min)) m(a,b,c) = (grad.min);
					if(m(a,b,c) > (grad.max)) m(a,b,c) = (grad.max);
				}
		
		printf("Model is calculated\n");
		printf("----------------\n");

	#undef m

	return(1);
}


/*-------------------------------------------------------------*/
/* Add the boundaries to the models*/
/* Parameter: *grid       := pointer on the grid structure*/
/*            *vel1,vel2  := Pointer on the velocity model before and after adding the boundaries*/
/*            *border_index1,*border_index2 := Highlighting the cells that are not belonging to the boundaries or the air*/

int AddBoundaries(GRID_STRUCT *grid, double *vel1, int *border_index1, double *vel2, int *border_index2)
{  
	int a,b,c;
	long nx,ny,nz,nyz,nx2,ny2,nz2,nyz2,nborder;

	#define v(x,y,z) vel1[nyz*(x) + nz*(y) + (z)]
	#define b_index(x,y,z) border_index1[nyz*(x) + nz*(y) + (z)]
	#define v2(x,y,z) vel2[nyz2*(x) + nz2*(y) + (z)]
	#define b_index2(x,y,z) border_index2[nyz2*(x) + nz2*(y) + (z)]

	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nborder = grid->nborder;
	nyz= ny*nz;

	nx2=nx+2*nborder;
	ny2=ny+2*nborder;
	nz2=nz+2*nborder;
	nyz2=ny2*nz2;

	/*Initialize velocity cube with boundaries*/
	for(a=0; a<nx2; a++)
		for(b=0; b<ny2;b++)
			for(c=0; c<nz2;c++)
			{
						b_index2(a,b,c) = 0;
						      v2(a,b,c) = (grid->vborder);
			}

	/*Introduce the former calculated velocity cube in the velocitiy cube with boundaries*/
	for(a=0; a<nx; a++)
		for(b=0; b<ny;b++)
			for(c=0; c<nz;c++)
			{
						b_index2(a+nborder,b+nborder,c+nborder) = b_index(a,b,c);
						      v2(a+nborder,b+nborder,c+nborder) = v(a,b,c);
			}

	

	printf("The boundaries are added to the velocity model\n");
	printf("----------------\n");

	#undef v
	#undef b_index
	#undef v2
	#undef b_index2

	return(1);
}


