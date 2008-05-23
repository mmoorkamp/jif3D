#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "dcomplex.h"
#include "inv3d.h"
#include "fkt.h"


/*------------------------------------------------------------*/
/* Calculate the slowness dependent derivatives of the density/resistivity*/
/*Parameter:    *grid   := Pointer on the grid structure*/
/*			    Y2 -Y1  := Perturbation of the density OR resistivity for all forward cells*/
/*				rel_vdr	:= Structure that organizes the relationship of velocity and density and velocity and resistivity*/
/* Output: Perturbation of the  slowness*/

double CalcDerivDensSlow(GRID_STRUCT *grid, double *Y2, double *Y1, REL_VDR_STRUCT rel_vdr)
{
	long a,b,c,a1,b1,c1,nx,ny,nz,nyz;
	long ny3,nz3,nyz3;
	long ny2,nz2,nyz2;

	double s1,s2,f1,f2;
	double deltaX; /*Perturbation of the slowness in [s/m]*/
	double eps = 0.00001;

	/*corresponds to the chosen ds in s/m*/ 
	#define seis(x,y,z) grid->slow[nyz3*(x) + nz3*(y) + (z)] 
	#define bindex(x,y,z) grid->border_index[nyz2*(x) + nz2*(y) + (z)] 

	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nyz = ny*nz;

	ny2 = ny + (2*grid->nborder);
	nz2 = nz + (2*grid->nborder);
	nyz2 = ny2*nz2;

	ny3 = ny2 + 1;
	nz3 = nz2 + 1;
	nyz3 = ny3*nz3;

	/*Loop over all cells in the model*/
	for(a=0; a<nx; a++)
		for(b=0; b<ny;b++)
			for(c=0;c<nz;c++)
			{
				a1 = a + grid->nborder;
				b1 = b + grid->nborder;
				c1 = c + grid->nborder;

					/*slowness values s1 and s2*/
				s1 = (seis(a1,b1,c1)/grid->h) - (eps/2);
				s2 = (seis(a1,b1,c1)/grid->h) + (eps/2);
					
				if(bindex(a1,b1,c1) != 0) /*Not considered for cells in the air*/
				{
					/*gravity or resistivity values f(s1) and f(s2) in [mgal] or [ohmm]*/
					f1 = VeltoDensRes(1/s1, rel_vdr);
					f2 = VeltoDensRes(1/s2, rel_vdr);

					/*Check if the resistivity/density values are close to 0*/
					if(f1 < 0 || f2 < 0)
					{
						printf("WARNING !!:\n");
						printf("After determing the resistivity/density perturbations from the velocity\n");
						printf("perturbations, some restivity OR density values occur that are smaller than 0 !!\n");
					}

					/*Calculating the perturbations*/
					Y1[nyz*a + nz*b + c] = f1;
					Y2[nyz*a + nz*b + c] = f2;
				}

			}

	deltaX = eps;

	#undef seis
	#undef bindex

	return(deltaX);
}


/*------------------------------------------------------------*/
/* Make a velocity, density or resistivty field with constant perturbations:*/ 
/*Required to calculate the derivatives of the gravity/MT data with the density/resistivity*/
/* ("simple" gravity or MT inversions)*/
/*Parameter:    *grid   := Pointer on the grid structure*/
/*             *values  := Density or resistivity field + constant perturbation*/
/*				   eps	:= constant perturbation*/
/*		 kind_of_model	:= 1 == slowness; 2 == density; 3 == resistivity*/

int FillwithValues(GRID_STRUCT *grid, double *values, double eps, int kind_of_model)
{
	long a,b,c,nx,ny,nz,nyz;
	long ny2,nz2,nyz2;
	long ny3,nz3,nyz3;

	nx = grid->nx;
	ny = grid->ny;
	nz = grid->nz;
	nyz = ny*nz;

	ny2 = ny + 2*grid->nborder;
	nz2 = nz + 2*grid->nborder;
	nyz2 = ny2*nz2;

	nz3 = nz + 2*grid->nborder + 1;
	ny3 = ny + 2*grid->nborder + 1;
	nyz3 = ny3*nz3;

	/*Loop over all cells in the model*/
	for(a=0; a<nx; a++)
		for(b=0; b<ny;b++)
			for(c=0;c<nz;c++)
			{
				if(kind_of_model == 1) /*slowness model(The border cells remain unaffected)*/
				{
					if(grid->border_index[nyz2*(a + grid->nborder) + nz2*(b + grid->nborder) + (c + grid->nborder)] != 0)
						values[nyz3*(a + grid->nborder) + nz3*(b + grid->nborder) + (c + grid->nborder)] = eps + grid->slow[nyz3*(a + grid->nborder) + nz3*(b + grid->nborder) + (c + grid->nborder)];
					else
						values[nyz3*(a + grid->nborder) + nz3*(b + grid->nborder) + (c + grid->nborder)] = grid->slow[nyz3*(a + grid->nborder) + nz3*(b + grid->nborder) + (c + grid->nborder)];
				}
				else if(kind_of_model == 2) /*density model*/
				{
					if(grid->border_index[nyz2*(a + grid->nborder) + nz2*(b + grid->nborder) + (c + grid->nborder)] != 0)
						values[nyz*a + nz*b + c] = eps + grid->dens[nyz*a + nz*b + c];
					else
						values[nyz*a + nz*b + c] = grid->dens[nyz*a + nz*b + c];
				}
				else	/*resistivity model*/
				{
					if(grid->border_index[nyz2*(a + grid->nborder) + nz2*(b + grid->nborder) + (c + grid->nborder)] != 0)
						values[nyz*a + nz*b + c] = eps + grid->res[nyz*a + nz*b + c];
					else
						values[nyz*a + nz*b + c] = grid->res[nyz*a + nz*b + c];
				}
			}

	return(0);
}

/*-------------------------------------------------------------*/
/*Calculating the slowness derivatives for the gravity following I.Trinks and R.Hobbs, 200?, Gravity modelling based on small cells, Earth Sciences*/
/*Parameter:	geo  := Geometry structure  */
/*              grid  := Grid structure */
/*              *data  := Pointer on data structure*/
/*				*Rho2 -Rho1 := Perturbation of the densities in g/cm^3*/
/*				deltaSlow := Perturbation of the slowness in s/m*/

int DerivativesGrav(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, double *Rho2, double *Rho1, double deltaSlow)
{
	long   i,a,b,c,ix,iy,iz,ixiy,n1;
	long   nx,ny,nz,nyz,ny2,nz2,nyz2,nborder;
	double org1,eps;
	double xs,ys,zs,xys;								/*position of the station in m*/
	double x,x1,x2,y,y1,y2,z,xy,xy1,xy2;				/*position of the cell centers in m*/
	double dx,dy,dz,dx1,dx2,dy1,dy2,dxy,dxz,dyz,dxy1,dxy2,dxyz;	/*distance between the cell centers and the gravity stations in m*/

	double const1;
	double Grav2, Grav1, *tmp_delta_grav;
	double tmp_grav;
	double tmp_grav11, tmp_grav12, tmp_grav21, tmp_grav22;
	double tmp_grav31, tmp_grav32, tmp_grav41, tmp_grav42;

	char fname[40];
	FILE *ouf;

	#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/
	#define rho13(x,y,z) Rho1[((grid.ny)*(grid.nz))*(x) + (grid.nz)*(y) + (z)]
	#define rho23(x,y,z) Rho2[((grid.ny)*(grid.nz))*(x) + (grid.nz)*(y) + (z)]
	#define tmp_deltagrav3(x,y,z) tmp_delta_grav[(grid.ny*grid.nz)*(x) + (grid.nz)*(y) + (z)]
	#define rho1(x_y,z) Rho1[(grid.nz)*(x_y) + (z)]
	#define rho2(x_y,z) Rho2[(grid.nz)*(x_y) + (z)]
	#define tmp_deltagrav(x_y,z) tmp_delta_grav[grid.nz*(x_y) + (z)]
	#define FACTOR 100 /*Specify the distance to the cell centers (unit: number of cell length), for which the response will be calculated exactly*/

	nx = grid.nx;
	ny = grid.ny;
	nz = grid.nz;
	nyz = ny*nz;

	ny2 = ny + 2*grid.nborder;
	nz2 = nz + 2*grid.nborder;
	nyz2 = ny2*nz2;

	eps = 0.0000000001;

	if(data->ndata_grav == 0)
		printf("!!WARNING!! NO gravity data exists but the calcultion of frechet derivatives is activated\n\n");


	tmp_delta_grav = (double *)memory(NULL,(nx*ny*nz),sizeof(double),"DerivativeGrav");
	for(i=0;i<(grid.nx*grid.ny*grid.nz);i++)
		tmp_delta_grav[i] = 0;

	/*Decide if the calculation is 1D, 2D or 3D by considering the number of cells in the horizontal directions:*/

	/******************************************************/
					/*1D-calculation*/
	/******************************************************/
	/* ATTENTION!!! 1-D WAS NOT TESTED !!!*/

	if(grid.nx == 1 && grid.ny == 1)
	{
		/*Loop over all stations*/
		for(i=0;i<geo.nstat_grav;i++)
		{
			/*Positions of the gravity stations:*/
			xs = (double)geo.x[(data->gravs[i]-1)];
			ys = (double)geo.y[(data->gravs[i]-1)];
			zs = (double)geo.z[(data->gravs[i]-1)];

			/**********************/
			for(iz=0;iz<grid.nz;iz++)
			{
				z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
				dz = z - zs;	/*z-component of distance between grid cell and station*/

				/*Calculate the gravitation effect from the cell at the station*/
				Grav2 = 2*PI*G_KONST*Rho2[iz]*(grid.h); /*double-infinite sheet*/
				Grav1 = 2*PI*G_KONST*Rho1[iz]*(grid.h);
				tmp_delta_grav[iz] = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/
			}

			sprintf(fname,"tmp_grav_deriv%d.dat",i);
			ouf = fopen(fname,"wb");

			fwrite(&(nx),sizeof(long),1,ouf);
			fwrite(&(ny),sizeof(long),1,ouf);
			fwrite(&(nz),sizeof(long),1,ouf);
			nborder = grid.nborder;
			fwrite(&(nborder),sizeof(long),1,ouf);

			/*Write the derivatives in a temporary file*/
			for(a=0;a<nx;a++)
				for(b=0;b<ny;b++)
					for(c=0;c<nz;c++)
					{
						tmp_grav = tmp_delta_grav[nyz*a + nz*b +c];
						fwrite(&(tmp_grav),sizeof(double),1,ouf);					/*Write out the derivatives*/
						
					}

			fclose(ouf);

		if(i%10 == 0)
			printf("Derivatives are determined for %d of %d stations\n", i+1, geo.nstat_grav);

		}
	}
	
	/******************************************************/
					/*2D-calculation*/
	/******************************************************/

	/*calculate the product of cell size and gravitational const.: dx*dy*G_konst in m^5/(Kg*s^2)*/
    const1 = (grid.h)*(grid.h)*G_KONST;

	if(grid.nx == 1 || grid.ny == 1)
	{

		/*Loop over all stations*/
		for(i=0;i<geo.nstat_grav;i++)
		{

			/*Positions of the gravity stations:*/
			xs = (double)geo.x[(data->gravs[i]-1)];
			ys = (double)geo.y[(data->gravs[i]-1)];
			zs = (double)geo.z[(data->gravs[i]-1)];

			/*Distinguish between the cases nx==1 and ny==1 */
			if(grid.nx == 1)
			{
				n1 = grid.ny;		/*First running index*/
				org1 = (double)grid.org[1]; /*choose the x-component of the origin for further calculations*/
				xys = ys;			/*y coordinate of position of the gravity station will be used*/
			}
			else
			{
				n1 = grid.nx;		/*First running index*/
				org1 = (double)grid.org[0]; /*choose the x-component of the origin for further calculations*/
				xys = xs;			/*x coordinate of position of the gravity station will be used*/
			}

			/**********************/
			/*1.Case: Cells at the borders of the 2-D grid*/

			/*left side*/
			ixiy = 0;
			xy1 = org1 + (ixiy+0.5)*(grid.h);		 /*x resp.y component of the cell centers in m*/
			dxy1 = fabs(xy1 - xys);				 /*xy-component of distance between grid cell and station*/
			/*right side*/
			ixiy = n1-1;
			xy2 = org1 + (ixiy-0.5)*(grid.h);		 /*x resp.y component of the cell centers in m*/
			dxy2 = fabs(xy2 - xys);				 /*xy-component of distance between grid cell and station*/

			/*Loop over all cell at the left (ix,iy = 0) and (ix,iy = nx-1,ny-1) right border*/
			for(iz=0;iz<grid.nz;iz++)
				{
					z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					dz = z - zs;				/*z-component of distance between grid cell and station*/
					if(dz == 0)
						dz = eps;

					/*Calculate the gravitation effect from the cell at the station*/
					/*left border*/
					if(dz >= 0)
					{
						Grav1 = (2*G_KONST*rho1(0,iz)*grid.h)*((PI/2) - atan(dxy1/dz)); /*semi-infinite horizontal sheet*/
						Grav2 = (2*G_KONST*rho2(0,iz)*grid.h)*((PI/2) - atan(dxy1/dz));
					}
					else
					{
						Grav1 = (2*G_KONST*rho1(0,iz)*grid.h)*((-PI/2) - atan(dxy1/dz)); /*semi-infinite horizontal sheet*/
						Grav2 = (2*G_KONST*rho2(0,iz)*grid.h)*((-PI/2) - atan(dxy1/dz));
					}
					tmp_deltagrav(0,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/				
					/*right border*/
					if(dz >= 0)
					{
						Grav1 = (2*G_KONST*rho1(n1-1,iz)*grid.h)*((PI/2) - atan(dxy2/dz)); /*semi-infinite horizontal sheet*/
						Grav2 = (2*G_KONST*rho2(n1-1,iz)*grid.h)*((PI/2) - atan(dxy2/dz));
					}
					else
					{
						Grav1 = (2*G_KONST*rho1(n1-1,iz)*grid.h)*((-PI/2) - atan(dxy2/dz)); /*semi-infinite horizontal sheet*/
						Grav2 = (2*G_KONST*rho2(n1-1,iz)*grid.h)*((-PI/2) - atan(dxy2/dz));
					}
					tmp_deltagrav(n1-1,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/
					
				}

			/**********************/
			/*2.Case: Internal cells of the 2-D grid*/
			/*Loop over all remaining cells*/
			for(ixiy=1;ixiy<(n1-1);ixiy++)
			{
				xy = org1 + ixiy*(grid.h);		 /*x resp.y component of the cell centers in m*/
				dxy = xy -xys;					 /*xy-component of distance between grid cell and station*/

				for(iz=0;iz<grid.nz;iz++)
				{
					z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					dz = z - zs;					/*z-component of distance between grid cell and station*/
					if(dz == 0)
						dz = eps;

					dxyz = sqrt(dxy*dxy + dz*dz);

					/*Calculate the gravitation from the cell at the station*/
					if(dxy < (grid.h * FACTOR))
					{
						/*exact calculation of attraction of two-dimensional rectangular blocks for cells close to the corresponding station following Heiland (1946)*/
						Grav1 = Response_2D_Rectangle(rho1(ixiy,iz),dxy-(grid.h/2),dxy+(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
						Grav2 = Response_2D_Rectangle(rho2(ixiy,iz),dxy-(grid.h/2),dxy+(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
					}
					else
					{
						/*infinitely long horizontal rod*/
						Grav1 = (2*const1*rho1(ixiy,iz))/(dz*(1 +((dxy*dxy)/(dz*dz)))); 
						Grav2 = (2*const1*rho2(ixiy,iz))/(dz*(1 +((dxy*dxy)/(dz*dz))));
					}

					tmp_deltagrav(ixiy,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/

				}
			}

			sprintf(fname,"tmp_grav_deriv%d.dat",i);
			ouf = fopen(fname,"wb");

			fwrite(&(nx),sizeof(long),1,ouf);
			fwrite(&(ny),sizeof(long),1,ouf);
			fwrite(&(nz),sizeof(long),1,ouf);
			nborder = grid.nborder;
			fwrite(&(nborder),sizeof(long),1,ouf);
	
			/*Write the derivatives in a temporary file*/
			for(a=0;a<nx;a++)
				for(b=0;b<ny;b++)
					for(c=0;c<nz;c++)
					{
						tmp_grav = tmp_delta_grav[nyz*a + nz*b +c];
						fwrite(&(tmp_grav),sizeof(double),1,ouf);					/*Write out the derivatives*/
					}

			fclose(ouf);

		if(i%10 == 0)
			printf("Derivatives are determined for %d of %d stations\n", i+1, geo.nstat_grav);

		}
	}


	/******************************************************/
					/*3D-calculation*/
	/******************************************************/


	else
	{

		/*Loop over all stations*/
		for(i=0;i<geo.nstat_grav;i++)
		{
			/*Positions of the gravity stations:*/
			xs = (double)geo.x[(data->gravs[i]-1)];
			ys = (double)geo.y[(data->gravs[i]-1)];
			zs = (double)geo.z[(data->gravs[i]-1)];

			/*left left (ix==0)*/
			ix = 0;
			x1 = (double)grid.org[0] + (ix+0.5)*(grid.h);	 /*x component of the cell centers in m*/
			dx1 = fabs(x1 - xs);				 /*x-component of distance between grid cell and station*/
			/*right side (ix = nx-1)*/
			ix = grid.nx-1;
			x2 = (double)grid.org[0] + (ix-0.5)*(grid.h);	 /*x component of the cell centers in m*/
			dx2 = fabs(x2 - xs);				 /*x-component of distance between grid cell and station*/
			/*front side (iy==0)*/
			iy = 0;
			y1 = (double)grid.org[1] + (iy+0.5)*(grid.h);	 /*y component of the cell centers in m*/
			dy1 = fabs(y1 - ys);				 /*y-component of distance between grid cell and station*/
			/*back side (iy = ny-1)*/
			iy = grid.ny-1;
			y2 = (double)grid.org[1] + (iy-0.5)*(grid.h);	 /*y component of the cell centers in m*/
			dy2 = fabs(y2 - ys);				 /*y-component of distance between grid cell and station*/


			/**********************/
			/*1.Case: Cells at the edges of the 3-D grid*/

			/*Loop over all cell at the edges (ix=0,iy=0; ix=0,iy=ny-1; ix=nx-1,iy=0; ix=nx-1,iy=ny-1)*/
			for(iz=0;iz<grid.nz;iz++)
			{
				z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
				dz = z - zs;					/*z-component of distance between grid cell and station*/
				if(dz == 0)
						dz = eps;

				/*Calculate the gravitation effect from the cell at the station*/
				/*left front edge*/
				dxyz = sqrt(dx1*dx1 + dy1*dy1 + dz*dz);
				if(dxyz == 0)
					dxyz = eps;
				if(dz >= 0)
				{
					Grav1 = (G_KONST*rho13(0,0,iz)*grid.h)*((PI/2) - atan(dx1/dz) - atan(dy1/dz) + atan((dx1*dy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					Grav2 = (G_KONST*rho23(0,0,iz)*grid.h)*((PI/2) - atan(dx1/dz) - atan(dy1/dz) + atan((dx1*dy1)/(dz*dxyz)));
				}
				else
				{
					Grav1 = (G_KONST*rho13(0,0,iz)*grid.h)*((-PI/2) - atan(dx1/dz) - atan(dy1/dz) + atan((dx1*dy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					Grav2 = (G_KONST*rho23(0,0,iz)*grid.h)*((-PI/2) - atan(dx1/dz) - atan(dy1/dz) + atan((dx1*dy1)/(dz*dxyz)));
				}
				tmp_deltagrav3(0,0,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/

				/*right front edge*/
				dxyz = sqrt(dx2*dx2 + dy1*dy1 + dz*dz);
				if(dxyz == 0)
					dxyz = eps;
				if(dz >= 0)
				{
					Grav1 = (G_KONST*rho13(grid.nx-1,0,iz)*grid.h)*((PI/2) - atan(dx2/dz) - atan(dy1/dz) + atan((dx2*dy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					Grav2 = (G_KONST*rho23(grid.nx-1,0,iz)*grid.h)*((PI/2) - atan(dx2/dz) - atan(dy1/dz) + atan((dx2*dy1)/(dz*dxyz)));
				}
				else
				{
					Grav1 = (G_KONST*rho13(grid.nx-1,0,iz)*grid.h)*((-PI/2) - atan(dx2/dz) - atan(dy1/dz) + atan((dx2*dy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					Grav2 = (G_KONST*rho23(grid.nx-1,0,iz)*grid.h)*((-PI/2) - atan(dx2/dz) - atan(dy1/dz) + atan((dx2*dy1)/(dz*dxyz)));
				}
				tmp_deltagrav3(grid.nx-1,0,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/				

				/*left back edge*/
				dxyz = sqrt(dx1*dx1 + dy2*dy2 + dz*dz);
				if(dxyz == 0)
					dxyz = eps;
				if(dz >= 0)
				{
					Grav1 = (G_KONST*rho13(0,grid.ny-1,iz)*grid.h)*((PI/2) - atan(dx1/dz) - atan(dy2/dz) + atan((dx1*dy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					Grav2 = (G_KONST*rho23(0,grid.ny-1,iz)*grid.h)*((PI/2) - atan(dx1/dz) - atan(dy2/dz) + atan((dx1*dy2)/(dz*dxyz)));
				}
				else
				{
					Grav1 = (G_KONST*rho13(0,grid.ny-1,iz)*grid.h)*((-PI/2) - atan(dx1/dz) - atan(dy2/dz) + atan((dx1*dy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					Grav2 = (G_KONST*rho23(0,grid.ny-1,iz)*grid.h)*((-PI/2) - atan(dx1/dz) - atan(dy2/dz) + atan((dx1*dy2)/(dz*dxyz)));
				}
				tmp_deltagrav3(0,grid.ny-1,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/				

				/*right back edge*/
				dxyz = sqrt(dx2*dx2 + dy2*dy2 + dz*dz);
				if(dxyz == 0)
					dxyz = eps;
				if(dz >= 0)
				{
					Grav1 = (G_KONST*rho13(grid.nx-1,grid.ny-1,iz)*grid.h)*((PI/2) - atan(dx2/dz) - atan(dy2/dz) + atan((dx2*dy2)/(dz*dxyz))) ; /*semi-semi-infinite horizontal sheet*/
					Grav2 = (G_KONST*rho23(grid.nx-1,grid.ny-1,iz)*grid.h)*((PI/2) - atan(dx2/dz) - atan(dy2/dz) + atan((dx2*dy2)/(dz*dxyz))) ; 
				}
				else
				{
					Grav1 = (G_KONST*rho13(grid.nx-1,grid.ny-1,iz)*grid.h)*((-PI/2) - atan(dx2/dz) - atan(dy2/dz) + atan((dx2*dy2)/(dz*dxyz))) ; /*semi-semi-infinite horizontal sheet*/
					Grav2 = (G_KONST*rho23(grid.nx-1,grid.ny-1,iz)*grid.h)*((-PI/2) - atan(dx2/dz) - atan(dy2/dz) + atan((dx2*dy2)/(dz*dxyz))) ; 
				}
				tmp_deltagrav3(grid.nx-1,grid.ny-1,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/				

			}

						
			/**********************/
			/*2.Case: Cells at the sides (except for the edges) of the 3-D grid*/

			/*Loop over all cell at the right and left sides (ix = 0 and ix = nx -1)*/
			for(iy=1;iy<(grid.ny-1);iy++)
			{
				y = grid.org[1] + iy*(grid.h);		 /*y component of the cell centers in m*/
				dy = fabs(y - ys);					 /*y-component of distance between grid cell and station*/

				for(iz=0;iz<grid.nz;iz++)
				{
					z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					dz = z - zs;					/*z-component of distance between grid cell and station*/

					/*Calculate the gravitation from the cell at the station*/
					/*left side*/
					dxyz = sqrt(dx1*dx1 + dy*dy + dz*dz);
					if(dxyz == 0)
						dxyz = eps;
					dyz = dz*dz + dy*dy;
					if(dyz == 0)
						dyz = eps;
					Grav1 = (G_KONST*rho13(0,iy,iz)*(grid.h)*(grid.h))*(dz/dyz)*(1-(dx1/dxyz)); /*semi-infinitely long horizontal rod*/
					Grav2 = (G_KONST*rho23(0,iy,iz)*(grid.h)*(grid.h))*(dz/dyz)*(1-(dx1/dxyz));
					tmp_deltagrav3(0,iy,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/

					/*right side*/
					dxyz = sqrt(dx2*dx2 + dy*dy + dz*dz);
					if(dxyz == 0)
						dxyz = eps;
					dyz = dz*dz + dy*dy;
					if(dyz == 0)
						dyz = eps;
					Grav1 = (G_KONST*rho13(grid.nx-1,iy,iz)*(grid.h)*(grid.h))*(dz/dyz)*(1-(dx2/dxyz)); /*semi-infinitely long horizontal rod*/
					Grav2 = (G_KONST*rho23(grid.nx-1,iy,iz)*(grid.h)*(grid.h))*(dz/dyz)*(1-(dx2/dxyz)); 
					tmp_deltagrav3(grid.nx-1,iy,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/				

				}
			}

			/*Loop over all cell at the front and back sides (iy = 0 and iy = ny -1)*/
			for(ix=1;ix<(grid.nx-1);ix++)
			{
				x = grid.org[0] + ix*(grid.h);		 /*x component of the cell centers in m*/
				dx = fabs(x - xs);					 /*x-component of distance between grid cell and station*/

				for(iz=0;iz<grid.nz;iz++)
				{
					z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					dz = z - zs;					/*z-component of distance between grid cell and station*/

					/*Calculate the gravitation from the cell at the station*/
					/*front side*/
					dxyz = sqrt(dx*dx + dy1*dy1 + dz*dz);
					if(dxyz == 0)
						dxyz = eps;
					dxz = dz*dz + dx*dx;
					if(dxz == 0)
						dxz = eps;
					Grav1 = (G_KONST*rho13(ix,0,iz)*(grid.h)*(grid.h))*(dz/dxz)*(1-(dy1/dxyz)); /*semi-infinitely long horizontal rod*/
					Grav2 = (G_KONST*rho23(ix,0,iz)*(grid.h)*(grid.h))*(dz/dxz)*(1-(dy1/dxyz));
					tmp_deltagrav3(ix,0,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/				

					/*back side*/
					dxyz = sqrt(dx*dx + dy2*dy2 + dz*dz);
					if(dxyz == 0)
						dxyz = eps;
					dxz = dz*dz + dx*dx;
					if(dxz == 0)
						dxz = eps;
					Grav1 = (G_KONST*rho13(ix,grid.ny-1,iz)*(grid.h)*(grid.h))*(dz/dxz)*(1-(dy2/dxyz)); /*semi-infinitely long horizontal rod*/
					Grav2 = (G_KONST*rho23(ix,grid.ny-1,iz)*(grid.h)*(grid.h))*(dz/dxz)*(1-(dy2/dxyz));
					tmp_deltagrav3(ix,grid.ny-1,iz) = 1.0E8 * (Grav2 -Grav1)/deltaSlow; /*Calculate the total derivatives*/									

				}

			}

			/**********************/
			/*3.Case: Internal cells of the 3-D grid*/
			/*Loop over all remaining cells*/
			for(ix=1;ix<grid.nx-1;ix++)
			{
				x = grid.org[0] + ix*(grid.h);		 /*x component of the cell centers in m*/
				dx = x - xs;						 /*x-component of distance between grid cell and station*/

				for(iy=1;iy<grid.ny-1;iy++)
				{
					y = grid.org[1] + iy*(grid.h);		 /*y component of the cell centers in m*/
					dy = y - ys;						 /*y-component of distance between grid cell and station*/

					for(iz=0;iz<grid.nz;iz++)
					{
						z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
						dz = z - zs;					/*z-component of distance between grid cell and station*/

						/*Calculate the gravitation from the cell at the station*/
						dxyz = sqrt(dx*dx + dy*dy + dz*dz);

						/*For cells nearby the stations the response will be calculated exactly (see Nagy, 1966, Geophysics)*/
						if(dxyz < (grid.h * FACTOR))
						{ 
							/********************************************************/
							/*Case I: The X and Y coordinates of the station are not located within the edges of the rectangle*/
							if(fabs(dx) >= grid.h/2 && fabs(dy) >= grid.h/2)
							{
								Grav1 = Response_3D_Rectangle(rho13(ix,iy,iz),dx-(grid.h/2),dx +(grid.h/2),dy-(grid.h/2),dy +(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
								Grav2 = Response_3D_Rectangle(rho23(ix,iy,iz),dx-(grid.h/2),dx +(grid.h/2),dy-(grid.h/2),dy +(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
							}

							/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
							else if(fabs(dx) < grid.h/2 && fabs(dy) >= grid.h/2)
							{
								/*The integral have to be splitted (see Nagy, 1966)*/
								tmp_grav11 = Response_3D_Rectangle(rho13(ix,iy,iz),dx-(grid.h/2),0,dy-(grid.h/2),dy +(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav21 = Response_3D_Rectangle(rho13(ix,iy,iz),0,dx+(grid.h/2),dy-(grid.h/2),dy +(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));

								tmp_grav12 = Response_3D_Rectangle(rho23(ix,iy,iz),dx-(grid.h/2),0,dy-(grid.h/2),dy +(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav22 = Response_3D_Rectangle(rho23(ix,iy,iz),0,dx+(grid.h/2),dy-(grid.h/2),dy +(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));

								Grav1 = tmp_grav11 + tmp_grav12;
								Grav2 = tmp_grav21 + tmp_grav22;
								
							}

							/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
							else if(fabs(dx) >= grid.h/2 && fabs(dy) < grid.h/2)
							{
								/*The integral have to be splitted (see Nagy, 1966)*/
								tmp_grav11 = Response_3D_Rectangle(rho13(ix,iy,iz),dx-(grid.h/2),dx+(grid.h/2),dy-(grid.h/2),0,dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav21 = Response_3D_Rectangle(rho13(ix,iy,iz),dx-(grid.h/2),dx+(grid.h/2),0,dy+(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));

								tmp_grav12 = Response_3D_Rectangle(rho23(ix,iy,iz),dx-(grid.h/2),dx+(grid.h/2),dy-(grid.h/2),0,dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav22 = Response_3D_Rectangle(rho23(ix,iy,iz),dx-(grid.h/2),dx+(grid.h/2),0,dy+(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
	
								Grav1 = tmp_grav11 + tmp_grav12;
								Grav2 = tmp_grav21 + tmp_grav22;
								
							}
							/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
							else
							{
								/*The integral have to be splitted (see Nagy, 1966)*/
								tmp_grav11 = Response_3D_Rectangle(rho13(ix,iy,iz),0,dx+(grid.h/2),dy-(grid.h/2),0,dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav21 = Response_3D_Rectangle(rho13(ix,iy,iz),0,dx+(grid.h/2),0,dy+(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav31 = Response_3D_Rectangle(rho13(ix,iy,iz),dx-(grid.h/2),0,dy-(grid.h/2),0,dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav41 = Response_3D_Rectangle(rho13(ix,iy,iz),dx-(grid.h/2),0,0,dy+(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));

								tmp_grav12 = Response_3D_Rectangle(rho23(ix,iy,iz),0,dx+(grid.h/2),dy-(grid.h/2),0,dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav22 = Response_3D_Rectangle(rho23(ix,iy,iz),0,dx+(grid.h/2),0,dy+(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav32 = Response_3D_Rectangle(rho23(ix,iy,iz),dx-(grid.h/2),0,dy-(grid.h/2),0,dz-(grid.h/2),dz +(grid.h/2));
								tmp_grav42 = Response_3D_Rectangle(rho23(ix,iy,iz),dx-(grid.h/2),0,0,dy+(grid.h/2),dz-(grid.h/2),dz +(grid.h/2));
	
								Grav1 = tmp_grav11 + tmp_grav21 + tmp_grav31 + tmp_grav41;
								Grav2 = tmp_grav12 + tmp_grav22 + tmp_grav32 + tmp_grav42;
							}

							/********************************************************/
		
						}
						else
						{	
							/*For cells that are for away from the station (approximated like suggested from Trinks and Hobbs)*/
							Grav1 = (G_KONST*rho13(ix,iy,iz)*(grid.h)*(grid.h)*(grid.h))*(dz/(dxyz*dxyz*dxyz)); /*3D-voxels*/
							Grav2 = (G_KONST*rho23(ix,iy,iz)*(grid.h)*(grid.h)*(grid.h))*(dz/(dxyz*dxyz*dxyz));
						}

						/*Calculate the total derivatives*/
						tmp_deltagrav3(ix,iy,iz) = 1.0E8 * (Grav2 - Grav1)/deltaSlow; 

					}
				}
			}

			sprintf(fname,"tmp_grav_deriv%d.dat",i);
			ouf = fopen(fname,"wb");

			fwrite(&(nx),sizeof(long),1,ouf);
			fwrite(&(ny),sizeof(long),1,ouf);
			fwrite(&(nz),sizeof(long),1,ouf);
			nborder = grid.nborder;
			fwrite(&(nborder),sizeof(long),1,ouf);

			/*Write the derivatives in a temporary file*/
			for(a=0;a<nx;a++)
				for(b=0;b<ny;b++)
					for(c=0;c<nz;c++)
					{
						tmp_grav = tmp_delta_grav[nyz*a + nz*b +c];
						fwrite(&(tmp_grav),sizeof(double),1,ouf);					/*Write out the derivatives*/
					}

			fclose(ouf);

		if(i%10 == 0)
			printf("Derivatives are determined for %d of %d stations\n", i+1, geo.nstat_grav);

		}
	}


	free(tmp_delta_grav);

	#undef G_KONST
	#undef rho1
	#undef rho2
	#undef rho13
	#undef rho23
	#undef tmp_deltagrav
	#undef tmp_deltagrav3
	#undef FACTOR

	return(1);
}

/*-------------------------------------------------------------*/
/*Calculating the slowness derivatives for the 1-D MT model following the script from M.Jegen for layered media*/
/*Parameter:	geo  := Geometry structure  */
/*              grid  := Grid structure */
/*              *data  := Pointer on data structure (the MT data will be (re-)determined in this routine*/
/*				*R2 -R1 := Perturbation of the densities in ohmm*/
/*				deltaSlow := Perturbation of the slowness in s/m*/
/*				*mt := MT structure including the slowness derivatives for the forward cells*/

/*Remark: the magnetic permeability is considered to correspond to the one of the free air*/

#define res(x,y,z) grid.res[((grid.ny)*(grid.nz))*(x) + (grid.nz)*(y) + (z)]
#define R1(x,y,z) R1[((grid.ny)*(grid.nz))*(x) + (grid.nz)*(y) + (z)]
#define R2(x,y,z) R2[((grid.ny)*(grid.nz))*(x) + (grid.nz)*(y) + (z)]
#define M_KONST 4*PI*1E-7 /*magnetic permeability of the air( considered for the earth) in [H/m]*/

int Derivatives1DMT(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, double *R2, double *R1, double deltaSlow, MT_STRUCT *mt)
{
	int nr_xcells, nr_ycells;
	long i,j,k,m,n,r,s,t,pos_count, nr_count;
	long int_tmp, nx_int[2], ny_int[2], nz_int;
	long ny2,nz2,nyz2;
	long *act[3];
	float xs,ys,zs;
	double ang_freq; /*angular frequency*/
	double nx_float, ny_float, nz_float, dz;
	double alpha_square_imag, app_res;
	dcomplex alpha,Q_value_all[2], Q_value[2];
	dcomplex cmplx_nu, imag_one, h_times_alpha, nu_div_alpha, alpha_div_nu, tanh_div_alpha, tanh_times_Qvalue;
	dcomplex  nominator, denominator;

//FILE *out;
//out = fopen("test.txt","w");

	ny2 = grid.ny +2*grid.nborder;
	nz2 = grid.nz +2*grid.nborder;
	nyz2 = ny2*nz2;

	if(data->ndata_mt == 0)
		printf("!!WARNING!! NO MT data exists but the calcultion of frechet derivatives is activated\n\n");

	cmplx_nu = Complex((M_KONST),0);
	imag_one = Complex(0,1);

	/******************************************************/
	/*1D-calculation (HORIZONTAL LAYERED MEDIA) in a 2D or 3D grid*/
	/******************************************************/

	/*Loop over all stations*/
	for(i=0;i<geo.nstat_mt;i++)
	{
			/*Positions of the MT stations:*/
			xs = (float)geo.x[(data->mts[i]-1)];
			ys = (float)geo.y[(data->mts[i]-1)];
			zs = (float)geo.z[(data->mts[i]-1)];

			/*Calculate the x,y,z positions of the stations in the grid */
				/*cell nr in x-direction*/
			nx_float = (xs - grid.org[0])/grid.h;
			nr_xcells = 0; /*Number of forward cells in x-direction that will be considered for the calculation*/

			for(j=0;j<2;j++)
				nx_int[j] = -9999;

			if(nx_float - (int)nx_float <= 0.5)
			{
				int_tmp = (int)floor(nx_float);
				if(int_tmp >= 0 && int_tmp < grid.nx)
				{
					nx_int[nr_xcells] = int_tmp;
					nr_xcells++;
				}
			}

			if(nx_float - (int)nx_float >= 0.5)
			{
				int_tmp = (int)ceil(nx_float);
				if(int_tmp >= 0 && int_tmp < grid.nx)
				{
					nx_int[nr_xcells] = int_tmp;
					nr_xcells++;
				}
			}

			if(nr_xcells == 0)
			{
				printf("The MT station %d is located outside of the grid\n",i);
				printf("x-position of station = %f m\n",xs);
				exit(0);
			}

				/*cell nr in y-direction*/
			ny_float = (ys - grid.org[1])/grid.h;
			nr_ycells = 0; /*Number of forward cells in y-direction that will be considered for the calculation*/

			for(j=0;j<2;j++)
				ny_int[j] = -9999;

			if(ny_float - (int)ny_float <= 0.5)
			{
				int_tmp = (int)floor(ny_float);
				if(int_tmp >= 0 && int_tmp < grid.ny)
				{
					ny_int[nr_ycells] = int_tmp;
					nr_ycells++;
				}
			}

			if(ny_float - (int)ny_float >= 0.5)
			{
				int_tmp = (int)ceil(ny_float);
				if(int_tmp >= 0 && int_tmp < grid.ny)
				{
					ny_int[nr_ycells] = int_tmp;
					nr_ycells++;
				}
			}

			if(nr_ycells == 0)
			{
				printf("The MT station %d is located outside of the grid\n",i);
				printf("y-position of station = %f m\n",ys);
				exit(0);
			}

				/*cell nr in z-direction*/
			nz_float = (zs - grid.org[2])/grid.h;
			
			if(nz_float - (int)nz_float <= 0.5)
			{
				nz_int = (int)floor(nz_float);
				dz = (0.5 - (nz_float - nz_int)) * grid.h; /*Thickness of the "layer" immediately below the MT stations*/
			}
			else
			{
				nz_int = (int)ceil(nz_float);
				dz = (0.5 + (nz_int - nz_float)) * grid.h; /*Thickness of the "layer" immediately below the MT stations*/
			}

			if(nz_int < 0 || nz_int >= grid.nz)
			{
				printf("The MT station %d is located outside of the grid\n",i);
				printf("z-position of station = %f m\n",zs);
				exit(0);
			}

			/********************************************/
			/*Calculate the MT response*/

			/*Relate the stations to the corresponding MT measurements*/
			for(m=0;m<data->ndata_mt;m++)
			{
				if(data->mts[i] == data->mno[m])
				{

					/*Determine the cells which will be perturbated*/
					nr_count = 0;

					for(j=0;j<3;j++)
						act[j] = (long *)memory(NULL,1,sizeof(long),"DerivativesMT");
						
					/*Loop over all used cells in xy direction*/
					for(j=0;j<nr_xcells;j++)
						for(k=0;k<nr_ycells;k++)
						{
							for(r=(grid.nz - 1);r >= nz_int;r--)
							{
								nr_count++;
								for(s=0;s<3;s++)
									act[s] = (long *)memory((char *)act[s],nr_count,sizeof(long),"DerivativesMT");

								act[0][nr_count -1] = nx_int[j];
								act[1][nr_count -1] = ny_int[k];
								act[2][nr_count -1] = r;

							}
						}

					/*Set parameters for the MT structure (for the forward cells)*/
					mt[m].ncell = nr_count;				/*Number of cells that are used for the MT measurement*/
					mt[m].nfreq = data->nfreq_mt[m];	/*Number of frequencies for the measurement*/
					mt[m].n = (long *)memory(NULL,1,sizeof(long),"DerivativesMT"); /*NOT used for the foward cells!!!*/
					
					if(mt[m].ncell != 0)
					{
						mt[m].ele = (long *)memory(NULL,mt[m].ncell,sizeof(long),"DerivativesMT");
					    mt[m].deriv = (double **)memory(NULL,mt[m].ncell,sizeof(double *),"DerivativesMT");
					}
					else
					{
						mt[m].ele = (long *)memory(NULL,1,sizeof(long),"DerivativesMT");
						mt[m].deriv = (double **)memory(NULL,mt[m].ncell,sizeof(double *),"DerivativesMT");
					}

					for(n=0;n<mt[m].ncell;n++)
					{
						if(mt[m].nfreq != 0)
							mt[m].deriv[n] = (double *)memory(NULL,2*mt[m].nfreq,sizeof(double),"DerivativesMT");
						else
							mt[m].deriv[n] = (double *)memory(NULL,1,sizeof(double),"DerivativesMT");

						mt[m].ele[n] = 0;
						
						for(s=0;s<(2*mt[m].nfreq);s++)
							mt[m].deriv[n][s] = 0.0;
					}


					/*Loop over all frequencies*/
					for(n=0;n<data->nfreq_mt[m];n++)
					{
						/*Use angular frequencies 2*PI*f*/
						ang_freq = 2*PI*(data->freq_mt[m][n]);

						/*Loop over all perturbations*/
						for(s=0;s<nr_count;s++)
						{
							/*Loop over the upper and lower perturbation*/
							for(t=0;t<2;t++)
							{
	
								pos_count = 0;

								Q_value_all[t].r = 0.0;
								Q_value_all[t].i = 0.0;
			
								/*Loop over all used cells in xy direction*/
								for(j=0;j<nr_xcells;j++)
									for(k=0;k<nr_ycells;k++)
									{
										/*calculate the response from the deepest layer (homogeneous half-space)*/
										if(act[0][s] == nx_int[j] && act[1][s] == ny_int[k] && act[2][s] == grid.nz-1) /*The layer is perturbated*/
										{
											if(t == 0) /*Pertubation dRho/d(-1/s)*/
												alpha_square_imag = (M_KONST)*ang_freq*(1/R1(nx_int[j],ny_int[k],grid.nz-1)); /*Imag.part of squared spatial wavenumber a*/
											else	/*Pertubation dRho/d(1/s)*/
												alpha_square_imag = (M_KONST)*ang_freq*(1/R2(nx_int[j],ny_int[k],grid.nz-1)); /*Imag.part of squared spatial wavenumber a*/
										}
										else
											alpha_square_imag = (M_KONST)*ang_freq*(1/res(nx_int[j],ny_int[k],grid.nz-1)); /*Imag.part of squared spatial wavenumber a*/
										alpha = Complex(0,alpha_square_imag); /*Make a complex number*/
										alpha = Csqrt(alpha);				  /*Calculate the complex spatial wavenumber a*/

										Q_value[t] = Cdiv(cmplx_nu,alpha);		/*Calculate the ratio: nu*E_y/i*w*B_x*/

										/**********************/
										/*calculate the recursively the response from the other layers*/
										/*(Backward)loop over the used cells in z-direction (ONLY cells benneath the station are considered)*/
										for(r=(grid.nz - 2);r >= nz_int;r--)
										{
											/*Recursive calculation of the Q-ratio*/
											if(act[0][s] == nx_int[j] && act[1][s] == ny_int[k] && act[2][s] == r)/*The layer is perturbated*/
											{
												if(t == 0) /*Pertubation dRho/d(-1/s)*/
													alpha_square_imag = (M_KONST)*ang_freq*(1/R1(nx_int[j],ny_int[k],r)); /*Imag.part of squared spatial wavenumber a*/
												else	/*Pertubation dRho/d(1/s)*/
													alpha_square_imag = (M_KONST)*ang_freq*(1/R2(nx_int[j],ny_int[k],r)); /*Imag.part of squared spatial wavenumber a*/
											}
											else
												alpha_square_imag = (M_KONST)*ang_freq*(1/res(nx_int[j],ny_int[k],r)); /*Imag.part of squared spatial wavenumber a*/
											alpha = Complex(0,alpha_square_imag); /*Make a complex number*/
											alpha = Csqrt(alpha);				  /*Calculate the complex spatial wavenumber a*/

											/*Determine the tanh*/
											if(r == nz_int)
												h_times_alpha = RCmul(dz,alpha);        /*Determine the product of slowness and the REMAINING thickness of the LAST CELL*/
											else
												h_times_alpha = RCmul(grid.h,alpha);	/*Determine the product of slowness and thickness of the layer/cells*/
											h_times_alpha = CBtanh(h_times_alpha);		/*Calculate the tangens hyperbolicus*/

											/*Calculate the nominator*/
											nu_div_alpha = Cdiv(cmplx_nu,alpha);				/*Determine the quotient of nu_0 and the spatial wavenumber a*/
											tanh_div_alpha = Cmul(h_times_alpha,nu_div_alpha);	/*Determine the product of tanh(...)*nu_0/a */
											nominator = Cadd(tanh_div_alpha,Q_value[t]);

											/*Calculate the denominator*/
											alpha_div_nu = Cdiv(alpha,cmplx_nu);				/*Determine the quotient of the spatial wavenumber a and nu_0*/
											tanh_times_Qvalue = Cmul(h_times_alpha,Q_value[t]);	/*Determine the product of tanh(...)*Q_value (of the former iteration)*/
											denominator = Cmul(tanh_times_Qvalue,alpha_div_nu);
											denominator.r = denominator.r + 1.0;

											/*Determine the new Q-ratio*/
											Q_value[t] = Cdiv(nominator,denominator);

										}

										/*If the station are located exactly at the border between to cells, the results are summed*/
										Q_value_all[t].r = Q_value_all[t].r + Q_value[t].r;
										Q_value_all[t].i = Q_value_all[t].i + Q_value[t].i;

										pos_count++;

									}

									Q_value_all[t].r = Q_value_all[t].r/pos_count;
									Q_value_all[t].i = Q_value_all[t].i/pos_count;

									/*Determine the impedance from the Q-ratio*/
									Q_value_all[t] = Cmul(Q_value_all[t], imag_one);
									Q_value_all[t].r = ((-1)*ang_freq*Q_value_all[t].r);
									Q_value_all[t].i = ((-1)*ang_freq*Q_value_all[t].i);

									/*Apparent resistivity*/
									app_res = (double)(((Q_value_all[t].r)*(Q_value_all[t].r)) + ((Q_value_all[t].i)*(Q_value_all[t].i)))*(1/(ang_freq*(M_KONST)));

								}
									
								/*Calculate the derivatives dZ/ds*/
								mt[m].deriv[s][2*n] = (double)(Q_value_all[1].r - Q_value_all[0].r)/deltaSlow;
								mt[m].deriv[s][2*n+1] = (double)(Q_value_all[1].i - Q_value_all[0].i)/deltaSlow;
								
								/*Positions of the forward cells*/
								mt[m].ele[s] = (nyz2)*(act[0][s] +  grid.nborder) + (nz2)*(act[1][s] + grid.nborder) + (act[2][s] + grid.nborder);

//fprintf(out,"%d %f %f %f %f %d\n",m,data->freq_mt[m][n],Q_value_all[0].r, Q_value_all[0].i, app_res, t);				
					
							}
					}

				for(j=0;j<3;j++)
					free(act[j]);

				}
			}

			if(i%10 == 0)
				printf("Derivatives are calculated for %d of %d stations\n", i+1, geo.nstat_mt);
	}

//fclose(out);

	return(1);
}


/*-------------------------------------------------------------*/
/*Calculating the slowness derivatives for the 2-D MT model*/
/*Method: Calculating analytic impression following Jegen, 1997*/
/*Parameter: *calc_2D_mt := 2D mt structure*/
/*     nr_of_slices_2D_mt:= Number of 2D slices*/
/*                  grid := grid structure*/
/*					data := data structure*/
/*				*R2 -R1 := Perturbation of the resistivities in ohmm*/
/*				deltaSlow := Perturbation of the slowness in s/m*/
/*					*mt := MT structure including the slowness derivatives for the forward cells*/
/*			    *flag   := Flag structure*/

/*Remark: the magnetic permeability is considered to correspond to the one of the free air*/

int Derivavtives2DMTAnaly(CALC_2D_MT *calc_2D_mt, long nr_of_slices_2D_mt, GRID_STRUCT grid, DATA_STRUCT *data, double *R2, double *R1, double deltaSlow, MT_STRUCT *mt, FLAG_STRUCT *flag)
{
	int *already_used_mt_structure; /*Specify, if values were already assigned to the corresponding mt (forward structure; output structure) yes == 1; no == 0*/
	long i,j,k,m,a,b,c;
	long ny2,nz2,nyz2;
	long nx,nz, index_y_cell;
	double *res_2D;
	double *perturb1, *perturb2;
	double *b_index, *tmp_3D_b_index;
	long nz_without_air; /*Nr of cells in z-direction without air*/
	double **d_Er1, **d_Ei1, **d_Hr1, **d_Hi1; /*Perturbances of the E- and the H-fields*/
	double **d_Er2, **d_Ei2, **d_Hr2, **d_Hi2; /*Perturbances of the E- and the H-fields*/
	double *d_Zr, *d_Zi;
	double dZr, dZi; /*Real and imaginary part of the derivative of the impedance*/
	double dEr, dEi ,dHr, dHi, Er, Ei, Hr, Hi; 
	
	char fname[40];

//FILE *out;

	/*Set the index that specify, if parameters has been already assigned to a mt (forward) structure*/
	if(data->ndata_mt != 0)
		already_used_mt_structure = (int *)memory(NULL,data->ndata_mt,sizeof(int),"Derivavtives2DMTAnaly");
	else
		already_used_mt_structure = (int *)memory(NULL,1,sizeof(int),"Derivavtives2DMTAnaly");

	for(i=0;i<data->ndata_mt;i++)
		already_used_mt_structure[i] = 0;

	ny2 = 2*grid.nborder + grid.ny;
	nz2 = 2*grid.nborder + grid.nz;
	nyz2 = ny2*nz2; 

	tmp_3D_b_index = (double *)memory(NULL,(grid.nx*grid.ny*grid.nz),sizeof(double),"Derivavtives2DMTAnaly");

	/*Adjusting the index that specify border cells*/
	for(a=0;a<grid.nx;a++)
		for(b=0;b<grid.ny;b++)
			for(c=0;c<grid.nz;c++)
			{
				tmp_3D_b_index[a*(grid.ny*grid.nz) + b*grid.nz + c] = (double)grid.border_index[(a+grid.nborder)*nyz2 + (b+grid.nborder)*nz2 +(c+grid.nborder)];
			}


	/*Specify the direction of the 2d slices*/
	if(flag->direc_2D_mt != 2) 
		nx = grid.nx;
	else
		nx = grid.ny;

	nz = grid.nz;

	index_y_cell = 0; /*First cell in the of the slice in the strike direction*/


	printf("\n---------------\n");
	printf("Start calculating the derivatives for the 2-D MT measurements\n");
	printf("---------------\n");

	/*Loop over all slices*/
	for(i=0;i<nr_of_slices_2D_mt;i++)
	{
		/*Loop over all frequencies*/
		for(j=0;j<calc_2D_mt[i].nfreq;j++)
		{
			nz_without_air = calc_2D_mt[i].nz[j]- calc_2D_mt[i].nz_atm[j] - calc_2D_mt[i].nz_ionos[j]; /*Nr of cells in z-direction without air*/

			perturb1 = (double *)memory(NULL,(nx*nz), sizeof(double),"Derivavtives2DMTAnaly");
			perturb2 = (double *)memory(NULL,(nx*nz), sizeof(double),"Derivavtives2DMTAnaly");
			/************************************************/
			res_2D = (double *)memory(NULL,(nx*nz), sizeof(double),"Derivavtives2DMTAnaly");
			b_index = (double *)memory(NULL,(nx*nz), sizeof(double),"Derivavtives2DMTAnaly");


			/*Determine the average resistivity perturbation along the strike direction of the slices*/
			CalcAveValues2DSlice(grid, index_y_cell, flag->direc_2D_mt, calc_2D_mt[i].nr_cells_2d_mt, nx, nz, perturb1, R1);
			CalcAveValues2DSlice(grid, index_y_cell, flag->direc_2D_mt, calc_2D_mt[i].nr_cells_2d_mt, nx, nz, perturb2, R2);
			CalcAveValues2DSlice(grid, index_y_cell, flag->direc_2D_mt, calc_2D_mt[i].nr_cells_2d_mt, nx, nz, res_2D, grid.res);
			CalcAveValues2DSlice(grid, index_y_cell, flag->direc_2D_mt, calc_2D_mt[i].nr_cells_2d_mt, nx, nz, b_index, tmp_3D_b_index);

			/*Determine the pertubations of the conductivity by means of the perturbations of the resistivity*/
            /*d_sigma = - d_r/(r)^2*/
			for(k=0;k<(nx*nz);k++)
			{
				perturb1[k] = (-1)*(perturb1[k] - res_2D[k])/(res_2D[k]*res_2D[k]);
				perturb2[k] = (-1)*(perturb2[k] - res_2D[k])/(res_2D[k]*res_2D[k]);
			}

			free(res_2D);

			/************************************************/
			/* Assign the perturbation to the refined grid  */
			/************************************************/
			perturb1 = (double *)memory((char *)perturb1,calc_2D_mt[i].nx[j]*(calc_2D_mt[i].nz[j]- calc_2D_mt[i].nz_atm[j] - calc_2D_mt[i].nz_ionos[j]), sizeof(double),"AssignRefinedGridMT");
			perturb2 = (double *)memory((char *)perturb2,calc_2D_mt[i].nx[j]*(calc_2D_mt[i].nz[j]- calc_2D_mt[i].nz_atm[j] - calc_2D_mt[i].nz_ionos[j]), sizeof(double),"AssignRefinedGridMT");

			AssignRefinedGridMT(calc_2D_mt[i], perturb1, grid.h, nx, nz, j);
			AssignRefinedGridMT(calc_2D_mt[i], perturb2, grid.h, nx, nz, j);

			/*Allocate memory for the perturbations of the field components*/
			if(calc_2D_mt[i].freq_ndata[j] != 0)
			{
				d_Er1 = (double **)memory(NULL,calc_2D_mt[i].freq_ndata[j],sizeof(double *),"Derivavtives2DMTAnaly");
				d_Ei1 = (double **)memory(NULL,calc_2D_mt[i].freq_ndata[j],sizeof(double *),"Derivavtives2DMTAnaly");
				d_Hr1 = (double **)memory(NULL,calc_2D_mt[i].freq_ndata[j],sizeof(double *),"Derivavtives2DMTAnaly");
				d_Hi1 = (double **)memory(NULL,calc_2D_mt[i].freq_ndata[j],sizeof(double *),"Derivavtives2DMTAnaly");

				d_Er2 = (double **)memory(NULL,calc_2D_mt[i].freq_ndata[j],sizeof(double *),"Derivavtives2DMTAnaly");
				d_Ei2 = (double **)memory(NULL,calc_2D_mt[i].freq_ndata[j],sizeof(double *),"Derivavtives2DMTAnaly");
				d_Hr2 = (double **)memory(NULL,calc_2D_mt[i].freq_ndata[j],sizeof(double *),"Derivavtives2DMTAnaly");
				d_Hi2 = (double **)memory(NULL,calc_2D_mt[i].freq_ndata[j],sizeof(double *),"Derivavtives2DMTAnaly");
			
				for(k=0;k<calc_2D_mt[i].freq_ndata[j];k++)
				{
					d_Er1[k] = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");
					d_Ei1[k] = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");
					d_Hr1[k] = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");
					d_Hi1[k] = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");

					d_Er2[k] = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");
					d_Ei2[k] = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");
					d_Hr2[k] = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");
					d_Hi2[k] = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");

				}

			}
			else
			{
				d_Er1 = (double **)memory(NULL,1,sizeof(double *),"Derivavtives2DMTAnaly");
				d_Ei1 = (double **)memory(NULL,1,sizeof(double *),"Derivavtives2DMTAnaly");
				d_Hr1 = (double **)memory(NULL,1,sizeof(double *),"Derivavtives2DMTAnaly");
				d_Hi1 = (double **)memory(NULL,1,sizeof(double *),"Derivavtives2DMTAnaly");

				d_Er2 = (double **)memory(NULL,1,sizeof(double *),"Derivavtives2DMTAnaly");
				d_Ei2 = (double **)memory(NULL,1,sizeof(double *),"Derivavtives2DMTAnaly");
				d_Hr2 = (double **)memory(NULL,1,sizeof(double *),"Derivavtives2DMTAnaly");
				d_Hi2 = (double **)memory(NULL,1,sizeof(double *),"Derivavtives2DMTAnaly");

			}

			d_Zr = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");
			d_Zi = (double *)memory(NULL,(calc_2D_mt[i].nx[j] * nz_without_air),sizeof(double),"Derivavtives2DMTAnaly");
	
			/***********************************************************************************/
			/*Calculation of the derivatives for the TE mode following the analytic impressions*/
			/*from Marion Jegen (described in her thesis,1997)*/
			if(flag->kind_of_data_mt != 2)
			{
				
				for(m=0;m<(calc_2D_mt[i].nx[j] * nz_without_air);m++)
				{
					for(k=0;k<calc_2D_mt[i].freq_ndata[j];k++)
					{
						d_Er1[k][m]=0.0;
						d_Ei1[k][m]=0.0;
						d_Hr1[k][m]=0.0;
						d_Hi1[k][m]=0.0;

						d_Er2[k][m]=0.0;
						d_Ei2[k][m]=0.0;
						d_Hr2[k][m]=0.0;
						d_Hi2[k][m]=0.0;
					}

					d_Zr[m]=0.0;
					d_Zi[m]=0.0;
				}

				/*Calculate the perturbations of the fields*/
				Deriv_MT_analytical(calc_2D_mt[i],perturb1,i,j,1, d_Er1, d_Ei1, d_Hr1, d_Hi1, grid.h);
				Deriv_MT_analytical(calc_2D_mt[i],perturb2,i,j,1, d_Er2, d_Ei2, d_Hr2, d_Hi2, grid.h);

				/*****************************************************/
					/*Delete the temporary files including the field component Ey*/
				sprintf(fname,"tmp_EyR_layer%d_freq%f.dat",i,calc_2D_mt[i].freq[j]);
				remove(fname);
				sprintf(fname,"tmp_EyI_layer%d_freq%f.dat",i,calc_2D_mt[i].freq[j]);
				remove(fname);
				/*****************************************************/

				for(k=0;k<calc_2D_mt[i].freq_ndata[j];k++)
				{
					for(m=0;m<(calc_2D_mt[i].nx[j] * nz_without_air);m++)
					{

						/*Determine the derivatives for the impedances*/
						/*Calculate dE/ds and dH/ds*/
						dEr = (d_Er2[k][m]- d_Er1[k][m])/deltaSlow;
						dEi = (d_Ei2[k][m]- d_Ei1[k][m])/deltaSlow;
						dHr = (d_Hr2[k][m]- d_Hr1[k][m])/deltaSlow;
						dHi = (d_Hi2[k][m]- d_Hi1[k][m])/deltaSlow;

						Er = calc_2D_mt[i].Ey_r[j][k];
						Ei = calc_2D_mt[i].Ey_i[j][k];
						Hr = calc_2D_mt[i].Hx_r[j][k];
						Hi = calc_2D_mt[i].Hx_i[j][k];

						/*calculate dZ/ds= ((dE/ds)*H - (dH/ds)*E)/H^2 */
						CDerivQuo(Er,Ei,Hr,Hi,dEr,dEi,dHr,dHi,&dZr,&dZi);

						d_Zr[m] = dZr;
						d_Zi[m] = dZi;

//d_Zr[m] = dEr;
//d_Zi[m] = dEi;


					}

					/*Transfer the sensitivities to the forward grid*/
					SumSensToForwardCells(calc_2D_mt[i], d_Zr, grid.h, nx, nz, j);
					SumSensToForwardCells(calc_2D_mt[i], d_Zi, grid.h, nx, nz, j);

					/*Build the mt-forward structure*/
					makeMTstructure2D(calc_2D_mt[i], d_Zr, d_Zi, b_index, grid, &mt[calc_2D_mt[i].freq_data[j][k]],j ,k, 1 , flag->direc_2D_mt, data, &already_used_mt_structure[calc_2D_mt[i].freq_data[j][k]], index_y_cell);

/********************/
/*out = fopen("test_Ey_16b.txt","w");
	for(m=0;m<(nx*nz);m++)
		fprintf(out,"%20.18f %20.18f\n",d_Zr[m],d_Zi[m]);
fclose(out);
/********************/

				}

				/*****************************************************/

				printf("Deriv. of the TE-Mode are calculated for the layer %d and frequency %fHz\n",i,calc_2D_mt[i].freq[j]);

			}
			/***********************************************************************************/
			/*Calculation of the derivatives for the TM mode following the analytic impressions*/
			/*from Marion Jegen (described in her thesis,1997)*/
			if(flag->kind_of_data_mt != 1)
			{

				for(m=0;m<(calc_2D_mt[i].nx[j] * nz_without_air);m++)
				{
					for(k=0;k<calc_2D_mt[i].freq_ndata[j];k++)
					{

						d_Er1[k][m]=0.0;
						d_Ei1[k][m]=0.0;
						d_Hr1[k][m]=0.0;
						d_Hi1[k][m]=0.0;

						d_Er2[k][m]=0.0;
						d_Ei2[k][m]=0.0;
						d_Hr2[k][m]=0.0;
						d_Hi2[k][m]=0.0;

						
					}
					d_Zr[m]=0.0;
					d_Zi[m]=0.0;
				}

				/*Calculate the perturbations of the fields*/
				Deriv_MT_analytical(calc_2D_mt[i],perturb1,i,j,2, d_Er1, d_Ei1, d_Hr1, d_Hi1, grid.h);
				Deriv_MT_analytical(calc_2D_mt[i],perturb2,i,j,2, d_Er2, d_Ei2, d_Hr2, d_Hi2, grid.h);


				/*****************************************************/
				/*Delete the temporary files including the field components Ex and Ez*/
				sprintf(fname,"tmp_ExR_layer%d_freq%f.dat",i,calc_2D_mt[i].freq[j]);
				remove(fname);
				sprintf(fname,"tmp_ExI_layer%d_freq%f.dat",i,calc_2D_mt[i].freq[j]);
				remove(fname);
				sprintf(fname,"tmp_EzR_layer%d_freq%f.dat",i,calc_2D_mt[i].freq[j]);
				remove(fname);
				sprintf(fname,"tmp_EzI_layer%d_freq%f.dat",i,calc_2D_mt[i].freq[j]);
				remove(fname);
				/*****************************************************/

				/*Determine the derivatives for the impedances*/
				for(k=0;k<calc_2D_mt[i].freq_ndata[j];k++)
				{
					for(m=0;m<(calc_2D_mt[i].nx[j] * nz_without_air);m++)
					{
						/*Calculate dE/ds and dH/ds*/
						dEr = (d_Er2[k][m]- d_Er1[k][m])/deltaSlow;
						dEi = (d_Ei2[k][m]- d_Ei1[k][m])/deltaSlow;
						dHr = (d_Hr2[k][m]- d_Hr1[k][m])/deltaSlow;
						dHi = (d_Hi2[k][m]- d_Hi1[k][m])/deltaSlow;

						Er = calc_2D_mt[i].Ex_r[j][k];
						Ei = calc_2D_mt[i].Ex_i[j][k];
						Hr = calc_2D_mt[i].Hy_r[j][k];
						Hi = calc_2D_mt[i].Hy_i[j][k];

						/*calculate dZ/ds= ((dE/ds)*H - (dH/ds)*E)/H^2 */
						CDerivQuo(Er,Ei,Hr,Hi,dEr,dEi,dHr,dHi,&dZr,&dZi);

						d_Zr[m] = dZr;
						d_Zi[m] = dZi;

//d_Zr[m] = dEr;
//d_Zi[m] = dEi;

					}

					/*Transfer the sensitivities to the forward grid*/
					SumSensToForwardCells(calc_2D_mt[i], d_Zr, grid.h, nx, nz, j);
					SumSensToForwardCells(calc_2D_mt[i], d_Zi, grid.h, nx, nz, j);
					
					/*Build the mt-forward structure*/
					makeMTstructure2D(calc_2D_mt[i], d_Zr, d_Zi, b_index, grid, &mt[calc_2D_mt[i].freq_data[j][k]],j ,k, 2, flag->direc_2D_mt ,data, &already_used_mt_structure[calc_2D_mt[i].freq_data[j][k]], index_y_cell);

/********************/
/*	out = fopen("test_Ex_17.txt","w");
	for(m=0;m<(nx*nz);m++)
		fprintf(out,"%20.18f %20.18f\n",d_Zr[m],d_Zi[m]);
	fclose(out);
/********************/

				}

				/*****************************************************/
				printf("Deriv. of the TM-Mode are calculated for the layer %d and frequency %fHz\n",i,calc_2D_mt[i].freq[j]);
			}

			for(k=0;k<calc_2D_mt[i].freq_ndata[j];k++)
			{
				free(d_Er1[k]);
				free(d_Ei1[k]);
				free(d_Hr1[k]);
				free(d_Hi1[k]);

				free(d_Er2[k]);
				free(d_Ei2[k]);
				free(d_Hr2[k]);
				free(d_Hi2[k]);
			}

			free(d_Er1);
			free(d_Ei1);
			free(d_Hr1);
			free(d_Hi1);

			free(d_Er2);
			free(d_Ei2);
			free(d_Hr2);
			free(d_Hi2);

			free(d_Zr);
			free(d_Zi);

			/*Delete the temporary files including the resistivity values (and indeces)*/
			sprintf(fname,"tmp_res_layer%d_freq%f.dat",i,calc_2D_mt[i].freq[j]);
			remove(fname);
			sprintf(fname,"tmp_indeces_layer%d_freq%f.dat",i,calc_2D_mt[i].freq[j]);
			remove(fname);
				
			free(perturb1);
			free(perturb2);

			free(b_index);

		}

		/*Set first cell in y-direction*/
		index_y_cell = calc_2D_mt[i].nr_cells_2d_mt + index_y_cell;

	}

	free(already_used_mt_structure);
	free(tmp_3D_b_index);

	return(0);
}

/*-------------------------------------------------------------*/
/*Calculate analytically the perturbations of the fields at the stations (following the PhD thesis from Marion Jegen, 1997, pages 60-70)*/
/*Parameter: calc_2D_mt := 2D mt structure*/
/*             *perturb := perturbations of the conductivities*/
/*              i_layer := Index of the slice*/
/*              j_freq  := Index of the frequency*/
/*           kind_of_mt := kind_of_mt (1 == TE mode; 2 == TM mode)*/
/*               grid_h := spatial interval of the original grid (without refinenent)*/
/*	  **d_Er,**d_Ei,**d_Hr, **d_Hi := Pertubation of the E- and the H-field (output of the routine); Remark: The first index specify the measurement following the sorting in "calc_2d_mt.freq_ndata"*/
/*                                     The second index specify the cell number of the grid*/

/*REMARK: Only for cells that are located "BELOW" the considered station, the derivatives are calculated; Otherwise, they remain =0*/
/*		  (This inaccuracy for strong topographies can be surely removed one time)*/

int Deriv_MT_analytical(CALC_2D_MT calc_2D_mt, double *perturb, long i_layer,long j_freq, int kind_of_mt, double **d_Er, double **d_Ei, double **d_Hr,double **d_Hi, double grid_h)
{

	int flag_ray_based; /*Specify if the calculations of the resistivities are ray-based (yes = 1; no= 0)*/
	int flag_reflect_coeff; /*Specify if the reflection coefficient will be considered in the calculation (yes = 1; no = 0)*/
	int **index_water_layer; /*Specify if there is a "water" layer between the air and the stations*/
	int potenz;
	long i,j,k,m,n,q,s, index_i;
	long cell_index, mult_cell; /*Number of rays that are used to calculate the average resistivity*/
	long nr_of_samples, nr_of_samples_Ex, i_sample, nhz;
	long nz_without_air, nz_air, nx_body, nr_of_cells_in_int;
	long index_stat, index_data; /*Indices of the stations and measurements*/
	long nr_of_zcoord; /*Number of different z-coordinates*/
	long *nr_of_stat_per_zcoord; /*Number of stations per z-coordinate*/
	long **l_index_stat, **l_index_data; /*List of station indices for a number of z-coordinates*/
	long index_interval, *z_int_stat_index; /*Index of the layers where the stations are placed*/
	long x_int_stat_index;
	double **x_index_stat; /*cell indeces in x-direction for the station positions*/
	double index_1,index_2;
	double nr_dsamples;
	double delta_x;
	double length_x, x_size_model;
	double *z_int_stat; /*z-position of the station*/
	double *Ey_r,*Ey_i; /*Components of the E-field*/
	double *first_cell_dz; 
	double res_station; /*resistivity exactly below the station*/
	double *Ex_r,*Ex_i,*Ez_r,*Ez_i;
	double *res_slice; /*resistivities*/
	double a_current_R, a_current_I; /*current from the perturbation*/
	double a_magmoment_R1, a_magmoment_I1; /*current dipole moment (x-component) from the perturbation*/
	double a_magmoment_R2, a_magmoment_I2; /*current dipole moment (z-component) from the perturbation*/
	double ave_Ey_r, ave_Ey_i, ave_Ex_r, ave_Ex_i, ave_Ez_r, ave_Ez_i;
	double  dx; /*z-distance between station and considered cell*/
	long idx; /*Number of samples in the x-direction*/
	double tmp_z_pos, z_pos, x_pos; /*position of the cell centers*/
	double *average_res, *hz;
	double tmp_sigma;
	double *x_index;
	double *p, *p_Ex; /*real part of the wavenumber*/
	double ****fft_input_B1, ****fft_input_B2, ****fft_input_E, ****fft_input_Ex1, ****fft_input_Ex2;
	double sampling_int; /*sampling interval*/
	double sampling_int_Ex; /*sampling interval for the Ex component*/
	double *waterdepths, **res_water; /*Thickness and resistivity of the water column*/
	dcomplex ave_Ex, ave_Ez;
	dcomplex Hx_stat, Ey_stat, Hy_stat1, Ex_stat1, Hy_stat2, Ex_stat2; /*fields at the stations*/
	dcomplex a_current; /*current from the pertubation*/
	dcomplex a_magmoment_x, a_magmoment_z; /*dipole moment from the perturbation*/

	/*****************************************/
	/*****************************************/
	/*DEFINITIONS*/
	/*****************************************/
	/*****************************************/
	#define EPS 1e-5		/*Interval that is related to the same z-coordinate for the seismic stations*/
	#define RAY_BASED 0		/*Specify, if the average resistivity of each layer (==0) or the resistivity along "rays" (!=0)*/
							/*is used for the calculation of the analytical derivatives*/
	#define MULT_CELL 10	/*(If RAY_BASED != 0) Specify, how many "rays" will be calculated to determine the average resistivity*/		

	/*Settings for the FFT in the p_x domain*/
	#define DELTA_X	100		/*Interval in x-direction in m required for the fourier transform into the p_x domain*/ 
    #define DELTA_X_EX 100	/*Interval in x-direction in m for the Ex component required for the fourier transform into the p_x domain; finer sampling is required because the fourier function is NOT continiously*/
	#define MULTI_X  1.3	/*Length of the line that are sampled by the fouriertranform corresponds to MULTI_X*(length_of_model in x-direction)*/

	#define REFLECT_COEFF 1 /*Specify, if the complex wavenumbers in the analytical formula are derived from the layered impedances (!=0)*/
							/* or by using an average values of the resistivities (==0)*/

	#define RES_AIR	1E10	/*considered resistivity of air in OHMm; (used for the calculation of the thickness of the water column)*/
	/*****************************************/
	/*****************************************/

	flag_ray_based = RAY_BASED;
	flag_reflect_coeff = REFLECT_COEFF;

	nx_body = calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq] - calc_2D_mt.nx_left_border[j_freq]; /*Nr of cells in x-direction in the body*/
	nz_without_air = calc_2D_mt.nz[j_freq]- calc_2D_mt.nz_atm[j_freq] - calc_2D_mt.nz_ionos[j_freq]; /*Nr of cells in z-direction without air*/
	nz_air = calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq]; /*Nr of cells in z-direction in the air*/

	/************************/
	/*Number of rays used to calculate the average resistivities*/
	mult_cell = MULT_CELL;

	if(mult_cell < 1)
		mult_cell = 1;

	if(mult_cell > nx_body)
		mult_cell = nx_body;

	nr_of_cells_in_int = (long) ceil(((double)nx_body)/(double)mult_cell);
    
	if(nx_body%nr_of_cells_in_int == 0)
		mult_cell = nx_body/nr_of_cells_in_int;

	/************************/
	/*Calculate the starting points for the "rays" in x-direction*/
	if(flag_ray_based != 0)
	{
		x_index = (double *)memory(NULL,mult_cell,sizeof(double),"Deriv_MT_analytical");
		index_i = calc_2D_mt.nx_left_border[j_freq];

		for(i=0;i<mult_cell;i++)
		{
			index_1 = (double)index_i; 

			for(j=0;j<nr_of_cells_in_int && index_i < (calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq]);j++)
				index_i++;

			index_2 = (double)(index_i -1.0); 

			/*Starting position of the ray*/
			x_index[i] = (index_1 + index_2)*0.5;
		}
	}
	else
	{
		x_index = (double *)memory(NULL,1,sizeof(double),"Deriv_MT_analytical");
		x_index[0] = 0.0;
	}

	/********************/
	res_slice = (double *)memory(NULL,calc_2D_mt.nx[j_freq]*nz_without_air,sizeof(double),"Deriv_MT_analytical");

	/*Read in the resistivity values from a temporary file*/
	ReadTResistivityIn(i_layer, calc_2D_mt.freq[j_freq], calc_2D_mt.nx[j_freq]*nz_without_air, res_slice, 0);

	/*****************************************************************************************************/
	/*****************************************************************************************************/
	/*Re-organize the indices for the stations and measurements by z-coordinates of the station location*/
	/*This can shorting strongly the runing-time of the code*/

	nr_of_stat_per_zcoord = (long *)memory(NULL,1,sizeof(long),"Deriv_MT_analytical");
	z_int_stat = (double *)memory(NULL,1,sizeof(double),"Deriv_MT_analytical");
	z_int_stat_index = (long *)memory(NULL,1,sizeof(long),"Deriv_MT_analytical");
	first_cell_dz = (double *)memory(NULL,1,sizeof(double),"Deriv_MT_analytical");
	l_index_stat = (long **)memory(NULL,1,sizeof(long *),"Deriv_MT_analytical");
	l_index_data = (long **)memory(NULL,1,sizeof(long *),"Deriv_MT_analytical");

	nr_of_zcoord = 0;

	/*Loop over all measurements*/
	for(k=0;k<calc_2D_mt.freq_ndata[j_freq];k++)
	{
		for(m=0;m<calc_2D_mt.nstat_mt;m++)
		{
			if(calc_2D_mt.freq_stat[j_freq][k] == calc_2D_mt.mts[m])
			{
				index_stat = m;
				index_data = k;
				goto dubidu;
			}
		}

		printf("NO station is found for the data measurement:\n");
		exit(0);

		dubidu:;

		index_interval = -1;
		/*Look if the station was already were used one time in the depth interval*/
		for(i=0;i<nr_of_zcoord;i++)
		{
			if(z_int_stat[i] - EPS <= calc_2D_mt.z[index_stat] && z_int_stat[i] + EPS >= calc_2D_mt.z[index_stat])
			{
				index_interval = i;
				break;
			}
		}

		/*z-coordinate NOT used by a station before*/
		if(index_interval == -1)
		{
			nr_of_zcoord++;

			nr_of_stat_per_zcoord = (long *)memory((char *)nr_of_stat_per_zcoord,nr_of_zcoord,sizeof(long),"Deriv_MT_analytical");
			z_int_stat = (double *)memory((char *)z_int_stat,nr_of_zcoord,sizeof(double),"Deriv_MT_analytical");
			z_int_stat_index = (long *)memory((char *)z_int_stat_index,nr_of_zcoord,sizeof(long),"Deriv_MT_analytical");
			first_cell_dz = (double *)memory((char *)first_cell_dz,nr_of_zcoord,sizeof(double),"Deriv_MT_analytical");
			l_index_stat = (long **)memory((char *)l_index_stat,nr_of_zcoord,sizeof(long *),"Deriv_MT_analytical");
			l_index_data = (long **)memory((char *)l_index_data,nr_of_zcoord,sizeof(long *),"Deriv_MT_analytical");

			l_index_stat[nr_of_zcoord -1] = (long *)memory(NULL,1,sizeof(long),"Deriv_MT_analytical");
			l_index_data[nr_of_zcoord -1] = (long *)memory(NULL,1,sizeof(long),"Deriv_MT_analytical");

			nr_of_stat_per_zcoord[nr_of_zcoord -1] = 1;
			z_int_stat[nr_of_zcoord -1] = calc_2D_mt.z[index_stat];
			z_int_stat_index[nr_of_zcoord -1] = 0;
			l_index_stat[nr_of_zcoord -1][0] = index_stat;
			l_index_data[nr_of_zcoord -1][0] = index_data;
		}
		else
		{
			nr_of_stat_per_zcoord[index_interval]++;

			l_index_stat[index_interval] = (long *)memory((char *)l_index_stat[index_interval],nr_of_stat_per_zcoord[index_interval],sizeof(long),"Deriv_MT_analytical");
			l_index_data[index_interval] = (long *)memory((char *)l_index_data[index_interval],nr_of_stat_per_zcoord[index_interval],sizeof(long),"Deriv_MT_analytical");

			l_index_stat[index_interval][nr_of_stat_per_zcoord[index_interval]-1] = index_stat;
			l_index_data[index_interval][nr_of_stat_per_zcoord[index_interval]-1] = index_data;
		}
	}

	/**********************************************************************/
	/*Determine the z-indices of coordinates where the stations are placed*/

	/*Loop over all depth intervals where stations are placed*/
	for(i=0;i<nr_of_zcoord;i++)
	{
		tmp_z_pos = calc_2D_mt.org[1][j_freq];
		z_pos = tmp_z_pos + calc_2D_mt.hz[j_freq][nz_air];

		/*Loop over all cells in z-direction*/
		j=0;

		while(z_pos < z_int_stat[i] && (j+1) < nz_without_air)
		{
			tmp_z_pos = z_pos;
			z_pos = z_pos + calc_2D_mt.hz[j_freq][nz_air + j + 1];
			j++;
		}
		
		first_cell_dz[i] = (z_int_stat[i] -  tmp_z_pos)/calc_2D_mt.hz[j_freq][nz_air + j];
		z_int_stat_index[i] = j;
	}

	/***********************************************************************/
	/*Determine the indeces of the cells in x-direction, where the stations are placed*/

	if(flag_ray_based != 0)
	{

		if(nr_of_zcoord != 0)
			x_index_stat = (double **)memory(NULL,nr_of_zcoord,sizeof(double *),"Deriv_MT_analytical");
		else
			x_index_stat = (double **)memory(NULL,1,sizeof(double *),"Deriv_MT_analytical");

		/*Loop over all layers, where stations are placed*/
		for(i=0;i<nr_of_zcoord;i++)
		{

			if(nr_of_stat_per_zcoord[i] != 0)
				x_index_stat[i] = (double *)memory(NULL,nr_of_stat_per_zcoord[i],sizeof(double),"Deriv_MT_analytical");
			else
				x_index_stat[i] = (double *)memory(NULL,1,sizeof(double),"Deriv_MT_analytical");

			/*Loop over all stations in a layer*/
			for(n=0;n<nr_of_stat_per_zcoord[i];n++)
			{
				x_pos = calc_2D_mt.org[0][j_freq];

				/*Loop over all cells in x-direction*/
				for(k=0;k<calc_2D_mt.nx[j_freq];k++)
				{
					x_pos = x_pos + calc_2D_mt.hx[j_freq][k];

					/*x-position is found*/
					if(x_pos >= calc_2D_mt.x[l_index_stat[i][n]])
					{
						delta_x = (x_pos - calc_2D_mt.x[l_index_stat[i][n]])/calc_2D_mt.hx[j_freq][k];

						x_index_stat[i][n] = k + delta_x;
						break;
					}
				}
			}
		}
	}
	else
	{
		x_index_stat = (double **)memory(NULL,1,sizeof(double *),"Deriv_MT_analytical");
	}
	/**********************************************************************/
	/**********************************************************************/
	/*Determine the sampling interval (using the horizontal interval between the cells)*/
	sampling_int = DELTA_X;
	/*... the sampling interval for Ex is different from the other components*/
	sampling_int_Ex = DELTA_X_EX; 

	/*The sampling interval should be never smaller than the grid sizes*/
	if(sampling_int > grid_h)
		sampling_int = grid_h;

	if(sampling_int_Ex > sampling_int)
		sampling_int_Ex = sampling_int;

	/*Determine the length of the sampled line (by means of the size of the modell and a scaling factor)*/
	x_size_model = 0.0;

	for(i=calc_2D_mt.nx_left_border[j_freq];i<(calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq]);i++)
		x_size_model = x_size_model + calc_2D_mt.hx[j_freq][i];

	length_x = MULTI_X * x_size_model;

	/*Determine the numbers of samples used for the spatial fourier transformation*/

		/******************/
		/*for the Ey, Bx and By component*/
	nr_dsamples = length_x/sampling_int;

	potenz = 0;
	while(pow(2,potenz) < nr_dsamples)
		potenz++;

	nr_of_samples = (long) pow(2,potenz);
	
	if(potenz == 0)
	{
		printf("NO samples exist for the fourier transform for determing the derivatives of the 2D-MT\n!!!");
		exit(0);
	}

	sampling_int = length_x/(double)nr_of_samples;

		/******************/
		/*for the Ex component*/
	nr_dsamples = length_x/sampling_int_Ex;

	potenz = 0;

	while(pow(2,potenz) < nr_dsamples)
		potenz++;

	nr_of_samples_Ex = (long) pow(2,potenz);

	if(potenz == 0)
	{
		printf("NO samples exist for the fourier transform for determing the derivatives of the 2D-MT\n!!!");
		exit(0);
	}

	sampling_int_Ex = length_x/(double)nr_of_samples_Ex;

	/**********************************************************************/
	/*Determine the real part of the wavenumber*/
	/*Allocate memory for the real part of the wavenumber*/
	p = (double *)memory(NULL,nr_of_samples+1,sizeof(double),"Deriv_MT_analytical");
	p_Ex = (double *)memory(NULL,nr_of_samples_Ex+1,sizeof(double),"Deriv_MT_analytical");

	if(nr_of_zcoord == 0)
	{
		fft_input_B1 = (double ****)memory(NULL,1, sizeof(double ***),"Deriv_MT_analytical");
		fft_input_B2 = (double ****)memory(NULL,1, sizeof(double ***),"Deriv_MT_analytical");
		fft_input_E = (double ****)memory(NULL,1, sizeof(double ***),"Deriv_MT_analytical");
		fft_input_Ex1 = (double ****)memory(NULL,1, sizeof(double ***),"Deriv_MT_analytical");
		fft_input_Ex2 = (double ****)memory(NULL,1, sizeof(double ***),"Deriv_MT_analytical");
	}
	else
	{
		
		fft_input_B1 = (double ****)memory(NULL,nr_of_zcoord, sizeof(double ***),"Deriv_MT_analytical");
		fft_input_B2 = (double ****)memory(NULL,nr_of_zcoord, sizeof(double ***),"Deriv_MT_analytical");
		fft_input_E = (double ****)memory(NULL,nr_of_zcoord, sizeof(double ***),"Deriv_MT_analytical");
		fft_input_Ex1 = (double ****)memory(NULL,nr_of_zcoord, sizeof(double ***),"Deriv_MT_analytical");
		fft_input_Ex2 = (double ****)memory(NULL,nr_of_zcoord, sizeof(double ***),"Deriv_MT_analytical");

		/*If the conductivity are calculated by averaging along the layers*/
		if(flag_ray_based == 0)
		{
			for(i=0;i<nr_of_zcoord;i++)
			{
				fft_input_B1[i] = (double ***)memory(NULL,1, sizeof(double **),"Deriv_MT_analytical");
				fft_input_B2[i] = (double ***)memory(NULL,1, sizeof(double **),"Deriv_MT_analytical");
				fft_input_E[i] = (double ***)memory(NULL,1, sizeof(double **),"Deriv_MT_analytical");
				fft_input_Ex1[i] = (double ***)memory(NULL,1, sizeof(double **),"Deriv_MT_analytical");
				fft_input_Ex2[i] = (double ***)memory(NULL,1, sizeof(double **),"Deriv_MT_analytical");

				fft_input_B1[i][0] = (double **)memory(NULL,1, sizeof(double *),"Deriv_MT_analytical");
				fft_input_B2[i][0] = (double **)memory(NULL,1, sizeof(double *),"Deriv_MT_analytical");
				fft_input_E[i][0] = (double **)memory(NULL,1, sizeof(double *),"Deriv_MT_analytical");
				fft_input_Ex1[i][0] = (double **)memory(NULL,1, sizeof(double *),"Deriv_MT_analytical");
				fft_input_Ex2[i][0] = (double **)memory(NULL,1, sizeof(double *),"Deriv_MT_analytical");

				fft_input_B1[i][0][0] = (double *)memory(NULL,(4*nr_of_samples), sizeof(double),"Deriv_MT_analytical");
				fft_input_B2[i][0][0] = (double *)memory(NULL,(4*nr_of_samples), sizeof(double),"Deriv_MT_analytical");
				fft_input_E[i][0][0] = (double *)memory(NULL,(4*nr_of_samples), sizeof(double),"Deriv_MT_analytical");
				fft_input_Ex1[i][0][0] = (double *)memory(NULL,(4*nr_of_samples_Ex), sizeof(double),"Deriv_MT_analytical");
				fft_input_Ex2[i][0][0] = (double *)memory(NULL,(4*nr_of_samples_Ex), sizeof(double),"Deriv_MT_analytical");
			}
		}
		/*If the conductivity is calculated along rays between the station and the perturbated cell*/
		else
		{
			for(i=0;i<nr_of_zcoord;i++)
			{

				fft_input_B1[i] = (double ***)memory(NULL,nr_of_stat_per_zcoord[i], sizeof(double **),"Deriv_MT_analytical");
				fft_input_B2[i] = (double ***)memory(NULL,nr_of_stat_per_zcoord[i], sizeof(double **),"Deriv_MT_analytical");
				fft_input_E[i] = (double ***)memory(NULL,nr_of_stat_per_zcoord[i], sizeof(double **),"Deriv_MT_analytical");
				fft_input_Ex1[i] = (double ***)memory(NULL,nr_of_stat_per_zcoord[i], sizeof(double **),"Deriv_MT_analytical");
				fft_input_Ex2[i] = (double ***)memory(NULL,nr_of_stat_per_zcoord[i], sizeof(double **),"Deriv_MT_analytical");

				for(j=0;j<nr_of_stat_per_zcoord[i];j++)
				{
					fft_input_B1[i][j] = (double **)memory(NULL,mult_cell, sizeof(double *),"Deriv_MT_analytical");
					fft_input_B2[i][j] = (double **)memory(NULL,mult_cell, sizeof(double *),"Deriv_MT_analytical");
					fft_input_E[i][j] = (double **)memory(NULL,mult_cell, sizeof(double *),"Deriv_MT_analytical");
					fft_input_Ex1[i][j] = (double **)memory(NULL,mult_cell, sizeof(double *),"Deriv_MT_analytical");
					fft_input_Ex2[i][j] = (double **)memory(NULL,mult_cell, sizeof(double *),"Deriv_MT_analytical");

					for(k=0;k<mult_cell;k++)
					{
						fft_input_B1[i][j][k] = (double *)memory(NULL,(4*nr_of_samples), sizeof(double),"Deriv_MT_analytical");
						fft_input_B2[i][j][k] = (double *)memory(NULL,(4*nr_of_samples), sizeof(double),"Deriv_MT_analytical");
						fft_input_E[i][j][k] = (double *)memory(NULL,(4*nr_of_samples), sizeof(double),"Deriv_MT_analytical");
						fft_input_Ex1[i][j][k] = (double *)memory(NULL,(4*nr_of_samples_Ex), sizeof(double),"Deriv_MT_analytical");
						fft_input_Ex2[i][j][k] = (double *)memory(NULL,(4*nr_of_samples_Ex), sizeof(double),"Deriv_MT_analytical");
					}
				}
			}
		}
	}


	/******************************/
	/*Determine p for the Ex, Hy and Hx components*/
	/*Loop over all samples*/
	for(i_sample = 0;i_sample<nr_of_samples+1; i_sample++)
	{
		/*Calculate p = sample_nr/(dx*2*Nr_of_samples)*/
		/*Attention: Because we have to consider the left and the right side of the symmetric function there is a factor 2*/
		p[i_sample] = i_sample/(2*nr_of_samples*sampling_int);
	}

	/*Determine p for the Ey component*/
	/*Loop over all samples*/
	for(i_sample = 0;i_sample<nr_of_samples_Ex+1; i_sample++)
	{
		/*Calculate p = sample_nr/(dx*2*Nr_of_samples)*/
		/*Attention: Because we have to consider the left and the right side of the symmetric function there is a factor 2*/
		p_Ex[i_sample] = i_sample/(2*nr_of_samples_Ex*sampling_int_Ex);
	}

	/******************************/
	/*Calculate the thickness and the resistivity of the water layer*/
	/*ATTENTION !!! The routine will only work reliably, if there is max. one horizontal boundary between the air and the water layer*/
	/*Internal layering within the watercolumn will be ignored*/

	if(nr_of_zcoord != 0)
	{
		waterdepths = (double *)memory(NULL,nr_of_zcoord,sizeof(double),"Deriv_MT_analytical");
		res_water = (double **)memory(NULL,nr_of_zcoord,sizeof(double *),"Deriv_MT_analytical");
		index_water_layer = (int **)memory(NULL,nr_of_zcoord,sizeof(int *),"Deriv_MT_analytical");

		for(i=0;i<nr_of_zcoord;i++)
		{
			if(flag_ray_based != 0)
			{
				res_water[i] = (double *)memory(NULL,nr_of_stat_per_zcoord[i],sizeof(double),"Deriv_MT_analytical");
				index_water_layer[i] = (int *)memory(NULL,nr_of_stat_per_zcoord[i],sizeof(int),"Deriv_MT_analytical");
			}
			else
			{
				res_water[i] = (double *)memory(NULL,1,sizeof(double),"Deriv_MT_analytical");
				index_water_layer[i] = (int *)memory(NULL,1,sizeof(int),"Deriv_MT_analytical");
			}
		}
	}
	else
	{
		waterdepths = (double *)memory(NULL,1,sizeof(double),"Deriv_MT_analytical");
		res_water = (double **)memory(NULL,1,sizeof(double *),"Deriv_MT_analytical");
		index_water_layer = (int **)memory(NULL,1,sizeof(int *),"Deriv_MT_analytical");
	}

	/*Loop over all depth intervals where stations are placed*/
	for(i=0;i<nr_of_zcoord;i++)
	{
		waterdepths[i] = 0.0;


		/*resistivities exactly above the stations are used*/
		if(flag_ray_based != 0)
		{
			for(j=0;j<nr_of_stat_per_zcoord[i];j++)
			{
				res_water[i][j] = 0.0;
				index_water_layer[i][j] = 0;

				x_int_stat_index = (long) x_index_stat[i][j];

				if(x_index_stat[i][j] - x_int_stat_index > 0.5)
					x_int_stat_index = x_int_stat_index + 1;

				/*Calculate the resistivity and the thickness of the water column*/
				DetWaterColumn(&waterdepths[i], &res_water[i][j], RES_AIR, z_int_stat_index[i], x_int_stat_index ,res_slice, calc_2D_mt, first_cell_dz[i], nz_without_air, j_freq, &index_water_layer[i][j], flag_ray_based);
			}
		}
		/*resistivities are averaged along a  layer*/
		else
		{
			res_water[i][0] = 0.0;
			index_water_layer[i][0] = 0;

			/*Calculate the resistivity and the thickness of the water column*/
			DetWaterColumn(&waterdepths[i], &res_water[i][0], RES_AIR, z_int_stat_index[i], 0,res_slice, calc_2D_mt, first_cell_dz[i], nz_without_air, j_freq, &index_water_layer[i][0], flag_ray_based);
		}
	}

	/***************************************************************************************************************/
	/***************************************************************************************************************/

	/************************/
	/************************/
		  /*TE-mode*/
	/************************/
	/************************/

	if(kind_of_mt == 1)
	{
		Ey_r = (double *)memory(NULL,calc_2D_mt.nx[j_freq]*nz_without_air,sizeof(double),"Deriv_MT_analytical");
		Ey_i = (double *)memory(NULL,calc_2D_mt.nx[j_freq]*nz_without_air,sizeof(double),"Deriv_MT_analytical");

		/*Read in the Ey_component from a temporary file*/
		ReadTResistivityIn(i_layer, calc_2D_mt.freq[j_freq], calc_2D_mt.nx[j_freq]*nz_without_air, Ey_r, 3);
		ReadTResistivityIn(i_layer, calc_2D_mt.freq[j_freq], calc_2D_mt.nx[j_freq]*nz_without_air, Ey_i, 4);

		/*z-position of the cell centers*/
		z_pos = calc_2D_mt.org[1][j_freq] + 0.5*calc_2D_mt.hz[j_freq][nz_air];

		/*Loop over all cells in z-direction*/
		for(j=0;j<nz_without_air;j++)
		{

			/*************************************************************************************/
			/*Loop over all depth intervals where stations are placed*/
			for(i=0;i<nr_of_zcoord;i++)
			{
				/*REMARK: Consider only cells BELOW the station*/
				if(z_int_stat_index[i]<=j)
				{
					/***********************************************/
					/*If the conductivity are calculated by averaging along the layers*/
					if(flag_ray_based == 0)
					{

						/***********************************************/ 
						/*Determine the average resistivity within the layers between the station and the perturbated cells*/
						average_res = (double *)memory(NULL,(j +1 -z_int_stat_index[i]),sizeof(double),"Deriv_MT_analytical");
						hz	= (double *)memory(NULL,(j +1 -z_int_stat_index[i]),sizeof(double),"Deriv_MT_analytical");

						nhz = j +1 -z_int_stat_index[i];

						for(q=0;q<(j +1 -z_int_stat_index[i]);q++)
						{
							average_res[q] = 0.0;
							hz[q] = 0.0;
						}

						/*The average resistivity for each layer will be calculated*/
						DetMeanResWithinLayer(calc_2D_mt,res_slice, average_res, hz, z_int_stat_index[i], first_cell_dz[i], j, nz_without_air, j_freq, &res_station, 0.0, 0.0, flag_ray_based);
						/*************************************************************************************/

						for(q=0;q<4*nr_of_samples;q++)
						{
							fft_input_B1[i][0][0][q] = 0.0;
							fft_input_E[i][0][0][q] = 0.0;
						}
		
						/*Calculate the Bx and Ey component Eq.4.30 and 4.31 in Marion thesis*/
						CalcBx_Ey(nr_of_samples, p, fft_input_B1[i][0][0], fft_input_E[i][0][0], average_res, hz, nhz, calc_2D_mt.freq[j_freq], waterdepths[i], res_water[i][0], index_water_layer[i][0], sampling_int, res_station, flag_reflect_coeff);

						/*Fourier Transformation*/
						four1(fft_input_B1[i][0][0]-1,2*nr_of_samples,-1);
						four1(fft_input_E[i][0][0]-1,2*nr_of_samples,-1);

						free(average_res);
						free(hz);
					}

					/***********************************************/
					/*If the conductivity is calculated along rays between the station and the perturbated cell*/
					else
					{

						/*Loop over all stations in a layer*/
						for(n=0;n<nr_of_stat_per_zcoord[i];n++)
						{

							/*Loop over all "rays"*/
							for(s=0;s<mult_cell;s++)
							{
							/***********************************************/ 
								/*Determine the average resistivity within the layers between the station and the perturbated cells*/
							
								average_res = (double *)memory(NULL,(j +1 -z_int_stat_index[i]),sizeof(double),"Deriv_MT_analytical");
								hz	= (double *)memory(NULL,(j +1 -z_int_stat_index[i]),sizeof(double),"Deriv_MT_analytical");

								nhz = j +1 -z_int_stat_index[i];
	
								for(q=0;q<(j +1 -z_int_stat_index[i]);q++)
								{
									average_res[q] = 0.0;
									hz[q] = 0.0;
								}
							
								/*The resistivity will be determined along the ray between the stations and the perturbated cells*/
								DetMeanResWithinLayer(calc_2D_mt,res_slice, average_res, hz, z_int_stat_index[i], first_cell_dz[i], j, nz_without_air, j_freq, &res_station, x_index[s], x_index_stat[i][n], flag_ray_based);

								for(q=0;q<4*nr_of_samples;q++)
								{
									fft_input_B1[i][n][s][q] = 0.0;
									fft_input_E[i][n][s][q] = 0.0;
								}

								/*Calculate the Bx and Ey component Eq.4.30 and 4.31 in Marion thesis*/
								CalcBx_Ey(nr_of_samples, p, fft_input_B1[i][n][s], fft_input_E[i][n][s], average_res, hz, nhz, calc_2D_mt.freq[j_freq], waterdepths[i], res_water[i][n], index_water_layer[i][n], sampling_int, res_station, flag_reflect_coeff);

								/*Fourier Transformation*/
								four1(fft_input_B1[i][n][s]-1,2*nr_of_samples,-1);
								four1(fft_input_E[i][n][s]-1,2*nr_of_samples,-1);

								free(average_res);
								free(hz);
							}
						}
					}
				}
			}

			/*x-position of the cell centers*/
			x_pos = calc_2D_mt.org[0][j_freq] + 0.5*calc_2D_mt.hx[j_freq][0];

			/*Loop over all cells in x-direction*/
			for(i=0;i<calc_2D_mt.nx[j_freq];i++)
			{
				
				/*Find the E-fields in the cells*/
				if(j < nz_without_air -1)
				{
					ave_Ey_r = 0.5*(Ey_r[i*nz_without_air + j] + Ey_r[i*nz_without_air + (j +1)]);
					ave_Ey_i = 0.5*(Ey_i[i*nz_without_air + j] + Ey_i[i*nz_without_air + (j +1)]);
 				}
				else
				{
					ave_Ey_r = Ey_r[i*nz_without_air + j];
					ave_Ey_i = Ey_i[i*nz_without_air + j];
				}

				/***********************************************/
				/*Calculate the currents from each perturbation (eq.4.29 in Marions thesis)*/
				a_current_R = perturb[i*nz_without_air +j]*ave_Ey_r*calc_2D_mt.hx[j_freq][i]*calc_2D_mt.hz[j_freq][j+nz_air]; /*real part*/
				a_current_I = perturb[i*nz_without_air +j]*ave_Ey_i*calc_2D_mt.hx[j_freq][i]*calc_2D_mt.hz[j_freq][j+nz_air]; /*imaginary part*/
				/***********************************************/

				/*Loop over all measurements*/
				for(k=0;k<nr_of_zcoord;k++)
				{
					/*REMARK: Consider only cells BELOW the station*/
					if(z_int_stat_index[k]<=j)
					{
						for(n=0;n<nr_of_stat_per_zcoord[k];n++)
						{
							/**************************************************/
							/*Determine the sample that corresponds to the horizontal distance between the station and the perturbated cell:*/
							dx = fabs(calc_2D_mt.x[l_index_stat[k][n]] - x_pos)/(sampling_int);
							idx = (long) dx;

							/*Find the closest sample*/
							if((dx-idx) > 0.5)
							  idx = idx + 1;

							/**************************************************/
							/*ATTENTION !!	DO NOT CONSIDER the sensitivities of CELLS that are have a LARGER horizontal DISTANCE from the*/
							/*				station position than the last calculated sample:*/

							if(idx <= nr_of_samples)
							{
								/*If the conductivity are calculated by averaging along the layers*/
								if(flag_ray_based == 0)
								{
									/*Determine the Bx and Ey fields at the station positions*/
									Hx_stat.r = fft_input_B1[k][0][0][2*idx];
									Hx_stat.i = fft_input_B1[k][0][0][2*idx + 1];

									Ey_stat.r = fft_input_E[k][0][0][2*idx];
									Ey_stat.i = fft_input_E[k][0][0][2*idx + 1];
								}
								/*If the conductivity is calculated along rays between the station and the perturbated cell*/
								else
								{

									/*For the cells at the left border*/
									if(i < calc_2D_mt.nx_left_border[j_freq])
									{
										/*Determine the Bx and Ey fields at the station positions*/
										Hx_stat.r = fft_input_B1[k][n][0][2*idx];
										Hx_stat.i = fft_input_B1[k][n][0][2*idx + 1];

										Ey_stat.r = fft_input_E[k][n][0][2*idx];
										Ey_stat.i = fft_input_E[k][n][0][2*idx + 1];
									}
									/*For the cells at the right border*/
									else if(i >= calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq])
									{
										/*Determine the Bx and Ey fields at the station positions*/
										Hx_stat.r = fft_input_B1[k][n][mult_cell-1][2*idx];
										Hx_stat.i = fft_input_B1[k][n][mult_cell-1][2*idx + 1];

										Ey_stat.r = fft_input_E[k][n][mult_cell-1][2*idx];
										Ey_stat.i = fft_input_E[k][n][mult_cell-1][2*idx + 1];
									}
									/*For the cells in the modell*/
									else
									{
										cell_index = (i - calc_2D_mt.nx_left_border[j_freq])/nr_of_cells_in_int;
									
										/*Determine the Bx and Ey fields at the station positions*/
										Hx_stat.r = fft_input_B1[k][n][cell_index][2*idx];
										Hx_stat.i = fft_input_B1[k][n][cell_index][2*idx + 1];
	
										Ey_stat.r = fft_input_E[k][n][cell_index][2*idx];
										Ey_stat.i = fft_input_E[k][n][cell_index][2*idx + 1];
									}
								}

								/*Multiply the disturbances with the output of the fourier transform*/
								Hx_stat = RCmul(0.5,Hx_stat);
								Ey_stat = RCmul(M_KONST*0.5,Ey_stat);

								a_current = Complex(a_current_R, a_current_I);
							
								Hx_stat = Cmul(Hx_stat,a_current);
								Ey_stat = Cmul(Ey_stat,a_current);
							}
							else
							{
								Hx_stat.r = 0.0;
								Hx_stat.i = 0.0;

								Ey_stat.r = 0.0;
								Ey_stat.i = 0.0;
							}

							/**************************************************/
							/*Assign the Ey and the Hx to the corresponding measurement:*/

							d_Er[l_index_data[k][n]][i*nz_without_air + j] = (double)Ey_stat.r;
							d_Ei[l_index_data[k][n]][i*nz_without_air + j] = (double)Ey_stat.i;
							d_Hr[l_index_data[k][n]][i*nz_without_air + j] = (double)Hx_stat.r;
							d_Hi[l_index_data[k][n]][i*nz_without_air + j] = (double)Hx_stat.i;
							/**************************************************/

						}
					}
				}
					
				if(i != calc_2D_mt.nx[j_freq] -1)
					x_pos = x_pos + 0.5*(calc_2D_mt.hx[j_freq][i] + calc_2D_mt.hx[j_freq][i + 1]);

			}

			/*z-position of the cell centers*/
			if(j != nz_without_air -1)
				z_pos = z_pos + 0.5*(calc_2D_mt.hz[j_freq][nz_air + j] + calc_2D_mt.hz[j_freq][nz_air + j + 1]);

		}

		free(Ey_r);
		free(Ey_i);
	}


	/************************/
	/************************/
		/*TM-mode*/
	/************************/
	/************************/

	if(kind_of_mt == 2)
	{
		Ex_r = (double *)memory(NULL,calc_2D_mt.nx[j_freq]*nz_without_air,sizeof(double),"Deriv_MT_analytical");
		Ex_i = (double *)memory(NULL,calc_2D_mt.nx[j_freq]*nz_without_air,sizeof(double),"Deriv_MT_analytical");
		Ez_r = (double *)memory(NULL,calc_2D_mt.nx[j_freq]*nz_without_air,sizeof(double),"Deriv_MT_analytical");
		Ez_i = (double *)memory(NULL,calc_2D_mt.nx[j_freq]*nz_without_air,sizeof(double),"Deriv_MT_analytical");

		/*Read in the Ey an Ez component from a temporary file*/
		ReadTResistivityIn(i_layer, calc_2D_mt.freq[j_freq], calc_2D_mt.nx[j_freq]*nz_without_air, Ex_r, 1);
		ReadTResistivityIn(i_layer, calc_2D_mt.freq[j_freq], calc_2D_mt.nx[j_freq]*nz_without_air, Ex_i, 2);
		ReadTResistivityIn(i_layer, calc_2D_mt.freq[j_freq], calc_2D_mt.nx[j_freq]*nz_without_air, Ez_r, 9);
		ReadTResistivityIn(i_layer, calc_2D_mt.freq[j_freq], calc_2D_mt.nx[j_freq]*nz_without_air, Ez_i, 10);

		/*z-position of the cell centers*/
		z_pos = calc_2D_mt.org[1][j_freq] + 0.5*calc_2D_mt.hz[j_freq][nz_air];

		/*Loop over all cells*/
		for(j=0;j<nz_without_air;j++)
		{

			/*************************************************************************************/			
			/*Loop over all depth intervals where stations are placed*/
			for(i=0;i<nr_of_zcoord;i++)
			{
				/*REMARK: Consider only cells BELOW the station*/
				if(z_int_stat_index[i]<=j)
				{
					/*If the conductivity are calculated by averaging along the layers*/
					if(flag_ray_based == 0)
					{
						/***********************************************/ 
						/*Determine the average resistivity within the layers between the station and the perturbated cells*/
						average_res = (double *)memory(NULL,(j +1 -z_int_stat_index[i]),sizeof(double),"Deriv_MT_analytical");
						hz = (double *)memory(NULL,(j +1 -z_int_stat_index[i]),sizeof(double),"Deriv_MT_analytical");

						nhz = j +1 -z_int_stat_index[i];

						for(q=0;q<(j +1 -z_int_stat_index[i]);q++)
						{
							average_res[q] = 0.0;
							hz[q] = 0.0;
						}

						/*The average resistivity for each layer will be calculated*/
						DetMeanResWithinLayer(calc_2D_mt,res_slice, average_res, hz, z_int_stat_index[i], first_cell_dz[i], j, nz_without_air, j_freq, &res_station, 0.0, 0.0, flag_ray_based);
						/*************************************************************************************/
				
						for(q=0;q<4*nr_of_samples;q++)
						{
							fft_input_B1[i][0][0][q] = 0.0;
							fft_input_B2[i][0][0][q] = 0.0;
						}
					
						for(q=0;q<4*nr_of_samples_Ex;q++)
						{
							fft_input_Ex1[i][0][0][q] = 0.0;
							fft_input_Ex2[i][0][0][q] = 0.0;
						}
	
						/*Calculate the By and Ex component Eq.4.42 and 4.43 in Marion thesis*/
						/*x-directed dipole component*/
						CalcBy_Ex_Ix(nr_of_samples, nr_of_samples_Ex, p, p_Ex, fft_input_B1[i][0][0], fft_input_Ex1[i][0][0], average_res, hz, nhz, calc_2D_mt.freq[j_freq],waterdepths[i], res_water[i][0], sampling_int, sampling_int_Ex, res_station, flag_reflect_coeff);
						/*z-directed dipole component*/
						CalcBy_Ex_Iz(nr_of_samples, nr_of_samples_Ex, p, p_Ex, fft_input_B2[i][0][0], fft_input_Ex2[i][0][0], average_res, hz, nhz, calc_2D_mt.freq[j_freq],waterdepths[i], res_water[i][0], sampling_int, sampling_int_Ex, res_station, flag_reflect_coeff);

							/*Fourier Transformation*/
							/*x-directed dipole component*/
						four1(fft_input_B1[i][0][0]-1,2*nr_of_samples,-1);
						four1(fft_input_Ex1[i][0][0]-1,2*nr_of_samples_Ex,-1);
							/*z-directed dipole component*/
						four1(fft_input_B2[i][0][0]-1,2*nr_of_samples,-1);
						four1(fft_input_Ex2[i][0][0]-1,2*nr_of_samples_Ex,-1);

						free(average_res);
						free(hz);
					}
					/*If the conductivity is calculated along rays between the station and the perturbated cell*/
					else
					{
						/*Loop over all stations in a layer*/
						for(n=0;n<nr_of_stat_per_zcoord[i];n++)
						{
							/*Loop over all "rays"*/
							for(s=0;s<mult_cell;s++)
							{
								/***********************************************/ 
								/*Determine the average resistivity within the layers between the station and the perturbated cells*/
								average_res = (double *)memory(NULL,(j +1 -z_int_stat_index[i]),sizeof(double),"Deriv_MT_analytical");
								hz	= (double *)memory(NULL,(j +1 -z_int_stat_index[i]),sizeof(double),"Deriv_MT_analytical");

								nhz = j +1 -z_int_stat_index[i];

								for(q=0;q<(j +1 -z_int_stat_index[i]);q++)
								{
									average_res[q] = 0.0;
									hz[q] = 0.0;
								}
							

								/*The resistivity will be determined along the ray between the stations and the perturbated cells*/
								DetMeanResWithinLayer(calc_2D_mt,res_slice, average_res, hz, z_int_stat_index[i], first_cell_dz[i], j, nz_without_air, j_freq, &res_station, x_index[s], x_index_stat[i][n], flag_ray_based);

								for(q=0;q<4*nr_of_samples;q++)
								{
									fft_input_B1[i][n][s][q] = 0.0;
									fft_input_B2[i][n][s][q] = 0.0;
								}
						
								for(q=0;q<4*nr_of_samples_Ex;q++)
								{
									fft_input_Ex1[i][n][s][q] = 0.0;
									fft_input_Ex2[i][n][s][q] = 0.0;
								}

								/*Calculate the By and Ex component Eq.4.42 and 4.43 in Marion thesis*/
								/*x-directed dipole component*/
								CalcBy_Ex_Ix(nr_of_samples, nr_of_samples_Ex, p, p_Ex, fft_input_B1[i][n][s], fft_input_Ex1[i][n][s], average_res, hz, nhz, calc_2D_mt.freq[j_freq],waterdepths[i], res_water[i][n], sampling_int, sampling_int_Ex, res_station, flag_reflect_coeff);
								/*z-directed dipole component*/
								CalcBy_Ex_Iz(nr_of_samples, nr_of_samples_Ex, p, p_Ex, fft_input_B2[i][n][s], fft_input_Ex2[i][n][s], average_res, hz, nhz, calc_2D_mt.freq[j_freq],waterdepths[i], res_water[i][n], sampling_int, sampling_int_Ex, res_station, flag_reflect_coeff);

								/*Fourier Transformation*/
								/*x-directed dipole component*/
								four1(fft_input_B1[i][n][s]-1,2*nr_of_samples,-1);
								four1(fft_input_B2[i][n][s]-1,2*nr_of_samples,-1);
								/*z-directed dipole component*/
								four1(fft_input_Ex1[i][n][s]-1,2*nr_of_samples_Ex,-1);
								four1(fft_input_Ex2[i][n][s]-1,2*nr_of_samples_Ex,-1);
	
								free(average_res);
								free(hz);
							}
						}
					}
				}
			}

			/*x-position of the cell centers*/
			x_pos = calc_2D_mt.org[0][j_freq] + 0.5*calc_2D_mt.hx[j_freq][0];

			for(i=0;i<calc_2D_mt.nx[j_freq];i++)
			{
				/*Find the E-fields in the cells*/
				if(j < nz_without_air -1)
				{
					ave_Ex_r = 0.5*(Ex_r[i*nz_without_air + j] + Ex_r[i*nz_without_air + (j +1)]);
					ave_Ex_i = 0.5*(Ex_i[i*nz_without_air + j] + Ex_i[i*nz_without_air + (j +1)]);
					ave_Ez_r = 0.5*(Ez_r[i*nz_without_air + j] + Ez_r[i*nz_without_air + (j +1)]);
					ave_Ez_i = 0.5*(Ez_i[i*nz_without_air + j] + Ez_i[i*nz_without_air + (j +1)]);
 				}
				else
				{
					ave_Ex_r = Ex_r[i*nz_without_air + j];
					ave_Ex_i = Ex_i[i*nz_without_air + j];
					ave_Ez_r = Ez_r[i*nz_without_air + j];
					ave_Ez_i = Ez_i[i*nz_without_air + j];
				}
	
				/***********************************************/
				/*Calculate the magnetic moments from each perturbation (eq.4.41 in Marions thesis)*/
				tmp_sigma = (perturb[i*nz_without_air +j]*(1/res_slice[i*nz_without_air +j]))/((2/res_slice[i*nz_without_air +j])+ perturb[i*nz_without_air +j]);
	
				ave_Ex = Complex(ave_Ex_r, ave_Ex_i);
				ave_Ez = Complex(ave_Ez_r, ave_Ez_i);
	
				/*x-component*/
				a_magmoment_R1 = 2*calc_2D_mt.hx[j_freq][i]*calc_2D_mt.hz[j_freq][j+nz_air]*(double)ave_Ex.r*tmp_sigma;
				a_magmoment_I1 = 2*calc_2D_mt.hx[j_freq][i]*calc_2D_mt.hz[j_freq][j+nz_air]*(double)ave_Ex.i*tmp_sigma;
				/*z-component*/
				a_magmoment_R2 = 2*calc_2D_mt.hx[j_freq][i]*calc_2D_mt.hz[j_freq][j+nz_air]*(double)ave_Ez.r*tmp_sigma;
				a_magmoment_I2 = 2*calc_2D_mt.hx[j_freq][i]*calc_2D_mt.hz[j_freq][j+nz_air]*(double)ave_Ez.i*tmp_sigma;
				/***********************************************/

				/*Loop over all measurements*/
				for(k=0;k<nr_of_zcoord;k++)
				{
					/*REMARK: Consider only cells BELOW the station*/
					if(z_int_stat_index[k]<=j)
					{
						for(n=0;n<nr_of_stat_per_zcoord[k];n++)
						{
							/**************************************************/
										/*For the By component:*/
							/**************************************************/
							/*Determine the sample that corresponds to the horizontal distance between the station and the perturbated cell:*/
							dx = fabs(calc_2D_mt.x[l_index_stat[k][n]] - x_pos)/(sampling_int);
							idx = (long) dx;

							/*Find the closest sample*/
							if((dx-idx) > 0.5)
							  idx = idx + 1;

							/**************************************************/
							/*ATTENTION !!	DO NOT CONSIDER the sensitivities of CELLS that are have a LARGER horizontal DISTANCE from the*/
							/*				station position than the last calculated sample:*/
						
							if(idx <= nr_of_samples)
							{
								/*Determine the Bx fields at the station positions*/
									/*If the conductivity are calculated by averaging along the layers*/
								if(flag_ray_based == 0)
								{
									Hy_stat1.r = fft_input_B1[k][0][0][2*idx];
									Hy_stat1.i = fft_input_B1[k][0][0][2*idx + 1];
							
									Hy_stat2.r = fft_input_B2[k][0][0][2*idx];
									Hy_stat2.i = fft_input_B2[k][0][0][2*idx + 1];
								}
								/*If the conductivity is calculated along rays between the station and the perturbated cell*/
								else
								{
									/*For the cells at the left border*/
									if(i < calc_2D_mt.nx_left_border[j_freq])
									{
										Hy_stat1.r = fft_input_B1[k][n][0][2*idx];
										Hy_stat1.i = fft_input_B1[k][n][0][2*idx + 1];
									
										Hy_stat2.r = fft_input_B2[k][n][0][2*idx];
										Hy_stat2.i = fft_input_B2[k][n][0][2*idx + 1];
									}
									/*For the cells at the right border*/
									else if(i >= calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq])
									{
										Hy_stat1.r = fft_input_B1[k][n][mult_cell-1][2*idx];
										Hy_stat1.i = fft_input_B1[k][n][mult_cell-1][2*idx + 1];
		
										Hy_stat2.r = fft_input_B2[k][n][mult_cell-1][2*idx];
										Hy_stat2.i = fft_input_B2[k][n][mult_cell-1][2*idx + 1];
									}
									/*For the cells in the modell*/
									else
									{
										cell_index = (i - calc_2D_mt.nx_left_border[j_freq])/nr_of_cells_in_int;
									
										/*Determine the Bx and Ey fields at the station positions*/
										Hy_stat1.r = fft_input_B1[k][n][cell_index][2*idx];
										Hy_stat1.i = fft_input_B1[k][n][cell_index][2*idx + 1];
						
										Hy_stat2.r = fft_input_B2[k][n][cell_index][2*idx];
										Hy_stat2.i = fft_input_B2[k][n][cell_index][2*idx + 1];
									}
								}

								/*Multiply the disturbances with the output of the fourier transform*/
									/*x-component*/
								Hy_stat1 = RCmul(1/2.0,Hy_stat1);
								a_magmoment_x = Complex(a_magmoment_R1, a_magmoment_I1);
								Hy_stat1 = Cmul(Hy_stat1,a_magmoment_x);
									/*z-component*/
								Hy_stat2 = RCmul(1/2.0,Hy_stat2);
								a_magmoment_z = Complex(a_magmoment_R2, a_magmoment_I2);
								Hy_stat2 = Cmul(Hy_stat2,a_magmoment_z);

								/*Sum up the effects of both components*/
								Hy_stat1 = Cadd(Hy_stat1,Hy_stat2);
	
							}
							else
							{
								Hy_stat1.r = 0.0;
								Hy_stat1.i = 0.0;
							}
	
							/**************************************************/
									/*For the Ex component:*/
							/**************************************************/
							/*Determine the sample that corresponds to the horizontal distance between the station and the perturbated cell:*/
							dx = fabs(calc_2D_mt.x[l_index_stat[k][n]] - x_pos)/(sampling_int_Ex);
							idx = (long) dx;	

							/*Find the closest sample*/
							if((dx-idx) > 0.5)
							  idx = idx + 1;

							/**************************************************/
							/*ATTENTION !!	DO NOT CONSIDER the sensitivities of CELLS that are have a LARGER horizontal DISTANCE from the*/
							/*				station position than the last calculated sample:*/
							if(idx <= nr_of_samples_Ex)
							{
								/*Determine the Ey fields at the station positions*/
									/*If the conductivity are calculated by averaging along the layers*/
								if(flag_ray_based == 0)
								{
									Ex_stat1.r = fft_input_Ex1[k][0][0][2*idx];
									Ex_stat1.i = fft_input_Ex1[k][0][0][2*idx + 1];	

									Ex_stat2.r = fft_input_Ex2[k][0][0][2*idx];
									Ex_stat2.i = fft_input_Ex2[k][0][0][2*idx + 1];
								}
								/*If the conductivity is calculated along rays between the station and the perturbated cell*/
								else
								{
									/*For the cells at the left border*/
									if(i < calc_2D_mt.nx_left_border[j_freq])
									{
										Ex_stat1.r = fft_input_Ex1[k][n][0][2*idx];
										Ex_stat1.i = fft_input_Ex1[k][n][0][2*idx + 1];	

										Ex_stat2.r = fft_input_Ex2[k][n][0][2*idx];
										Ex_stat2.i = fft_input_Ex2[k][n][0][2*idx + 1];
									}
									/*For the cells at the right border*/
									else if(i >= calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq])
									{
										Ex_stat1.r = fft_input_Ex1[k][n][mult_cell -1][2*idx];
										Ex_stat1.i = fft_input_Ex1[k][n][mult_cell -1][2*idx + 1];

										Ex_stat2.r = fft_input_Ex2[k][n][mult_cell -1][2*idx];
										Ex_stat2.i = fft_input_Ex2[k][n][mult_cell -1][2*idx + 1];
									}
									/*For the cells in the modell*/
									else
									{
										cell_index = (i - calc_2D_mt.nx_left_border[j_freq])/nr_of_cells_in_int;
										
										/*Determine the Bx and Ey fields at the station positions*/
										Ex_stat1.r = fft_input_Ex1[k][n][cell_index][2*idx];
										Ex_stat1.i = fft_input_Ex1[k][n][cell_index][2*idx + 1];	

										Ex_stat2.r = fft_input_Ex2[k][n][cell_index][2*idx];
										Ex_stat2.i = fft_input_Ex2[k][n][cell_index][2*idx + 1];
									}
								}

								/*Multiply the disturbances with the output of the fourier transform*/
									/*x-component*/
								Ex_stat1 = RCmul(M_KONST/2.0,Ex_stat1);
								a_magmoment_x = Complex(a_magmoment_R1, a_magmoment_I1);
								Ex_stat1 = Cmul(Ex_stat1,a_magmoment_x);
									/*z-component*/
								Ex_stat2 = RCmul(M_KONST/2.0,Ex_stat2);
								a_magmoment_z = Complex(a_magmoment_R2, a_magmoment_I2);
								Ex_stat2 = Cmul(Ex_stat2,a_magmoment_z);

								/*Sum up the effects of both components*/
								Ex_stat1 = Cadd(Ex_stat1,Ex_stat2);
							}
							else
							{
								Ex_stat1.r = 0.0;
								Ex_stat1.i = 0.0;
							}

							/**************************************************/
							/*Assign the Ex and the Hy to the corresponding measurement:*/
							d_Er[l_index_data[k][n]][i*nz_without_air + j] = (double)Ex_stat1.r;
							d_Ei[l_index_data[k][n]][i*nz_without_air + j] = (double)Ex_stat1.i;
							d_Hr[l_index_data[k][n]][i*nz_without_air + j] = (double)Hy_stat1.r;
							d_Hi[l_index_data[k][n]][i*nz_without_air + j] = (double)Hy_stat1.i;
							/**************************************************/

						}
					}
				}

				if(i != calc_2D_mt.nx[j_freq] -1)
					x_pos = x_pos + 0.5*(calc_2D_mt.hx[j_freq][i] + calc_2D_mt.hx[j_freq][i + 1]);
			}

			/*z-position of the cell centers*/
			if(j != nz_without_air -1)
				z_pos = z_pos + 0.5*(calc_2D_mt.hz[j_freq][nz_air + j] + calc_2D_mt.hz[j_freq][nz_air + j + 1]);

		}

		free(Ex_r);
		free(Ex_i);
		free(Ez_r);
		free(Ez_i);
	}


	/*********************************************/
	/*Free the memory*/	

	free(res_slice);

	for(i=0;i<nr_of_zcoord;i++)
	{
		free(res_water[i]);
		free(index_water_layer[i]);
		free(l_index_data[i]);
		free(l_index_stat[i]);
		if(flag_ray_based != 0)
			free(x_index_stat[i]);
	}

	for(i=0;i<nr_of_zcoord;i++)
	{
		/*If the conductivity are calculated by averaging along the layers*/
		if(flag_ray_based == 0)
		{
			free(fft_input_B1[i][0][0]);
			free(fft_input_B2[i][0][0]);
			free(fft_input_E[i][0][0]);
			free(fft_input_Ex1[i][0][0]);
			free(fft_input_Ex2[i][0][0]);

			free(fft_input_B1[i][0]);
			free(fft_input_B2[i][0]);
			free(fft_input_E[i][0]);
			free(fft_input_Ex1[i][0]);
			free(fft_input_Ex2[i][0]);

			free(fft_input_B1[i]);
			free(fft_input_B2[i]);
			free(fft_input_E[i]);
			free(fft_input_Ex1[i]);
			free(fft_input_Ex2[i]);
		}
		/*If the conductivity is calculated along rays between the station and the perturbated cell*/

		else
		{
			for(j=0;j<nr_of_stat_per_zcoord[i];j++)
			{
				for(k=0;k<mult_cell;k++)
				{
					free(fft_input_B1[i][j][k]);
					free(fft_input_B2[i][j][k]);
					free(fft_input_E[i][j][k]);
					free(fft_input_Ex1[i][j][k]);
					free(fft_input_Ex2[i][j][k]);
				}

				free(fft_input_B1[i][j]);
				free(fft_input_B2[i][j]);
				free(fft_input_E[i][j]);
				free(fft_input_Ex1[i][j]);
				free(fft_input_Ex2[i][j]);
			}
			
			free(fft_input_B1[i]);
			free(fft_input_B2[i]);
			free(fft_input_E[i]);
			free(fft_input_Ex1[i]);
			free(fft_input_Ex2[i]);

		}
	}

	free(nr_of_stat_per_zcoord);
	free(z_int_stat_index);
	free(z_int_stat);
	free(first_cell_dz);
	free(p);
	free(p_Ex);

	free(l_index_data);
	free(l_index_stat);
	free(x_index_stat);
	free(fft_input_B1);
	free(fft_input_B2);
	free(fft_input_E);
	free(fft_input_Ex1);
	free(fft_input_Ex2);

	free(waterdepths);
	free(res_water);
	free(index_water_layer);

	free(x_index);

	return(0);
}

#undef EPS 
#undef RES_AIR
#undef DELTA_X
#undef DELTA_X_EX
#undef MULTI_X
#undef MULT_CELL
#undef RAY_BASED
#undef REFLECT_COEFF

/*-------------------------------------------------------------*/
/*Assign the perturbations to the refined grid*/
/*Parameter: calc_2D_mt := 2D MT structure*/
/*           *perturb := perturbations	    (input:	 parameters in the original constant grid*/
/*											(output: parameters in the refined grid for the MT calculations)*/
/*             grid_h := grid size of the original grid*/
/*              nx,nz := number of cells in x and z direction in the original grid*/
/*                  j := frequency index*/

int AssignRefinedGridMT(CALC_2D_MT calc_2D_mt, double *perturb, double grid_h, long nx, long nz, long j)
{
	long k,m,n,p;
	long nz_index, nx_index;
	long nz_without_base, nz_without_air; /*nr of cells in z-direction without basement and air, respectively*/
	long l_scal_factor_x ,l_scal_factor_z;
	double *tmp_perturb;
	double scal_factor_x ,scal_factor_z ,resd;

	/*Number of cells in the z-direction for the refined grid */
	nz_without_base = calc_2D_mt.nz[j] - calc_2D_mt.nz_basement[j]; /*(without basement)*/
	nz_without_air = calc_2D_mt.nz[j] - calc_2D_mt.nz_atm[j] - calc_2D_mt.nz_ionos[j];  /*(without air and ionosphere)*/
			
	tmp_perturb = (double *)memory(NULL,(nx*nz), sizeof(double),"AssignRefinedGridMT");

	for(k=0;k<(nx*nz);k++)
		tmp_perturb[k] = perturb[k];

	nz_index = calc_2D_mt.nz_atm[j] + calc_2D_mt.nz_ionos[j];

	/*Loop over cells in the z-direction for the unrefined grid*/
	for(k=0;k<nz;k++)
	{
		/**************/
		/*ATTENTION!! The routine is only working, if the size of the forward cells is a multiple of the size of the regridded cells*/
		/*within wach forward cell*/
		/*determine the scaling factor in z-direction*/
		resd = modf((grid_h/calc_2D_mt.hz[j][nz_index]), &scal_factor_z);

		if(resd > 0.5)
			scal_factor_z = scal_factor_z + 1.0;

		l_scal_factor_z = (long)scal_factor_z;
		/**************/

		/*Loop over all cells (in z-direction) of the refined grid which belongs to the corresponding unrefined cell*/
		for(m=nz_index; m<(l_scal_factor_z + nz_index) ;m++)
		{
			/*Add the perturbations perturb1 and perturb2 at the left and right border of the refined grid*/
				/*left*/
			for(n=0;n<calc_2D_mt.nx_left_border[j];n++)
				perturb[nz_without_air*n + (m - calc_2D_mt.nz_atm[j] - calc_2D_mt.nz_ionos[j])] = tmp_perturb[k];

				/*right*/
			for(n=(calc_2D_mt.nx[j] - calc_2D_mt.nx_right_border[j]); n<calc_2D_mt.nx[j]; n++)
				perturb[nz_without_air*n + (m - calc_2D_mt.nz_atm[j] - calc_2D_mt.nz_ionos[j])] = tmp_perturb[(nx-1)*nz + k];

			nx_index = calc_2D_mt.nx_left_border[j];


			/*Loop over cells in the x-direction for the unrefined grid*/
			for(n=0;n<nx;n++)
			{
				/**************/
				/*determine the scaling factor in x-direction*/
				resd = modf((grid_h/calc_2D_mt.hx[j][nx_index]), &scal_factor_x);

				if(resd > 0.5)
					scal_factor_x = scal_factor_x + 1;

				l_scal_factor_x = (long)scal_factor_x;
				/**************/

				/*Loop over all cells (in x-direction) of the refined grid which belongs to the corresponding unrefined cell*/
				for(p=nx_index; p<(l_scal_factor_x + nx_index) ;p++)
					/*Add the perturbances in the central part of the refined grid*/
					perturb[nz_without_air*p + (m - calc_2D_mt.nz_atm[j] - calc_2D_mt.nz_ionos[j])] = tmp_perturb[n*nz + k];

				nx_index = p;
			}

			if(nx_index != (calc_2D_mt.nx[j] - calc_2D_mt.nx_right_border[j]))
			{
				printf("The number of cells in x-direction do NOT correspond!!\n");
				printf("during determing the perturbances in the refined cells\n");
				printf("Originally: %d, Lateron: %d\n", (calc_2D_mt.nx[j] - calc_2D_mt.nx_right_border[j]), nx_index);
				exit(0);
			}

		}

		nz_index = m;
	}

	if(nz_index != nz_without_base)
	{
		printf("The number of cells in z-direction do NOT correspond!!\n");
		printf("during determing the perturbances in the refined cells\n");
		printf("Originally: %d, Lateron: %d\n", nz_without_base, nz_index);
		exit(0);
	}

	/*****************/
	/*Fill the cells at the basement with the perturbations*/
	for(k=(calc_2D_mt.nz[j]-calc_2D_mt.nz_basement[j]);k<calc_2D_mt.nz[j] ;k++)
	{
		for(n=0;n<calc_2D_mt.nx[j];n++)
			perturb[nz_without_air*n + k - calc_2D_mt.nz_atm[j] - calc_2D_mt.nz_ionos[j]] =
			perturb[nz_without_air*n + nz_without_base - 1 - calc_2D_mt.nz_atm[j] - calc_2D_mt.nz_ionos[j]];
	}

	free(tmp_perturb);

	return(0);
}


#undef res
#undef R1
#undef R2

/*-----------------------------------------------*/
/*Transfer the sensitivities to the cells of the forward grid*/
/*Parameter: calc_2D_mt := 2D MT structure*/
/*           *dsens := perturbations	    (input:	 sensitivities in the original refined grid*/
/*											(output: sensitivities in the forward grid)*/
/*             grid_h := grid size of the forward grid*/
/*              nx,nz := number of cells in x and z direction in the forward grid*/
/*             j_freq := frequency index*/

int SumSensToForwardCells(CALC_2D_MT calc_2D_mt, double *dSens, double grid_h, long nx, long nz, long j_freq)
{
	long i,j,m,n,p;
	long nx_index, nz_index, l_scal_factor_x, l_scal_factor_z;
	long nz_without_base, nz_without_air, nz_air;
	double scal_factor_x, scal_factor_z, resd;
	double *tmp_dSens;

		/*Number of cells in the z-direction for the refined grid */
	nz_without_base = calc_2D_mt.nz[j_freq] - calc_2D_mt.nz_basement[j_freq]; /*(without basement)*/
	nz_without_air = calc_2D_mt.nz[j_freq] - calc_2D_mt.nz_atm[j_freq] - calc_2D_mt.nz_ionos[j_freq];  /*(without air and ionosphere)*/
	nz_air = calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq];

	nz_index = calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq];

	tmp_dSens = (double *)memory(NULL,(nx*nz),sizeof(double),"SumSensToForwardCells");

	for(i=0;i<(nx*nz);i++)
		tmp_dSens[i] = 0.0;

	/*Loop over all forward cells in z-direction*/
	for(i=0;i<nz;i++)
	{

		/**************/
		/*ATTENTION!! The routine is only working, if the size of the forward cells is a multiple of the size of the regridded cells*/
		/*within wach forward cell*/
		/*determine the scaling factor in z-direction*/
		resd = modf((grid_h/calc_2D_mt.hz[j_freq][nz_index]), &scal_factor_z);

		if(resd > 0.5)
			scal_factor_z = scal_factor_z + 1.0;

		l_scal_factor_z = (long)scal_factor_z;
		/**************/

		/*Loop over all cells (in z-direction) of the refined grid which belongs to the corresponding unrefined cell*/
		for(m=nz_index; m<(l_scal_factor_z + nz_index) ;m++)
		{

			/*Sum sensitivities of the cells at the left border*/
			for(n=0;n<calc_2D_mt.nx_left_border[j_freq];n++)
				tmp_dSens[i] = tmp_dSens[i] + dSens[n*nz_without_air + (m - nz_air)]; 
				
			nx_index = calc_2D_mt.nx_left_border[j_freq];

			/*Sum sensitivities of the cells in the central part*/
			for(j=0;j<nx;j++)
			{
				/**************/
				/*determine the scaling factor in x-direction*/
				resd = modf((grid_h/calc_2D_mt.hx[j_freq][nx_index]), &scal_factor_x);

				if(resd > 0.5)
					scal_factor_x = scal_factor_x + 1.0;

				l_scal_factor_x = (long)scal_factor_x;
				/**************/

				/*Loop over all cells (in x-direction) of the refined grid which belongs to the corresponding unrefined cell*/
				for(p=nx_index; p<(l_scal_factor_x + nx_index) ;p++)
					/*Add the sensitivities in the central part of the refined grid*/
					tmp_dSens[j*nz + i] = tmp_dSens[j*nz + i] + dSens[p*nz_without_air + (m - nz_air)]; 

				nx_index = p;

			}

			if(nx_index != (calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq]))
			{
				printf("The number of cells in x-direction do NOT correspond!!\n");
				printf("during determing the perturbances in the refined cells\n");
				printf("Originally: %d, Lateron: %d\n", (calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq]), nx_index);
				exit(0);
			}


			/*Sum sensitivities of the cells at the right border*/
			for(n=calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq];n<calc_2D_mt.nx[j_freq];n++)
				tmp_dSens[(nx-1)*nz + i] = tmp_dSens[(nx-1)*nz + i] + dSens[n*nz_without_air + (m - nz_air)];

		}

		nz_index = m;
	}

	if(nz_index != nz_without_base)
	{
		printf("The number of cells in z-direction do NOT correspond!!\n");
		printf("during determing the perturbances in the refined cells\n");
		printf("Originally: %d, Lateron: %d\n", nz_without_base, nz_index);
		exit(0);
	}

	/*Sum the sensitivities of the cells at the basement*/
	/*Loop over all cells (in z-direction) of the refined grid which belongs to the corresponding unrefined cell*/
	for(m=nz_without_base; m<calc_2D_mt.nz[j_freq] ;m++)
	{
		nx_index = calc_2D_mt.nx_left_border[j_freq];

		for(j=0;j<nx;j++)
		{
			/**************/
			/*determine the scaling factor in x-direction*/
			resd = modf((grid_h/calc_2D_mt.hx[j_freq][nx_index]), &scal_factor_x);

			if(resd > 0.5)
				scal_factor_x = scal_factor_x + 1.0;

			l_scal_factor_x = (long)scal_factor_x;
			/**************/

			/*Loop over all cells (in x-direction) of the refined grid which belongs to the corresponding unrefined cell*/
			for(p=nx_index; p<(l_scal_factor_x + nx_index) ;p++)
				tmp_dSens[j*nz + (nz-1)] = tmp_dSens[j*nz + (nz-1)] + dSens[p*nz_without_air + (m - nz_air)];

			nx_index = p;
		}
	}

	for(i=0;i<(nx*nz);i++)
		dSens[i] = tmp_dSens[i];


	free(tmp_dSens);

	return(0);
}


/*-----------------------------------------------*/
/*Build the mt (forward) structure for the 2D mt measurements*/
/*REMARK: The derivatives will be written into a temporary file*/
/*Parameter:  calc_2D_mt = 2D MT structure*/
/*            d_Zr, d_Zi = all derivatives of the impedances (real, and imaginary part) for the station k*/
/*               b_index = border_index for the considered slice*/
/*                  grid = grid structure*/
/*                   *mt = mt (forward) structure required for the inversion (will be modified by this routine)*/
/*				  j_freq = frequency index*/
/*                k_stat = station index*/
/*             index_mode = specify the mode (1== TE; 2 == TM)*/
/*			  direc_2D_mt = specify the direction of the slices (1 == parallel to the x-direction, 2 == parallel to the y-direction)*/
/*				   *data  = data structure*/
/*           already_used = Specify, if this is the first time that parameters are assigned to this structure (yes == 0, no == 1)*/
/*      first_cell_strike = first cell of the slice in the strike direction*/

int makeMTstructure2D(CALC_2D_MT calc_2D_mt, double *d_Zr, double *d_Zi, double *b_index, GRID_STRUCT grid, MT_STRUCT *mt, long j_freq, long k_stat, int index_mode, int direc_2D_mt, DATA_STRUCT *data, int *already_used, long first_cell_strike)
{
	int index_found;
	long i,j,k,freq_index;
	long nx,ny,nz,count;
	long ny2,nz2,nyz2;
	double **tmp_deriv;
	double help_deriv;
	
	char fname[40];
	FILE *ouf;

	/********************/
	/*Determine the number of considered cells for the stations*/
	/*Number of cells parallel to the slice*/
	nz = grid.nz;

	if(direc_2D_mt != 2)
		nx = grid.nx;
	else
		nx = grid.ny;

	/*Number of cells vertical to th slice*/
	ny = calc_2D_mt.nr_cells_2d_mt;

	/********************/
	nz2 = (2*grid.nborder) + grid.nz;
	ny2 = (2*grid.nborder) + grid.ny;
	nyz2 = ny2*nz2;
	/********************/

	/*Check, if the structure is already predefined*/
	if(*already_used == 0)
	{
		/********************/
		mt->ncell = nx*ny*nz;
		/********************/

		/*Allocate the memory*/
		mt->n = (long *)memory(NULL,1,sizeof(long),"makeMTstructure2D"); /*NOT used for the forward grid !!!!*/

		/*For the mt->deriv will not be used for the forward cells, because it requires too much memory*/
		/*Instead: the derivatives will be written in a temporary file*/
		mt->deriv = (double **)memory(NULL,1,sizeof(double *),"makeMTstructure2D");

		if(mt->ncell != 0)
			mt->ele = (long *)memory(NULL,mt->ncell,sizeof(long),"makeMTstructure2D");
		else
			mt->ele = (long *)memory(NULL,1,sizeof(long),"makeMTstructure2D");

		*already_used = 1;
	}


	if(mt->ncell != 0)
		tmp_deriv = (double **)memory(NULL,mt->ncell,sizeof(double *),"makeMTstructure2D");
	else
		tmp_deriv = (double **)memory(NULL,1,sizeof(double *),"makeMTstructure2D");

	for(i=0;i<mt->ncell;i++)
	{
		tmp_deriv[i] = (double *)memory(NULL,2,sizeof(double),"makeMTstructure2D");

		for(j=0;j<2;j++)	
			tmp_deriv[i][j] = 0.0;
	}

	/**********************************************/
	/**********************************************/

	/*Find the number of frequencies for this station*/
	mt->nfreq = data->nfreq_mt[calc_2D_mt.freq_data[j_freq][k_stat]];

	/*Find the right frequency index*/
	freq_index = -1;
	index_found = 0;
	for(i=0;i<data->nfreq_mt[calc_2D_mt.freq_data[j_freq][k_stat]];i++)
	{
		if(data->freq_mt[calc_2D_mt.freq_data[j_freq][k_stat]][i] == calc_2D_mt.freq[j_freq])
		{
			freq_index = i;
			index_found = 1;
			break;
		}
	}

	if(index_found == 0)
	{
		printf("The frequency %f Hz was NOT found for the station %d\n", calc_2D_mt.freq[j_freq], data->mno[calc_2D_mt.freq_data[j_freq][k_stat]]);
		exit(0);
	}

	/**********************************************/
	/**********************************************/
	count = 0;

		/*slice is parallel to the x-direction*/
	if(direc_2D_mt != 2)
	{
		for(i=0;i<nx;i++)
			for(j=first_cell_strike;j<first_cell_strike + ny;j++)
				for(k=0;k<nz;k++)
				{
					/*Determine the cell numbers*/
					mt->ele[count] = (nyz2 *(grid.nborder + i)) + (nz2 * (grid.nborder + j)) + (grid.nborder + k);

					/*Fill the derivatives*/

					/*Do NOT fill cells of the air*/
					if(b_index[i*nz + k] > 0.0)
					{
						tmp_deriv[count][0] =	d_Zr[i*nz + k]/((double)ny);				/*Real part*/
						tmp_deriv[count][1] =   d_Zi[i*nz + k]/((double)ny);			/*Imaginary part*/
					}

					count++;

					if(count > mt->ncell)
					{
						printf("The number of cells (%d) in the mt forward structure is specified too small !!!\n",  mt->ncell);
						exit(0);
					}
				}
	}
	else
	{	
		/*slice is parallel to the y-direction*/
		for(i=first_cell_strike;i<first_cell_strike + ny;i++)
			for(j=0;j<nx;j++)
				for(k=0;k<nz;k++)
				{
					/*Determine the cell numbers*/
					mt->ele[count] = (nyz2 *(grid.nborder + i)) + (nz2 * (grid.nborder + j)) + (grid.nborder + k);

					/*Fill the derivatives*/
					/*Do NOT fill cells of the air*/
					if(b_index[j*nz + k] > 0.0)
					{
						tmp_deriv[count][0] =	d_Zr[j*nz + k]/((double)ny);				/*Real part*/
						tmp_deriv[count][1] =   d_Zi[j*nz + k]/((double)ny);			/*Imaginary part*/
					}

					count++;

					if(count > mt->ncell)
					{
						printf("The number of cells (%d) in the mt forward structure is specified too small !!!\n", mt->ncell);
						exit(0);
					}
				}
	}


	/*Write the derivatives into a temporary file*/
	if(index_mode == 1) /*TE-mode*/
		sprintf(fname,"tmp_mt_TE_deriv%d_freq%d.dat",calc_2D_mt.freq_data[j_freq][k_stat],freq_index);
	else				/*TM-mode*/
		sprintf(fname,"tmp_mt_TM_deriv%d_freq%d.dat",calc_2D_mt.freq_data[j_freq][k_stat],freq_index);

	ouf = fopen(fname,"wb");

	for(i=0;i<mt->ncell;i++)
		for(j=0;j<2;j++)
		{
			help_deriv = tmp_deriv[i][j];
			fwrite(&(help_deriv),sizeof(double),1,ouf);					/*Write out the derivatives*/
		}


	fclose(ouf);

	for(i=0;i<mt->ncell;i++)
		free(tmp_deriv[i]);

	free(tmp_deriv);
	
	return(1);
}

/*-----------------------------------------------*/
/*Calculate the resistivity and the thickness of the water column*/
/*Parameter:    waterdepths = thickness of the water column*/
/*                res_water = resistivity of the water in ohmm*/
/*               res_slice  = resistivities in the slice*/
/*				 calc_2D_mt = 2D MT structure*/
/*               x_int_stat_index, z_int_stat_index := location index of the considered station in x and z-direction*/
/*               nz			= Nr of cells in z-direction*/
/*               j_freq     = frequency index*/
/*               *index_no_water = Index that indicates if a water layer exists (=1) or not (=0)*/
/*			 flag_ray_based = Specify, if the average resistivity within a layer (= 0) or only the resistivities exactly above the stations (=1) are used*/
/*Remark: The program considers the average resistivity of every row above the stations; if the resistivity changes in a row the program brings out a warning*/
/*If the average resistvity of a row is equal or higher than the resistivity of the air, the iteration breaks and the thickness and resistivity of the water layer */
/*will be written out!!*/


int DetWaterColumn(double *waterdepths, double *res_water, double res_air,  long z_int_stat_index, long x_int_stat_index, double *res_slice, CALC_2D_MT calc_2D_mt, double first_cell_dz, long nz, long j_freq, int *index_no_water, int flag_ray_based)
{

	long i,m;
	double res_row, tmp_res_row;
	double all_res; /*Average resistivity of the complete water layer*/
	double *dz, sum_dz, all_sum_dz, sum_dz2;

	all_sum_dz = 0.0;
	all_res = 0.0;
	sum_dz2 = 0.0;

	dz = (double *)memory(NULL,(z_int_stat_index+1), sizeof(double),"DetWaterColumn");

	dz[0] = first_cell_dz * calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + z_int_stat_index];

	for(i=1;i<z_int_stat_index+1;i++)
		dz[i] = calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + (z_int_stat_index - i)];

	/*Loop over all cells in z-direction*/
	for(i=0;i<z_int_stat_index+1;i++)
	{

		res_row = 0.0;

		if(flag_ray_based == 0)
		{
			sum_dz =  0.0;

			/*Loop over all cells in x-direction*/
			for(m=calc_2D_mt.nx_left_border[j_freq];m<(calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq]);m++)
			{
				res_row = res_row + (res_slice[m*nz + (z_int_stat_index - i)]*dz[i]);
				sum_dz = dz[i] + sum_dz;
			}

			/*Determine the average resistivity in the row*/
			tmp_res_row = (res_row/sum_dz);

		}
		else
		{
			res_row = res_slice[x_int_stat_index*nz + (z_int_stat_index - i)]*dz[i];
			tmp_res_row = res_row;
			sum_dz = dz[i];
		}


		/*Identify if the resistivity belongs to the air layer*/
		if(tmp_res_row < res_air)
		{
			all_sum_dz = all_sum_dz + sum_dz;
			all_res = all_res + res_row;
			sum_dz2 = sum_dz2 + dz[i];
		}
		else
		{
			/*The layer was found that corresponds to the air*/
			goto found_air;
		}

	}

	found_air:;

	/*NO water layer exists*/
	if(i==0)
	{
		(*index_no_water) = 0;
	}
	/*There is a boundary between the water and the air layer above the stations*/
	else
	{
		/*Determine the resistivity of the water*/
		(*res_water) = (all_res/all_sum_dz);
		/*Determine the depths of the water*/
		(*waterdepths) = sum_dz2;
		/*Specify that there is a "water" layer between the stations and the air*/
		(*index_no_water) = 1;
	}


	free(dz);

	return(0);
}


/*-----------------------------------------------*/
/*Determine the average resistivity in an depth region surrounded by the z-axis from the stations and the cells*/
/*Parameter:   calc_2D_mt := 2D MT structure*/
/*            res_slice  := resistivity slice*/
/*           average_res := average resistivity for each layer in z-direction (will be calculated in this routine)*/
/*					  dz := cell thicknesses for each layer in z-direction (will be calculated in this routine)*/
/*         first_cell_dz := relative z-position of the station in the first cell*/
/*           z_int_stat_index := Depths for the station indices*/
/*             iz        := Index of the cell in the z-direction*/
/*            nz         := Nr of cells in z-direction */
/*           j_freq      := Frequency index*/
/*        *station_res   := resisistivity directly below (ONE cell below) the station (will be calculated in this routine)*/
/*			    x_index  := Index of the pertubated cell in x-direction*/
/*         x_index_stat  := Index of the cell in x-direction where the station is placed*/
/*         kind_of_calc  := Specify, if the average resistivity of the complete layer (==0) or the resistivity along a ray will be calculated (==1)*/
/*Remark: weighting of the cells by means of the thickness in z-direction*/

int DetMeanResWithinLayer(CALC_2D_MT calc_2D_mt,double *res_slice, double *average_res, double *dz, long z_int_stat_index, double first_cell_dz, long iz, long nz, long j_freq, double *station_res, double x_index, double x_index_stat, int kind_of_calc)
{

	#define EPS 1e-5 

	long k,m;
	long nr_cells, x1_int, x2_int;
	long x_index_stat_int;
	double sum_x, sum_z, sum_dx, sum_dz, x1, x2;
	double factor, sum_factor;

	/*The perturbated cells are located below the stations*/
	if(z_int_stat_index <= iz)
	{
		dz[0] = (1.0 - first_cell_dz)*calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + z_int_stat_index];

		/*station is in the perturbated cell*/
		if(iz == z_int_stat_index)
			dz[0] = fabs((0.5 - first_cell_dz)*calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + z_int_stat_index]);

		if(dz[0] == 0.0)
			dz[0] = EPS;

		for(k=1;k<(iz - z_int_stat_index);k++)
			dz[k]= calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + z_int_stat_index + k];

		dz[iz - z_int_stat_index] = (calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + z_int_stat_index + k])/2.0;

		/****************************************/
		/*Determine the average resistivity of the layer*/
		if(kind_of_calc == 0)
		{

			for(k=z_int_stat_index;k<iz+1;k++)
			{
				nr_cells = 0;

				for(m=calc_2D_mt.nx_left_border[j_freq];m<(calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq]);m++)
				{
					average_res[k - z_int_stat_index] = res_slice[m*nz + k] + average_res[k - z_int_stat_index];
					nr_cells++;
				}

				average_res[k - z_int_stat_index] = average_res[k - z_int_stat_index]/nr_cells;
			}
		}
		/****************************************/
		/*Determine the resistivity along a ray of station and perturbated cell*/
		else
		{
			 /*horizontal (x-)distance between the cells (in number of cell indeces)*/
			 sum_x = x_index - x_index_stat;

			 sum_z = 0.0;
			 for(k=0;k<(iz + 1 - z_int_stat_index);k++)
				 sum_z = dz[k] + sum_z;
			
			 sum_dz = 0.0;
			 x1 = x_index_stat + 0.5;

			 /*Determine the intersections of the "rays" with the edges of the cells*/
			 for(k=z_int_stat_index;k<(iz+1);k++)
			 {
				 sum_dz = dz[k - z_int_stat_index] + sum_dz;

				 sum_dx = (sum_dz*sum_x)/sum_z;
				 x2 = sum_dx + x_index_stat + 0.5;
				 
				 x1_int = (long) x1;
				 x2_int = (long) x2;

				 if(x1_int == x2_int)
					 average_res[k - z_int_stat_index] = res_slice[x1_int*nz + k];
				 else if(x1_int < x2_int)
				 {

					 factor = (1.0 + x1_int) - x1;
					 average_res[k - z_int_stat_index] = res_slice[x1_int*nz + k]*factor;

					 sum_factor = factor;

					 for(m=(x1_int+1);m<(x2_int-1);m++)
					 {
						 average_res[k - z_int_stat_index] = average_res[k - z_int_stat_index] + res_slice[m*nz + k];
						 sum_factor = sum_factor + 1.0;
					 }

					 factor = x2 - x2_int;
					 average_res[k - z_int_stat_index] = average_res[k - z_int_stat_index] + res_slice[x2_int*nz + k]*factor;
					 sum_factor = sum_factor + factor;

					 average_res[k - z_int_stat_index] = average_res[k - z_int_stat_index]/sum_factor;

				 }
				 else
				 {
					 factor = (1.0 + x2_int) - x2;
					 average_res[k - z_int_stat_index] = res_slice[x2_int*nz + k]*factor;

					 sum_factor = factor;

					 for(m=(x2_int+1);m<(x1_int-1);m++)
					 {
						 average_res[k - z_int_stat_index] = average_res[k - z_int_stat_index] + res_slice[m*nz + k];
						 sum_factor = sum_factor + 1.0;
					 }

					 factor = x1 - x1_int;		 
					 average_res[k - z_int_stat_index] = average_res[k - z_int_stat_index] + res_slice[x1_int*nz + k]*factor;
					 sum_factor = sum_factor + factor;

					 average_res[k - z_int_stat_index] = average_res[k - z_int_stat_index]/sum_factor;
				 }

				 x1 = x2;
			 }
		}
	}


	/****************************************/
	/****************************************/
	/*The perturbated cells are located above the stations*/
	else
	{

		/*Center of the first cell has to be used*/
		dz[0] = (calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + iz])/2.0;

		for(k=1;k<(z_int_stat_index-iz);k++)
			dz[k] = calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + iz + k];

		dz[z_int_stat_index-iz] = first_cell_dz * calc_2D_mt.hz[j_freq][calc_2D_mt.nz_atm[j_freq] + calc_2D_mt.nz_ionos[j_freq] + z_int_stat_index];

		if(dz[z_int_stat_index-iz] == 0.0)
			dz[z_int_stat_index-iz] = EPS;

		/****************************************/
		/*Determine the average resistivity of the layer*/
		if(kind_of_calc == 0)
		{

			for(k=iz;k<z_int_stat_index+1;k++)
			{
				nr_cells = 0;

				for(m=calc_2D_mt.nx_left_border[j_freq];m<(calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq]);m++)
				{
					average_res[k-iz] = res_slice[m*nz + k] + average_res[k-iz];
					nr_cells++;
				}
			
				average_res[k-iz] = average_res[k-iz]/nr_cells;
			}
		}
		/****************************************/
		/*Determine the resistivity along a ray connecting station and perturbated cell*/
		else
		{
			
			 /*horizontal (x-)distance between the cells (in number of cell indeces)*/
			 sum_x = x_index_stat - x_index;

			 sum_z = 0.0;
			 for(k=0;k<(z_int_stat_index - iz + 1);k++)
				 sum_z = dz[k] + sum_z;
			
			 sum_dz = 0.0;
			 x1 = x_index + 0.5;

			 /*Determine the intersections of the "rays" with the edges of the cells*/
			 for(k=iz;k<(z_int_stat_index + 1);k++)
			 {
				 sum_dz = dz[k-iz] + sum_dz;

				 sum_dx = (sum_dz*sum_x)/sum_z;
				 x2 = sum_dx + x_index + 0.5;
				 
				 x1_int = (long) x1;
				 x2_int = (long) x2;

				 if(x1_int == x2_int)
					 average_res[k - iz] = res_slice[x1_int*nz + k];
				 else if(x1_int < x2_int)
				 {

					 factor = (1.0 + x1_int) - x1;
					 average_res[k - iz] = res_slice[x1_int*nz + k]*factor;

					 sum_factor = factor;

					 for(m=(x1_int+1);m<(x2_int-1);m++)
					 {
						 average_res[k - iz] = average_res[k - iz] + res_slice[m*nz + k];
						 sum_factor = sum_factor + 1.0;
					 }

					 factor = x2 - x2_int;
					 average_res[k - iz] = average_res[k - iz] + res_slice[x2_int*nz + k]*factor;
					 sum_factor = sum_factor + factor;

					 average_res[k - iz] = average_res[k - iz]/sum_factor;

				 }
				 else
				 {
					 factor = (1.0 + x2_int) - x2;
					 average_res[k - iz] = res_slice[x2_int*nz + k]*factor;

					 sum_factor = factor;

					 for(m=(x2_int+1);m<(x1_int-1);m++)
					 {
						 average_res[k - iz] = average_res[k - iz] + res_slice[m*nz + k];
						 sum_factor = sum_factor + 1.0;
					 }

					 factor = x1 - x1_int;		 
					 average_res[k - iz] = average_res[k - iz] + res_slice[x1_int*nz + k]*factor;
					 sum_factor = sum_factor + factor;

					 average_res[k - iz] = average_res[k - iz]/sum_factor;
				 }

				 x1 = x2;
			 }
		}

	}

	/**********************************************/
	/*Determine the resistivity directly below the stations (Remark: One cell BELOW the station to make sure that the cell has the resistivity of the subsurface and not of the water)*/
	/****************************************/
	/*Determine the average resistivity of the layer directly below the station*/
	if(kind_of_calc == 0)
	{
		nr_cells = 0;
		(*station_res) = 0.0;

		for(m=calc_2D_mt.nx_left_border[j_freq];m<(calc_2D_mt.nx[j_freq] - calc_2D_mt.nx_right_border[j_freq]);m++)
		{
			if((z_int_stat_index+1) != nz)
				(*station_res) = res_slice[m*nz + z_int_stat_index+1] + (*station_res);
			else
				(*station_res) = res_slice[m*nz + z_int_stat_index] + (*station_res);

			nr_cells++;
		}

		(*station_res) = (*station_res)/nr_cells;

	}
	/****************************************/
	/*Use the resistivity directly below the station*/
	else
	{
		x_index_stat_int = (long)x_index_stat;

		if((x_index_stat - x_index_stat_int) != 0.5)
		{
			if((x_index_stat - x_index_stat_int) > 0.5)
				x_index_stat_int = x_index_stat_int + 1; 

			if((z_int_stat_index+1) != nz)
				(*station_res) = res_slice[x_index_stat_int*nz + z_int_stat_index+1];
			else
				(*station_res) = res_slice[x_index_stat_int*nz + z_int_stat_index];
		}
		else
		{
			if((z_int_stat_index+1) != nz)
				(*station_res) = 0.5*(res_slice[x_index_stat_int*nz + z_int_stat_index+1] + res_slice[(x_index_stat_int + 1)*nz + z_int_stat_index+1]);
			else
				(*station_res) = 0.5*(res_slice[x_index_stat_int*nz + z_int_stat_index] + res_slice[(x_index_stat_int + 1)*nz + z_int_stat_index]);
		}

	}


	return(0);

	#undef EPS 
}

/*-------------------------------------------------------------*/
/*Calculate the expressions for the Bx and Ey components*/
/*that are used for the ifft (formulas 4.30 and 4.31) in Marion Jegens Thesis*/
/*Parameters: nr_of_samples := Nr of samples in the fft*/
/*                        *p := Real horizontal wavemumber in x-direction*/
/*                       *dz := vertical distance between the station and the perturbated cell*/
/*            **fft_input_Ey := Entries of the Ey components in the wavenumber domain (even samples == real part; odd samples == imaginary part)*/
/*            **fft_input_Bx := Entries of the Bx components in the wavenumber domain (even samples == real part; odd samples == imaginary part)*/
/*                 *int_res  :=  Average resistivity (in ohmm) of each layer*/
/*						*hz := Average thickness of each layer (in m)*/
/*                      nhz := Number of cells in z direction between the station and the perturbated cell*/
/*                     freq := frequency in Hz*/
/*              waterdepths := Water depths above the stations*/
/*         res_water_layer  := resistivity of the water layer*/
/*         index_water_layer:= specify, if a "water layer" (with resistivities larger than the one of air) exists (yes==1, no==0)*/
/*               sample_int := Sampling interval*/
/*				res_station := resistivity at the stations*/
/*				index_imped := determine, if the complex wavenumbers are determined by the 1-D impedances (==1) or not (==0)*/

#define E_KONST 8.8543E-12

int CalcBx_Ey(long nr_of_samples, double *p, double *fft_input_Bx, double *fft_input_Ey, double *int_res, double *hz, long nhz, double freq, double waterdepths, double res_water_layer, int index_water_layer, double sample_int, double res_station, int index_imped)
{

	long i_sample,i;
	double hz_all, ang_freq;
	double tmp1, tmp2, tmp3, tmp4, tmp5;
	double z_real_part, z_imag_part;
	double R_TE_r, R_TE_i; /*reflection coefficient of the air-water layer*/
	dcomplex complex_p, complex_p1, theta, tmp_cplx;
	dcomplex theta_1, theta_0, theta_air; /*horizontal spatial wavenumbers at the ocean bottom, in the "water" and in the "air"*/
	dcomplex R_TE;	/*Complex reflection coefficient of the air water layer*/
	dcomplex tmp_complex_p;
	dcomplex layered_z; /*layered impedance*/

	hz_all = 0.0;

	for(i=0;i<nhz;i++)
		hz_all = hz_all + hz[i];

	ang_freq = 2*PI*freq;

	tmp1 = 4*PI*PI;
	tmp2 = 2*PI*M_KONST*freq;
	tmp3 = (2*PI*freq)*(2*PI*freq)*M_KONST*E_KONST;
	tmp4 = 2*nr_of_samples*sample_int;

		/*Loop over all samples*/
	for(i_sample = 0;i_sample<(nr_of_samples+1); i_sample++)
	{
		tmp5 = p[i_sample]*p[i_sample]*tmp1;

		/*****************************************************************/
		/*Calculate the real part of the complex wavenumber k= sqrt(i*w*nu*sigma + p^2) with p = sample_nr/(dx*Nr_of_samples)*/
		/*Eq.4.32 in M.Jegen PhD-thesis*/

		complex_p.r = 0.0;
		complex_p.i = 0.0;

		for(i=0;i<nhz;i++)
		{
			tmp_cplx = Complex(tmp5, tmp2/int_res[i]);
			tmp_cplx = Csqrt(tmp_cplx);

			/*Calculate exp(-h*phi) Eq.4.30 in M.Jegen PhD-thesis*/
			tmp_complex_p = RCmul((-1)*hz[i],tmp_cplx);
			complex_p = Cadd(complex_p, tmp_complex_p);
		}

		theta.r = 0.0;
		theta.i = 0.0;

		/*Using an average complex wavenumber*/
		if(index_imped == 0)
			theta = RCmul((-1/hz_all),complex_p);
		/*Determine the average complex wavenumber from the 1-D layered impedance*/
		else
		{
			z_real_part = 0.0;
			z_imag_part = 0.0;

			MT_1D_CALC(ang_freq, &z_real_part, &z_imag_part, hz, nhz, int_res, p[i_sample]);
			
			/*theta = -iw*nu/Z (see Marions script)*/
			theta = Complex(0,-tmp2);
			layered_z = Complex(z_real_part,z_imag_part);
			theta = Cdiv(theta,layered_z);
		}

		if(complex_p.r < -300)
		{
			complex_p.r = 0.0;
			complex_p.i = 0.0;
		}
		else
			complex_p = Cexp(complex_p);

		/*****************************************************************/
		/*Add the reflection coefficients*/
		 /*Complex wavenumber in the air (see below 4.34)*/
		theta_air = Complex(tmp5 + tmp3, 0);
		theta_air = Csqrt(theta_air);

		 /*Complex wavenumber in the "water"*/
		/*Like eq. 4.32*/
		theta_0 = Complex(tmp5, tmp2/res_water_layer);
		theta_0 = Csqrt(theta_0);

		 /*Complex wavenumber in the "ocean bottom"*/
		/*Like eq. 4.32*/
		theta_1 = Complex(tmp5, tmp2/res_station);
		theta_1 = Csqrt(theta_1);
        
		/*Calculate the reflection coefficient*/
		R_TE_r = 0.0; 
		R_TE_i = 0.0;

		Reflex_TE(waterdepths, (double)theta.r, (double)theta.i, (double)theta_1.r, (double)theta_1.i, (double)theta_air.r, (double)theta_air.i, &R_TE_r, &R_TE_i, index_water_layer);

		/*Add the terms from the reflection coefficient to the function for the fourier transform*/
		R_TE = Complex(R_TE_r, R_TE_i);
		tmp_cplx.r = 1.0 + R_TE.r;
		tmp_cplx.i = R_TE.i;

		complex_p1 =Cmul(complex_p,tmp_cplx);
		complex_p1 = RCmul(-1.0, complex_p1);
		/*****************************************************************/
		/*Prepare the input for the ifft of Bx*/
		if(i_sample == 0)
		{
			/*Real part*/
			fft_input_Bx[0] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_Bx[1] = (double)complex_p1.i/tmp4;
		}
		else
		{
			/*Real part*/
			fft_input_Bx[2*i_sample] = (double)complex_p1.r/tmp4;
			fft_input_Bx[4*nr_of_samples - 2*i_sample] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_Bx[2*i_sample + 1] = (double)complex_p1.i/tmp4;
			fft_input_Bx[4*nr_of_samples - 2*i_sample + 1] = (double)complex_p1.i/tmp4;
		}

		/*****************************************************************/
		/*Calculate (-iw/phi)exp(-h*phi)*/
		tmp_cplx = Complex(0,-1);
		tmp_cplx = RCmul(2*PI*freq,tmp_cplx);
		tmp_cplx = Cdiv(tmp_cplx, theta);
		complex_p = Cmul(tmp_cplx, complex_p);

		/*Add the terms from the reflection coefficient to the function for the fourier transform*/
		tmp_cplx.r = 1.0 - R_TE.r;
		tmp_cplx.i = (-1.0)*R_TE.i;
		complex_p1 = Cmul(tmp_cplx, complex_p);
		/*****************************************************************/
		/*Prepare the input for the ifft of Ey*/
		if(i_sample == 0)
		{
			fft_input_Ey[0] = (double)complex_p1.r/tmp4;
			fft_input_Ey[1] = (double)complex_p1.i/tmp4;
		}
		else
		{
			/*Real part*/
			fft_input_Ey[2*i_sample] = (double)complex_p1.r/tmp4;
			fft_input_Ey[4*nr_of_samples - 2*i_sample] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_Ey[2*i_sample + 1] = (double)complex_p1.i/tmp4;
			fft_input_Ey[4*nr_of_samples - 2*i_sample + 1] = (double)complex_p1.i/tmp4;
		}
		/*****************************************************************/
	}
	
	return(0);
}

/*-------------------------------------------------------------*/
/*Calculate the expressions for the By and Ex components for a current perturbation in x-direction*/
/*that are used for the ifft (formulas 4.42 and 4.43) in Marion Jegens Thesis*/
/*Parameters: nr_of_samples := Nr of samples in the fft for the By component*/
/*			  nr_of_samples_Ex :=  Nr of samples in the fft for the Ex component (it is different from the one in the By mode, because a finer sampling is required for the Ex component)*/
/*                        *p := Real horizontal wavemumber in x-direction (By component)*/
/*					  *p_Ex  := Real horizontal wavemumber in x-direction (Ex component)*/
/*                        dz := vertical distance between the station and the perturbated cell*/
/*            **fft_input_Ey := Entries of the Ey components in the wavenumber domain (even samples == real part; odd samples == imaginary part)*/
/*            **fft_input_Bx := Entries of the Bx components in the wavenumber domain (even samples == real part; odd samples == imaginary part)*/
/*             *average_res  :=  Average resistivity (in ohmm) of each layer*/
/*						*hz	:= Thickness of each layer in m*/
/*                      nhz := Number of cells in z direction between the station and the perturbated cell*/
/*                     freq := frequency in Hz*/
/*              waterdepths := Water depths above the stations*/
/*         res_water_layer  := resistivity of the water layer*/
/*                         i:= Index nr of the considered layer*/
/*               sample_int := Sampling interval for the By component*/
/*			  sample_int_Ex := Sampling interval for the Ex component*/
/*				res_station := resistivity at the stations*/
/*				index_imped := determine, if the complex wavenumbers are determined by the 1-D impedances (==1) or not (==0)*/
/*										*/
/* Remark: For the Ex component are more samples used, because the Ex is a discontinious function in the x-domain*/

int CalcBy_Ex_Ix(long nr_of_samples, long nr_of_samples_Ex, double *p, double *p_Ex, double *fft_input_By, double *fft_input_Ex, double *int_res, double *hz, long nhz, double freq, double waterdepths, double res_water_layer, double sample_int, double sample_int_Ex, double res_station, int index_imped)
{

	long i_sample,i;
	double average_res, hz_all, ang_freq, z_real_part, z_imag_part;
	double tmp1, tmp2, tmp3, tmp4, tmp5;
	double R_TM_r, R_TM_i; /*reflection coefficient of the air-water layer*/
	dcomplex theta_0, theta_air; /*horizontal spatial wavenumbers in the "water" and in the "air"*/
	dcomplex  theta, complex_p, complex_p1, tmp_cplx;
	dcomplex R_TM; /*Complex reflection coefficient*/
	dcomplex tmp_complex_p, layered_z; /*Layered impedances*/

	#define MINIMUM_H 20 /*Minimum thickness of the most upper cell in m to avoid instabilities of the fourier transform of the E-field*/ 

	/*Calculate the average resistivity*/
	average_res = 0.0;
	hz_all = 0.0;

	ang_freq = 2*PI*freq;
	
	for(i=0;i<nhz;i++)
	{
		average_res = average_res + int_res[i]*hz[i];
		hz_all = hz_all + hz[i];
	}

	average_res = average_res/hz_all;

	/*Since the E-field is not continous at h=0, h, smaller a specific threshold, is set equal to the threshold:*/
	//EVENTUELL ABAENDERN !!!!
	if(hz_all < MINIMUM_H)
	{
		for(i=0;i<nhz;i++)
			hz[i] = hz[i]*(MINIMUM_H/hz_all);

		hz_all = MINIMUM_H;
	}


	tmp1 = 4*PI*PI;
	tmp2 = 2*PI*M_KONST*freq;
	tmp3 = (2*PI*freq)*(2*PI*freq)*M_KONST*E_KONST;
	tmp4 = 2*nr_of_samples*sample_int;

	/*Loop over all samples (By-component)*/
	for(i_sample = 0;i_sample<(nr_of_samples+1); i_sample++)
	{

		tmp5 = p[i_sample]*p[i_sample]*tmp1;

		/*****************************************************************/
		/*Calculate the real part of the complex wavenumber k= sqrt(i*w*nu*sigma + p^2) with p = sample_nr/(dx*Nr_of_samples)*/
		/*Eq.4.32 in M.Jegen PhD-thesis*/

		complex_p.r = 0.0;
		complex_p.i = 0.0;

		for(i=0;i<nhz;i++)
		{
			tmp_cplx = Complex(tmp5, tmp2/int_res[i]);
			tmp_cplx = Csqrt(tmp_cplx);
			/*Calculate exp(-h*phi) Eq.4.42 in M.Jegen PhD-thesis*/
			tmp_complex_p = RCmul((-1)*hz[i],tmp_cplx);
			complex_p = Cadd(complex_p,tmp_complex_p);
		}

		theta.r = 0.0;
		theta.i = 0.0;

		/*Using an average complex wavenumber*/
		if(index_imped == 0)
			theta = RCmul((-1/hz_all),complex_p);
		/*Determine the average complex wavenumber from the 1-D layered impedance*/
		else
		{
			z_real_part = 0.0;
			z_imag_part = 0.0;

			MT_1D_CALC(ang_freq, &z_real_part, &z_imag_part, hz, nhz, int_res, p[i_sample]);
			
			/*theta = -iw*nu/Z (see Marions script)*/
			theta = Complex(0,-tmp2);
			layered_z = Complex(z_real_part,z_imag_part);
			theta = Cdiv(theta,layered_z);
		}

		if(complex_p.r < -300)
		{
			complex_p.r = 0.0;
			complex_p.i = 0.0;
		}
		else
			complex_p = Cexp(complex_p);
		/*****************************************************************/
		/*Add the reflection coefficients*/

		 /*Complex wavenumber in the air (see below 4.34)*/
		theta_air = Complex(tmp5 + tmp3, 0);
		theta_air = Csqrt(theta_air);

		 /*Complex wavenumber in the "water"*/
		/*Like eq. 4.32*/
		theta_0 = Complex(tmp5, tmp2/res_water_layer);
		theta_0 = Csqrt(theta_0);
        
		/*Calculate the reflection coefficient*/
		R_TM_r = 0.0; 
		R_TM_i = 0.0;

		Reflex_TM(waterdepths, (1/res_station), (1/res_water_layer), (double)theta.r, (double)theta.i, (double)theta_0.r, (double)theta_0.i, &R_TM_r, &R_TM_i);

		/*Add the terms from the reflection coefficient to the function for the fourier transform*/
		R_TM = Complex(R_TM_r, R_TM_i);
		tmp_cplx.r = 1.0 + R_TM.r;
		tmp_cplx.i = R_TM.i;

		complex_p1 =Cmul(complex_p,tmp_cplx);
		/*****************************************************************/

		/*Prepare the input for the ifft of Bx*/
		
		if(i_sample == 0)
		{
			/*Real part*/
			fft_input_By[0] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_By[1] = (double)complex_p1.i/tmp4;
		}
		else
		{
			/*Real part*/
			fft_input_By[2*i_sample] = (double)complex_p1.r/tmp4;
			fft_input_By[4*nr_of_samples - 2*i_sample] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_By[2*i_sample + 1] = (double)complex_p1.i/tmp4;
			fft_input_By[4*nr_of_samples - 2*i_sample +1] = (double)complex_p1.i/tmp4;
		}
		/*****************************************************************/
	}

	tmp4 = 2*nr_of_samples_Ex*sample_int_Ex;

	/*Loop over all samples (Ex-component)*/
	for(i_sample = 0;i_sample<(nr_of_samples_Ex+1); i_sample++)
	{
		tmp5 = p_Ex[i_sample]*p_Ex[i_sample]*tmp1;

		/*****************************************************************/
		/*Calculate the real part of the complex wavenumber k= sqrt(i*w*nu*sigma + p^2) with p = sample_nr/(dx*Nr_of_samples)*/
		/*Eq.4.32 in M.Jegen PhD-thesis*/

		complex_p.r = 0.0;
		complex_p.i = 0.0;

		for(i=0;i<nhz;i++)
		{
			tmp_cplx = Complex(tmp5, tmp2/int_res[i]);
			tmp_cplx = Csqrt(tmp_cplx);
			/*Calculate exp(-h*phi) Eq.4.42 in M.Jegen PhD-thesis*/
			tmp_complex_p = RCmul((-1)*hz[i],tmp_cplx);
			complex_p = Cadd(complex_p,tmp_complex_p);
		}

		theta.r = 0.0;
		theta.i = 0.0;

		theta = RCmul((-1/hz_all),complex_p);

		if(complex_p.r < -300)
		{
			complex_p.r = 0.0;
			complex_p.i = 0.0;
		}
		else
			complex_p = Cexp(complex_p);
		/*****************************************************************/
		/*Add the reflection coefficients*/

		 /*Complex wavenumber in the air (see below 4.34)*/
		theta_air = Complex(tmp5 + tmp3, 0);
		theta_air = Csqrt(theta_air);

		 /*Complex wavenumber in the "water"*/
		/*Like eq. 4.32*/
		theta_0 = Complex(tmp5, tmp2/res_water_layer);
		theta_0 = Csqrt(theta_0);
        
		/*Calculate the reflection coefficient*/
		R_TM_r = 0.0; 
		R_TM_i = 0.0;

		Reflex_TM(waterdepths, (1/res_station), (1/res_water_layer), (double)theta.r, (double)theta.i, (double)theta_0.r, (double)theta_0.i, &R_TM_r, &R_TM_i);
		/*****************************************************************/

		/*Calculate (-theta/(mu*sigma))exp(-h*theta)*/
		tmp_cplx = Complex(((-1)*(average_res))/(M_KONST),0);
		tmp_cplx = Cmul(tmp_cplx, theta);
		complex_p = Cmul(tmp_cplx,complex_p);

		/*Add the terms from the reflection coefficient to the function for the fourier transform*/
		R_TM = Complex(R_TM_r, R_TM_i);
		tmp_cplx.r = 1.0 - R_TM.r;
		tmp_cplx.i = (-1.0)*R_TM.i;
		complex_p1 = Cmul(tmp_cplx, complex_p);
		/*****************************************************************/
		/*Prepare the input for the ifft of Ex*/
		
		if(i_sample == 0)
		{
			/*Real part*/
			fft_input_Ex[0] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_Ex[1] = (double)complex_p1.i/tmp4;
		}
		else
		{
			/*Real part*/
			fft_input_Ex[2*i_sample] = (double)complex_p1.r/tmp4;
			fft_input_Ex[4*nr_of_samples_Ex - 2*i_sample] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_Ex[2*i_sample + 1] = (double)complex_p1.i/tmp4;
			fft_input_Ex[4*nr_of_samples_Ex - 2*i_sample + 1] = (double)complex_p1.i/tmp4;
		}
		/*****************************************************************/
	}

	return(0);
}


/*-------------------------------------------------------------*/
/*Calculate the expressions for the By and Ex components for a current perturbation in z-direction*/
/*Parameters: nr_of_samples := Nr of samples in the fft for the By component*/
/*			  nr_of_samples_Ex :=  Nr of samples in the fft for the Ex component (it is different from the one in the By mode, because a finer sampling is required for the Ex component)*/
/*                        *p := Real horizontal wavemumber in x-direction (By component)*/
/*					  *p_Ex  := Real horizontal wavemumber in x-direction (Ex component)*/
/*                        dz := vertical distance between the station and the perturbated cell*/
/*            **fft_input_Ey := Entries of the Ey components in the wavenumber domain (even samples == real part; odd samples == imaginary part)*/
/*            **fft_input_Bx := Entries of the Bx components in the wavenumber domain (even samples == real part; odd samples == imaginary part)*/
/*             *average_res  :=  Average resistivity (in ohmm) of each layer*/
/*						*hz	:= Thickness of each layer in m*/
/*                      nhz := Number of cells in z direction between the station and the perturbated cell*/
/*                     freq := frequency in Hz*/
/*              waterdepths := Water depths above the stations*/
/*         res_water_layer  := resistivity of the water layer*/
/*                         i:= Index nr of the considered layer*/
/*               sample_int := Sampling interval for the By component*/
/*			  sample_int_Ex := Sampling interval for the Ex component*/
/*				res_station := resistivity at the stations*/
/*				index_imped := determine, if the complex wavenumbers are determined by the 1-D impedances (==1) or not (==0)*/
/*										*/
/* Remark: For the Ex component are usually more samples used, because the Ex is a discontinious function in the x-domain*/

int CalcBy_Ex_Iz(long nr_of_samples, long nr_of_samples_Ex, double *p, double *p_Ex, double *fft_input_By, double *fft_input_Ex, double *int_res, double *hz, long nhz, double freq, double waterdepths, double res_water_layer, double sample_int, double sample_int_Ex, double res_station, int index_imped)
{
	
	long i_sample,i;
	double average_res, hz_all, ang_freq, z_real_part, z_imag_part;
	double tmp1, tmp2, tmp3, tmp4, tmp5;
	double R_TM_r, R_TM_i; /*reflection coefficient of the air-water layer*/
	dcomplex theta_0, theta_air; /*horizontal spatial wavenumbers in the "water" and in the "air"*/
	dcomplex  theta, complex_p, complex_p1, tmp_cplx;
	dcomplex R_TM; /*Complex reflection coefficient*/
	dcomplex tmp_complex_p, layered_z; /*Layered impedances*/

	#define MINIMUM_H 20 /*Minimum thickness of the most upper cell in m to avoid instabilities of the fourier transform of the E-field*/ 

	/*Calculate the average resistivity*/
	average_res = 0.0;
	hz_all = 0.0;
	
	ang_freq = 2*PI*freq;

	for(i=0;i<nhz;i++)
	{
		average_res = average_res + int_res[i]*hz[i];
		hz_all = hz_all + hz[i];
	}

	average_res = average_res/hz_all;

	/*Since the E-field is not continous at h=0, h, smaller a specific threshold, is set equal to the threshold:*/
	if(hz_all < MINIMUM_H)
	{
		for(i=0;i<nhz;i++)
			hz[i] = hz[i]*(MINIMUM_H/hz_all);

		hz_all = MINIMUM_H;
	}

	tmp1 = 4*PI*PI;
	tmp2 = 2*PI*M_KONST*freq;
	tmp3 = (2*PI*freq)*(2*PI*freq)*M_KONST*E_KONST;
	tmp4 = 2*nr_of_samples*sample_int;

	/*Loop over all samples (By-component)*/
	for(i_sample = 0;i_sample<(nr_of_samples+1); i_sample++)
	{

		tmp5 = p[i_sample]*p[i_sample]*tmp1;

		/*****************************************************************/
		/*Calculate the real part of the complex wavenumber k= sqrt(i*w*nu*sigma + p^2) with p = sample_nr/(dx*Nr_of_samples)*/
		/*Eq.4.32 in M.Jegen PhD-thesis*/

		complex_p.r = 0.0;
		complex_p.i = 0.0;

		for(i=0;i<nhz;i++)
		{
			tmp_cplx = Complex(tmp5, tmp2/int_res[i]);
			tmp_cplx = Csqrt(tmp_cplx);
			/*Calculate exp(-h*phi) Eq.4.42 in M.Jegen PhD-thesis*/
			tmp_complex_p = RCmul((-1)*hz[i],tmp_cplx);
			complex_p = Cadd(complex_p,tmp_complex_p);
		}

		theta.r = 0.0;
		theta.i = 0.0;

		/*Using an average complex wavenumber*/
		if(index_imped == 0)
			theta = RCmul((-1/hz_all),complex_p);
		/*Determine the average complex wavenumber from the 1-D layered impedance*/
		else
		{
			z_real_part = 0.0;
			z_imag_part = 0.0;

			MT_1D_CALC(ang_freq, &z_real_part, &z_imag_part, hz, nhz, int_res, p[i_sample]);
			
			/*theta = -iw*nu/Z (see Marions script)*/
			theta = Complex(0,-tmp2);
			layered_z = Complex(z_real_part,z_imag_part);
			theta = Cdiv(theta,layered_z);
		}

		if(complex_p.r < -300)
		{
			complex_p.r = 0.0;
			complex_p.i = 0.0;
		}
		else
			complex_p = Cexp(complex_p);
		/*****************************************************************/
		/*Add the reflection coefficients*/

		 /*Complex wavenumber in the air (see below 4.34)*/
		theta_air = Complex(tmp5 + tmp3, 0);
		theta_air = Csqrt(theta_air);

		 /*Complex wavenumber in the "water"*/
		/*Like eq. 4.32*/
		theta_0 = Complex(tmp5, tmp2/res_water_layer);
		theta_0 = Csqrt(theta_0);
        
		/*Calculate the reflection coefficient*/
		R_TM_r = 0.0; 
		R_TM_i = 0.0;

		Reflex_TM(waterdepths, (1/res_station), (1/res_water_layer), (double)theta.r, (double)theta.i, (double)theta_0.r, (double)theta_0.i, &R_TM_r, &R_TM_i);

		/*Add the terms from the reflection coefficient to the function for the fourier transform*/
		R_TM = Complex(R_TM_r, R_TM_i);
		tmp_cplx.r = 1.0 + R_TM.r;
		tmp_cplx.i = R_TM.i;

		complex_p1 =Cmul(complex_p,tmp_cplx);

		/********************************************/
		/*Calculate ip/theta*/
		tmp_cplx = Complex(0,p[i_sample]);
		tmp_cplx = Cdiv(tmp_cplx,theta);

		complex_p1 = Cmul(tmp_cplx ,complex_p1);
		/*****************************************************************/

		/*Prepare the input for the ifft of Bx*/
		
		if(i_sample == 0)
		{
			/*Real part*/
			fft_input_By[0] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_By[1] = (double)complex_p1.i/tmp4;
		}
		else
		{
			/*Real part*/
			fft_input_By[2*i_sample] = (double)complex_p1.r/tmp4;
			fft_input_By[4*nr_of_samples - 2*i_sample] = (double)(-1)*complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_By[2*i_sample + 1] = (double)complex_p1.i/tmp4;
			fft_input_By[4*nr_of_samples - 2*i_sample +1] = (double)(-1)*complex_p1.i/tmp4;
		}
		/*****************************************************************/
	}

	tmp4 = 2*nr_of_samples_Ex*sample_int_Ex;

	/*Loop over all samples (Ex-component)*/
	for(i_sample = 0;i_sample<(nr_of_samples_Ex+1); i_sample++)
	{
		tmp5 = p_Ex[i_sample]*p_Ex[i_sample]*tmp1;

		/*****************************************************************/
		/*Calculate the real part of the complex wavenumber k= sqrt(i*w*nu*sigma + p^2) with p = sample_nr/(dx*Nr_of_samples)*/
		/*Eq.4.32 in M.Jegen PhD-thesis*/

		complex_p.r = 0.0;
		complex_p.i = 0.0;

		for(i=0;i<nhz;i++)
		{
			tmp_cplx = Complex(tmp5, tmp2/int_res[i]);
			tmp_cplx = Csqrt(tmp_cplx);
			/*Calculate exp(-h*phi) Eq.4.42 in M.Jegen PhD-thesis*/
			tmp_complex_p = RCmul((-1)*hz[i],tmp_cplx);
			complex_p = Cadd(complex_p,tmp_complex_p);
		}

		theta.r = 0.0;
		theta.i = 0.0;

		theta = RCmul((-1/hz_all),complex_p);

		if(complex_p.r < -300)
		{
			complex_p.r = 0.0;
			complex_p.i = 0.0;
		}
		else
			complex_p = Cexp(complex_p);
		/*****************************************************************/
		/*Add the reflection coefficients*/

		 /*Complex wavenumber in the air (see below 4.34)*/
		theta_air = Complex(tmp5 + tmp3, 0);
		theta_air = Csqrt(theta_air);

		 /*Complex wavenumber in the "water"*/
		/*Like eq. 4.32*/
		theta_0 = Complex(tmp5, tmp2/res_water_layer);
		theta_0 = Csqrt(theta_0);
        
		/*Calculate the reflection coefficient*/
		R_TM_r = 0.0; 
		R_TM_i = 0.0;

		Reflex_TM(waterdepths, (1/res_station), (1/res_water_layer), (double)theta.r, (double)theta.i, (double)theta_0.r, (double)theta_0.i, &R_TM_r, &R_TM_i);
		/*****************************************************************/

		/*Calculate (-ip/(mu*sigma))exp(-h*theta)*/
		tmp_cplx.i = ((-1.0)*p[i_sample]*average_res)/M_KONST;
		tmp_cplx.r = 0.0;
		complex_p = Cmul(tmp_cplx,complex_p);

		/*Add the terms from the reflection coefficient to the function for the fourier transform*/
		R_TM = Complex(R_TM_r, R_TM_i);
		tmp_cplx.r = 1.0 - R_TM.r;
		tmp_cplx.i = (-1.0)*R_TM.i;
		complex_p1 = Cmul(tmp_cplx, complex_p);
		/*****************************************************************/
		/*Prepare the input for the ifft of Ex*/
		
		if(i_sample == 0)
		{
			/*Real part*/
			fft_input_Ex[0] = (double)complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_Ex[1] = (double)complex_p1.i/tmp4;
		}
		else
		{
			/*Real part*/
			fft_input_Ex[2*i_sample] = (double)complex_p1.r/tmp4;
			fft_input_Ex[4*nr_of_samples_Ex - 2*i_sample] = (double)(-1)*complex_p1.r/tmp4;
			/*Imaginary part*/
			fft_input_Ex[2*i_sample + 1] = (double)complex_p1.i/tmp4;
			fft_input_Ex[4*nr_of_samples_Ex - 2*i_sample + 1] = (double)(-1)*complex_p1.i/tmp4;
		}
		/*****************************************************************/
	}

	return(0);
}


#undef M_KONST
#undef E_KONST
#undef MINIMUM_H 


/*-------------------------------------------------------*/
/*Calculate the reflection coefficient for the TE-mode*/
/*to consider the effects of the water layer and the air layer above*/
/*(formula 4.33 and 4.34 in Marion Jegens thesis)*/
/*Parameter:	waterdepths := waterdepths above the stations in m*/
/*				theta_r, theta_i := complex wavenumber of the underground media*/
/*				theta_0_r,theta_0_i := complex wavenumber of the watercolumn*/
/*				theta_air_r, theta_air_i := complex wavenumber of the air*/
/*				R_TE_r, R_TE_i := Real and imaginary part of the reflection coefficient*/
/*              index_water_layer := Index that specify if a "water" layer with resistivities lower than the one of air exists above the stations*/

int Reflex_TE(double waterdepths, double theta_r, double theta_i, double theta_0_r, double theta_0_i, double theta_air_r, double theta_air_i, double *R_TE_r, double *R_TE_i, int index_water_layer)
{
	dcomplex theta, theta_0, theta_air;
	dcomplex tmp_cplx1, tmp_cplx2, tmp_cplx3;
	dcomplex nominator, denominator;
	dcomplex R_TE;

    theta = Complex(theta_r, theta_i);
	theta_0 = Complex(theta_0_r, theta_0_i);
	theta_air = Complex(theta_air_r, theta_air_i);

	/*Equation 4.34*/
	/*A water layer exists*/
	if(index_water_layer != 0)
	{
		tmp_cplx1 = RCmul(2*waterdepths,theta_0);

		if(tmp_cplx1.r < -300)
		{
			tmp_cplx1.r = 0.0;
			tmp_cplx1.i = 0.0;
		}
		else if(tmp_cplx1.r > 700)
		{
			tmp_cplx1.r = 700;
			tmp_cplx1 = Cexp(tmp_cplx1);
		}
		else
			tmp_cplx1 = Cexp(tmp_cplx1);

		tmp_cplx2 = Cdiv(theta_air, theta_0);

		tmp_cplx3.r = tmp_cplx2.r - 1.0;
		tmp_cplx3.i = tmp_cplx2.i;

		tmp_cplx2.r = tmp_cplx2.r + 1.0;

		/*Calculate the nominator*/
		nominator = Cdiv(tmp_cplx3,tmp_cplx2);
		nominator = Csub(tmp_cplx1,nominator);

		/*Calculate the denominator*/
		denominator = Cdiv(tmp_cplx3,tmp_cplx2);
		denominator = Cadd(tmp_cplx1,denominator);

		R_TE = Cdiv(nominator, denominator);
	}
	/*A water layer does not exist*/
	else
	{
		R_TE = Complex(1,0);
	}

	/*Equation 4.33*/
	/*Calculate the nominator and denominator*/
	tmp_cplx1 = Cdiv(theta_0, theta);
	nominator = Csub(tmp_cplx1, R_TE);
	denominator = Cadd(tmp_cplx1, R_TE);

	R_TE = Cdiv(nominator,denominator);

	*R_TE_r = (double)R_TE.r;
	*R_TE_i = (double)R_TE.i;

	return(0);
}

/*-------------------------------------------------------*/
/*Calculate the reflection coefficient for the TM-mode*/
/*to consider the effects of the water layer and the air layer above*/
/*(formula 4.44 and 4.45 in Marion Jegens thesis)*/
/*Parameter:	waterdepths := waterdepths above the stations in m*/
/*				sigma, sigma_0 := conductivity of the layer 0 (ocean layer) and the subsurface*/
/*				theta_r, theta_i := complex wavenumber of the underground media*/
/*				theta_0_r,theta_0_i := complex wavenumber of the watercolumn*/
/*				R_TM_r, R_TM_i := Real and imaginary part of the reflection coefficient*/

int Reflex_TM(double waterdepths, double sigma, double sigma_0, double theta_r, double theta_i, double theta_0_r, double theta_0_i, double *R_TM_r, double *R_TM_i)
{
	dcomplex theta_0, theta;
	dcomplex tmp_cplx, tmp_cplx1, R_TM;
	dcomplex nominator, denominator, zaehler, nenner;

	theta = Complex(theta_r, theta_i);
	theta_0 = Complex(theta_0_r, theta_0_i);

	/*Equation 4.45*/
	tmp_cplx = RCmul(2*waterdepths,theta_0);

	if(tmp_cplx.r < -300)
	{
		tmp_cplx.r = 0.0;
		tmp_cplx.i = 0.0;
	}
	else if(tmp_cplx.r > 700)
	{
		tmp_cplx.r = 700.0;
		tmp_cplx = Cexp(tmp_cplx);
	}
	else
		tmp_cplx = Cexp(tmp_cplx);

	zaehler.r = tmp_cplx.r + 1.0;
	zaehler.i = tmp_cplx.i;

	nenner.r = tmp_cplx.r - 1.0;
	nenner.i = tmp_cplx.i;

	R_TM = Cdiv(zaehler,nenner);

	/*Equation 4.44*/
	tmp_cplx = RCmul(sigma,theta_0);
	tmp_cplx1 = RCmul(sigma_0, theta);

	tmp_cplx = Cdiv(tmp_cplx1, tmp_cplx);
	/*nominator*/
	nominator = Csub(tmp_cplx, R_TM);
	denominator = Cadd(tmp_cplx, R_TM);

	R_TM = Cdiv(nominator, denominator);

	*R_TM_r = (double)R_TM.r;
	*R_TM_i = (double)R_TM.i;

	return(0);
}