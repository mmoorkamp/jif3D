#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "inv3d.h"
#include "fkt.h"


/*-------------------------------------------------------------*/
/*Performing the gravity modeling described by: I.Trinks and R.Hobbs, 200?, Gravity modelling based on small cells, Earth Sciences*/
/*Parameter:	geo  := Geometry structure  */
/*              grid  := Grid structure */
/*              *data  := Pointer on data structure (the gravity will be (re-)determined in this routine)*/
/*				inv		:= Inversion structure */

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/

#define FACTOR 10 /*Specify the distance to the cell centers (unit: number of cell length), for which the response will be calculated exactly*/
				   /*REMARK: If you re-define this parameter, the same parameter have to be changed in the "derivative" routine, too !!! */

#define refine 3  /*3D refinement factor: in x-y plane one cell is refined by refine*refine smaller cells*/


int ForwardModGrav(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, TOPO_STRUCT topo, INV inv) /*Jin,12_03*/
{
	int		specify_direc_index; /*2D case: specify 2D in x_direction (1); y_direction (0).*/
	long   i,j,l,k,m,n,ix,iy,iz,ixiy,n1,p;
	long   corner,corner1,corner2,corner3,corner4; /*The corners index of the cells*/
	double eps;
	double org1,sum,half_h,half_h_refine;
	double xs,ys,zs,xys;										/*position of the station in m*/
	double x,x1,x2,y,y1,y2,z,xy,xy1,xy2;						/*position of the cell centers in m*/
	/*distance between the cell centers and the gravity stations in m*/
	double dx,dy,dz,dx1,dx2,dy1,dy2,dx_1,dx_2,dy_1,dy_2,dxy,dxy1,dxy2,dxyz,dz1,dz2,dxyz1,dxz,dyz,dxx1,dxx2,dyy1,dyy2;
	double zcut1,zcut2,zcut3,zcut4,aver_z1,aver_z2; /*The cutting points in 3D on the boder and corners, and the average of the cutting points.*/
    double *aver,*aver_dis; /*The average value of the cutting point in smaller cells.*/
	double grav, tmp_grav, tmp_grav1, tmp_grav2, tmp_grav3, tmp_grav4;
	/* Calculate the cutting points of the triangles with the cells*/
	double const1;

	/*double *cutting_points_x,*cutting_points_y,*cutting_points_z;  */
	double *xgrid_edges,*ygrid_edges,*zgrid_edges;
    double min_tmp,max_tmp; /*The temporary min and max of the cutting points.*/

	FILE *out;


	for(i=0;i<data->ndata_grav;i++)
			data->calc_grav[i] = -1.0;

	if(data->ndata_grav == 0)
		printf("!!WARNING!! NO gravity data exists but the gravity forward modelling is activated\n\n");

	eps = 0.0000000001;
	half_h = grid.h/2.0;
	half_h_refine = grid.h/refine/2.0;

	/*Decide if the calculation is 1D, 2D or 3D by considering the number of cells in the horizontal directions:*/

	/******************************************************/
					/*1D-calculation*/
	/******************************************************/
	/* ATTENTION!!! 1-D WAS NOT TESTED !!!*/

	if(grid.nx == 1 && grid.ny == 1)
	{
		printf("The gravity model is identified as 1-D\n");
		printf("Start with the gravity modeling\n");
		printf("-------------------------------\n\n");


		/*Loop over all stations*/
		for(i=0;i<geo.nstat_grav;i++)
		{
			/*Positions of the gravity stations:*/
			xs = (double)geo.x[(data->gravs[i]-1)];
			ys = (double)geo.y[(data->gravs[i]-1)];
			zs = (double)geo.z[(data->gravs[i]-1)];

			grav = 0.0;

			/**********************/
			for(iz=0;iz<grid.nz;iz++)
			{
				z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
				dz = z - zs;					/*z-component of distance between grid cell and station*/

				/*Calculate the gravitation from the cell at the station*/
				tmp_grav = 2*PI*G_KONST*(grid.dens[iz])*(grid.h); /*double-infinite sheet*/
				tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

				grav = grav + tmp_grav;

			}

			/*Assign the calculated density value to the corresponding measurements*/
			for(j=0;j<data->ndata_grav;j++)
			{
				if(data->gravs[i] == data->gno[j])
					data->calc_grav[j] = grav;
			}

			if(i%10 == 0)
				printf("Modeling for %d of %d stations is completed\n", i+1, geo.nstat_grav);

		}
	}
	
	/******************************************************/
					/*2D-calculation*/
	/******************************************************/

	/*calculate the product of cell size and gravitational const.: dx*dy*G_konst in m^5/(Kg*s^2)*/
    const1 = (grid.h)*(grid.h)*G_KONST;

	#define dens(x_y,z) grid.dens[(grid.nz)*(x_y) + (z)]


	if(grid.nx == 1 || grid.ny == 1)
	{
		printf("The gravity model is identified as 2-D\n");
		printf("Start with the gravity modeling\n");
		printf("-------------------------------\n\n");

		/*Determine the x-y coordinates of the edges of the grid cells*/
		if(grid.nx == 1)
		{
	        xgrid_edges = (double *)memory(NULL,(grid.ny+1),sizeof(double),"ForwardModGrav");
	        ygrid_edges = (double *)memory(NULL,(grid.ny+1),sizeof(double),"ForwardModGrav");
	        zgrid_edges = (double *)memory(NULL,(grid.ny+1),sizeof(double),"ForwardModGrav");
		}
		else
		{
	        xgrid_edges = (double *)memory(NULL,(grid.nx+1),sizeof(double),"ForwardModGrav");
	        ygrid_edges = (double *)memory(NULL,(grid.nx+1),sizeof(double),"ForwardModGrav");
	        zgrid_edges = (double *)memory(NULL,(grid.nx+1),sizeof(double),"ForwardModGrav");
		}

        if(topo.nr_of_triangles != 0)
	    {		
	       if(grid.nx == 1)
		   {
			   	for(j=0;j<(grid.ny+1);j++)
				{
					xgrid_edges[j] = grid.org[0];                       /*The x-coordinate of the midpoint of the cells*/
					ygrid_edges[j] = grid.org[1] +(j+0.5)*(grid.h); /*The y-coordinate of the midpoints of the cells*/
				}

			   specify_direc_index = 0;

			   /*Interpolation by means of triangulation/ using barycentric coordinates (see in: Contouring (1992) from D.F. Watson; Pergamon Press; pages 76-78)*/
    
			   TriangleInterpol(xgrid_edges,ygrid_edges,zgrid_edges,(grid.ny+1),topo.x,topo.y,topo.z,topo.nr_of_topo_points,topo.index,topo.nr_of_triangles);
		   }
		   else
		   {
			  
			   for(j=0;j<(grid.nx+1);j++)
			   {
				   xgrid_edges[j] = grid.org[1];                       /*The x-coordinate of the midpoint of the cells*/
				   ygrid_edges[j] = grid.org[0] +(j+0.5)*(grid.h); /*The y-coordinate of the midpoints of the cells*/
			   }

			   specify_direc_index = 1;

			   /*Interpolation by means of triangulation/ using barycentric coordinates (see in: Contouring (1992) from D.F. Watson; Pergamon Press; pages 76-78)*/

	           TriangleInterpol(xgrid_edges,ygrid_edges,zgrid_edges,(grid.nx+1),topo.y,topo.x,topo.z,topo.nr_of_topo_points,topo.index,topo.nr_of_triangles);
		   }
	    }

		/*Loop over all stations*/
		for(i=0;i<geo.nstat_grav;i++)
		{

			/*Positions of the gravity stations:*/
			xs = (double)geo.x[(data->gravs[i]-1)];
			ys = (double)geo.y[(data->gravs[i]-1)];
			zs = (double)geo.z[(data->gravs[i]-1)];

			grav = 0.0;

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

			/**********************/
			
		    /*Loop over all cell at the left (ix,iy = 0) border*/
			if((topo.nr_of_triangles != 0) && (zgrid_edges[0] != -99999.9))
			{						
			    for(iz=0;iz<grid.nz;iz++)/*Consider the topography information.*/
				{
					z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					dz = z - zs;				/*z-component of distance between grid cell and station*/
					if(dz == 0)
						dz = eps;

					/*Calculate the gravitation from the cell at the station*/
					/*left border*/
					if(dz >= 0)
					{
						if(zgrid_edges[0]>=z-half_h && zgrid_edges[0]<=z+half_h) /*Calculate the gravity of the cell cut by topography*/
						{
							/*semi-infinite horizontal sheet divided by two parts: the upper part is water.*/
						    tmp_grav = (2*G_KONST*grid.dens_air_water*(zgrid_edges[0]-(z-half_h)))*((PI/2) - atan(dxy1/(grid.org[2]+(iz-0.5)*(grid.h)+(zgrid_edges[0]-(z-half_h)/2)-zs)));
							tmp_grav = tmp_grav+(2*G_KONST*dens(0,iz)*(z+half_h-zgrid_edges[0]))*((PI/2) - atan(dxy1/(grid.org[2]+(iz)*(grid.h)+(zgrid_edges[0]-(z-half_h)/2)-zs))); 
						}
						else

						    tmp_grav = (2*G_KONST*dens(0,iz)*grid.h)*((PI/2) - atan(dxy1/dz)); /*semi-infinite horizontal sheet*/
					}
					else
					{
                        if(zgrid_edges[0]>=z-half_h && zgrid_edges[0]<=z+half_h)
						{
							/*semi-infinite horizontal sheet divided by two parts: the upper part is water.*/
						    tmp_grav = (2*G_KONST*grid.dens_air_water*(zgrid_edges[0]-(z-half_h)))*((-PI/2) - atan(dxy1/(grid.org[2]+(iz-0.5)*(grid.h)+(zgrid_edges[0]-(z-half_h)/2)-zs)));
							tmp_grav = tmp_grav+(2*G_KONST*dens(0,iz)*(z+half_h-zgrid_edges[0]))*((-PI/2) - atan(dxy1/(grid.org[2]+(iz)*(grid.h)+(zgrid_edges[0]-(z-half_h)/2)-zs))); 
						}
						else

     						tmp_grav = (2*G_KONST*dens(0,iz)*grid.h)*((-PI/2) - atan(dxy1/dz)); /*semi-infinite horizontal sheet*/
					}
					
					tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					grav = grav + tmp_grav;
				}
			}
			else
			{            /*Loop over all cell at the left (ix,iy = 0) border*/

			    for(iz=0;iz<grid.nz;iz++)/*There is no topography information*/
				{
					z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					dz = z - zs;				/*z-component of distance between grid cell and station*/
					if(dz == 0)
						dz = eps;

					/*Calculate the gravitation from the cell at the station*/
					/*left border*/
					if(dz >= 0)
						tmp_grav = (2*G_KONST*dens(0,iz)*grid.h)*((PI/2) - atan(dxy1/dz)); /*semi-infinite horizontal sheet*/
					else
						tmp_grav = (2*G_KONST*dens(0,iz)*grid.h)*((-PI/2) - atan(dxy1/dz)); /*semi-infinite horizontal sheet*/
					
					tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					grav = grav + tmp_grav;
				}
			}

		   /*Loop over all cell at (ix,iy = nx-1,ny-1) right border*/
			/*Consider the topography and the triangulations.*/
			if((topo.nr_of_triangles != 0) && (zgrid_edges[n1-1] != -99999.9))
			{						
			    for(iz=0;iz<grid.nz;iz++)/*Consider the topography information.*/
				{
					z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					dz = z - zs;				/*z-component of distance between grid cell and station*/
					if(dz == 0)
						dz = eps;

					/*right border*/
					if(dz >= 0)
					{/*If the cell is cut by topography.*/
                        if(zgrid_edges[n1-1]>=z-half_h && zgrid_edges[n1-1]<=z+half_h)
						{
							/*semi-infinite horizontal sheet divided by two parts*//*The upper part is water.*/
						    tmp_grav = (2*G_KONST*grid.dens_air_water*(zgrid_edges[n1-1]-(z-half_h)))*((PI/2) - atan(dxy2/(grid.org[2]+(iz-0.5)*(grid.h)+(zgrid_edges[n1-1]-(z-half_h)/2)-zs)));
							tmp_grav = tmp_grav+(2*G_KONST*dens(n1-1,iz)*(z+half_h-zgrid_edges[n1-1]))*((PI/2) - atan(dxy2/(grid.org[2]+(iz)*(grid.h)+(zgrid_edges[n1-1]-(z-half_h)/2)-zs))); 
						}
						else
						    tmp_grav = (2*G_KONST*dens(n1-1,iz)*grid.h)*((PI/2) - atan(dxy2/dz)); /*semi-infinite horizontal sheet*/
					}
					else
					{
                        if(zgrid_edges[n1-1]>=z-half_h && zgrid_edges[n1-1]<=z+half_h)
						{
							/*semi-infinite horizontal sheet divided by two parts*//*The upper part is water.*/
						    tmp_grav = (2*G_KONST*grid.dens_air_water*(zgrid_edges[n1-1]-(z-half_h)))*((-PI/2) - atan(dxy2/(grid.org[2]+(iz-0.5)*(grid.h)+(zgrid_edges[n1-1]-(z-half_h)/2)-zs)));
							tmp_grav = tmp_grav+(2*G_KONST*dens(n1-1,iz)*(z+half_h-zgrid_edges[n1-1]))*((-PI/2) - atan(dxy2/(grid.org[2]+(iz)*(grid.h)+(zgrid_edges[n1-1]-(z-half_h)/2)-zs))); 
						}
						else
						tmp_grav = (2*G_KONST*dens(n1-1,iz)*grid.h)*((-PI/2) - atan(dxy2/dz)); /*semi-infinite horizontal sheet*/
					}
					
					tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					grav = grav + tmp_grav;
				}
			}
			else
			{            /*Loop over all cell at the (ix,iy = nx-1,ny-1) right border*/
 
			    for(iz=0;iz<grid.nz;iz++)/*There is no topography information*/
				{
					z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					dz = z - zs;				/*z-component of distance between grid cell and station*/
					if(dz == 0)
						dz = eps;

					/*right border*/
					if(dz >= 0)
						tmp_grav = (2*G_KONST*dens(n1-1,iz)*grid.h)*((PI/2) - atan(dxy2/dz)); /*semi-infinite horizontal sheet*/
					else
						tmp_grav = (2*G_KONST*dens(n1-1,iz)*grid.h)*((-PI/2) - atan(dxy2/dz)); /*semi-infinite horizontal sheet*/
					tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					grav = grav + tmp_grav;
				}
			}

			/**********************/
			/*2.Case: Internal cells of the 2-D grid*/
			/*Loop over all remaining cells*/

			for(ixiy=1;ixiy<(n1-1);ixiy++)
			{
				xy = org1 + ixiy*(grid.h);		 /*x resp.y component of the cell centers in m*/
				dxy = xy -xys;					 /*xy-component of distance between grid cell and station*/

				/*Consider the topography and the triangulations.*/
			    if((topo.nr_of_triangles != 0) && (zgrid_edges[ixiy-1] != -99999.9) && (zgrid_edges[ixiy] != -99999.9))
				{ /*Temporary min and max of the cutting points.*/
					min_tmp = min(zgrid_edges[ixiy-1],zgrid_edges[ixiy]);
					max_tmp = max(zgrid_edges[ixiy-1],zgrid_edges[ixiy]);

					l = 0; /*The first cell is cut.*/
					k = 0; /*The last cell is cut.*/

					for(iz=0;iz<grid.nz;iz++)
					{						
						if(grid.org[2]+(iz-0.5)*(grid.h) <= min_tmp && grid.org[2]+(iz+0.5)*(grid.h) > min_tmp)
					        l = iz;
						if(grid.org[2]+(iz-0.5)*(grid.h) <= max_tmp && grid.org[2]+(iz+0.5)*(grid.h) > max_tmp)
						{
						    k = iz;
							break;
						}
					}

					/*Loop before cutting.*/
                    for(iz=0;iz<l;iz++)
				    {
					    z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					    dz = z - zs;					/*z-component of distance between grid cell and station*/
					    if(dz == 0)
						    dz = eps;
 
						dxyz = sqrt(dxy*dxy + dz*dz);

					    if(dxyz < (grid.h * FACTOR))
					    {
					     	/*exact calculation of attraction of two-dimensional rectangular blocks for cells close to the corresponding station following Heiland (1946)*/
						    tmp_grav = Response_2D_Rectangle(dens(ixiy,iz),dxy-half_h,dxy+half_h,dz-half_h,dz +half_h);
					    }
					    else
						    /*infinitely long horizontal rod*/
						    tmp_grav = (2*const1*dens(ixiy,iz))/(dz*(1 +((dxy*dxy)/(dz*dz)))); 

					    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					    grav = grav + tmp_grav;
				    }

					/*Loop the cells being cut.*/
					/*Cut just one cell.*/
  			        if(l == k)
				    {
					    z  = grid.org[2] + l*(grid.h);  /*z component of the cell centers in m*/
					    dz = z - zs;					/*z-component of distance between grid cell and station*/
					    if(dz == 0)
						    dz = eps;

					    
						
						if(zgrid_edges[ixiy-1] == min_tmp)
						{ /*The cell divided in 4 parts: two rectangulars and two triangulars.*//*The upper part is water.*/
							tmp_grav = Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy+half_h,dz-half_h,zgrid_edges[ixiy-1]-zs);
							tmp_grav = tmp_grav+Response_2D_upper_triangle(grid.dens_air_water,dxy-half_h,dxy+half_h,zgrid_edges[ixiy-1]-zs,zgrid_edges[ixiy]-zs);
			   	            tmp_grav = tmp_grav+Response_2D_lower_triangle(dens(ixiy,l),dxy-half_h,dxy+half_h,zgrid_edges[ixiy-1]-zs,zgrid_edges[ixiy]-zs);
			   	            tmp_grav = tmp_grav+Response_2D_Rectangle(dens(ixiy,l),dxy-half_h,dxy+half_h,zgrid_edges[ixiy]-zs,dz +half_h);
						}
						else
						{  /*The cell divided in 4 parts: two rectangulars and two triangulars.*//*The upper part is water.*/
							tmp_grav = Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy+half_h,dz-half_h,zgrid_edges[ixiy]-zs);
							tmp_grav = tmp_grav+Response_2D_upper_triangle1(grid.dens_air_water,dxy-half_h,dxy+half_h,zgrid_edges[ixiy]-zs,zgrid_edges[ixiy-1]-zs);
			   	            tmp_grav = tmp_grav+Response_2D_lower_triangle1(dens(ixiy,l),dxy-half_h,dxy+half_h,zgrid_edges[ixiy]-zs,zgrid_edges[ixiy-1]-zs);
			   	            tmp_grav = tmp_grav+Response_2D_Rectangle(dens(ixiy,l),dxy-half_h,dxy+half_h,zgrid_edges[ixiy-1]-zs,dz +half_h);						
						}
					    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					    grav = grav + tmp_grav;
					}						
					else 
				    {  /*Cut two cells.*/
						if(k-l == 1)
						{
				            if(zgrid_edges[ixiy-1] == min_tmp)
				            { /*Two cells are divided into 8 parts.*//*The upper part is water.*/
                                tmp_grav = Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy+half_h,grid.org[2]+(l-0.5)*(grid.h)-zs,zgrid_edges[ixiy-1]-zs);
								tmp_grav = tmp_grav+Response_2D_upper_triangle(grid.dens_air_water,dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),zgrid_edges[ixiy-1]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
								tmp_grav = tmp_grav+Response_2D_lower_triangle(dens(ixiy,l),dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),zgrid_edges[ixiy-1]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
                                tmp_grav = tmp_grav+Response_2D_Rectangle(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,zgrid_edges[ixiy-1]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
                                tmp_grav = tmp_grav+Response_2D_Rectangle(dens(ixiy,k),dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy]-zs);
                                tmp_grav = tmp_grav+Response_2D_upper_triangle(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy]-zs);
                                tmp_grav = tmp_grav+Response_2D_lower_triangle(dens(ixiy,k),dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy]-zs);
                                tmp_grav = tmp_grav+Response_2D_Rectangle(dens(ixiy,k),dxy-half_h,dxy+half_h,zgrid_edges[ixiy]-zs,grid.org[2]+(k+0.5)*(grid.h)-zs);								
							}
							else
							{   /*Two cells are divided into 8 parts.*//*The upper part is water.*/
                                tmp_grav = Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy+half_h,grid.org[2]+(l-0.5)*(grid.h)-zs,zgrid_edges[ixiy]-zs);
								tmp_grav = tmp_grav+Response_2D_upper_triangle1(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,zgrid_edges[ixiy]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
								tmp_grav = tmp_grav+Response_2D_lower_triangle1(dens(ixiy,l),dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,zgrid_edges[ixiy]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
                                tmp_grav = tmp_grav+Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),zgrid_edges[ixiy]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
                                tmp_grav = tmp_grav+Response_2D_Rectangle(dens(ixiy,k),dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy-1]-zs);
                                tmp_grav = tmp_grav+Response_2D_upper_triangle1(grid.dens_air_water,dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy-1]-zs);
                                tmp_grav = tmp_grav+Response_2D_lower_triangle1(dens(ixiy,k),dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy-1]-zs);
                                tmp_grav = tmp_grav+Response_2D_Rectangle(dens(ixiy,k),dxy-half_h,dxy+half_h,zgrid_edges[ixiy-1]-zs,grid.org[2]+(k+0.5)*(grid.h)-zs);								
							}

					    }
						else
						{ /*More than 2 cells are cut.*/
							if(zgrid_edges[ixiy-1] == min_tmp)
							{ /*The first and the last cells are divided into 8 parts.*//*The upper part is water.*/
								tmp_grav = Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy+half_h,grid.org[2]+(l-0.5)*(grid.h)-zs,zgrid_edges[ixiy-1]-zs);
								tmp_grav =tmp_grav+ Response_2D_upper_triangle(grid.dens_air_water,dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(l+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),zgrid_edges[ixiy-1]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
								tmp_grav =tmp_grav+ Response_2D_lower_triangle(dens(ixiy,l),dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(l+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),zgrid_edges[ixiy-1]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
								tmp_grav =tmp_grav+ Response_2D_Rectangle(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(l+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,zgrid_edges[ixiy-1]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);

								tmp_grav =tmp_grav+ Response_2D_Rectangle(dens(ixiy,k),dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy]-zs);
								tmp_grav =tmp_grav+ Response_2D_upper_triangle(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy]-zs);
								tmp_grav =tmp_grav+ Response_2D_lower_triangle(dens(ixiy,k),dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy]-zs);
								tmp_grav =tmp_grav+ Response_2D_Rectangle(dens(ixiy,k),dxy-half_h,dxy+half_h,zgrid_edges[ixiy]-zs,grid.org[2]+(k+0.5)*(grid.h)-zs);

								for(m=l+1;m<k;m++)
								{/*Each cell in the middle are devided into 4 parts.*//*The upper part is water.*/
		     						tmp_grav =tmp_grav+ Response_2D_Rectangle(dens(ixiy,m),dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(m-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(m-0.5)*(grid.h)-zs,grid.org[2]+(m+0.5)*(grid.h)-zs);
		     						tmp_grav =tmp_grav+ Response_2D_lower_triangle(dens(ixiy,m),dxy-half_h+(grid.h)*(grid.org[2]+(m-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy-half_h+(grid.h)*(grid.org[2]+(m+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(m-0.5)*(grid.h)-zs,grid.org[2]+(m+0.5)*(grid.h)-zs);
		     						tmp_grav =tmp_grav+ Response_2D_upper_triangle(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(m-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy-half_h+(grid.h)*(grid.org[2]+(m+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(m-0.5)*(grid.h)-zs,grid.org[2]+(m+0.5)*(grid.h)-zs);
		     						tmp_grav =tmp_grav+ Response_2D_Rectangle(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(m+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,grid.org[2]+(m-0.5)*(grid.h)-zs,grid.org[2]+(m+0.5)*(grid.h)-zs);
								}
							}
							else
							{/*The upper part is water.*/
                                tmp_grav = Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy+half_h,grid.org[2]+(l-0.5)*(grid.h)-zs,zgrid_edges[ixiy]-zs);
								tmp_grav =tmp_grav+ Response_2D_upper_triangle1(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(l+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,zgrid_edges[ixiy]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
								tmp_grav =tmp_grav+ Response_2D_lower_triangle1(dens(ixiy,l),dxy-half_h+(grid.h)*(grid.org[2]+(l+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,zgrid_edges[ixiy]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);
								tmp_grav =tmp_grav+ Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(l+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),zgrid_edges[ixiy]-zs,grid.org[2]+(l+0.5)*(grid.h)-zs);

								tmp_grav =tmp_grav+ Response_2D_Rectangle(dens(ixiy,k),dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy-1]-zs);
								tmp_grav =tmp_grav+ Response_2D_upper_triangle1(grid.dens_air_water,dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy-1]-zs);
								tmp_grav =tmp_grav+ Response_2D_lower_triangle1(dens(ixiy,k),dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(k-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(k-0.5)*(grid.h)-zs,zgrid_edges[ixiy-1]-zs);
								tmp_grav =tmp_grav+ Response_2D_Rectangle(dens(ixiy,k),dxy-half_h,dxy+half_h,zgrid_edges[ixiy-1]-zs,grid.org[2]+(k+0.5)*(grid.h)-zs);

								for(m=l+1;m<k;m++)
								{/*The upper part is water.*/
		     						tmp_grav =tmp_grav+ Response_2D_Rectangle(grid.dens_air_water,dxy-half_h,dxy-half_h+(grid.h)*(grid.org[2]+(m+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(m-0.5)*(grid.h)-zs,grid.org[2]+(m+0.5)*(grid.h)-zs);
		     						tmp_grav =tmp_grav+ Response_2D_lower_triangle1(dens(ixiy,m),dxy-half_h+(grid.h)*(grid.org[2]+(m+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy-half_h+(grid.h)*(grid.org[2]+(m-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(m-0.5)*(grid.h)-zs,grid.org[2]+(m+0.5)*(grid.h)-zs);
		     						tmp_grav =tmp_grav+ Response_2D_upper_triangle1(grid.dens_air_water,dxy-half_h+(grid.h)*(grid.org[2]+(m+0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy-half_h+(grid.h)*(grid.org[2]+(m-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),grid.org[2]+(m-0.5)*(grid.h)-zs,grid.org[2]+(m+0.5)*(grid.h)-zs);
		     						tmp_grav =tmp_grav+ Response_2D_Rectangle(dens(ixiy,m),dxy-half_h+(grid.h)*(grid.org[2]+(m-0.5)*(grid.h)-zgrid_edges[ixiy-1])/(zgrid_edges[ixiy]-zgrid_edges[ixiy-1]),dxy+half_h,grid.org[2]+(m-0.5)*(grid.h)-zs,grid.org[2]+(m+0.5)*(grid.h)-zs);
								}
							}

						}

					    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					    grav = grav + tmp_grav;					
					
				    }
					
					for(iz=k+1;iz<grid.nz;iz++)
					{
                        z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					    dz = z - zs;					/*z-component of distance between grid cell and station*/
					    if(dz == 0)
						    dz = eps;

						dxyz = sqrt(dxy*dxy + dz*dz);
						
					    if(dxyz < (grid.h * FACTOR))
					    {
						    /*exact calculation of attraction of two-dimensional rectangular blocks for cells close to the corresponding station following Heiland (1946)*/
						    tmp_grav = Response_2D_Rectangle(dens(ixiy,iz),dxy-half_h,dxy+half_h,dz-half_h,dz +half_h);
					    }
					    else
						    /*infinitely long horizontal rod*/
						    tmp_grav = (2*const1*dens(ixiy,iz))/(dz*(1 +((dxy*dxy)/(dz*dz)))); 

					    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					    grav = grav + tmp_grav;
					}
				}
				else
				{
                    for(iz=0;iz<grid.nz;iz++)
				    {
					    z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
					    dz = z - zs;					/*z-component of distance between grid cell and station*/
					    if(dz == 0)
						    dz = eps;

						dxyz = sqrt(dxy*dxy + dz*dz);

					    if(dxyz < (grid.h * FACTOR))
					    {
						    /*exact calculation of attraction of two-dimensional rectangular blocks for cells close to the corresponding station following Heiland (1946)*/
						    tmp_grav = Response_2D_Rectangle(dens(ixiy,iz),dxy-half_h,dxy+half_h,dz-half_h,dz +half_h);
					    }
					    else
						    /*infinitely long horizontal rod*/
						    tmp_grav = (2*const1*dens(ixiy,iz))/(dz*(1 +((dxy*dxy)/(dz*dz)))); 

					    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

					    grav = grav + tmp_grav;
				    }
				}
			}

			/*Assign the calculated density value to the corresponding measurements*/
			for(j=0;j<data->ndata_grav;j++)
			{
				if(data->gravs[i] == data->gno[j])
					data->calc_grav[j] = grav;
			}

			if(i%10 == 0)
				printf("Modeling for %d of %d stations is completed\n", i+1, geo.nstat_grav);
		}

	    free(xgrid_edges);
	    free(ygrid_edges);
	    free(zgrid_edges);

	}

	#undef dens

	/******************************************************/
					/*3D-calculation*/
	/******************************************************/
   
  	#define dens(x,y,z) grid.dens[((grid.ny)*(grid.nz))*(x) + (grid.nz)*(y) + (z)]

    #define Xgrid(x,y) xgrid_edges[(refine*(grid.ny-2)+1)*(x)+(y)] /*x coordinates of the grid points*/
    #define Ygrid(x,y) ygrid_edges[(refine*(grid.ny-2)+1)*(x)+(y)] /*y coordinates of the grid points*/
    #define Zgrid(x,y) zgrid_edges[(refine*(grid.ny-2)+1)*(x)+(y)] /*z coordinates of the cutting points of the (x,y) grid edge*/

	else
	{

		printf("The gravity model is identified as 3-D\n");
		printf("Start with the gravity modeling\n");
		printf("-------------------------------\n\n");

        xgrid_edges = (double *)memory(NULL,(refine*(grid.nx-2)+1)*(refine*(grid.ny-2)+1),sizeof(double),"ForwardModGrav");
        ygrid_edges = (double *)memory(NULL,(refine*(grid.nx-2)+1)*(refine*(grid.ny-2)+1),sizeof(double),"ForwardModGrav");
        zgrid_edges = (double *)memory(NULL,(refine*(grid.nx-2)+1)*(refine*(grid.ny-2)+1),sizeof(double),"ForwardModGrav");

        for(i=0;i<(refine*(grid.nx-2)+1);i++)
			for(j=0;j<(refine*(grid.ny-2)+1);j++)
			{
				Xgrid(i,j) = grid.org[0] + half_h + (grid.h/refine)*i;
				Ygrid(i,j) = grid.org[1] + half_h + (grid.h/refine)*j;
			}
		/*calculate the z coordinate of the cutting points*/
		if(topo.nr_of_triangles != 0) /*BJOERN_MOD5*/
			TriangleInterpol(xgrid_edges,ygrid_edges,zgrid_edges,(refine*(grid.nx-2)+1)*(refine*(grid.ny-2)+1),topo.x,topo.y,topo.z,topo.nr_of_topo_points,topo.index,topo.nr_of_triangles);

		/*Loop over all stations*/
		for(i=0;i<geo.nstat_grav;i++)
		{
			/*Positions of the gravity stations:*/
			xs = (double)geo.x[(data->gravs[i]-1)];
			ys = (double)geo.y[(data->gravs[i]-1)];
			zs = (double)geo.z[(data->gravs[i]-1)];

			grav = 0.0;

			/*left left (ix==0)*/
			ix = 0;
			x1 = (double)grid.org[0] + (ix+0.5)*(grid.h);	 /*x component of the right edge of the cell in m*/
			dx1 = x1 - xs;				 /*x-component of distance between grid cell and station*/
			dxx1 = fabs(dx1);
			/*right side (ix = nx-1)*/
			ix = grid.nx-1;
			x2 = (double)grid.org[0] + (ix-0.5)*(grid.h);	 /*x component of the left edge of the cell in m*/
			dx2 = x2 - xs;				 /*x-component of distance between grid cell and station*/
			dxx2 = fabs(dx2);
			/*front side (iy==0)*/
			iy = 0;
			y1 = (double)grid.org[1] + (iy+0.5)*(grid.h);	 /*y component of the front edge of the cell in m*/
			dy1 = y1 - ys;				 /*y-component of distance between grid cell and station*/
			dyy1 = fabs(dy1);
			/*back side (iy = ny-1)*/
			iy = grid.ny-1;
			y2 = (double)grid.org[1] + (iy-0.5)*(grid.h);	 /*y component of the back edge of the cell in m*/
			dy2 = y2 - ys;				 /*y-component of distance between grid cell and station*/
			dyy2 = fabs(dy2);


			/**********************/
			/*1.Case: Cells at the edges of the 3-D grid*/

			/*Loop over all cell at the edges (ix=0,iy=0; ix=0,iy=ny-1; ix=nx-1,iy=0; ix=nx-1,iy=ny-1)*/

			if(topo.nr_of_triangles != 0)
			{
				zcut1=zgrid_edges[0];
				zcut2=zgrid_edges[(refine*(grid.nx-2))*(refine*(grid.ny-2)+1)];
                zcut3=zgrid_edges[refine*(grid.ny-2)];
				zcut4=zgrid_edges[(refine*(grid.nx-2))*(refine*(grid.ny-2)+1)+refine*(grid.ny-2)];

			    for(iz=0;iz<grid.nz;iz++)
			    {
				    z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
				    dz = z - zs;					/*z-component of distance between grid cell and station*/
				    if(dz == 0)
					    dz = eps;

					dz1=dz-half_h; /*z coordinate of the upper border in m*/
					dz2=dz+half_h;/*z coordinate of the lower border in m*/

				    dxyz = sqrt(dx1*dx1 + dy1*dy1 + dz*dz);
				    if(dxyz == 0)
					    dxyz = eps;


				    /*Calculate the gravitation from the cell at the station*/

					/*left front edge*/
					if((zcut1 != -99999.9) && (zcut1 >= z - half_h) && (zcut1 < z + half_h))
					{/* The cell is divided in two part*/

					    if(dx1 <= 0 && dy1 <= 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_2(grid.dens_air_water,-dx1,-dy1,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),-dx1,-dy1,zcut1-zs,dz2);
						}
						else if(dx1 > 0 && dy1 <= 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,0,dx1,-dy1,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,0,iz),0,dx1,-dy1,zcut1-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,0,-dy1,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),0,-dy1,zcut1-zs,dz2);
						}
						else if(dx1 <= 0 && dy1 > 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,0,dy1,-dx1,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,0,iz),0,dy1,-dx1,zcut1-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,-dx1,0,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),-dx1,0,zcut1-zs,dz2);
						}
				        else
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle(grid.dens_air_water,0,dx1,0,dy1,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle(dens(0,0,iz),0,dx1,0,dy1,zcut1-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dx1,0,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,0,iz),0,dx1,0,zcut1-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dy1,0,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,0,iz),0,dy1,0,zcut1-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,0,0,dz1,zcut1-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),0,0,zcut1-zs,dz2);						
						}
					}
					else
					{
						if(dxyz < (grid.h * FACTOR))
						{
						    if(dx1 <= 0 && dy1 <= 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_2(dens(0,0,iz),-dx1,-dy1,dz1,dz2);
							}
							else if(dx1 > 0 && dy1 <= 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_1(dens(0,0,iz),0,dx1,-dy1,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),0,-dy1,dz1,dz2);
							}
							else if(dx1 <= 0 && dy1 > 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_1(dens(0,0,iz),0,dy1,-dx1,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),-dx1,0,dz1,dz2);
							}
							else
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle(dens(0,0,iz),0,dx1,0,dy1,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,0,iz),0,dx1,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,0,iz),0,dy1,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),0,0,dz1,dz2);						
							}
						}
						else
 					    {
					        if(dz >= 0)
								tmp_grav = (G_KONST*dens(0,0,iz)*grid.h)*((PI/2) - atan(dxx1/dz) - atan(dyy1/dz) + atan((dxx1*dyy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
							else
								tmp_grav = (G_KONST*dens(0,0,iz)*grid.h)*((-PI/2) - atan(dxx1/dz) - atan(dyy1/dz) + atan((dxx1*dyy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
						}
					}
					tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

				    grav = grav + tmp_grav;

				    /*right front edge*/
				    dxyz = sqrt(dx2*dx2 + dy1*dy1 + dz*dz);
				    if(dxyz == 0)
					    dxyz = eps;

				    if(( zcut2!= -99999.9) && (zcut2 >= z - half_h) && (zcut2 < z + half_h))
					{ /*The cell is divided into 2 parts.*/
					    if(dx2 >= 0 && dy1 <= 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_2(grid.dens_air_water,dx2,-dy1,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),dx2,-dy1,zcut2-zs,dz2);
						}
						else if(dx2 < 0 && dy1 <= 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dx2,0,-dy1,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),dx2,0,-dy1,zcut2-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,0,-dy1,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),0,-dy1,zcut2-zs,dz2);
						}
						else if(dx2 >= 0 && dy1 > 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,0,dy1,dx2,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),0,dy1,dx2,zcut2-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,dx2,0,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),dx2,0,zcut2-zs,dz2);
						}
				        else
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx2,0,0,dy1,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle(dens(grid.nx-1,0,iz),dx2,0,0,dy1,zcut2-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dx2,0,0,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),dx2,0,0,zcut2-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dy1,0,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),0,dy1,0,zcut2-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,0,0,dz1,zcut2-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),0,0,zcut2-zs,dz2);						
						}
					}
					else
					{
						if(dxyz < (grid.h * FACTOR))
						{
						   if(dx2 >= 0 && dy1 <= 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),dx2,-dy1,dz1,dz2);
							}
							else if(dx2 < 0 && dy1 <= 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),dx2,0,-dy1,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),0,-dy1,dz1,dz2);
							}
							else if(dx2 >= 0 && dy1 > 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),0,dy1,dx2,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),dx2,0,dz1,dz2);
							}
							else
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle(dens(grid.nx-1,0,iz),dx2,0,0,dy1,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),dx2,0,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),0,dy1,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),0,0,dz1,dz2);						
							}
						}
						else
						{
					        if(dz >= 0)
						        tmp_grav = (G_KONST*dens(grid.nx-1,0,iz)*grid.h)*((PI/2) - atan(dxx2/dz) - atan(dyy1/dz) + atan((dxx2*dyy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					        else
						        tmp_grav = (G_KONST*dens(grid.nx-1,0,iz)*grid.h)*((-PI/2) - atan(dxx2/dz) - atan(dyy1/dz) + atan((dxx2*dyy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
						}
					}
				    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

				    grav = grav + tmp_grav;

				    /*left back edge*/
				    dxyz = sqrt(dx1*dx1 + dy2*dy2 + dz*dz);
				    if(dxyz == 0)
					    dxyz = eps;

				    if(( zcut3!= -99999.9) && (zcut3 >= z - half_h) && (zcut3 < z + half_h))
					{/*The cell is divided into 2 parts.*/
					    if(dx1 <= 0 && dy2 >= 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_2(grid.dens_air_water,-dx1,dy2,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),-dx1,dy2,zcut3-zs,dz2);
						}
						else if(dx1 > 0 && dy2 >= 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,0,dx1,dy2,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),0,dx1,dy2,zcut3-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,0,dy2,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),0,dy2,zcut3-zs,dz2);
						}
						else if(dx1 <= 0 && dy2 < 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dy2,0,-dx1,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),dy2,0,-dx1,zcut3-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,-dx1,0,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),-dx1,0,zcut3-zs,dz2);
						}
				        else
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle(grid.dens_air_water,0,dx1,dy2,0,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle(dens(0,grid.ny-1,iz),0,dx1,dy2,0,zcut3-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dx1,0,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),0,dx1,0,zcut3-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dy2,0,0,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),dy2,0,0,zcut3-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,0,0,dz1,zcut3-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),0,0,zcut3-zs,dz2);						
						}
					}
					else
					{
						if(dxyz < grid.h * FACTOR)
						{
							if(dx1 <= 0 && dy2 >= 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),-dx1,dy2,dz1,dz2);
							}	
							else if(dx1 > 0 && dy2 >= 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),0,dx1,dy2,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),0,dy2,dz1,dz2);
							}
							else if(dx1 <= 0 && dy2 < 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),dy2,0,-dx1,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),-dx1,0,dz1,dz2);
							}
							else
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle(dens(0,grid.ny-1,iz),0,dx1,dy2,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),0,dx1,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),dy2,0,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),0,0,dz1,dz2);						
							}
						}
						else
						{
					        if(dz >= 0)
						        tmp_grav = (G_KONST*dens(0,grid.ny-1,iz)*grid.h)*((PI/2) - atan(dxx1/dz) - atan(dyy2/dz) + atan((dxx1*dyy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					        else
						        tmp_grav = (G_KONST*dens(0,grid.ny-1,iz)*grid.h)*((-PI/2) - atan(dxx1/dz) - atan(dyy2/dz) + atan((dxx1*dyy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
						}
					}
				    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
				    grav = grav + tmp_grav;

				    /*right back edge*/
				    dxyz = sqrt(dx2*dx2 + dy2*dy2 + dz*dz);
				    if(dxyz == 0)
					    dxyz = eps;

				    if(( zcut4!= -99999.9) && (zcut4 >= z - half_h) && (zcut4 < z + half_h))
					{/*The cell is divided into 2 parts.*/
					    if(dx2 >= 0 && dy2 >= 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_2(grid.dens_air_water,dx2,dy2,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),dx2,dy2,zcut4-zs,dz2);
						}
						else if(dx2 < 0 && dy2 >= 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dx2,0,dy2,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dx2,0,dy2,zcut4-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,0,dy2,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),0,dy2,zcut4-zs,dz2);
						}
						else if(dx2 >= 0 && dy2 < 0)
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dy2,0,dx2,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dy2,0,dx2,zcut4-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,dx2,0,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),dx2,0,zcut4-zs,dz2);
						}
				        else
						{/*Two semi-semi-infinite horizontal sheets*//*The upper part is water.*/
							tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx2,0,dy2,0,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle(dens(grid.nx-1,grid.ny-1,iz),dx2,0,dy2,0,zcut4-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dx2,0,0,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dx2,0,0,zcut4-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dy2,0,0,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dy2,0,0,zcut4-zs,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(grid.dens_air_water,0,0,dz1,zcut4-zs);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),0,0,zcut4-zs,dz2);						
						}

					}
					else
					{
						if(dxyz < grid.h * FACTOR)
						{
							if(dx2 >= 0 && dy2 >= 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),dx2,dy2,dz1,dz2);
							}
							else if(dx2 < 0 && dy2 >= 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dx2,0,dy2,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),0,dy2,dz1,dz2);
							}
							else if(dx2 >= 0 && dy2 < 0)
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dy2,0,dx2,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),dx2,0,dz1,dz2);
							}
							else
							{/*Two semi-semi-infinite horizontal sheets*/
								tmp_grav = Response_3D_Rectangle(dens(grid.nx-1,grid.ny-1,iz),dx2,0,dy2,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dx2,0,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dy2,0,0,dz1,dz2);
								tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),0,0,dz1,dz2);						
							}
						}
						else
	 					{
							if(dz >= 0)
						        tmp_grav = (G_KONST*dens(grid.nx-1,grid.ny-1,iz)*grid.h)*((PI/2) - atan(dxx2/dz) - atan(dyy2/dz) + atan((dxx2*dyy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					        else
						        tmp_grav = (G_KONST*dens(grid.nx-1,grid.ny-1,iz)*grid.h)*((-PI/2) - atan(dxx2/dz) - atan(dyy2/dz) + atan((dxx2*dyy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
						}
					}
				    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
				    grav = grav + tmp_grav;
			    }
			}
			else
			{
			    for(iz=0;iz<grid.nz;iz++)
			    {
				    z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
				    dz = z - zs;					/*z-component of distance between grid cell and station*/
				    if(dz == 0)
					    dz = eps;

					dz1=dz-half_h; /*z coordinate of the upper border in m*/
					dz2=dz+half_h;/*z coordinate of the lower border in m*/

				    /*Calculate the gravitation from the cell at the station*/
				    /*left front edge*/
				    dxyz = sqrt(dx1*dx1 + dy1*dy1 + dz*dz);
				    if(dxyz == 0)
					    dxyz = eps;

					if(dxyz < grid.h * FACTOR)
					{
						if(dx1 <= 0 && dy1 <= 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_2(dens(0,0,iz),-dx1,-dy1,dz1,dz2);
						}
						else if(dx1 > 0 && dy1 <= 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_1(dens(0,0,iz),0,dx1,-dy1,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),0,-dy1,dz1,dz2);
						}
						else if(dx1 <= 0 && dy1 > 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_1(dens(0,0,iz),0,dy1,-dx1,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),dx1,0,dz1,dz2);
						}
						else
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle(dens(0,0,iz),0,dx1,0,dy1,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,0,iz),0,dx1,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,0,iz),0,dy1,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,0,iz),0,0,dz1,dz2);						
						}
					}
					else
					{
						if(dz >= 0)
						    tmp_grav = (G_KONST*dens(0,0,iz)*grid.h)*((PI/2) - atan(dxx1/dz) - atan(dyy1/dz) + atan((dxx1*dyy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
						else
						    tmp_grav = (G_KONST*dens(0,0,iz)*grid.h)*((-PI/2) - atan(dxx1/dz) - atan(dyy1/dz) + atan((dxx1*dyy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					}
				    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

				    grav = grav + tmp_grav;

				    /*right front edge*/
				    dxyz = sqrt(dx2*dx2 + dy1*dy1 + dz*dz);
				    if(dxyz == 0)
					    dxyz = eps;

					if(dxyz < grid.h * FACTOR)
					{
					    if(dx2 >= 0 && dy1 <= 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),dx2,-dy1,dz1,dz2);
						}
						else if(dx2 < 0 && dy1 <= 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),dx2,0,-dy1,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),0,-dy1,dz1,dz2);
						}
						else if(dx2 >= 0 && dy1 > 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),0,dy1,dx2,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),dx2,0,dz1,dz2);
						}
						else
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle(dens(grid.nx-1,0,iz),dx2,0,0,dy1,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),dx2,0,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,0,iz),0,dy1,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,0,iz),0,0,dz1,dz2);						
						}
					}
					else
					{
					    if(dz >= 0)
						    tmp_grav = (G_KONST*dens(grid.nx-1,0,iz)*grid.h)*((PI/2) - atan(dxx2/dz) - atan(dyy1/dz) + atan((dxx2*dyy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					    else
						    tmp_grav = (G_KONST*dens(grid.nx-1,0,iz)*grid.h)*((-PI/2) - atan(dxx2/dz) - atan(dyy1/dz) + atan((dxx2*dyy1)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					}
				    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

				    grav = grav + tmp_grav;

				    /*left back edge*/
				    dxyz = sqrt(dx1*dx1 + dy2*dy2 + dz*dz);
				    if(dxyz == 0)
					    dxyz = eps;

					if(dxyz < grid.h * FACTOR)
					{
					   if(dx1 <= 0 && dy2 >= 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),-dx1,dy2,dz1,dz2);
						}
						else if(dx1 > 0 && dy2 >= 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),0,dx1,dy2,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),0,dy2,dz1,dz2);
						}
						else if(dx1 <= 0 && dy2 < 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),dy2,0,-dx1,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),-dx1,0,dz1,dz2);
						}
						else
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle(dens(0,grid.ny-1,iz),0,dx1,dy2,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),0,dx1,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,grid.ny-1,iz),dy2,0,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(0,grid.ny-1,iz),0,0,dz1,dz2);						
						}
					}
					else
					{
					    if(dz >= 0)
						    tmp_grav = (G_KONST*dens(0,grid.ny-1,iz)*grid.h)*((PI/2) - atan(dxx1/dz) - atan(dyy2/dz) + atan((dxx1*dyy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					    else
						    tmp_grav = (G_KONST*dens(0,grid.ny-1,iz)*grid.h)*((-PI/2) - atan(dxx1/dz) - atan(dyy2/dz) + atan((dxx1*dyy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					}
				    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
				    grav = grav + tmp_grav;

				    /*right back edge*/
				    dxyz = sqrt(dx2*dx2 + dy2*dy2 + dz*dz);
				    if(dxyz == 0)
					    dxyz = eps;

					if(dxyz < grid.h * FACTOR)
					{
					  if(dx2 >= 0 && dy2 >= 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),dx2,dy2,dz1,dz2);
						}
						else if(dx2 < 0 && dy2 >= 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dx2,0,dy2,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),0,dy2,dz1,dz2);
						}
						else if(dx2 >= 0 && dy2 < 0)
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dy2,0,dx2,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),dx2,0,dz1,dz2);
						}
						else
						{/*Two semi-semi-infinite horizontal sheets*/
							tmp_grav = Response_3D_Rectangle(dens(grid.nx-1,grid.ny-1,iz),dx2,0,dy2,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dx2,0,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,grid.ny-1,iz),dy2,0,0,dz1,dz2);
							tmp_grav = tmp_grav + Response_3D_Rectangle_2(dens(grid.nx-1,grid.ny-1,iz),0,0,dz1,dz2);						
						}
					}
					else
					{
					    if(dz >= 0)
						    tmp_grav = (G_KONST*dens(grid.nx-1,grid.ny-1,iz)*grid.h)*((PI/2) - atan(dxx2/dz) - atan(dyy2/dz) + atan((dxx2*dyy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					    else
						    tmp_grav = (G_KONST*dens(grid.nx-1,grid.ny-1,iz)*grid.h)*((-PI/2) - atan(dxx2/dz) - atan(dyy2/dz) + atan((dxx2*dyy2)/(dz*dxyz))); /*semi-semi-infinite horizontal sheet*/
					}
				    tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
				    grav = grav + tmp_grav;
			    }
			}
						
			/**********************/
			/*2.Case: Cells at the sides (except for the edges) of the 3-D grid*/

			/*Loop over all cell at the right and left sides (ix = 0 and ix = nx -1)*/
			if(topo.nr_of_triangles != 0)
			{
			    for(iy=1;iy<(grid.ny-1);iy++)
				{
					y = grid.org[1] + iy*(grid.h);		 /*y component of the cell centers in m*/
					dy = y - ys;					 /*y-component of distance between grid cell and station*/

					dy_1 = y - ys - half_h;	/*y-component of distance between front edge of the grid cell and station*/
					dy_2 = y - ys + half_h;  /*y-component of distance between back edge of the grid cell and station*/

					/*Left side, two z cutting value.*/
					zcut1 = zgrid_edges[refine*(iy-1)];
					zcut2 = zgrid_edges[refine*iy];
					aver_z1 = 0.5*(zcut1+zcut2); /*The average of the values.*/
					/*Right side, two z cutting value.*/
					zcut3 = zgrid_edges[(refine*(grid.nx-2))*(refine*(grid.ny-2)+1)+refine*(iy-1)];
					zcut4 = zgrid_edges[(refine*(grid.nx-2))*(refine*(grid.ny-2)+1)+refine*iy];
					aver_z2 = 0.5*(zcut3+zcut4);/*The average of the values.*/

					for(iz=0;iz<grid.nz;iz++)
					{
						z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
						dz = z - zs;					/*z-component of distance between grid cell and station*/

						dz1=dz-half_h;
						dz2=dz+half_h;							

						/*Calculate the gravitation from the cell at the station*/
						/*left side*/
						dxyz = sqrt(dx1*dx1 + dy*dy + dz*dz);
						if(dxyz == 0)
							dxyz = eps;
						dyz = dz*dz + dy*dy;
						if(dyz == 0)
							dyz = eps;

						if((zcut1 != -99999.9) && (zcut2 != -99999.9) && (aver_z1 >= z-half_h) && (aver_z1 < z+half_h))
						{/*The cell is divided into 2 parts.*/

							/*The upper part is water.*/
							if(dx1 <= 0)
							{
								if(dy_1 <= 0 && dy_2 > 0)
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dy_1,0,-dx1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dy_2,-dx1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,0,-dx1,aver_z1-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),0,dy_2,-dx1,aver_z1-zs,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dy_1,dy_2,-dx1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,dy_2,-dx1,aver_z1-zs,dz2);
								}
							}
							else
							{
								if(dy_1 <= 0 && dy_2 > 0)
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle(grid.dens_air_water,0,dx1,dy_1,0,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(grid.dens_air_water,0,dx1,0,dy_2,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(0,iy,iz),0,dx1,dy_1,0,aver_z1-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle(dens(0,iy,iz),0,dx1,0,dy_2,aver_z1-zs,dz2);

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dy_1,0,0,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dy_2,0,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,0,0,aver_z1-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),0,dy_2,0,aver_z1-zs,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle(grid.dens_air_water,0,dx1,dy_1,dy_2,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(0,iy,iz),0,dx1,dy_1,dy_2,aver_z1-zs,dz2);

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dy_1,dy_2,0,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,dy_2,0,aver_z1-zs,dz2);
								}
							}
						}
						else
						{
							if(dxyz < grid.h * FACTOR)
							{
								if(dx1 <= 0)
								{
									if(dy_1 <= 0 && dy_2 > 0)
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,0,-dx1,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),0,dy_2,-dx1,dz1,dz2);
									}
									else
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,dy_2,-dx1,dz1,dz2);
									}
								}
								else
								{
									if(dy_1 <= 0 && dy_2 > 0)
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle(dens(0,iy,iz),0,dx1,dy_1,0,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle(dens(0,iy,iz),0,dx1,0,dy_2,dz1,dz2);

										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,0,0,dz1,dz2);
                						tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),0,dy_2,0,dz1,dz2);
									}
									else
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle(dens(0,iy,iz),0,dx1,dy_1,dy_2,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,dy_2,0,dz1,dz2);
									}
								}
							}
							else
		                        tmp_grav = (G_KONST*dens(0,iy,iz)*(grid.h)*(grid.h))*(dz/dyz)*(1-(dxx1/dxyz)); /*semi-infinitely long horizontal rod*/
						}
						tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
						grav = grav + tmp_grav;

						/*right side*/
						dxyz = sqrt(dx2*dx2 + dy*dy + dz*dz);
						if(dxyz == 0)
							dxyz = eps;

						if((zcut3 != -99999.9) && (zcut4 != -99999.9) && (aver_z2 >= z-half_h) && (aver_z2 < z+half_h))
						{
							/*The upper part is water.*/
							if(dx2 >= 0)
							{
								if(dy_1 <= 0 && dy_2 > 0)
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dy_1,0,dx2,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dy_2,dx2,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,0,dx2,aver_z2-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),0,dy_2,dx2,aver_z2-zs,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dy_1,dy_2,dx2,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,dy_2,dx2,aver_z2-zs,dz2);
								}
							}
							else
							{
								if(dy_1 <= 0 && dy_2 > 0)
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx2,0,dy_1,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(grid.dens_air_water,dx2,0,0,dy_2,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,dy_1,0,aver_z2-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,0,dy_2,aver_z2-zs,dz2);

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dy_1,0,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dy_2,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,0,0,aver_z2-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),0,dy_2,0,aver_z2-zs,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx2,0,dy_1,dy_2,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,dy_1,dy_2,aver_z2-zs,dz2);

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dy_1,dy_2,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,dy_2,0,aver_z2-zs,dz2);
								}
							}
						}
						else
						{
							if(dxyz < grid.h * FACTOR)
							{
								if(dx2 >= 0)
								{
									if(dy_1 <= 0 && dy_2 > 0)
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,0,dx2,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),0,dy_2,dx2,dz1,dz2);
									}
									else
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,dy_2,dx2,dz1,dz2);
									}
								}
								else
								{
									if(dy_1 <= 0 && dy_2 > 0)
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,dy_1,0,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,0,dy_2,dz1,dz2);

										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,0,0,dz1,dz2);
                						tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),0,dy_2,0,dz1,dz2);
									}
									else
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,dy_1,dy_2,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,dy_2,0,dz1,dz2);
									}
								}
							}
							else
							    tmp_grav = (G_KONST*dens(grid.nx-1,iy,iz)*(grid.h)*(grid.h))*(dz/dyz)*(1-(dxx2/dxyz)); /*semi-infinitely long horizontal rod*/
						}
						tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
						grav = grav + tmp_grav;
					}
				}
			}
			else
			{
			    for(iy=1;iy<(grid.ny-1);iy++)
				{
					y = grid.org[1] + iy*(grid.h);		 /*y component of the cell centers in m*/
					dy = y - ys;					 /*y-component of distance between grid cell and station*/

					dy_1 = y - ys - half_h;	/*y-component of distance between front edge of the grid cell and station*/
					dy_2 = y - ys + half_h;  /*y-component of distance between back edge of the grid cell and station*/

					for(iz=0;iz<grid.nz;iz++)
					{
						z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
						dz = z - zs;					/*z-component of distance between grid cell and station*/

						dz1=dz-half_h;
						dz2=dz+half_h;							

						/*Calculate the gravitation from the cell at the station*/
						/*left side*/
						dxyz = sqrt(dx1*dx1 + dy*dy + dz*dz);
						if(dxyz == 0)
							dxyz = eps;
						dyz = dz*dz + dy*dy;
						if(dyz == 0)
							dyz = eps;

						if(dxyz < grid.h * FACTOR)
						{
							if(dx1 <= 0)
							{
								if(dy_1 <= 0 && dy_2 > 0)
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,0,-dx1,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),0,dy_2,-dx1,dz1,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,dy_2,-dx1,dz1,dz2);
								}
							}
							else
							{
								if(dy_1 <= 0 && dy_2 > 0)
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle(dens(0,iy,iz),0,dx1,dy_1,0,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(0,iy,iz),0,dx1,0,dy_2,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,0,0,dz1,dz2);
               						tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),0,dy_2,0,dz1,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle(dens(0,iy,iz),0,dx1,dy_1,dy_2,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(0,iy,iz),dy_1,dy_2,0,dz1,dz2);
								}
							}
						}
						else
							tmp_grav = (G_KONST*dens(0,iy,iz)*(grid.h)*(grid.h))*(dz/dyz)*(1-(dxx1/dxyz)); /*semi-infinitely long horizontal rod*/
						tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
						grav = grav + tmp_grav;

						/*right side*/
						dxyz = sqrt(dx2*dx2 + dy*dy + dz*dz);
						if(dxyz == 0)
							dxyz = eps;

						if(dxyz < grid.h * FACTOR)
						{
							if(dx2 >= 0)
							{
								if(dy_1 <= 0 && dy_2 > 0)
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,0,dx2,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),0,dy_2,dx2,dz1,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,dy_2,dx2,dz1,dz2);
								}
							}
							else
							{
								if(dy_1 <= 0 && dy_2 > 0)
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,dy_1,0,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,0,dy_2,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,0,0,dz1,dz2);
               						tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),0,dy_2,0,dz1,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle(dens(grid.nx-1,iy,iz),dx2,0,dy_1,dy_2,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(grid.nx-1,iy,iz),dy_1,dy_2,0,dz1,dz2);
								}
							}
						}
						else
							tmp_grav = (G_KONST*dens(grid.nx-1,iy,iz)*(grid.h)*(grid.h))*(dz/dyz)*(1-(dxx2/dxyz)); /*semi-infinitely long horizontal rod*/
						tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
						grav = grav + tmp_grav;
					}
				}
			}
			/*Loop over all cell at the front and back sides (iy = 0 and iy = ny -1)*/
			if(topo.nr_of_triangles != 0)
			{
				for(ix=1;ix<(grid.nx-1);ix++)
				{
					x = grid.org[0] + ix*(grid.h);		 /*x component of the cell centers in m*/
					dx = x - xs;					 /*x-component of distance between grid cell and station*/

					dx_1 = x - xs - half_h;	/*x-component of distance between front edge of the grid cell and station*/
					dx_2 = x - xs + half_h;  /*x-component of distance between back edge of the grid cell and station*/

					/*Front side, two z cutting value.*/
					zcut1 = zgrid_edges[(refine*(ix-1))*(refine*(grid.ny-2)+1)];
					zcut2 = zgrid_edges[(refine*ix)*(refine*(grid.ny-2)+1)];
					aver_z1 = 0.5*(zcut1+zcut2);/*The average of the values.*/
					/*Back side, two z cutting value.*/
					zcut3 = zgrid_edges[(refine*(ix-1))*(refine*(grid.ny-2)+1)+refine*(grid.ny-2)];
					zcut4 = zgrid_edges[(refine*ix)*(refine*(grid.ny-2)+1)+refine*(grid.ny-2)];
					aver_z2 = 0.5*(zcut3+zcut4);/*The average of the values.*/

					for(iz=0;iz<grid.nz;iz++)
					{
						z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
						dz = z - zs;					/*z-component of distance between grid cell and station*/

						dz1 = dz - half_h;
						dz2 = dz + half_h;
						/*Calculate the gravitation from the cell at the station*/
						/*front side*/
						dxyz = sqrt(dx*dx + dy1*dy1 + dz*dz);
						if(dxyz == 0)
							dxyz = eps;
						dxz = dz*dz + dx*dx;
						if(dxz == 0)
							dxz = eps;

						if((zcut1 != -99999.9) && (zcut2 != -99999.9) && (aver_z1 >= z-half_h) && (aver_z1 < z+half_h))
						{
							/*The upper part is water.*/
							if(dy1 <= 0)
							{
								if(dx_1 <= 0 && dx_2 > 0)
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dx_1,0,-dy1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dx_2,-dy1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,0,-dy1,aver_z1-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),0,dx_2,-dy1,aver_z1-zs,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dx_1,dx_2,-dy1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,dx_2,-dy1,aver_z1-zs,dz2);
								}
							}
							else
							{
								if(dx_1 <= 0 && dx_2 > 0)
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx_1,0,0,dy1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(grid.dens_air_water,0,dx_2,0,dy1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,0,iz),dx_1,0,0,dy1,aver_z1-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,0,iz),0,dx_2,0,dy1,aver_z1-zs,dz2);

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dx_1,0,0,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dx_2,0,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,0,0,aver_z1-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),0,dx_2,0,aver_z1-zs,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx_1,dx_2,0,dy1,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,0,iz),dx_1,dx_2,0,dy1,aver_z1-zs,dz2);

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dx_1,dx_2,0,dz1,aver_z1-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,dx_2,0,aver_z1-zs,dz2);
								}
							}
						}
						else
						{
							if(dxyz < grid.h * FACTOR)
							{
								if(dy1 <= 0)
								{
									if(dx_1 <= 0 && dx_2 > 0)
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,0,-dy1,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),0,dx_2,-dy1,dz1,dz2);
									}
									else
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,dx_2,-dy1,dz1,dz2);
									}
								}
								else
								{
									if(dx_1 <= 0 && dx_2 > 0)
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle(dens(ix,0,iz),dx_1,0,0,dy1,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,0,iz),0,dx_2,0,dy1,dz1,dz2);

										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,0,0,dz1,dz2);
                						tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),0,dx_2,0,dz1,dz2);
									}
									else
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle(dens(ix,0,iz),dx_1,dx_2,0,dy1,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,dx_2,0,dz1,dz2);
									}
								}
							}
							else
								tmp_grav = (G_KONST*dens(ix,0,iz)*(grid.h)*(grid.h))*(dz/dxz)*(1-(dyy1/dxyz)); /*semi-infinitely long horizontal rod*/
						}
						tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
						grav = grav + tmp_grav;

						/*back side*/
						dxyz = sqrt(dx*dx + dy2*dy2 + dz*dz);
						if(dxyz == 0)
							dxyz = eps;

						if((zcut3 != -99999.9) && (zcut4 != -99999.9) && (aver_z2 >= z-half_h) && (aver_z2 < z+half_h))
						{ /*The upper part is water.*/
							if(dy2 >= 0)
							{
								if(dx_1 <= 0 && dx_2 > 0)
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dx_1,0,dy2,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dx_2,dy2,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,0,dy2,aver_z2-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),0,dx_2,dy2,aver_z2-zs,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle_1(grid.dens_air_water,dx_1,dx_2,dy2,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,dx_2,dy2,aver_z2-zs,dz2);
								}
							}
							else
							{
								if(dx_1 <= 0 && dx_2 > 0)
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx_1,0,dy2,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(grid.dens_air_water,0,dx_2,dy2,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,grid.ny-1,iz),dx_1,0,dy2,0,aver_z2-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,grid.ny-1,iz),0,dx_2,dy2,0,aver_z2-zs,dz2);

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dx_1,0,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,0,dx_2,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,0,0,aver_z2-zs,dz2);
                					tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),0,dx_2,0,aver_z2-zs,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*//*The upper part is water.*/
									tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx_1,dx_2,dy2,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,grid.ny-1,iz),dx_1,dx_2,dy2,0,aver_z2-zs,dz2);

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(grid.dens_air_water,dx_1,dx_2,0,dz1,aver_z2-zs);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,dx_2,0,aver_z2-zs,dz2);
								}
							}
						}
						else
						{
							if(dxyz < grid.h * FACTOR)
							{
								if(dy2 >= 0)
								{
									if(dx_1 <= 0 && dx_2 > 0)
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,0,dy2,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),0,dx_2,dy2,dz1,dz2);
									}
									else
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,dx_2,dy2,dz1,dz2);
									}
								}
								else
								{
									if(dx_1 <= 0 && dx_2 > 0)
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle(dens(ix,grid.ny-1,iz),dx_1,0,dy2,0,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,grid.ny-1,iz),0,dx_2,dy2,0,dz1,dz2);

										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,0,0,dz1,dz2);
                						tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),0,dx_2,0,dz1,dz2);
									}
									else
									{/*Two semi-infinite horizontal sheets*/
										tmp_grav = Response_3D_Rectangle(dens(ix,grid.ny-1,iz),dx_1,dx_2,dy2,0,dz1,dz2);
										tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,dx_2,0,dz1,dz2);
									}
								}
							}
							else
								tmp_grav = (G_KONST*dens(ix,grid.ny-1,iz)*(grid.h)*(grid.h))*(dz/dxz)*(1-(dyy2/dxyz)); /*semi-infinitely long horizontal rod*/
						}
						tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
						grav = grav + tmp_grav;
					}
				}
			}
			else
			{
				for(ix=1;ix<(grid.nx-1);ix++)
				{
					x = grid.org[0] + ix*(grid.h);		 /*x component of the cell centers in m*/
					dx = x - xs;					 /*x-component of distance between grid cell and station*/

					dx_1 = x - xs - half_h;	/*x-component of distance between front edge of the grid cell and station*/
					dx_2 = x - xs + half_h;  /*x-component of distance between back edge of the grid cell and station*/

					for(iz=0;iz<grid.nz;iz++)
					{
						z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
						dz = z - zs;					/*z-component of distance between grid cell and station*/

						dz1=dz-half_h; /*z coordinate of the upper border in m*/
						dz2=dz+half_h;/*z coordinate of the lower border in m*/

						/*Calculate the gravitation from the cell at the station*/
						/*front side*/
						dxyz = sqrt(dx*dx + dy1*dy1 + dz*dz);
						if(dxyz == 0)
							dxyz = eps;
						dxz = dz*dz + dx*dx;
						if(dxz == 0)
							dxz = eps;

						if(dxyz < grid.h * FACTOR)
						{
							if(dy1 <= 0)
							{
								if(dx_1 <= 0 && dx_2 > 0)
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,0,-dy1,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),0,dx_2,-dy1,dz1,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,dx_2,-dy1,dz1,dz2);
								}
							}
							else
							{
								if(dx_1 <= 0 && dx_2 > 0)
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle(dens(ix,0,iz),dx_1,0,0,dy1,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,0,iz),0,dx_2,0,dy1,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,0,0,dz1,dz2);
               						tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),0,dx_2,0,dz1,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle(dens(ix,0,iz),dx_1,dx_2,0,dy1,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,0,iz),dx_1,dx_2,0,dz1,dz2);
								}
							}
						}
						else
							tmp_grav = (G_KONST*dens(ix,0,iz)*(grid.h)*(grid.h))*(dz/dxz)*(1-(dyy1/dxyz)); /*semi-infinitely long horizontal rod*/

						tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
						grav = grav + tmp_grav;

						/*back side*/
						dxyz = sqrt(dx*dx + dy2*dy2 + dz*dz);
						if(dxyz == 0)
							dxyz = eps;

						if(dxyz < grid.h * FACTOR)
						{
							if(dy2 >= 0)
							{
								if(dx_1 <= 0 && dx_2 > 0)
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,0,dy2,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),0,dx_2,dy2,dz1,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,dx_2,dy2,dz1,dz2);
								}
							}
							else
							{
								if(dx_1 <= 0 && dx_2 > 0)
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle(dens(ix,grid.ny-1,iz),dx_1,0,dy2,0,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle(dens(ix,grid.ny-1,iz),0,dx_2,dy2,0,dz1,dz2);	

									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,0,0,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),0,dx_2,0,dz1,dz2);
								}
								else
								{/*Two semi-infinite horizontal sheets*/
									tmp_grav = Response_3D_Rectangle(dens(ix,grid.ny-1,iz),dx_1,dx_2,dy2,0,dz1,dz2);
									tmp_grav = tmp_grav + Response_3D_Rectangle_1(dens(ix,grid.ny-1,iz),dx_1,dx_2,0,dz1,dz2);
								}
							}
						}
						else
							tmp_grav = (G_KONST*dens(ix,grid.ny-1,iz)*(grid.h)*(grid.h))*(dz/dxz)*(1-(dyy2/dxyz)); /*semi-infinitely long horizontal rod*/

						tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/
				
						grav = grav + tmp_grav;
					}
				}
			}

			/**********************/
			/*3.Case: Internal cells of the 3-D grid*/
			/*Loop over all remaining cells*/

			if(topo.nr_of_triangles != 0)
			{
                /*The average cutting values of the small cells in one big cell, indexed by y_direction then x_direction.*/
				aver = (double *)memory(NULL,refine*refine,sizeof(double),"ForwardModGrav");
				/*The distance of the averages to the station, indexed by y_direction then x_direction.*/
				aver_dis = (double *)memory(NULL,refine*refine,sizeof(double),"ForwardModGrav");

				for(ix=1;ix<grid.nx-1;ix++)
				{
					x = grid.org[0] + ix*(grid.h);		 /*x component of the cell centers in m*/
					dx = x - xs;						 /*x-component of distance between grid cell and station*/

					for(iy=1;iy<grid.ny-1;iy++)
					{
						y = grid.org[1] + iy*(grid.h);		 /*y component of the cell centers in m*/
						dy = y - ys;						 /*y-component of distance between grid cell and station*/
           
						/*The left-front corner index of the big cell.*/
						corner = (ix-1)*refine*(refine*(grid.ny-2)+1)+(iy-1)*refine;

						for(m=0;m<refine;m++)
						{
							for(n=0;n<refine;n++)
							{ /*The 4 corners index of the smaller cell.*/
								corner1 = corner + m*(refine*(grid.ny-2)+1)+n;
								corner2 = corner1 + 1;
								corner3 = corner1 + (refine*(grid.ny-2)+1);
								corner4 = corner3 + 1;
								if(zgrid_edges[corner1]!=-99999.9 && zgrid_edges[corner2]!=-99999.9 && zgrid_edges[corner3]!=-99999.9 && zgrid_edges[corner4]!=-99999.9)
								{
									/*The distance of the average of the cutting values of the smaller cell to the station.*/
									aver[m*refine+n]=(1.0/4.0)*(zgrid_edges[corner1]+zgrid_edges[corner2]+zgrid_edges[corner3]+zgrid_edges[corner4]);
                                    aver_dis[m*refine+n] = aver[m*refine+n] -zs;
								}
								else
								{
									aver[m*refine+n]=-99999.9;
									aver_dis[m*refine+n]=-99999.9;
								}
							}
						}

						/*Determine the maximun and minimum value of the average values in one big cell.*/
						max_tmp = aver[0];
						for(m=1;m<=(refine*refine-1);m++)
						{
							if(aver[m]>=max_tmp)
								max_tmp = aver[m];
						}
						min_tmp = max_tmp;
						for(m=0;m<=(refine*refine-1);m++)
						{
							if(aver[m]<=min_tmp && aver[m]!=-99999.9)
								min_tmp = aver[m];
						}

						l = 0; /*The frist cut cell.*/
						k = 0; /*The last cut cell.*/
						for(iz=0;iz<grid.nz;iz++)
						{						
							if(grid.org[2]+(iz-0.5)*(grid.h) <= min_tmp && grid.org[2]+(iz+0.5)*(grid.h) > min_tmp)
						     l = iz;
							if(grid.org[2]+(iz-0.5)*(grid.h) <= max_tmp && grid.org[2]+(iz+0.5)*(grid.h) > max_tmp)
							{
							   k = iz;
								break;
							}
						}

						/*Loop before cutting.*/
						for(iz=0;iz<l;iz++)
						{
							z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
							dz = z - zs;					/*z-component of distance between grid cell and station*/

							/*Calculate the gravitation from the cell at the station*/
							dxyz = sqrt(dx*dx + dy*dy + dz*dz);

							/*For cells nearby the stations the response will be calculated exactly (see Nagy, 1966, Geophysics)*/
							if(dxyz < (grid.h * FACTOR))
							{ 

								/********************************************************/
								/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
								if(fabs(dx) >= half_h && fabs(dy) >= half_h)
									tmp_grav = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx +half_h,dy-half_h,dy +half_h,dz-half_h,dz +half_h);
							
								/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
								else if(fabs(dx) < half_h && fabs(dy) >= half_h)
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,dy-half_h,dy +half_h,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,dy-half_h,dy +half_h,dz-half_h,dz +half_h);

									tmp_grav = tmp_grav1 + tmp_grav2;
								}

								/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
								else if(fabs(dx) >= half_h && fabs(dy) < half_h)
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx+half_h,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx+half_h,0,dy+half_h,dz-half_h,dz +half_h);
	
									tmp_grav = tmp_grav1 + tmp_grav2;
								}

								/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
								else
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,0,dy+half_h,dz-half_h,dz +half_h);

									tmp_grav3 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav4 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,0,dy+half_h,dz-half_h,dz +half_h);
	
									tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
								}
								/********************************************************/

							}
							else
							/*For cells that are for away from the station (approximated like suggested from Trinks and Hobbs)*/
								tmp_grav = (G_KONST*dens(ix,iy,iz)*(grid.h)*(grid.h)*(grid.h))*(dz/(dxyz*dxyz*dxyz)); /*3D-voxels*/

							tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

							grav = grav + tmp_grav;

						}

						/*Loop the cells being cut.*/
						/*Just one cell is cut.*/
						if(l == k)
						{
							z  = grid.org[2] + l*(grid.h);  /*z component of the cell centers in m*/
							dz = z - zs;					/*z-component of distance between grid cell and station*/

							for(m=0;m<refine;m++)
							{
								for(n=0;n<refine;n++)
								{ /*The smaller cell is divided into 2 parts.*/
									dx1=dx-half_h+half_h_refine+m*(grid.h)/refine;
									dy1=dy-half_h+half_h_refine+n*(grid.h)/refine;
									/*Calculate the gravitation from the cell at the station*/
									dxyz1 = sqrt(dx1*dx1 + dy1*dy1 + dz*dz);

									if(aver[m*refine+n] != -99999.9)
									{
										/********************************************************/
										/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
										if(fabs(dx1) >= half_h_refine && fabs(dy1) >= half_h_refine)
										{
											tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,dx1 +half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav = tmp_grav+Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,dx1 +half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,aver_dis[m*refine+n],dz +half_h);
										}
										/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
										else if(fabs(dx1) < half_h_refine && fabs(dy1) >= half_h_refine)
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											/*The upper part is water.*/
											tmp_grav1 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,0,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav1 = tmp_grav1+Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,0,dy1-half_h_refine,dy1 +half_h_refine,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav2 = Response_3D_Rectangle(grid.dens_air_water,0,dx1+half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav2 = tmp_grav2+Response_3D_Rectangle(dens(ix,iy,l),0,dx1+half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,aver_dis[m*refine+n],dz +half_h);

										    tmp_grav = tmp_grav1 + tmp_grav2;
										}

										/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
										else if(fabs(dx1) >= half_h_refine && fabs(dy1) < half_h_refine)
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											/*The upper part is water.*/
											tmp_grav1 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav1 = tmp_grav1+Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,dx1+half_h_refine,dy1-half_h_refine,0,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav2 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav2 = tmp_grav2+Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,dx1+half_h_refine,0,dy1+half_h_refine,aver_dis[m*refine+n],dz +half_h);
		
											tmp_grav = tmp_grav1 + tmp_grav2;
										}

										/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
										else
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											/*The upper part is water.*/
											tmp_grav1 = Response_3D_Rectangle(grid.dens_air_water,0,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav1 = tmp_grav1+Response_3D_Rectangle(dens(ix,iy,l),0,dx1+half_h_refine,dy1-half_h_refine,0,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav2 = Response_3D_Rectangle(grid.dens_air_water,0,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav2 = tmp_grav2+Response_3D_Rectangle(dens(ix,iy,l),0,dx1+half_h_refine,0,dy1+half_h_refine,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav3 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,0,dy-half_h_refine,0,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav3 = tmp_grav3+Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,0,dy-half_h_refine,0,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav4 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,0,0,dy1+half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav4 = tmp_grav4+Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,0,0,dy1+half_h_refine,aver_dis[m*refine+n],dz +half_h);
	
											tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
										}
									}
									else
									{
										/********************************************************/
										/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
										if(fabs(dx1) >= half_h_refine && fabs(dy1) >= half_h_refine)
											tmp_grav = Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,dx1 +half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);
						
										/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
										else if(fabs(dx1) < half_h_refine && fabs(dy1) >= half_h_refine)
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,0,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);
											tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,l),0,dx1+half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);

										    tmp_grav = tmp_grav1 + tmp_grav2;
										}


										/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
										else if(fabs(dx1) >= half_h_refine && fabs(dy1) < half_h_refine)
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,dz +half_h);
											tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,dz +half_h);
		
											tmp_grav = tmp_grav1 + tmp_grav2;
										}

										/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
										else
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,l),0,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,dz +half_h);
											tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,l),0,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,dz +half_h);

											tmp_grav3 = Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,0,dy-half_h_refine,0,dz-half_h,dz +half_h);
											tmp_grav4 = Response_3D_Rectangle(dens(ix,iy,l),dx1-half_h_refine,0,0,dy1+half_h_refine,dz-half_h,dz +half_h);
	
											tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
										}
									}
									/********************************************************/

									tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

									grav = grav + tmp_grav;
								}
							}
						}
						else
						{
							for(m=0;m<refine;m++)
							{
								for(n=0;n<refine;n++)
								{
									dx1=dx-half_h+half_h_refine+m*(grid.h)/refine;
									dy1=dy-half_h+half_h_refine+n*(grid.h)/refine;
									/*Calculate the gravitation from the cell at the station*/

									if(aver[m*refine+n] != -99999.9)
									{
										p = 0;
										for(iz=l;iz<=k;iz++)
										{
											if(aver[m*refine+n]>=grid.org[2]+(iz-0.5)*(grid.h) && aver[m*refine+n]<grid.org[2]+(iz+0.5)*(grid.h))
											{
												p=iz;
												break;
											}
										}
										for(iz=l;iz<p;iz++)
										{
	  										z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
											dz = z - zs;					/*z-component of distance between grid cell and station*/

											/********************************************************/
											/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
											if(fabs(dx1) >= half_h_refine && fabs(dy1) >= half_h_refine)
												tmp_grav = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1 +half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);
						
											/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
											else if(fabs(dx1) < half_h_refine && fabs(dy1) >= half_h_refine)
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);

											 tmp_grav = tmp_grav1 + tmp_grav2;
											}


											/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
											else if(fabs(dx1) >= half_h_refine && fabs(dy1) < half_h_refine)
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,dz +half_h);
		
												tmp_grav = tmp_grav1 + tmp_grav2;
											}

											/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
											else
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,dz +half_h);

												tmp_grav3 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,dy-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav4 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,0,dy1+half_h_refine,dz-half_h,dz +half_h);
	
												tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
											}
											tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

											grav = grav + tmp_grav;
										}

										z  = grid.org[2] + p*(grid.h);  /*z component of the cell centers in m*/
										dz = z - zs;					/*z-component of distance between grid cell and station*/

										/********************************************************/
										/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
										if(fabs(dx1) >= half_h_refine && fabs(dy1) >= half_h_refine)
										{/*The upper part is water.*/
											tmp_grav = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,dx1 +half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav = tmp_grav+Response_3D_Rectangle(dens(ix,iy,p),dx1-half_h_refine,dx1 +half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,aver_dis[m*refine+n],dz +half_h);
										}
										/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
										else if(fabs(dx1) < half_h_refine && fabs(dy1) >= half_h_refine)
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											/*The upper part is water.*/
											tmp_grav1 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,0,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav1 = tmp_grav1+Response_3D_Rectangle(dens(ix,iy,p),dx1-half_h_refine,0,dy1-half_h_refine,dy1 +half_h_refine,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav2 = Response_3D_Rectangle(grid.dens_air_water,0,dx1+half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav2 = tmp_grav2+Response_3D_Rectangle(dens(ix,iy,p),0,dx1+half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,aver_dis[m*refine+n],dz +half_h);

										    tmp_grav = tmp_grav1 + tmp_grav2;
										}


										/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
										else if(fabs(dx1) >= half_h_refine && fabs(dy1) < half_h_refine)
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											/*The upper part is water.*/
											tmp_grav1 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav1 = tmp_grav1+Response_3D_Rectangle(dens(ix,iy,p),dx1-half_h_refine,dx1+half_h_refine,dy1-half_h_refine,0,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav2 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav2 = tmp_grav2+Response_3D_Rectangle(dens(ix,iy,p),dx1-half_h_refine,dx1+half_h_refine,0,dy1+half_h_refine,aver_dis[m*refine+n],dz +half_h);
		
											tmp_grav = tmp_grav1 + tmp_grav2;
										}

											/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
										else
										{
											/*The integral have to be splitted (see Nagy, 1966)*/
											/*The upper part is water.*/
											tmp_grav1 = Response_3D_Rectangle(grid.dens_air_water,0,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav1 = tmp_grav1+Response_3D_Rectangle(dens(ix,iy,p),0,dx1+half_h_refine,dy1-half_h_refine,0,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav2 = Response_3D_Rectangle(grid.dens_air_water,0,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav2 = tmp_grav2+Response_3D_Rectangle(dens(ix,iy,p),0,dx1+half_h_refine,0,dy1+half_h_refine,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav3 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,0,dy-half_h_refine,0,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav3 = tmp_grav3+Response_3D_Rectangle(dens(ix,iy,p),dx1-half_h_refine,0,dy-half_h_refine,0,aver_dis[m*refine+n],dz +half_h);
											/*The upper part is water.*/
											tmp_grav4 = Response_3D_Rectangle(grid.dens_air_water,dx1-half_h_refine,0,0,dy1+half_h_refine,dz-half_h,aver_dis[m*refine+n]);
											tmp_grav4 = tmp_grav4+Response_3D_Rectangle(dens(ix,iy,p),dx1-half_h_refine,0,0,dy1+half_h_refine,aver_dis[m*refine+n],dz +half_h);

											tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
										}
										tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

										grav = grav + tmp_grav;

										for(iz=p+1;iz<=k;iz++)
										{
	  										z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
											dz = z - zs;					/*z-component of distance between grid cell and station*/

											/********************************************************/
											/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
											if(fabs(dx1) >= half_h_refine && fabs(dy1) >= half_h_refine)
												tmp_grav = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1 +half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);
						
											/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
											else if(fabs(dx1) < half_h_refine && fabs(dy1) >= half_h_refine)
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);

											 tmp_grav = tmp_grav1 + tmp_grav2;
											}


											/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
											else if(fabs(dx1) >= half_h_refine && fabs(dy1) < half_h_refine)
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,dz +half_h);
		
												tmp_grav = tmp_grav1 + tmp_grav2;
											}

											/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
											else
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,dz +half_h);

												tmp_grav3 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,dy-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav4 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,0,dy1+half_h_refine,dz-half_h,dz +half_h);
	
												tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
											}
											tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

											grav = grav + tmp_grav;
										}
									}
									else
									{ /*Loop the cells being cut.*/
										for(iz=l;iz<=k;iz++)
										{
	  										z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
											dz = z - zs;					/*z-component of distance between grid cell and station*/

											/********************************************************/
											/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
											if(fabs(dx1) >= half_h_refine && fabs(dy1) >= half_h_refine)
												tmp_grav = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1 +half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);
						
											/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
											else if(fabs(dx1) < half_h_refine && fabs(dy1) >= half_h_refine)
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,dy1-half_h_refine,dy1 +half_h_refine,dz-half_h,dz +half_h);

											 tmp_grav = tmp_grav1 + tmp_grav2;
											}


											/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
											else if(fabs(dx1) >= half_h_refine && fabs(dy1) < half_h_refine)
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,dz +half_h);
		
												tmp_grav = tmp_grav1 + tmp_grav2;
											}

											/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
											else
											{
												/*The integral have to be splitted (see Nagy, 1966)*/
												tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,dy1-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx1+half_h_refine,0,dy1+half_h_refine,dz-half_h,dz +half_h);

												tmp_grav3 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,dy-half_h_refine,0,dz-half_h,dz +half_h);
												tmp_grav4 = Response_3D_Rectangle(dens(ix,iy,iz),dx1-half_h_refine,0,0,dy1+half_h_refine,dz-half_h,dz +half_h);
	
												tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
											}
											tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

											grav = grav + tmp_grav;
										}
									}
								}
							}
						}

						/*Loop the cells after cut.*/
						for(iz=k+1;iz<grid.nz;iz++)
						{
							z  = grid.org[2] + iz*(grid.h);  /*z component of the cell centers in m*/
							dz = z - zs;					/*z-component of distance between grid cell and station*/

							/*Calculate the gravitation from the cell at the station*/
							dxyz = sqrt(dx*dx + dy*dy + dz*dz);

							/*For cells nearby the stations the response will be calculated exactly (see Nagy, 1966, Geophysics)*/
							if(dxyz < (grid.h * FACTOR))
							{ 

								/********************************************************/
								/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
								if(fabs(dx) >= half_h && fabs(dy) >= half_h)
									tmp_grav = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx +half_h,dy-half_h,dy +half_h,dz-half_h,dz +half_h);
							
								/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
								else if(fabs(dx) < half_h && fabs(dy) >= half_h)
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,dy-half_h,dy +half_h,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,dy-half_h,dy +half_h,dz-half_h,dz +half_h);

									tmp_grav = tmp_grav1 + tmp_grav2;
								}


								/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
								else if(fabs(dx) >= half_h && fabs(dy) < half_h)
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx+half_h,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx+half_h,0,dy+half_h,dz-half_h,dz +half_h);
	
									tmp_grav = tmp_grav1 + tmp_grav2;
								}

								/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
								else
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,0,dy+half_h,dz-half_h,dz +half_h);

									tmp_grav3 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav4 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,0,dy+half_h,dz-half_h,dz +half_h);
	
									tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
								}
								/********************************************************/


							}
							else
							/*For cells that are for away from the station (approximated like suggested from Trinks and Hobbs)*/
								tmp_grav = (G_KONST*dens(ix,iy,iz)*(grid.h)*(grid.h)*(grid.h))*(dz/(dxyz*dxyz*dxyz)); /*3D-voxels*/

							tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

							grav = grav + tmp_grav;

						}
					}
				}
				free(aver);
				free(aver_dis);
			}
			else
			{
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
								/*Case I: The X and Y coordinates of the station location are not located within the edges of the rectangle*/
								if(fabs(dx) >= half_h && fabs(dy) >= half_h)
									tmp_grav = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx +half_h,dy-half_h,dy +half_h,dz-half_h,dz +half_h);
							
								/*Case II: Only the X coordinates of the station location are located within the edges of the rectangle*/
								else if(fabs(dx) < half_h && fabs(dy) >= half_h)
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,dy-half_h,dy +half_h,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,dy-half_h,dy +half_h,dz-half_h,dz +half_h);

									tmp_grav = tmp_grav1 + tmp_grav2;
								}


								/*Case III: Only the Y coordinates of the station location are located within the edges of the rectangle*/
								else if(fabs(dx) >= half_h && fabs(dy) < half_h)
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx+half_h,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,dx+half_h,0,dy+half_h,dz-half_h,dz +half_h);
	
									tmp_grav = tmp_grav1 + tmp_grav2;
								}

								/*Case IV: Both the X AND Y coordinates of the station location are located within the edges of the rectangle*/
								else
								{
									/*The integral have to be splitted (see Nagy, 1966)*/
									tmp_grav1 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav2 = Response_3D_Rectangle(dens(ix,iy,iz),0,dx+half_h,0,dy+half_h,dz-half_h,dz +half_h);

									tmp_grav3 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,dy-half_h,0,dz-half_h,dz +half_h);
									tmp_grav4 = Response_3D_Rectangle(dens(ix,iy,iz),dx-half_h,0,0,dy+half_h,dz-half_h,dz +half_h);
	
									tmp_grav = tmp_grav1 + tmp_grav2 + tmp_grav3 + tmp_grav4;
								}
								/********************************************************/


							}
							else
							/*For cells that are for away from the station (approximated like suggested from Trinks and Hobbs)*/
								tmp_grav = (G_KONST*dens(ix,iy,iz)*(grid.h)*(grid.h)*(grid.h))*(dz/(dxyz*dxyz*dxyz)); /*3D-voxels*/

							tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

							grav = grav + tmp_grav;

						}
					}
				}
			}


			/*Assign the calculated density value to the corresponding measurements*/
			for(j=0;j<data->ndata_grav;j++)
			{
				if(data->gravs[i] == data->gno[j])
					data->calc_grav[j] = grav;
			}

			if(i%10 == 0)
				printf("Modeling for %d of %d stations is completed\n", i+1, geo.nstat_grav);


		}
	    
		free(xgrid_edges);
	    free(ygrid_edges);
	    free(zgrid_edges);


	}
	#undef dens 

    #undef Xgrid
    #undef Ygrid
    #undef Zgrid



	/*************************************/
	/*Remove a mean value from the calculated traveltimes*/
	//if(inv.ninv == 0)
	//	RemoveMeanGrav(data);

	/*Add the effect of the extra cell to the data; BJOERN_MOD5*/
	AddEffectExtraCell(inv.dens_extra_cell, data, grid);

	/****************************************************************************/
	/*Write out the calculated data:*/
	out = fopen("calc_grav.txt","w");
	for(i=0;i<data->ndata_grav;i++)
		fprintf(out,"%d %f\n",data->gno[i],data->calc_grav[i]);
	fclose(out);

	/*************************************/

	/*Determine the the RMS-value in mgal*/

	sum = 0.0;

	for(i=0;i<data->ndata_grav;i++)
		sum = sum + (data->obs_grav[i] - data->calc_grav[i])*(data->obs_grav[i] - data->calc_grav[i]);

	if(data->ndata_grav != 0)
		data->rms_grav = sqrt(sum/data->ndata_grav);
	else
		data->rms_grav = -99999.9;

	printf("Gravity forward modeling is finished:\n");
	printf("RMS-values: %10.5f mgal\n",data->rms_grav);
	printf("----------------\n\n\n");

	return(1);
}



#undef G_KONST
#undef FACTOR
#undef refine
/*-------------------------------------------------------------*/
/*Calculate accurate gravity response for a rectangle (2-D) block (see Heiland,1946, Geophysical Exploration)*/
/*Parameter:	dens       := Density of the cell*/
/*              x1,z1   := Distance from the (left,top) planes of the prsim to the considered gravity station*/
/*              x2,z2   := Distance from the (right,bottom) planes of the prsim to the considered gravity station*/
/* Remark:	- The output is the gravity response for the prism at the considered station*/

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/
#define EPS 1.0E-15

double Response_2D_Rectangle(double dens, double x1, double x2, double z1, double z2)
{
	double term_sq[4];
	double term_ln, term_arctan;
	double response;

	if((x1*x1 + z1*z1) != 0)
		term_sq[0] =(x1*x1 + z1*z1);
	else
		term_sq[0] =(x1*x1 + z1*z1) + EPS;

	if((x1*x1 + z2*z2) != 0)
		term_sq[1] =(x1*x1 + z2*z2);
	else
		term_sq[1] =(x1*x1 + z2*z2) + EPS;

	if((x2*x2 + z1*z1) != 0)
		term_sq[2] =(x2*x2 + z1*z1);
	else
		term_sq[2] =(x2*x2 + z1*z1) + EPS;

	if((x2*x2 + z2*z2) != 0)
		term_sq[3] =(x2*x2 + z2*z2);
	else
		term_sq[3] =(x2*x2 + z2*z2) + EPS;

	/*Calculate the response of the rectangular block in 2-D*/
	term_ln = x2*log(sqrt(term_sq[3]/term_sq[2])) - x1*log(sqrt(term_sq[1]/term_sq[0]));
	term_arctan = z2*(atan2(x2,z2) - atan2(x1,z2)) + z1*(atan2(x1,z1) - atan2(x2,z1));

    response = 2*G_KONST*dens*(term_ln + term_arctan);

	return(response);
}

#undef G_KONST
#undef EPS


/*-------------------------------------------------------------*/
/*Calculate accurate gravity response for a rectangle (3-D) prisms (see Nagy,1966, Geophysics)*/
/*Parameter:	dens       := Density of the cell*/
/*              x1,y1,z1   := Distance from the (left,front,top) planes of the prism to the considered gravity station*/
/*              x2,y2,z2   := Distance from the (right,back,bottom) planes of the prism to the considered gravity station*/
/* Remark:	- The output is the gravity response for the prism at the considered station*/
/*			- The x1 and x2(y1 and y2) have to have both >=0.0 or <=0.0 !!! */

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/
#define EPS 1.0E-15
#define EPS_2 1.0E-10

double Response_3D_Rectangle(double dens, double x1, double x2, double y1, double y2, double z1, double z2)
{
	long i,j,k,sign;
	double response, tmp_value;
	double x[2],y[2],z[2];
	double dist_r[2][2][2]; /*Distance from the station to the corners of the prism*/
	double term_ln;			/*The logarithmic terms*/
	double term_arcsin;	/*The arcsin terms*/

	/*Set the edges of the cells*/
	/*********************************************************/
	/*Remark: all x and y values have to be MADE positive!!*/

	if(EPS_2 > fabs(x1))
		x1 = 0.0;
	if(EPS_2 > fabs(x2))
		x2 = 0.0;
	if(EPS_2 > fabs(y1))
		y1 = 0.0;
	if(EPS_2 > fabs(y2))
		y2 = 0.0;


	if(x1 >= 0 && x2 >= 0)
	{
		x[0] = x1; /*left*/
		x[1] = x2; /*right*/
	}
	else if(x1 <= 0 && x2 <= 0)
	{
		x[0] = -x2; /*left*/
		x[1] = -x1; /*right*/
	}
	else
	{
		printf("There is an error for calculating the x-borders of prisms in the gravity forward modelling\n");
		exit(0);
	}

	if(y1 >= 0 && y2 >= 0)
	{
		y[0] = y1; /*front*/
		y[1] = y2; /*back*/
	}
	else if(y1 <= 0 && y2 <= 0)
	{
		y[0] = -y2; /*front*/
		y[1] = -y1; /*back*/
	}
	else
	{
		printf("There is an error for calculating the y-borders of prisms in the gravity forward modelling\n");
		exit(0);
	}

	z[0] = z1; /*top*/
	z[1] = z2; /*bottom*/
	/*********************************************************/

	/*Calculate the distances from the gravity station to the corners of the prism:*/
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
			{
				dist_r[i][j][k] = sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k]);
			}

	term_ln = 0.0;
	term_arcsin = 0.0;

	for(i=0;i<2;i++)		/*Loop over the x-index*/
		for(j=0;j<2;j++)	/*loop over the y-index*/
			for(k=0;k<2;k++) /*Loop over the z-index*/
			{
				if((i==j && k==1)||(i!=j && k==0)) /*Determine the sign of the terms*/
					sign = 1;
				else
					sign = -1;

				/*Calculate the logarithmic terms*/
				if((dist_r[i][j][k] + y[j]) != 0)
					term_ln = (-sign*x[i]*log(dist_r[i][j][k] + y[j])) + term_ln;
				else
					term_ln = (-sign*x[i]*log(EPS)) + term_ln;

				if((dist_r[i][j][k] + x[i]) != 0)
					term_ln = (-sign*y[j]*log(dist_r[i][j][k] + x[i])) + term_ln;
				else
					term_ln = (-sign*y[j]*log(EPS)) + term_ln;

				/*Calculate the arcsin terms*/
				if(((y[j]*y[j])+(z[k]*z[k])) != 0)
					tmp_value = ((y[j]*y[j])+(z[k]*z[k])+ (y[j]*dist_r[i][j][k]))/((y[j] + dist_r[i][j][k])*sqrt((y[j]*y[j])+(z[k]*z[k])));
				else
					tmp_value = ((y[j]*y[j])+(z[k]*z[k])+ (y[j]*dist_r[i][j][k]))/EPS;


				if((fabs(tmp_value)) < 1)
					term_arcsin = (sign*z[k]*asin(tmp_value)) + term_arcsin;
				else if(tmp_value >= 1)
					term_arcsin = (sign*z[k]*(PI/2)) + term_arcsin;
				else 
					term_arcsin = (sign*z[k]*(-PI/2)) + term_arcsin;
			}


	/*Calculate the response for the prism*/
	response = G_KONST*dens*(term_ln + term_arcsin);

	return(response);
}

#undef G_KONST
#undef EPS
#undef EPS_2


/////not tested, 31. 08. 2006!!!!!
/*-------------------------------------------------------------*/
/*Calculate accurate gravity response for the border cells in 3D.*/
/*In y- (or x-)direction, the cells are infinity long. i.e., y (or x) belong to [y1, +infinity] (or [x1, +infinity]). */
/*Parameter:	dens       := Density of the cell*/
/*              x1,y1,z1   := Distance from the (left,front,top) planes of the prism to the considered gravity station*/
/*              x2   ,z2   := Distance from the (right     ,bottom) planes of the prism to the considered gravity station*/
/* Remark:	- The output is the gravity response for the prism at the considered station*/
/*			- The x1 and x2 have to have both >=0.0 or <=0.0 !!! */

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/
#define EPS 1.0E-15
#define EPS_2 1.0E-10

double Response_3D_Rectangle_1(double dens, double x1, double x2, double y1, double z1, double z2)
{
	long i,k,sign;
	double response, tmp_value;
	double x[2],y,z[2];
	double dist_r[2][2]; /*Distance from the station to the corners of the prism*/
	double term_ln;			/*The logarithmic terms*/
	double term_arcsin;	/*The arcsin terms*/

	/*Set the edges of the cells*/
	/*********************************************************/
	/*Remark: all x and y values have to be MADE positive!!*/

	if(EPS_2 > fabs(x1))
		x1 = 0.0;
	if(EPS_2 > fabs(x2))
		x2 = 0.0;
	if(EPS_2 > fabs(y1))
		y1 = 0.0;
	
	if(x1 >= 0 && x2 >= 0)
	{
		x[0] = x1; /*left*/
		x[1] = x2; /*right*/
	}
	else if(x1 <= 0 && x2 <= 0)
	{
		x[0] = -x2; /*left*/
		x[1] = -x1; /*right*/
	}
	else
	{
		printf("There is an error for calculating the x-borders of prisms in the gravity forward modelling\n");
		exit(0);
	}

	if(y1 >= 0)
		y = y1; /*front*/
	else
		y = -y1; /*front*/

	z[0] = z1; /*top*/
	z[1] = z2; /*bottom*/
	/*********************************************************/

	/*Calculate the distances from the gravity station to the corners of the prism:*/
	for(i=0;i<2;i++)
			for(k=0;k<2;k++)
			{
				dist_r[i][k] = sqrt(x[i]*x[i] + y*y + z[k]*z[k]);
			}

	term_ln = 0.0;
	term_arcsin = 0.0;

	for(i=0;i<2;i++)		/*Loop over the x-index*/
			for(k=0;k<2;k++) /*Loop over the z-index*/
			{
				if((i==0 && k==1)||(i==1 && k==0)) /*Determine the sign of the terms*/
					sign = -1;
				else
					sign = 1;

				/*Calculate the logarithmic terms*/
				if((dist_r[i][k] + y) != 0)
					term_ln = (sign*x[i]*log(dist_r[i][k] + y)) + term_ln;
				else
					term_ln = (sign*x[i]*log(EPS)) + term_ln;

				if((dist_r[i][k] + x[i]) != 0)
					term_ln = (sign*y*log(dist_r[i][k] + x[i])) + term_ln;
				else
					term_ln = (sign*y*log(EPS)) + term_ln;

				/*Calculate the arcsin terms*/
				if(((y*y)+(z[k]*z[k])) != 0)
					tmp_value = ((y*y)+(z[k]*z[k])+ (y*dist_r[i][k]))/((y + dist_r[i][k])*sqrt((y*y)+(z[k]*z[k])));
				else
					tmp_value = ((y*y)+(z[k]*z[k])+ (y*dist_r[i][k]))/EPS;


				if((fabs(tmp_value)) < 1)
					term_arcsin = (-sign*z[k]*asin(tmp_value)) + term_arcsin;
				else if(tmp_value >= 1)
					term_arcsin = (-sign*z[k]*(PI/2)) + term_arcsin;
				else 
					term_arcsin = (-sign*z[k]*(-PI/2)) + term_arcsin;
			}


	/*Calculate the response for the prism*/
	response = G_KONST*dens*(term_ln + term_arcsin);

	return(response);
}

#undef G_KONST
#undef EPS
#undef EPS_2

/*-------------------------------------------------------------*/
/*Changed by Jin at 31.08.2006*/
/*Calculate accurate gravity response for the corner cells in 3D.*/
/*In y- and x-direction, the cells are infinity long. i.e., y and x belong to [y1, +infinity] and [x1, +infinity]. */
/*Parameter:	dens       := Density of the cell*/
/*              x1,y1,z1   := Distance from the (left,front,top) planes of the prism to the considered gravity station*/
/*                    z2   := Distance from the (          ,bottom) planes of the prism to the considered gravity station*/
/* Remark:	- The output is the gravity response for the prism at the considered station*/
			

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/
#define EPS 1.0E-15
#define EPS_2 1.0E-10

double Response_3D_Rectangle_2(double dens, double x1, double y1, double z1, double z2)
{
	long k,sign;
	double response, tmp_value,tmp_value_1;
	double x,y,z[2];
	double dist_r[2]; /*Distance from the station to the corners of the prism*/
	double term_ln;			/*The logarithmic terms*/
	double term_arcsin;	/*The arcsin terms*/

	/*Set the edges of the cells*/
	/*********************************************************/
	/*Remark: all x and y values have to be MADE positive!!*/

	if(EPS_2 > fabs(x1))
		x1 = 0.0;
	if(EPS_2 > fabs(y1))
		y1 = 0.0;
	
	if(x1 >= 0)
		x = x1; /*left*/
	else 
		x = -x1; /*left*/

	if(y1 >= 0)
		y = y1; /*front*/
	else
		y = -y1; /*front*/

	z[0] = z1; /*top*/
	z[1] = z2; /*bottom*/
	/*********************************************************/

	/*Calculate the distances from the gravity station to the corners of the prism:*/
	for(k=0;k<2;k++)
		dist_r[k] = sqrt(x*x + y*y + z[k]*z[k]);

	term_ln = 0.0;
	term_arcsin = 0.0;

	for(k=0;k<2;k++) /*Loop over the z-index*/
		{
			if(k==1) /*Determine the sign of the terms*/
				sign = -1;
			else
				sign = 1;

			/*Calculate the logarithmic terms*/
			if((dist_r[k] + y) != 0)
				term_ln = (sign*x*log(dist_r[k] + y)) + term_ln;
			else
				term_ln = (sign*x*log(EPS)) + term_ln;
			if((dist_r[k] + x) != 0)
				term_ln = (sign*y*log(dist_r[k] + x)) + term_ln;
			else
				term_ln = (sign*y*log(EPS)) + term_ln;

			/*Calculate the arcsin terms*/
			if(((y*y)+(z[k]*z[k])) != 0)
			{
				tmp_value = ((y*y)+(z[k]*z[k])+ (y*dist_r[k]))/((y + dist_r[k])*sqrt((y*y)+(z[k]*z[k])));
				tmp_value_1 = y/sqrt(y*y+z[k]*z[k]);
			}
			else
			{
				tmp_value = ((y*y)+(z[k]*z[k])+ (y*dist_r[k]))/EPS;
				tmp_value_1 = y/EPS;
			}
			if((fabs(tmp_value)) < 1 && (fabs(tmp_value_1)) < 1)
			{
				term_arcsin = (-sign*z[k]*asin(tmp_value)) + term_arcsin;
				term_arcsin = (sign*z[k]*asin(tmp_value_1)) + term_arcsin;
			}
			else if((fabs(tmp_value)) < 1 && tmp_value_1 >= 1)
			{
				term_arcsin = (-sign*z[k]*asin(tmp_value)) + term_arcsin;
				term_arcsin = (sign*z[k]*(PI/2)) + term_arcsin;
			}
			else if(tmp_value >= 1 && (fabs(tmp_value_1)) < 1)
			{
				term_arcsin = (-sign*z[k]*(PI/2)) + term_arcsin;
				term_arcsin = (sign*z[k]*asin(tmp_value_1)) + term_arcsin;
				
			}
			else if(tmp_value >= 1 && tmp_value_1 >= 1)
			{
				term_arcsin = (-sign*z[k]*(PI/2)) + term_arcsin;
				term_arcsin = (sign*z[k]*(PI/2)) + term_arcsin;
			}
			else 
			{
				term_arcsin = (-sign*z[k]*(-PI/2)) + term_arcsin;
				term_arcsin = (sign*z[k]*(-PI/2)) + term_arcsin;
			}
		}


	/*Calculate the response for the prism*/
	response = G_KONST*dens*(term_ln + term_arcsin);

	return(response);
}

#undef G_KONST
#undef EPS
#undef EPS_2

/*-------------------------------------------------------------*/
/*Remove a mean value (D_calc - D_obs) from the calculated gravity data*/
/*Parameter:	data = Data structure*/

int RemoveMeanGrav(DATA_STRUCT *data)
{
	long i;
	double diff_grav;

	diff_grav = 0.0;

	for(i=0;i<data->ndata_grav;i++)
	{
		diff_grav = diff_grav + (data->calc_grav[i] - data->obs_grav[i]);
	}

	diff_grav = diff_grav/data->ndata_grav;

	for(i=0;i<data->ndata_grav;i++)
		data->calc_grav[i] = data->calc_grav[i] - diff_grav;

	return(0);
}


/*-------------------------------------------------------------*/
/*BJOERN_MOD5*/
/*Add the effect from the extra cell to the calculated attraction*/
/*Parameter:	data  := Data structure*/
/*    dens_extra_cell := density of the extra cell*/
/*				grid  := Grid structure*/

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/

int AddEffectExtraCell(double dens_extra_cell, DATA_STRUCT *data, GRID_STRUCT grid)
{
	long i;
	double tmp_grav;

	/*Calculate the attraction of the extra cell*/
	tmp_grav = 2*PI*G_KONST*dens_extra_cell*(grid.h * grid.nz); /*infinite sheet*/
	tmp_grav = 1.0E8 * tmp_grav; /*Transform the calculations to mGal*/

	/*Add the effect of the rubbish bin cell to the gravitational attraction*/
	for(i=0;i<data->ndata_grav;i++)
		data->calc_grav[i] = data->calc_grav[i] + tmp_grav;

	return(1);
}

#undef G_KONST

/*-------------------------------------------------------------*/
/*Calculate the gravity response for a lower triangle (2-D) block with corners(x1,z1),(x1,z2),(x2,z2)*/
/*Parameter:	dens       := Density of the cell*/
/*              x1,z1   := Distance from the (left,top) corner of the triangle to the considered gravity station*/
/*              x2,z2   := Distance from the (right,bottom) corner of the triangle to the considered gravity station*/
/* Remark:	- The output is the gravity response for the triangle at the considered station*/

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/

double Response_2D_lower_triangle(double dens, double x1, double x2, double z1, double z2)
{
	double a,b,c,d,e,f,g,first_integ,second_integ;
	double response;

	if((x1 == x2) || (z1 == z2))
		response = 0;
	else
	{
    	/*Calculate the response of the rectangular block in 2-D*/
	    a=((x2-x1)*(x2-x1)+(z2-z1)*(z2-z1))/((x2-x1)*(x2-x1));
	    b=2*(z2-z1)*(x2*z1-x1*z2)/((x2-x1)*(x2-x1));
	    c=(z1-x1*(z2-z1)/(x2-x1))*(z1-x1*(z2-z1)/(x2-x1));

	    d=x2*log(a*x2*x2+b*x2+c)-x1*log(a*x1*x1+b*x1+c)-2*(x2-x1);
	    e=b/(2*a)*(log(a*x2*x2+b*x2+c)-log(a*x1*x1+b*x1+c));
	    f=(4*a*c-b*b);


	    first_integ=x2*log(x2*x2+z2*z2)-x1*log(x1*x1+z2*z2)-2*(x2-x1)+2*z2*atan2(x2,z2)-2*z2*atan2(x1,z2);
	

	    if(f > 0)
	     	g=sqrt(f)/a*(atan2(2*a*x2,sqrt(f))-atan2(2*a*x1,sqrt(f)));
	    else
		    if(f<0)
		        g=-sqrt(-f)/a*(log((2*a*x2+b-sqrt(-f))/(2*a*x2+b+sqrt(-f)))-log((2*a*x1+b-sqrt(-f))/(2*a*x1+b+sqrt(-f))));
		    else
		    	g=0;

	    second_integ=d+e+g;
	
	
	    response = G_KONST*dens*(first_integ-second_integ);
	}

	return(response);

#undef G_KONST

}



/*-------------------------------------------------------------*/
/*Calculate the gravity response for a upper triangle (2-D) block with ocrners(x1,z1),(x2,z1),(x2,z2)*/
/*Parameter:	dens       := Density of the cell*/
/*              x1,z1   := Distance from the (left,top) corner of the triangle to the considered gravity station*/
/*              x2,z2   := Distance from the (right,bottom) corner of the triangle to the considered gravity station*/
/* Remark:	- The output is the gravity response for the triangle at the considered station*/

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/

double Response_2D_upper_triangle(double dens, double x1, double x2, double z1, double z2)
{
	double a,b,c,d,e,f,g,first_integ,second_integ;
	double response;

	if((x1 == x2) || (z1 == z2))
		response = 0;
	else
	{
	
	    /*Calculate the response of the rectangular block in 2-D*/
	    a=((x2-x1)*(x2-x1)+(z2-z1)*(z2-z1))/((x2-x1)*(x2-x1));
	    b=2*(z2-z1)*(x2*z1-x1*z2)/((x2-x1)*(x2-x1));
	    c=(z1-x1*(z2-z1)/(x2-x1))*(z1-x1*(z2-z1)/(x2-x1));

	    d=x2*log(a*x2*x2+b*x2+c)-x1*log(a*x1*x1+b*x1+c)-2*(x2-x1);
	    e=b/(2*a)*(log(a*x2*x2+b*x2+c)-log(a*x1*x1+b*x1+c));
	    f=(4*a*c-b*b);


	    first_integ=x2*log(x2*x2+z1*z1)-x1*log(x1*x1+z1*z1)-2*(x2-x1)+2*z1*atan2(x2,z1)-2*z1*atan2(x1,z1);
	

	    if(f > 0)
		    g=sqrt(f)/a*(atan2(2*a*x2,sqrt(f))-atan2(2*a*x1,sqrt(f)));
	    else
		    if(f<0)
		         g=-sqrt(-f)/a*(log((2*a*x2+b-sqrt(-f))/(2*a*x2+b+sqrt(-f)))-log((2*a*x1+b-sqrt(-f))/(2*a*x1+b+sqrt(-f))));
		    else
			     g=0;

	    second_integ=d+e+g;
	
	
	    response = G_KONST*dens*(-first_integ+second_integ);
	}
	return(response);

#undef G_KONST

}

/*-------------------------------------------------------------*/
/*Calculate the gravity response for a lower triangle (2-D) block with corners(x1,z2),(x2,z1),(x2,z2)*/
/*Parameter:	dens       := Density of the cell*/
/*              x1,z1   := Distance from the (left,top) corner of the triangle to the considered gravity station*/
/*              x2,z2   := Distance from the (right,bottom) corner of the triangle to the considered gravity station*/
/* Remark:	- The output is the gravity response for the triangle at the considered station*/

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/

double Response_2D_lower_triangle1(double dens, double x1, double x2, double z1, double z2)
{
	double a,b,c,d,e,f,g,first_integ,second_integ;
	double response;

	if((x1 == x2) || (z1 == z2))
		response = 0;
	else
	{
	
     	/*Calculate the response of the rectangular block in 2-D*/
     	a=((x2-x1)*(x2-x1)+(z2-z1)*(z2-z1))/((x2-x1)*(x2-x1));
    	b=2*(z1-z2)*(x2*z2-x1*z1)/((x2-x1)*(x2-x1));
    	c=(z2-x1*(z1-z2)/(x2-x1))*(z2-x1*(z1-z2)/(x2-x1));

    	d=x2*log(a*x2*x2+b*x2+c)-x1*log(a*x1*x1+b*x1+c)-2*(x2-x1);
    	e=b/(2*a)*(log(a*x2*x2+b*x2+c)-log(a*x1*x1+b*x1+c));
	    f=(4*a*c-b*b);


	    first_integ=x2*log(x2*x2+z2*z2)-x1*log(x1*x1+z2*z2)-2*(x2-x1)+2*z2*atan2(x2,z2)-2*z2*atan2(x1,z2);
	

	    if(f > 0)
	    	g=sqrt(f)/a*(atan2(2*a*x2,sqrt(f))-atan2(2*a*x1,sqrt(f)));
	    else
	    	if(f<0)
	        	g=-sqrt(-f)/a*(log((2*a*x2+b-sqrt(-f))/(2*a*x2+b+sqrt(-f)))-log((2*a*x1+b-sqrt(-f))/(2*a*x1+b+sqrt(-f))));
		    else
		    	g=0;

	    second_integ=d+e+g;
	
	
	    response = G_KONST*dens*(first_integ-second_integ);
	}
	return(response);

#undef G_KONST

}



/*-------------------------------------------------------------*/
/*Calculate the gravity response for a upper triangle (2-D) block with corners(x1,z1),(x2,z1),(x1,z2)*/
/*Parameter:	dens       := Density of the cell*/
/*              x1,z1   := Distance from the (left,top) corner of the triangle to the considered gravity station*/
/*              x2,z2   := Distance from the (right,bottom) corner of the triangle to the considered gravity station*/
/* Remark:	- The output is the gravity response for the triangle at the considered station*/

#define G_KONST 6.672E-11 /*Gravitational constant in (m^3)/(Kg*s^2)*/

double Response_2D_upper_triangle1(double dens, double x1, double x2, double z1, double z2)
{
	double a,b,c,d,e,f,g,first_integ,second_integ;
	double response;

	if((x1 == x2) || (z1 == z2))
		response = 0;
	else
	{
	
    	/*Calculate the response of the rectangular block in 2-D*/
    	a=((x2-x1)*(x2-x1)+(z2-z1)*(z2-z1))/((x2-x1)*(x2-x1));
	    b=2*(z1-z2)*(x2*z2-x1*z1)/((x2-x1)*(x2-x1));
	    c=(z2-x1*(z1-z2)/(x2-x1))*(z2-x1*(z1-z2)/(x2-x1));

	    d=x2*log(a*x2*x2+b*x2+c)-x1*log(a*x1*x1+b*x1+c)-2*(x2-x1);
	    e=b/(2*a)*(log(a*x2*x2+b*x2+c)-log(a*x1*x1+b*x1+c));
	    f=(4*a*c-b*b);


	    first_integ=x2*log(x2*x2+z1*z1)-x1*log(x1*x1+z1*z1)-2*(x2-x1)+2*z1*atan2(x2,z1)-2*z1*atan2(x1,z1);
	

	    if(f > 0)
	    	g=sqrt(f)/a*(atan2(2*a*x2,sqrt(f))-atan2(2*a*x1,sqrt(f)));
	    else
		    if(f<0)
		        g=-sqrt(-f)/a*(log((2*a*x2+b-sqrt(-f))/(2*a*x2+b+sqrt(-f)))-log((2*a*x1+b-sqrt(-f))/(2*a*x1+b+sqrt(-f))));
		    else
			    g=0;

	    second_integ=d+e+g;
	
	
	    response = G_KONST*dens*(-first_integ+second_integ);
	}
	return(response);

#undef G_KONST

}