#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "inv3d.h"
#include "fkt.h"



/*-------------------------------------------------------------*/
/* Check if the positions of the shot and receiver locations are inside the model (if NOT: the program */
/* terminates) */
/* Parameter:	geo := geometry structure*/
/*              data:= data structure */
/*              grid := grid structure*/


int CheckCoordinates(GEOMETRY geo, DATA_STRUCT data, GRID_STRUCT grid, FLAG_STRUCT flag)
{
	int i,j,k,found[6];
	int a[6], max_outlier_i[6];
	double max_outlier[6];

	double *x_below[6] ; /*Coordinates of shots/receiver outside of the grid*/
	int *x_below_i[6] ; /*The indeces of the shots/receiver outside of the grid*/

	#define xs(a) geo.x[(data.shots[a])-1] /*x-coordinate of the shot locations*/
	#define ys(a) geo.y[(data.shots[a])-1] /*y-coordinate of the shot locations*/
	#define zs(a) geo.z[(data.shots[a])-1] /*z-coordinate of the shot locations*/
	#define xr(a) geo.x[(data.recs[a])-1] /*x-coordinate of the receiver locations*/
	#define yr(a) geo.y[(data.recs[a])-1] /*y-coordinate of the receiver locations*/
	#define zr(a) geo.z[(data.recs[a])-1] /*z-coordinate of the receiver locations*/
    #define xg(a) geo.x[(data.gravs[a])-1] /*x-coordinate of the gravity stations*/
	#define yg(a) geo.y[(data.gravs[a])-1] /*x-coordinate of the gravity stations*/
	#define zg(a) geo.z[(data.gravs[a])-1] /*x-coordinate of the gravity stations*/
    #define xm(a) geo.x[(data.mts[a])-1] /*x-coordinate of the MT stations*/
	#define ym(a) geo.y[(data.mts[a])-1] /*x-coordinate of the MT stations*/
	#define zm(a) geo.z[(data.mts[a])-1] /*x-coordinate of the MT stations*/

	
	for(k=0;k<6;k++)
	{
		x_below_i[k] = (int *)memory(NULL,1,sizeof(int),"CheckCoordinates");
		x_below[k]   = (double *)memory(NULL,1,sizeof(double),"CheckCoordinates");

		found[k]=0;
		a[k]=0;
	}

	printf("Start checking if all shot/receiver locations, gravity and MT stations are inside the model:\n");
	printf("----------------\n");


	if(flag.index_tseis != 0)
	{
			/*-----------------------------------------------------*/
			/*Control shots*/
			for(i=0;i<geo.nshot;i++)
			{
					/*Check if the x-coordinates of the shots are larger
					than the left grid x-boundary coordinate*/
				if(xs(i) < grid.org[0])
				{
					a[0]++;
					x_below_i[0] = (int *)memory((char *)x_below_i[0],a[0],sizeof(int),"CheckCoordinates");
					x_below[0]   = (double *)memory((char *)x_below[0],a[0],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and shot station numbers which are outside of the grid*/
					x_below_i[0][a[0]-1] = data.shots[i];
					x_below[0][a[0]-1] = xs(i);

					found[0] =1;

				}
			
				/*Check if the x-coordinates of the shots are smaller
				than the right grid x-boundary coordinate*/
				if(xs(i) > grid.org[0] + (grid.nx-1)*grid.h)
				{
					a[1]++;
					x_below_i[1] = (int *)memory((char *)x_below_i[1],a[1],sizeof(int),"CheckCoordinates");
					x_below[1]   = (double *)memory((char *)x_below[1],a[1],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and shot station numbers which are outside of the grid*/
					x_below_i[1][a[1]-1] = data.shots[i];
					x_below[1][a[1]-1] = xs(i);

					found[1] =1;
				}


				/*Check if the y-coordinates of the shots are larger
				  than the left grid y-boundary coordinate*/
				if(ys(i) < grid.org[1])
				{
					a[2]++;
					x_below_i[2] = (int *)memory((char *)x_below_i[2],a[2],sizeof(int),"CheckCoordinates");
					x_below[2]   = (double *)memory((char *)x_below[2],a[2],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and shot station numbers which are outside of the grid*/
					x_below_i[2][a[2]-1] = data.shots[i];
					x_below[2][a[2]-1] = ys(i);

					found[2] =1;
				}

				/*Check if the x-coordinates of the shots are smaller
				  than the right grid x-boundary coordinate*/
				if(ys(i) > grid.org[1] + (grid.ny-1)*grid.h)
				{
					a[3]++;
					x_below_i[3] = (int *)memory((char *)x_below_i[3],a[3],sizeof(int),"CheckCoordinates");
					x_below[3]   = (double *)memory((char *)x_below[3],a[3],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and shot station numbers which are outside of the grid*/
					x_below_i[3][a[3]-1] = data.shots[i];
					x_below[3][a[3]-1] = ys(i);

					found[3] =1;
				}

				/*Check if the z-coordinates of the shots are larger
				than the left grid z-boundary coordinate*/
				if(zs(i) < grid.org[2])
				{
					a[4]++;
					x_below_i[4] = (int *)memory((char *)x_below_i[4],a[4],sizeof(int),"CheckCoordinates");
					x_below[4]   = (double *)memory((char *)x_below[4],a[4],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and shot station numbers which are outside of the grid*/
					x_below_i[4][a[4]-1] = data.shots[i];
					x_below[4][a[4]-1] = zs(i);


					found[4] =1;
				}

				/*Check if the x-coordinates of the shots are smaller
				  than the right grid x-boundary coordinate*/
				if(zs(i) > grid.org[2] + (grid.nz-1)*grid.h)
				{
					a[5]++;
					x_below_i[5] = (int *)memory((char *)x_below_i[5],a[5],sizeof(int),"CheckCoordinates");
					x_below[5]   = (double *)memory((char *)x_below[5],a[5],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and shot station numbers which are outside of the grid*/
					x_below_i[5][a[5]-1] = data.shots[i];
					x_below[5][a[5]-1] = zs(i);

					found[5] =1;
				}

			}

		/*-----------------------------------------------------*/
		/*Control receivers*/
		for(i=0;i<geo.nrec;i++)
		{
				/*Check if the x-coordinates of the receivers are larger
				  than the left grid x-boundary coordinate*/
				if(xr(i) < grid.org[0])
				{
					a[0]++;
					x_below_i[0] = (int *)memory((char *)x_below_i[0],a[0],sizeof(int),"CheckCoordinates");
					x_below[0]   = (double *)memory((char *)x_below[0],a[0],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and receiver station numbers which are outside of the grid*/
					x_below_i[0][a[0]-1] = data.recs[i];
					x_below[0][a[0]-1] = xr(i);

					found[0] =1;

				}
			
				/*Check if the x-coordinates of the recs are smaller
				  than the right grid x-boundary coordinate*/
				if(xr(i) > grid.org[0] + (grid.nx-1)*grid.h)
				{
					a[1]++;
					x_below_i[1] = (int *)memory((char *)x_below_i[1],a[1],sizeof(int),"CheckCoordinates");
					x_below[1]   = (double *)memory((char *)x_below[1],a[1],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and receiver station numbers which are outside of the grid*/
					x_below_i[1][a[1]-1] = data.recs[i];
					x_below[1][a[1]-1] = xr(i);

					found[1] =1;
				}


				/*Check if the y-coordinates of the receivers are larger
				  than the left grid y-boundary coordinate*/
				if(yr(i) < grid.org[1])
				{
					a[2]++;
					x_below_i[2] = (int *)memory((char *)x_below_i[2],a[2],sizeof(int),"CheckCoordinates");
					x_below[2]   = (double *)memory((char *)x_below[2],a[2],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and receiver station numbers which are outside of the grid*/
					x_below_i[2][a[2]-1] = data.recs[i];
					x_below[2][a[2]-1] = yr(i);

					found[2] =1;
				}

				/*Check if the x-coordinates of the receivers are smaller
				than the right grid x-boundary coordinate*/
				if(yr(i) > grid.org[1] + (grid.ny-1)*grid.h)
				{
					a[3]++;
					x_below_i[3] = (int *)memory((char *)x_below_i[3],a[3],sizeof(int),"CheckCoordinates");
					x_below[3]   = (double *)memory((char *)x_below[3],a[3],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and receiver station numbers which are outside of the grid*/
					x_below_i[3][a[3]-1] = data.recs[i];
					x_below[3][a[3]-1] = yr(i);

					found[3] =1;
				}

				/*Check if the z-coordinates of the receivers are larger
				  than the left grid z-boundary coordinate*/
				if(zr(i) < grid.org[2])
				{
					a[4]++;
					x_below_i[4] = (int *)memory((char *)x_below_i[4],a[4],sizeof(int),"CheckCoordinates");
					x_below[4]   = (double *)memory((char *)x_below[4],a[4],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and receiver station numbers which are outside of the grid*/
					x_below_i[4][a[4]-1] = data.recs[i];
					x_below[4][a[4]-1] = zr(i);


					found[4] =1;
				}

				/*Check if the x-coordinates of the receivers are smaller
				  than the right grid x-boundary coordinate*/
				if(zr(i) > grid.org[2] + (grid.nz-1)*grid.h)
				{
					a[5]++;
					x_below_i[5] = (int *)memory((char *)x_below_i[5],a[5],sizeof(int),"CheckCoordinates");
					x_below[5]   = (double *)memory((char *)x_below[5],a[5],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and receiver station numbers which are outside of the grid*/
					x_below_i[5][a[5]-1] = data.recs[i];
					x_below[5][a[5]-1] = zr(i);

					found[5] =1;
				}

		}
	}

	/*-----------------------------------------------------*/

	if(flag.index_grav != 0)
	{
		/*Control gravity stations*/
		/*REMARK: In contrast to the seismic, where the shot and receiver have to be located within the region specified*/
		/*		  by the outer grid cell centers, the gravity stations are allowed to be located exactly at the border of the grid!!*/
		
		for(i=0;i<geo.nstat_grav;i++)
		{
				/*Check if the x-coordinates of the gravity stations are larger
				  than the left grid x-boundary coordinate*/
				if(xg(i) < (grid.org[0] - grid.h/2))
				{
					a[0]++;
					x_below_i[0] = (int *)memory((char *)x_below_i[0],a[0],sizeof(int),"CheckCoordinates");
					x_below[0]   = (double *)memory((char *)x_below[0],a[0],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and gravity station numbers which are outside of the grid*/
					x_below_i[0][a[0]-1] = data.gravs[i];
					x_below[0][a[0]-1] = xg(i);

					found[0] =1;

				}
			
				/*Check if the x-coordinates of the gravity stations are smaller
				  than the right grid x-boundary coordinate*/
				if(xg(i) > grid.org[0] + (grid.nx-0.5)*grid.h)
				{
					a[1]++;
					x_below_i[1] = (int *)memory((char *)x_below_i[1],a[1],sizeof(int),"CheckCoordinates");
					x_below[1]   = (double *)memory((char *)x_below[1],a[1],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and gravity station numbers which are outside of the grid*/
					x_below_i[1][a[1]-1] = data.gravs[i];
					x_below[1][a[1]-1] = xg(i);

					found[1] =1;
				}


				/*Check if the y-coordinates of the gravity stations are larger
				  than the left grid y-boundary coordinate*/
				if(yg(i) < (grid.org[1] - grid.h/2))
				{
					a[2]++;
					x_below_i[2] = (int *)memory((char *)x_below_i[2],a[2],sizeof(int),"CheckCoordinates");
					x_below[2]   = (double *)memory((char *)x_below[2],a[2],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and gravity station numbers which are outside of the grid*/
					x_below_i[2][a[2]-1] = data.gravs[i];
					x_below[2][a[2]-1] = yg(i);

					found[2] =1;
				}

				/*Check if the x-coordinates of the gravity stations are smaller
				than the right grid x-boundary coordinate*/
				if(yg(i) > grid.org[1] + (grid.ny-0.5)*grid.h)
				{
					a[3]++;
					x_below_i[3] = (int *)memory((char *)x_below_i[3],a[3],sizeof(int),"CheckCoordinates");
					x_below[3]   = (double *)memory((char *)x_below[3],a[3],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and gravity station numbers which are outside of the grid*/
					x_below_i[3][a[3]-1] = data.gravs[i];
					x_below[3][a[3]-1] = yg(i);

					found[3] =1;
				}

				/*Check if the z-coordinates of the gravity stations are larger
				  than the left grid z-boundary coordinate*/
				if(zg(i) < (grid.org[2] - grid.h/2))
				{
					a[4]++;
					x_below_i[4] = (int *)memory((char *)x_below_i[4],a[4],sizeof(int),"CheckCoordinates");
					x_below[4]   = (double *)memory((char *)x_below[4],a[4],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and gravity station numbers which are outside of the grid*/
					x_below_i[4][a[4]-1] = data.gravs[i];
					x_below[4][a[4]-1] = zg(i);

					found[4] =1;
				}

				/*Check if the x-coordinates of the gravity stations are smaller
				  than the right grid x-boundary coordinate*/
				if(zg(i) > grid.org[2] + (grid.nz-0.5)*grid.h)
				{
					a[5]++;
					x_below_i[5] = (int *)memory((char *)x_below_i[5],a[5],sizeof(int),"CheckCoordinates");
					x_below[5]   = (double *)memory((char *)x_below[5],a[5],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and gravity station numbers which are outside of the grid*/
					x_below_i[5][a[5]-1] = data.gravs[i];
					x_below[5][a[5]-1] = zg(i);

					found[5] =1;
				}

		}

	}
	
	/*-----------------------------------------------------*/

	if(flag.index_mt != 0)
	{
		/*Control MT stations*/
		for(i=0;i<geo.nstat_mt;i++)
		{
				/*Check if the x-coordinates of the MT stations are larger
				  than the left grid x-boundary coordinate*/
				if(xm(i) < grid.org[0])
				{
					a[0]++;
					x_below_i[0] = (int *)memory((char *)x_below_i[0],a[0],sizeof(int),"CheckCoordinates");
					x_below[0]   = (double *)memory((char *)x_below[0],a[0],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and MT station numbers which are outside of the grid*/
					x_below_i[0][a[0]-1] = data.mts[i];
					x_below[0][a[0]-1] = xm(i);

					found[0] =1;

				}
			
				/*Check if the x-coordinates of the MT stations are smaller
				  than the right grid x-boundary coordinate*/
				if(xm(i) > grid.org[0] + (grid.nx-1)*grid.h)
				{
					a[1]++;
					x_below_i[1] = (int *)memory((char *)x_below_i[1],a[1],sizeof(int),"CheckCoordinates");
					x_below[1]   = (double *)memory((char *)x_below[1],a[1],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and MT station numbers which are outside of the grid*/
					x_below_i[1][a[1]-1] = data.mts[i];
					x_below[1][a[1]-1] = xm(i);

					found[1] =1;
				}


				/*Check if the y-coordinates of the MT stations are larger
				  than the left grid y-boundary coordinate*/
				if(ym(i) < grid.org[1])
				{
					a[2]++;
					x_below_i[2] = (int *)memory((char *)x_below_i[2],a[2],sizeof(int),"CheckCoordinates");
					x_below[2]   = (double *)memory((char *)x_below[2],a[2],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and MT station numbers which are outside of the grid*/
					x_below_i[2][a[2]-1] = data.mts[i];
					x_below[2][a[2]-1] = ym(i);

					found[2] =1;
				}

				/*Check if the x-coordinates of the MT stations are smaller
				than the right grid x-boundary coordinate*/
				if(ym(i) > grid.org[1] + (grid.ny-1)*grid.h)
				{
					a[3]++;
					x_below_i[3] = (int *)memory((char *)x_below_i[3],a[3],sizeof(int),"CheckCoordinates");
					x_below[3]   = (double *)memory((char *)x_below[3],a[3],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and MT station numbers which are outside of the grid*/
					x_below_i[3][a[3]-1] = data.mts[i];
					x_below[3][a[3]-1] = ym(i);

					found[3] =1;
				}

				/*Check if the z-coordinates of the MT stations are larger
				  than the left grid z-boundary coordinate*/
				if(zm(i) < grid.org[2])
				{
					a[4]++;
					x_below_i[4] = (int *)memory((char *)x_below_i[4],a[4],sizeof(int),"CheckCoordinates");
					x_below[4]   = (double *)memory((char *)x_below[4],a[4],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and MT station numbers which are outside of the grid*/
					x_below_i[4][a[4]-1] = data.mts[i];
					x_below[4][a[4]-1] = zm(i);

					found[4] =1;
				}

				/*Check if the x-coordinates of the MT stations are smaller
				  than the right grid x-boundary coordinate*/
				if(zm(i) > grid.org[2] + (grid.nz-1)*grid.h)
				{
					a[5]++;
					x_below_i[5] = (int *)memory((char *)x_below_i[5],a[5],sizeof(int),"CheckCoordinates");
					x_below[5]   = (double *)memory((char *)x_below[5],a[5],sizeof(double),"CheckCoordinates");

					/*Collect the coordinates and MT station numbers which are outside of the grid*/
					x_below_i[5][a[5]-1] = data.mts[i];
					x_below[5][a[5]-1] = zm(i);

					found[5] =1;
				}

		}

	}

	/*-----------------------------------------------------*/
	/*Look for the largest outliers*/
	max_outlier[0] = grid.org[0];
	max_outlier[2] = grid.org[1];
	max_outlier[4] = grid.org[2];

	max_outlier[1] = grid.org[0] + (grid.nx - 1)*grid.h;
	max_outlier[3] = grid.org[1] + (grid.ny - 1)*grid.h;
	max_outlier[5] = grid.org[2] + (grid.nz - 1)*grid.h;


	for(i=0;i<6;i++)
	{
		if(found[i] == 1)
		{
			for(j=0;j<a[i];j++)
			{
				if(x_below[i][j] <= max_outlier[i] && i%2 == 0)
				{
					max_outlier[i] = x_below[i][j];
					max_outlier_i[i] = x_below_i[i][j];
				}
				else if(x_below[i][j] >= max_outlier[i] && (i+1)%2 == 0)
				{
					max_outlier[i] = x_below[i][j];
					max_outlier_i[i] = x_below_i[i][j];
				}
		
			}
		}
	}



		if(found[0] == 0)
		{
			printf("No receiver/shot/gravity/MT stations have coordinates that are\nsmaller than those of the lower x-border of the grid\n");
			printf("----------------\n");
		}
		else
		{
			printf("Altogether %d shot/receiver/gravity/MT stations have x-coordinates that are\nsmaller than those of the lower x-border of the grid:\n",a[0]);
			printf("The lower x-boundary: x= %fm \nThe receiver/shot/gravity/MT station %d has the smallest x-coordinate:\n x= %fm \n", grid.org[0], max_outlier_i[0], max_outlier[0]);
			printf("Re-define lower border in x-direction!!\n");
			printf("----------------\n");
		}

		if(found[1] == 0)
		{
			printf("No receiver/shot/gravity/MT stations have coordinates that are\nlarger than those of the upper x-border of the grid\n");
			printf("----------------\n");
		}
		else
		{
			printf("Altogether %d shot/receiver/gravity/MT stations have x-coordinates that are\nlarger than those of the upper x-border of the grid:\n",a[1]);
			printf("The upper x-boundary: x= %fm \nThe receiver/shot/gravity/MT station %d has the largest x-coordinate:\n x= %fm \n", grid.org[0] + (grid.nx-1)*grid.h, max_outlier_i[1], max_outlier[1]);
			printf("Re-define upper border in x-direction!!\n");
			printf("----------------\n");
		}

		if(found[2] == 0)
		{
			printf("No receiver/shot/gravity/MT stations have coordinates that are\nsmaller than those of the lower y-border of the grid\n");
			printf("----------------\n");
		}
		else
		{
			printf("Altogether %d shot/receiver/gravity/MT stations have y-coordinates that are\nsmaller than those of the lower x-border of the grid:\n",a[2]);
			printf("The lower y-boundary: y= %fm \nThe receiver/shot/gravity/MT station %d has the smallest y-coordinate:\n y= %fm \n", grid.org[1], max_outlier_i[2], max_outlier[2]);
			printf("Re-define lower border in y-direction!!\n");
			printf("----------------\n");
		}

		if(found[3] == 0)
		{
			printf("No receiver/shot/gravity/MT stations have coordinates that are\nlarger than those of the upper y-border of the grid\n");
			printf("----------------\n");
		}
		else
		{
			printf("Altogether %d shot/receiver/gravity/MT stations have y-coordinates that are\nlarger than those of the upper y-border of the grid:\n",a[3]);
			printf("The upper y-boundary: y= %fm \nThe receiver/shot/gravity/MT station %d has the largest y-coordinate:\n y= %fm \n", grid.org[1] + (grid.ny-1)*grid.h, max_outlier_i[3], max_outlier[3]);
			printf("Re-define upper border in y-direction!!\n");
			printf("----------------\n");
		}

		if(found[4] == 0)
		{
			printf("No receiver/shot/gravity/MT stations have coordinates that are\nsmaller than those of the lower z-border of the grid\n");
			printf("----------------\n");
		}
		else
		{
			printf("Altogether %d shot/receiver/gravity/MT stations have z-coordinates that are\nsmaller than those of the lower z-border of the grid:\n",a[4]);
			printf("The lower z-boundary: z= %fm \nThe receiver/shot/gravity/MT station %d has the smallest z-coordinate:\n z= %fm \n", grid.org[2], max_outlier_i[4], max_outlier[4]);
			printf("Re-define lower border in z-direction!!\n");
			printf("----------------\n");
		}

		if(found[5] == 0)
		{
			printf("No receiver/shot/gravity/MT stations have coordinates that are\nlarger than those of the upper z-border of the grid\n");
			printf("----------------\n");
		}
		else
		{
			printf("Altogether %d shot/receiver/gravity/MT stations have z-coordinates that are\nlarger than those of the upper z-border of the grid:\n",a[5]);
			printf("The upper z-boundary: z= %fm \nThe receiver/shot/gravity/MT station %d has the largest z-coordinate:\n z= %fm \n", grid.org[2] + (grid.nz-1)*grid.h, max_outlier_i[5], max_outlier[5]);
			printf("Re-define upper border in z-direction!!\n");
			printf("----------------\n");
		}

		for(i=0;i<6;i++)
			if(found[i]==1)
			{ 
				printf("Change the grid sizes in the parameter file in a proper way\nso that all receiver/shot/gravity/MT stations are inside the cube !!!\n");
	 			exit(0);
			}

		printf("All used shot/receiver/gravity/MT stations are located inside the model -> OK!!\n");
		printf("----------------\n\n");

		for(k=0;k<6;k++)
		{
			free(x_below_i[k]); 
			free(x_below[k]); 
		}

		#undef xs
		#undef ys
		#undef zs
		#undef xr
		#undef yr
		#undef zr
		#undef xg
		#undef yg
		#undef zg
		#undef xm
		#undef ym
		#undef zm

	return(1);
}



/*-------------------------------------------------------------*/
/* Check if the positions of the shot/receiver locations and gravity/MT stations are not in the air */
/* (the program give out warnings but will not terminate) */
/* Parameter:	geo := geometry structure*/
/*              data:= data structure */
/*              grid := grid structure*/


int CheckCoordinatesTopo(GEOMETRY geo,DATA_STRUCT data,GRID_STRUCT grid)
{
	int ix,iy,iz;
	int index_x, index_y, index_z;
	int i;

	long nyz,nz,ny,ny2,nz2,nyz2,nborder;

	double X,Y,Z;
	double deltaX,deltaY,deltaZ;

	#define xs(a) geo.x[(data.shots[a])-1] /*x-coordinate of the shot locations*/
	#define ys(a) geo.y[(data.shots[a])-1] /*y-coordinate of the shot locations*/
	#define zs(a) geo.z[(data.shots[a])-1] /*z-coordinate of the shot locations*/
	#define xr(a) geo.x[(data.recs[a])-1] /*x-coordinate of the receiver locations*/
	#define yr(a) geo.y[(data.recs[a])-1] /*y-coordinate of the receiver locations*/
	#define zr(a) geo.z[(data.recs[a])-1] /*z-coordinate of the receiver locations*/
    #define xg(a) geo.x[(data.gravs[a])-1] /*x-coordinate of the gravity stations*/
	#define yg(a) geo.y[(data.gravs[a])-1] /*x-coordinate of the gravity stations*/
	#define zg(a) geo.z[(data.gravs[a])-1] /*x-coordinate of the gravity stations*/
    #define xm(a) geo.x[(data.mts[a])-1] /*x-coordinate of the MT stations*/
	#define ym(a) geo.y[(data.mts[a])-1] /*x-coordinate of the MT stations*/
	#define zm(a) geo.z[(data.mts[a])-1] /*x-coordinate of the MT stations*/

	#define b_index(x,y,z) grid.border_index[nyz2*((x)+nborder) + nz2*((y)+nborder) + ((z)+nborder)]


	ny = grid.ny;
	nz = grid.nz;
	nyz = nz*ny;

	nborder = grid.nborder;
	ny2 = ny+2*nborder;
	nz2 = nz+2*nborder;
	nyz2 = ny2*nz2;

	printf("Start checking if all shot and receiver locations are NOT in the air:\n");
	 printf("----------------\n");

	/*Control shots*/
	for(i=0;i<geo.nshot;i++)
	{
		/*Determine the next cell center from the shot coordinates*/

		ix = (int)((xs(i)-grid.org[0])/grid.h); /* index of grid cell in x-direction (which has the next lower x-value compared to the shot coordinates)*/
		iy = (int)((ys(i)-grid.org[1])/grid.h); /* index of grid cell in y-direction */
		iz = (int)((zs(i)-grid.org[2])/grid.h); /* index of grid cell in z-direction */
		
		X = (double)((xs(i)-grid.org[0])/grid.h); /*normalized distances between Grid-origins and shot coordinates*/
		Y = (double)((ys(i)-grid.org[1])/grid.h);
		Z = (double)((zs(i)-grid.org[2])/grid.h);

		deltaX = X-ix; /*distance from the shot coordintes to next lower grid cell center*/
		deltaY = Y-iy;
		deltaZ = Z-iz;

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

		if(b_index(index_x,index_y,index_z)==0)
			printf("WARNING!! The shot %d (x=%f, y=%f, z=%f) is located\nin the air -> check this point!!\n",data.shots[i], xs(i),ys(i),zs(i));
	}

	/*Control receivers*/
	for(i=0;i<geo.nrec;i++)
	{
		/*Determine the next cell center from the shot coordinates*/

		ix = (int)((xr(i)-grid.org[0])/grid.h); /* index of grid cell in x-direction (which has the next lower x-value compared to the shot coordinates)*/
		iy = (int)((yr(i)-grid.org[1])/grid.h); /* index of grid cell in y-direction */
		iz = (int)((zr(i)-grid.org[2])/grid.h); /* index of grid cell in z-direction */
		
		X = (double)((xr(i)-grid.org[0])/grid.h); /*normalized distances between Grid-origins and rec. coordinates*/
		Y = (double)((yr(i)-grid.org[1])/grid.h);
		Z = (double)((zr(i)-grid.org[2])/grid.h);

		deltaX = X-ix; /*distance from the shot coordintes to next lower grid cell center*/
		deltaY = Y-iy;
		deltaZ = Z-iz;

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

		if(b_index(index_x,index_y,index_z)==0)
			printf("WARNING!! The rec. %d (x=%f, y=%f, z=%f) is located\nin the air -> check this point!!\n",data.recs[i], xr(i),yr(i),zr(i));
	}


	/*Control gravity stations*/
	for(i=0;i<geo.nstat_grav;i++)
	{
		/*Determine the next cell center from the gravity station coordinates*/

		ix = (int)((xg(i)-grid.org[0])/grid.h); /* index of grid cell in x-direction (which has the next lower x-value compared to the shot coordinates)*/
		iy = (int)((yg(i)-grid.org[1])/grid.h); /* index of grid cell in y-direction */
		iz = (int)((zg(i)-grid.org[2])/grid.h); /* index of grid cell in z-direction */
		
		X = (double)((xg(i)-grid.org[0])/grid.h); /*normalized distances between Grid-origins and gravity station coordinates*/
		Y = (double)((yg(i)-grid.org[1])/grid.h);
		Z = (double)((zg(i)-grid.org[2])/grid.h);

		deltaX = X-ix; /*distance from the shot coordintes to next lower grid cell center*/
		deltaY = Y-iy;
		deltaZ = Z-iz;

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

		if(b_index(index_x,index_y,index_z)==0)
			printf("WARNING!! The gravity st. %d (x=%f, y=%f, z=%f) is located\nin the air -> check this point!!\n",data.gravs[i], xg(i),yg(i),zg(i));
	}

	/*Control MT stations*/
	for(i=0;i<geo.nstat_mt;i++)
	{
		/*Determine the next cell center from the MT station coordinates*/

		ix = (int)((xm(i)-grid.org[0])/grid.h); /* index of grid cell in x-direction (which has the next lower x-value compared to the shot coordinates)*/
		iy = (int)((ym(i)-grid.org[1])/grid.h); /* index of grid cell in y-direction */
		iz = (int)((zm(i)-grid.org[2])/grid.h); /* index of grid cell in z-direction */
		
		X = (double)((xm(i)-grid.org[0])/grid.h); /*normalized distances between Grid-origins and MT station coordinates*/
		Y = (double)((ym(i)-grid.org[1])/grid.h);
		Z = (double)((zm(i)-grid.org[2])/grid.h);

		deltaX = X-ix; /*distance from the shot coordintes to next lower grid cell center*/
		deltaY = Y-iy;
		deltaZ = Z-iz;

		if(deltaX <= 0.5) 
			index_x = ix;	/*index of next grid cell center in x-direction*/
		else
			index_x = ix+1;

		if(deltaY <= 0.5) 
			index_y = iy;
		else
			index_y = iy+1;

		if(deltaZ < 0.5) 
			index_z = iz;
		else
			index_z = iz+1;

		if(b_index(index_x,index_y,index_z)==0)
			printf("WARNING!! The MT st. %d (x=%f, y=%f, z=%f) is located\nin the air -> check this point!!\n",data.mts[i], xm(i),ym(i),zm(i));
	}

	 printf("----------------\n");
	 printf("The coordinate check is completed:\n");
	 printf("----------------\n\n");

	 #undef xs
	 #undef ys
	 #undef zs
	 #undef xr
	 #undef yr
	 #undef zr
	 #undef xg
	 #undef yg
	 #undef zg
	 #undef xm
	 #undef ym
	 #undef zm

	 #undef b_index

	return(1);
}