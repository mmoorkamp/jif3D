#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "nrutil.h"
#include "inv3d.h"
#include "fkt.h"

/*-------------------------------------------------------------*/
/* Read in the relationsgip of velocity and density(resistivity) and check if the function is implicit*/
/* Parameter:  *fname	:= Filename of the parameter file including the relationsship*/
/*				rel_vdr	:= Structure that organizes the relationship of velocity and density and velocity and resistivity*/

int ReadRelVelDensRes( char *fname, REL_VDR_STRUCT *rel_vdr)
{
	int nr_col;
	long i; 
	FILE *inf;
	char line[128];

	inf = fopen(fname,"rt");
		if (inf == NULL)
		{
			fprintf(stderr,"Unable to open %s\n", fname);
			exit(0);
		}

	rel_vdr->tab_v = (double *)memory(NULL,1,sizeof(double),"REL_VDR_STRUCT");
	rel_vdr->tab_dr = (double *)memory(NULL,1,sizeof(double),"REL_VDR_STRUCT");

	/*************************************/
	/*Read in the parameter file*/
     GetNext(line,inf);
	 sscanf(line,"%d\n", &rel_vdr->index);

	 if(rel_vdr->index == 1) /*linear relationship of v-res(v-dens)*/
	 {
		 fgets(line,128,inf);
		 nr_col = sscanf(line,"%f %f\n", &rel_vdr->lin_para_a, &rel_vdr->lin_para_b);

		if(nr_col != 2)
			 printf("WARNING !!! In the file %s %d parameters instead of 2\n parameters are found for the linear relationship !!\n\n",fname,nr_col);
	 }
	 else if(rel_vdr->index == 2) /*v-res rel. of Dell'Aversana(2001)*/
	 {
		 fgets(line,128,inf);
		 nr_col = sscanf(line,"%f %f\n", &rel_vdr->dell_para_a, &rel_vdr->dell_para_b);

		 if(nr_col != 2)
			 printf("WARNING !!! In the file %s %d parameters instead of 2\n parameters are found for the relationship from Dell'Aversana !!\n\n",fname,nr_col);
	 }
	 else if(rel_vdr->index == 9) /*tabular relationship*/
	 {
		 i = 0;

		 fgets(line,128,inf);

		 while('#' != line[0])
		 {
		    i++;
			rel_vdr->tab_v = (double *)memory((char *)rel_vdr->tab_v,i+1,sizeof(double),"REL_VDR_STRUCT");	/*ATTENTION!!! Arrays start at 1 not at 0 !*/
			rel_vdr->tab_dr = (double *)memory((char *)rel_vdr->tab_dr,i+1,sizeof(double),"REL_VDR_STRUCT");	/*ATTENTION!!! Arrays start at 1 not at 0 !*/

			nr_col = sscanf(line,"%lf %lf\n", &rel_vdr->tab_v[i], &rel_vdr->tab_dr[i]);

			if(nr_col != 2)
				printf("WARNING !!! In the file %s %d columns instead of 2\n columns are used !!\n\n", fname,nr_col);


			fgets(line,128,inf);
		 }

		 rel_vdr->tab_nr_pairs = i;
	 }
	 else
	 {
		 printf("The relationship of velocities and densities \n or velocities and resistivities is not correctly specified !!\n\n");
		 exit(0);
	 }

	fclose(inf);

	/***********************************************************************/
	/*Check for the tabular if the function is implicit*/
	if(rel_vdr->index == 9)
	{
		/*Sort the velocities in ascending order*/
		sort2(rel_vdr->tab_nr_pairs,rel_vdr->tab_v,rel_vdr->tab_dr);

		/*Check if the all densities or resistivities have also all an ascending or descending order with increasing velocities*/
		CheckMonotom(rel_vdr->tab_nr_pairs,rel_vdr->tab_v,rel_vdr->tab_dr);
	}

	return(0);
}


/*-------------------------------------------------------------*/
/* Check if the function is strongly monoton*/
/* Parameter:  nr_pairs	:= number of pairs*/
/*				*x		:= x-values*/
/*				*y		:= y-values*/
/*Remark: - The vectors start at 1 NOT at 0*/
/*        - The x-values have to have an ascending order (sort them with "sort2.c" before)*/

int CheckMonotom(long nr_pairs,double *x, double *y)
{
	int trend, tmp_trend;
	long i;
	
	# define UPWARDS 1
	# define DOWNWARDS 0		

		tmp_trend = 2;
		trend = 2;

		for(i=1;i<nr_pairs; i++)
		{	
			trend = 2;

			if(y[i] < y[i+1])
				trend = UPWARDS;
			if(y[i] > y[i+1])
				trend = DOWNWARDS;

			if(i == 1)
				tmp_trend = trend;
			
			if(tmp_trend != trend || trend == 2 || (x[i] == x[i+1] && y[i] != y[i+1]))
			{
				printf("The function f(x) is not monotonically increasing or decreasing\nwith increasing x-values !!!\n");
				exit(0);
			}

			tmp_trend = trend;
		}

	#undef UPWARDS
	#undef DOWNWARDS

	return(0);
}

/*-------------------------------------------------------------*/
/* Determine the density[mgal] or resistivity[ohmm] from the velocity[m/s]*/
/* Parameter:   vel_value := Velocity value for which the density/resistivity value should be calculated */
/*				rel_vdr	  := Structure that organizes the relationship of velocity and density and velocity and resistivity*/
/*   Output: Density[mgal] or resistivity value[ohmm]*/
/*   Remark: For the tabular link the velocities have to be sorted in an ascending order*/

double VeltoDensRes(double vel_value, REL_VDR_STRUCT rel_vdr)
{
	long i;
	double dens_res_value;
	double e = 2.718281828;

	dens_res_value = 0.0; /*Density or resistivity value*/

	/*Choose the fomular*/
	switch(rel_vdr.index)
	{
		case 1: /*linear relationship of density/resistivity and velocity*/

			if(rel_vdr.lin_para_a == 0.0)
			{
				printf("The first parameter relating the velocity to the density/resistivity is badly chosen\n");
				printf("Change this parameter!!!\n");
				exit(0);
			}

			dens_res_value = (vel_value - rel_vdr.lin_para_b)/rel_vdr.lin_para_a;
			break;
	
		case 2:	/*v-res rel. of Dell'Aversana(2001)*/

			if(((vel_value - rel_vdr.dell_para_b)/rel_vdr.dell_para_a) > 6.5)
			{
				printf("The parameters for the relationships of resistivity/density and velocity from\nDell'Aversana are badly chosen!!\n");
				printf("Use parameters such that the resistivities values are in a useful range!!!\n");
				exit(0);
			}

			dens_res_value = pow(e, pow(e,((vel_value - rel_vdr.dell_para_b)/rel_vdr.dell_para_a)));
			break;

		case 9: /*Tabular relationship*/

			/*If only one data pair exists the value of the density/resistivity will NOT be calculated*/
			if(rel_vdr.tab_nr_pairs == 1) 
			{
				printf("At least TWO data pairs of densities/resistivities and velocities are required in the table\nto build up an implicite function\n!!");
				exit(0);
			}

			/*Remark: If the input data are smaller/larger than the smallest/largest density/resistivity value specified in the data file the data will be extrapolated*/
			/*		  Otherwise the data will be interpolated*/

			/*Extrapolation*/
			if(vel_value < rel_vdr.tab_v[1])
			{
				dens_res_value = (((rel_vdr.tab_dr[2] - rel_vdr.tab_dr[1])*(vel_value - rel_vdr.tab_v[1]))/(rel_vdr.tab_v[2] - rel_vdr.tab_v[1])) + rel_vdr.tab_dr[1];
				break;
			}

			/*Loop over all data points*/
			/*Interpolation*/
			for(i=1;i<rel_vdr.tab_nr_pairs + 1 ;i++)
			{
				if(vel_value < rel_vdr.tab_v[i+1])
				{
					dens_res_value = (((rel_vdr.tab_dr[i+1] - rel_vdr.tab_dr[i])*(vel_value - rel_vdr.tab_v[i]))/(rel_vdr.tab_v[i+1] - rel_vdr.tab_v[i])) + rel_vdr.tab_dr[i];
					goto beenden;
				}
			}

			/*Extrapolation*/
			dens_res_value = (((rel_vdr.tab_dr[rel_vdr.tab_nr_pairs] - rel_vdr.tab_dr[rel_vdr.tab_nr_pairs-1])*(vel_value - rel_vdr.tab_v[rel_vdr.tab_nr_pairs-1]))/(rel_vdr.tab_v[rel_vdr.tab_nr_pairs] - rel_vdr.tab_v[rel_vdr.tab_nr_pairs-1]))  + rel_vdr.tab_dr[rel_vdr.tab_nr_pairs-1];
			break;

			beenden:;
			break;

		default:
			printf("The relationship between the velocity and the density/resistivity is NOT correct specified\n");
			exit(0);
	}

	return(dens_res_value);
}

/*-------------------------------------------------------------*/
/* Determine velocity [m/s] from the density[mgal] or resistivity[ohmm]*/
/* Parameter:   dens_res_value := Density/resistivity value for which the velocity value should be calculated */
/*				rel_vdr	  := Structure that organizes the relationship of velocity and density and velocity and resistivity*/
/*   Output: Velocity[m/s]*/
/*Remark:-Inverse formulars of the ones in "VeltoDensRes"*/
/*       -For the tabular link the densities/resistivities have to be sorted in an ascending order*/

double DensRestoVel(double dens_res_value, REL_VDR_STRUCT rel_vdr)
{
	long i;
	double vel_value;

	vel_value = 0.0; /*Density or resistivity value*/

	/*Choose the fomular*/
	switch(rel_vdr.index)
	{
			case 1: /*linear relationship of density/resistivity and velocity*/
				vel_value = rel_vdr.lin_para_a * dens_res_value + rel_vdr.lin_para_b;
				break;

			case 2:	/*v-res rel. of Dell'Aversana(2001)*/
				vel_value = rel_vdr.dell_para_a*log(log(dens_res_value))+ rel_vdr.dell_para_b;
				break;

			case 9: /*Tabular relationship*/

			/*If only one data pair exists the value of the velocities will NOT be calculated*/
			if(rel_vdr.tab_nr_pairs == 1) 
			{
				printf("At least TWO data pairs of densities/resistivities and velocities are required in the table\nto build up an implicite function\n!!");
				exit(0);
			}

			/*Remark: If the input data are smaller/larger than the smallest/largest density/resistivity value specified in the data file the data will be extrapolated*/
			/*		  Otherwise the data will be interpolated*/

			/*Extrapolation*/
			if(dens_res_value <rel_vdr.tab_dr[1])
			{
				vel_value = (((rel_vdr.tab_v[2] - rel_vdr.tab_v[1])*(dens_res_value - rel_vdr.tab_dr[1]))/(rel_vdr.tab_dr[2] - rel_vdr.tab_dr[1])) + rel_vdr.tab_v[1];
				break;
			}

			/*Loop over all data points*/
			/*Interpolation*/
			for(i=1;i<rel_vdr.tab_nr_pairs + 1 ;i++)
			{
				if(dens_res_value < rel_vdr.tab_dr[i+1])
				{
					vel_value = (((rel_vdr.tab_v[i+1] - rel_vdr.tab_v[i])*(dens_res_value - rel_vdr.tab_dr[i]))/(rel_vdr.tab_dr[i+1] - rel_vdr.tab_dr[i])) + rel_vdr.tab_v[i];
					goto beenden;
				}
			}

			/*Extrapolation*/
			vel_value = (((rel_vdr.tab_v[rel_vdr.tab_nr_pairs] - rel_vdr.tab_v[rel_vdr.tab_nr_pairs-1])*(dens_res_value - rel_vdr.tab_dr[rel_vdr.tab_nr_pairs-1]))/(rel_vdr.tab_dr[rel_vdr.tab_nr_pairs] - rel_vdr.tab_dr[rel_vdr.tab_nr_pairs-1]))  + rel_vdr.tab_v[rel_vdr.tab_nr_pairs-1];
			break;

			beenden:;
			break;


			default:
				printf("The relationship between the velocity and the density/resistivity is NOT correct specified\n");
				exit(0);

	}

	return(vel_value);
}