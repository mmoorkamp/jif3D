#include "inv3d.h"
#include "modeling_seismic.h"
#include <cmath>
/*-------------------------------------------------------------*/
/*Performing the forward modeling(using the Podvin&Lecomte eikonal solver) and calculating fat-rays */
/*Parameter:	geo  := Geometry structure  */
/*              grid  := Grid structure */
/*              *data  := Pointer on data structure (the calculated traveltimes will be (re-)determined in this routine)*/
/*				*fat_rays := Vector of fat-ray structures (number of elements = number of shot-receiver combinations)*/
/*              time_start := Time when inv3d starts*/

/*REMARK: The Podvin&Lecomte forward modeling calculates the traveltimes at the grid cell edges; The number of grid cell edges is in each direction one higher than the number of grid cell centers*/
/*		Therefore the grid have to be readjusted*/
namespace jiba
  {
int ForwardModFatRay(GEOMETRY geo,GRID_STRUCT grid, DATA_STRUCT *data, F_RP_STRUCT *fat_rays, time_t time_start)
{
	int count							/*Number of active receivers for a shot*/;
	int i,j,k,nborder,exist;
	int *nact_rec;						/*Active receiver-numbers for the used shot*/
	long *nact_datapos;					/*Position of the active receivers in the data structure*/
	long counter_index,*rec_counter;	/*Vector of the location numbers of receivers for which the traveltimes models are already calculated*/
	long *rec_num_counter;				/* Counter of shots for each receiver, for which the calculation is already performed*/
	long nx,ny,nz,nx3,ny3,nz3,nyz3;
	long max_num;

	float *tt_shot,*tt_rec,*tmp_slow;
	float Xs,Ys,Zs,*Xr,*Yr,*Zr;			/*Normalized positions of the shots and receivers (referring to the grid cell nodes and NOT of the grid cell centers)*/
	float org[3];
	float delta_num = (float)0.001;
	double sum2;

	char s1[32],s2[32],s3[32],s1_tmp[32],s2_tmp[32],s3_tmp[32];

	time_t time_relative;
	double zeit_angabe_sekunden;
	/*FILE *out;*/

	nx= grid.nx+1;
	ny= grid.ny+1;
	nz= grid.nz+1;
	nborder =grid.nborder;

	org[0]= (float)(grid.org[0] - grid.h*0.5);
	org[1]= (float)(grid.org[1] - grid.h*0.5);
	org[2]= (float)(grid.org[2] - grid.h*0.5);

	nx3 = (nx + 2*nborder);
	ny3 = (ny + 2*nborder);
	nz3 = (nz + 2*nborder);
	nyz3 = ny3*nz3;

	sprintf(s1,"tmp_traveltimes");
	sprintf(s3,".dat");
	/*******************************************************************************************/
	/*******************************************************************************************/
	/*--Fat-rays--*/

		for(i=0;i<data->ndata_seis;i++)
			data->tcalc[i] = -1.0;

		/*Allocate memory for the travel-times (coming from the shots/receivers) that will be calculated by the forward algorithm*/
		tt_shot = (float *)memory(NULL,nx3*ny3*nz3,sizeof(float),"ForwardFatRayMod");
		/*Allocate memory for temporary slowness field in FLOAT*/
        tmp_slow = (float *)memory(NULL,nx3*ny3*nz3,sizeof(float),"ForwardFatRayMod");
		/*Allocate the memory for a vector of all receiver positions for which the traveltimes are already calculated*/
		rec_counter = (long *)memory(NULL, 1, sizeof(long),"ForwardModFatRay");
		counter_index = 0;

		/*Allocate the memory for the number of shots (for each receiver) for which the calculation is already performed*/
		rec_num_counter = (long *)memory(NULL, geo.nstat+1, sizeof(long),"ForwardModFatRay");
		for(i=0;i<geo.nstat+1;i++)
			rec_num_counter[i] = 0;

		for(j=0;j<nx3*ny3*nz3;j++)
			tmp_slow[j] = (float)grid.slow[j];

		if(data->ndata_seis == 0)
		{
			printf("!!WARNING!! NO seismic traveltime data exists although the seismic forward modelling is activated\n");
			printf("Seismic forward modeling is therefore finished\n");
			printf("-------------------------------\n\n");
			data->rms_seis = -99999.9;
			return(1);
		}
		/****************************************************************************************/
		/*Start the loop over all shots*/
		for (i=0;i<geo.nshot;i++)
		{
			if(data->lshots[i] != 0)
			{
				for(j=0;j<nx3*ny3*nz3;j++)
					tt_shot[j] = 0.0;

				Xs = ((geo.x[(data->shots[i])-1] - org[0])/grid.h + (float)nborder); /*normalized x-coordinate of the shot locations according to grid cell nodes*/
				Ys = ((geo.y[(data->shots[i])-1] - org[1])/grid.h + (float)nborder); /*normalized y-coordinate of the shot locations according to grid cell nodes*/
				Zs = ((geo.z[(data->shots[i])-1] - org[2])/grid.h + (float)nborder); /*normalized z-coordinate of the shot locations according to grid cell nodes*/

				/***************************************************************************************/
				/*Podvin&Lecomte forward algorithm*/
				/*tt_shot is the calculated traveltime for each grid cell node*/
				time_3d(tmp_slow,tt_shot, nx3, ny3, nz3 ,Xs,Ys,Zs,delta_num,0);

				/***************************************************************************************/

				/*Determine the receivers that are activate for the corresponding shot:*/
				count = 0;
				nact_rec =(int *)memory(NULL,1,sizeof(int),"ForwardModFatRay");
				nact_datapos = (long *)memory(NULL,1,sizeof(long),"ForwardModFatRay");
				Xr = (float *) memory(NULL,1,sizeof(float),"ForwardModFatRay");
				Yr = (float *) memory(NULL,1,sizeof(float),"ForwardModFatRay");
				Zr = (float *) memory(NULL,1,sizeof(float),"ForwardModFatRay");


				for(j=0; j<data->ndata_seis; j++)
				{
					if(data->shots[i] == data->sno[j])
					{
						nact_rec =(int *)memory((char *)nact_rec,count+1,sizeof(int),"ForwardModFatRay");
						nact_datapos = (long *)memory((char *)nact_datapos,count+1,sizeof(long),"ForwardModFatRay");
						Xr = (float *) memory((char *)Xr,count+1,sizeof(float),"ForwardModFatRay");
						Yr = (float *) memory((char *)Yr,count+1,sizeof(float),"ForwardModFatRay");
						Zr = (float *) memory((char *)Zr,count+1,sizeof(float),"ForwardModFatRay");

						nact_rec[count]		= data->rno[j];
						nact_datapos[count]	= j;

						count++;
					}
				}

				if(count != data->lshots[i])
				{
					printf("The number of active receivers of shot number %d is not correct\n", data->shots[i]);
					exit(0);
				}

				/***************************************************************************************/
				/*Determine the accurate traveltimes at the receiver-locations (by trilinear interpolation of the traveltimes at the grid cell edges)*/
				for(j=0;j<count;j++)
				{
					Xr[j] = ((geo.x[nact_rec[j]-1] - org[0])/grid.h + (float)nborder); /*normalized x-coordinate of the receiver locations according to grid cell nodes*/
					Yr[j] = ((geo.y[nact_rec[j]-1] - org[1])/grid.h + (float)nborder); /*normalized y-coordinate of the receiver locations according to grid cell nodes*/
					Zr[j] = ((geo.z[nact_rec[j]-1] - org[2])/grid.h + (float)nborder); /*normalized z-coordinate of the receiver locations according to grid cell nodes*/

					data->tcalc[nact_datapos[j]] = (double)(1000.0 * interpolate(Xr[j],Yr[j],Zr[j],&grid,tt_shot));

					if(nact_datapos[j] >= data->ndata_seis)
					 {
							printf("NOT enough memory is allocated: used %d, allocated %d\n",nact_datapos[j]+1, data->ndata_seis);
							exit(0);
					 }
				}



				time(&time_relative);
				zeit_angabe_sekunden = difftime(time_relative, time_start);

				printf(" For shot-nr. %d all traveltimes are calculated\n",data->shots[i]);
				printf("   (x=%f,y=%f,z=%f)\n",geo.x[(data->shots[i])-1], geo.y[(data->shots[i])-1], geo.z[(data->shots[i])-1]);
				printf("   Number of found receiver positions for the shot: %d\n",count);
				printf("   Computing time:%fh\n\n",zeit_angabe_sekunden/3600);


				/****************************************************************************************/
				/****************************************************************************************/
				/*Start the loop over the receivers belonging to a shot*/
				for(j=0;j<count;j++)
				{
					/*Check, if the traveltimes were already calculated*/
					exist = 0;
					for(k=0;k<counter_index;k++)
					{
						if(nact_rec[j] == rec_counter[k])
						{
							exist = 1;
							break;
						}
					}

					if(exist == 1)	/*The velocity model for the receiver has been already calculated*/
					{
						strcpy(s1_tmp,s1);
						sprintf(s2_tmp,"%d",nact_rec[j]);
						strcpy(s3_tmp,s3);
						strncat(s1_tmp,s2_tmp,10);
						strncat(s1_tmp,s3_tmp,4);

						/*Read in the traveltimes of the rec.location from the corresponding (temporarly) binary file*/
						tt_rec = ReadTravelTimeModIn(s1_tmp);
					}
					else
					{
						tt_rec = (float *)memory(NULL,nx3*ny3*nz3,sizeof(float),"ForwardFatRayMod");
						/***************************************************************************************/
						/*Podvin&Lecomte forward algorithm*/
						/*tt_rec is the calculated traveltime for each grid cell edge*/
						 time_3d(tmp_slow,tt_rec, nx3, ny3, nz3 ,Xr[j],Yr[j],Zr[j],delta_num,0);

						/***************************************************************************************/


						 /*Add the receiver to the "list of calculated traveltime models"*/
						 counter_index++;
						 rec_counter = (long *)memory((char *)rec_counter, counter_index, sizeof(long),"ForwardModFatRay");
					     rec_counter[counter_index - 1] = nact_rec[j];


						 time(&time_relative);
						 zeit_angabe_sekunden = difftime(time_relative, time_start);

						 printf(" For rec-nr. %d all traveltimes are calculated\n",nact_rec[j]);
						 printf("   (x=%f,y=%f,z=%f)\n",geo.x[nact_rec[j]-1], geo.y[nact_rec[j]-1], geo.z[nact_rec[j]-1]);
						 printf("   Computing time:%fh\n\n",zeit_angabe_sekunden/3600);

						/*Write the traveltimes from the receiver-position in a binary file:*/
						 WriteTravelTimeModOut(tt_rec,nact_rec[j],nx3,ny3,nz3);
					}
					/*********************************************************************************/
					/*Calculation of the fat rays:*/
					FatRayCalc(tt_shot,tt_rec,nx3,ny3,nz3,(float)data->tcalc[nact_datapos[j]],&fat_rays[nact_datapos[j]], grid.slow);
					fat_rays[nact_datapos[j]].n = nact_datapos[j];

					/*Check, if the receiver-position is used for following shots*/
					rec_num_counter[nact_rec[j]]++;


					/*Take the number of active shots for the considered receiver from the data structure*/
					max_num = 0;
					for(k=0;k<geo.nrec;k++)
					{
						if(nact_rec[j] == data->recs[k])
						{
							max_num = data->lrecs[k];
							break;
						}
					}

					/*This is the last shot location for which the receiver location was active*/
					if(rec_num_counter[nact_rec[j]] == max_num)
					{
						strcpy(s1_tmp,s1);
						sprintf(s2,"%d",nact_rec[j]);
						strcpy(s3_tmp,s3);
						strncat(s1_tmp,s2,10);
						strncat(s1_tmp,s3_tmp,4);

						remove(s1_tmp); /*Remove the file which will not be used any more*/

						printf("The temporary file %s is removed again\n\n",s1_tmp);
					}
					/**********************************************************************************/


					free(tt_rec);

				}
				/*End of the loop over the receivers belonging to a shot*/
				/****************************************************************************************/
				/****************************************************************************************/

				free(Xr);
				free(Yr);
				free(Zr);
				free(nact_rec);
				free(nact_datapos);

				printf("All rays for the shot %d are calculated\n",data->shots[i]);
				printf("----------------\n\n");
			}
		}
		/*End of the loop over all shots*/
		/****************************************************************************************/

		free(rec_counter);
		free(rec_num_counter);
		free(tt_shot);
		free(tmp_slow);

	/*******************************************************************************************/
	/*******************************************************************************************/
	/*Check the modified structures*/
	for(i=0;i<data->ndata_seis;i++)
	{
		if(data->tcalc[i] == -1.0)
		{
			printf("For the shot-receiver combination %d no traveltime was calculated\n->Check the program\n",i+1);
			exit(0);
		}

		if(fat_rays[i].ncell%1 != 0)
		{
			printf("For the shot-receiver combination %d no raypath was calculated\n->Check the program\n",i+1);
			exit(0);
		}
	}

	/****************************************************************************/
	/*Write out the fat-rays*/
	WriteFatRayOut(fat_rays, *data, grid);

	/*Determine the RMS-value*/
	sum2 = 0.0;
	for(i=0;i<data->ndata_seis;i++)
			sum2 = sum2 + (data->tcalc[i] - data->tobs[i])*(data->tcalc[i] - data->tobs[i]);

	if(data->ndata_seis != 0)
		data->rms_seis = sqrt(sum2/data->ndata_seis);
	else
		data->rms_seis = -99999.9;

	printf("Seismic forward modeling is finished:\n");
	printf("RMS-values: %10.5f ms\n",data->rms_seis);
	printf("----------------\n\n\n\n\n");

	return(1);
}

/*-------------------------------------------------------------*/
/*Calculate the fat-rays*/
/*Parameter:	*tt1,*tt2	:= traveltime field, calculated from a shot/receiver position */
/*              nx,ny,nz	:= Size of the grid; number of nodes*/
/*              tcalc		:= Traveltimes from a shot to a receiver*/
/*				fr			:= Fat-ray structure */
/*              slow		:= slowness in the grid cells (Remark: The size of the array is extended in every direction with one dummy cell)*/

int FatRayCalc(float *tt1,float *tt2, long nx, long ny,long nz, float tcalc,F_RP_STRUCT *fr, double *slow)
{
	long a,b,c,nxyz,nyz,nyz1,nx1,ny1,nz1,nxyz1,count;
	double dval; /*traveltime differences in ms*/
	double *tmp_fr_weight,dsum, sum_all;

	nx1=nx-1;
	ny1=ny-1;
	nz1=nz-1;

	nyz=ny*nz;
	nyz1 = ny1*nz1;
	nxyz=nx*ny*nz;
	nxyz1=nx1*ny1*nz1;

	tmp_fr_weight = (double *)memory(NULL,nxyz,sizeof(double),"FatRayCalc");
	fr->weight = (double *)memory(NULL,1,sizeof(double),"FatRayCalc");
	fr->ele = (long *)memory(NULL,1,sizeof(long),"FatRayCalc");

	/*Calculate the difference between the first-arrivals and the traveltimes to the grid nodes from the shot and the receiver positions*/
	sum_all = 0.0;

	for(a=0;a<nxyz;a++)
	{
		dval = fabs(tt1[a] + tt2[a] - (tcalc/1000));

		/*Weights of the fat rays at the nodes*/
		if(dval < (fr->fatthres/1000))
		{
			tmp_fr_weight[a] = exp(-dval*(fr->fatb*1000));
			sum_all = sum_all + tmp_fr_weight[a];
		}
		else
			tmp_fr_weight[a] = 0.0;
	}

	if(sum_all == 0.0)
		printf("WARNING !!! NO cell contributes to the fat-ray!!!\nThe size of the fat-ray should be larger choosen!!!\n\n");

	/*"Interpolate" the data (mean value) from the nodes to determine the value in the center*/
	#define tmp_fr(x,y,z) tmp_fr_weight[nyz*(x) + nz*(y) + (z)]
	#define slowness(x,y,z) slow[nyz*(x) + nz*(y) + (z)]

	count= 0;
	dsum = 0;

	for(a=0;a<nx1;a++)
		for(b=0;b<ny1;b++)
			for(c=0;c<nz1;c++)
			{
				dval  = (tmp_fr(a,b,c) + tmp_fr(a,b,c+1) + tmp_fr(a,b+1,c) + tmp_fr(a,b+1,c+1) +
					     tmp_fr(a+1,b,c) + tmp_fr(a+1,b,c+1) + tmp_fr(a+1,b+1,c) + tmp_fr(a+1,b+1,c+1))/8.0;



				if(dval>0.0)
				{
					fr->weight = (double *)memory((char *) fr->weight,count+1,sizeof(double),"FatRayCalc");
					fr->ele = (long *)memory((char *) fr->ele,count+1,sizeof(long),"FatRayCalc");

					fr->weight[count] = dval;			 /*Value of the weight*/
					fr->ele[count] = nyz1*a + nz1*b + c;	 /*corresponding cell number*/


					dsum = dsum + (slowness(a,b,c)*dval);	/*Calculate the normalization*/

					count++;
				}
			}

	fr->ncell = count; /*Number of cells in the fat-ray*/

	/*Normalize the fat-ray:*/
	for(a=0;a<count;a++)
		(fr->weight[a]) = ((fr->weight[a])*(tcalc/1000))/dsum;

	free(tmp_fr_weight);

	#undef tmp_fr
	#undef slowness

	return(1);
}




/*-------------------------------------------------------------*/
/*Performing the forward modeling(using the Podvin&Lecomte eikonal solver) and calculating conventional rays */
/*Parameter:	 geo  := Geometry structure  */
/*               grid  := Grid structure*/
/*              *data  := Pointer on data structure (the calculated traveltimes will be (re-)determined in this routine)*/
/*				*raypaths := Vector of the raypath structures (number of elements = number of shot-receiver combinations)*/
/*              time_start := Time when inv3d starts*/

/*REMARK: The Podvin&Lecomte forward modeling calculates the traveltimes at the grid cell nodes; The number of grid cell nodes is in each*/
/*		direction one higher than the number of grid cell centers.*/
/*		Therefore the grid have to be readjusted*/

#define travel_t(x,y,z) tt[nyz3*(x) + nz3*(y) + (z)]

int ForwardModRay(GEOMETRY geo,GRID_STRUCT grid, DATA_STRUCT *data, RP_STRUCT *raypath , time_t time_start)
{
	int /*a,b,c,*/count/*Number of active receivers for a shot*/;
	int i,k,nborder;
	int *nact_rec; /*active receiver-numbers for the used shot*/
	long *nact_datapos; /*Position of the active receivers in the data structure*/
	long j;
	long nx,ny,nz,nx3,ny3,nz3,nyz3;
	float *tt,*tmp_slow;
	float Xs,Ys,Zs,*Xr,*Yr,*Zr;		/*Normalized positions of the shots and receivers (referring to the grid cell nodes and NOT of the grid cell centers)*/
	float org[3];
	float delta_num = (float)0.001;
	double sum2;

	time_t time_relative;
	double zeit_angabe_sekunden;
	FILE *out;

	RP_STRUCT *raypath_tmp;


	nx= grid.nx +1;
	ny= grid.ny +1;
	nz= grid.nz +1;
	nborder =grid.nborder;

	org[0]= (float)(grid.org[0] - grid.h*0.5);
	org[1]= (float)(grid.org[1] - grid.h*0.5);
	org[2]= (float)(grid.org[2] - grid.h*0.5);

	nx3 = (nx + 2*nborder);
	ny3 = (ny + 2*nborder);
	nz3 = (nz + 2*nborder);
	nyz3 = ny3*nz3;

	/*******************************************************************************************/
	/*******************************************************************************************/
	/*--Conventional-rays--*/

		for(i=0;i<data->ndata_seis;i++)
			data->tcalc[i] = -1.0;

		/*Allocate memory for the travel-times that will be calculated by the forward algorithm*/
		tt = (float *)memory(NULL,nx3*ny3*nz3,sizeof(float),"ForwardModRay");


		if(data->ndata_seis == 0)
		{
			data->ndata_seis_act = 0;

			printf("!!WARNING!! NO seismic traveltime data exists although the seismic forward modelling is activated\n\n");
			printf("Seismic forward modeling is therefore finished\n\n");
			printf("-------------------------------\n\n");
			data->rms_seis = -99999.9;
			return(1);
		}

		/*Allocate memory for temporary slowness field in FLOAT*/
		tmp_slow = (float *)memory(NULL,nx3*ny3*nz3,sizeof(float),"ForwardModRay");

		for(i=0;i<nx3*ny3*nz3;i++)
			tmp_slow[i] = 0.0;


		for(i=0;i<nx3-1;i++)
			for(j=0;j<ny3-1;j++)
				for(k=0;k<nz3-1;k++)
				{
					tmp_slow[i*nyz3 + j*nz3 +k]=(float)grid.slow[i*nyz3 + j*nz3 +k];
				}


		/*Start the loop over all shots*/
		for (i=0;i<geo.nshot;i++)
		{
			if(data->lshots[i] != 0)
			{
				for(j=0;j<nx3*ny3*nz3;j++)
					tt[j] = 0.0;

				Xs = ((geo.x[(data->shots[i])-1] - org[0])/grid.h + (float)nborder); /*normalized x-coordinate of the shot locations according to grid cell nodes*/
				Ys = ((geo.y[(data->shots[i])-1] - org[1])/grid.h + (float)nborder); /*normalized y-coordinate of the shot locations according to grid cell nodes*/
				Zs = ((geo.z[(data->shots[i])-1] - org[2])/grid.h + (float)nborder); /*normalized z-coordinate of the shot locations according to grid cell nodes*/

				/***************************************************************************************/
				/*Podvin&Lecomte forward algorithm*/
				/*tt is the calculated traveltime for each grid cell node*/
				 time_3d(tmp_slow,tt, nx3, ny3, nz3 ,Xs,Ys,Zs,delta_num,0);
				/***************************************************************************************/

				/*Determine the receivers that are activate for the corresponding shot:*/
				count = 0;
				nact_rec =(int *)memory(NULL,1,sizeof(int),"ForwardModRay");
				nact_datapos = (long *)memory(NULL,1,sizeof(long),"ForwardModRay");
				Xr = (float *) memory(NULL,1,sizeof(float),"ForwardModRay");
				Yr = (float *) memory(NULL,1,sizeof(float),"ForwardModRay");
				Zr = (float *) memory(NULL,1,sizeof(float),"ForwardModRay");


				for(j=0; j<data->ndata_seis; j++)
				{
					if(data->shots[i] == data->sno[j])
					{
						nact_rec =(int *)memory((char *)nact_rec,count+1,sizeof(int),"ForwardModRay");
						nact_datapos = (long *)memory((char *)nact_datapos,count+1,sizeof(long),"ForwardModRay");
						Xr = (float *) memory((char *)Xr,count+1,sizeof(float),"ForwardModRay");
						Yr = (float *) memory((char *)Yr,count+1,sizeof(float),"ForwardModRay");
						Zr = (float *) memory((char *)Zr,count+1,sizeof(float),"ForwardModRay");

						nact_rec[count]		= data->rno[j];
						nact_datapos[count]	= j;

						count++;
					}
				}

				if(count != data->lshots[i])
				{
					printf("The number of active receivers of shot number %d is not correct\n", data->shots[i]);
					exit(0);
				}

				/***************************************************************************************/
				/*Determine the accurate traveltimes at the receiver-locations (by trilinear interpolation of the traveltimes at the grid cell edges)*/
				for(j=0;j<count;j++)
				{
					Xr[j] = ((geo.x[nact_rec[j]-1] - org[0])/grid.h + (float)nborder); /*normalized x-coordinate of the receiver locations according to grid cell EDGES*/
					Yr[j] = ((geo.y[nact_rec[j]-1] - org[1])/grid.h + (float)nborder); /*normalized y-coordinate of the receiver locations according to grid cell EDGES*/
					Zr[j] = ((geo.z[nact_rec[j]-1] - org[2])/grid.h + (float)nborder); /*normalized z-coordinate of the receiver locations according to grid cell EDGES*/

					 data->tcalc[nact_datapos[j]] = (double)(1000.0 * interpolate(Xr[j],Yr[j],Zr[j],&grid,tt));

					 if(nact_datapos[j] >= data->ndata_seis)
					 {
							printf("NOT enough memory is allocated: used %d, allocated %d\n",nact_datapos[j]+1, data->ndata_seis);
							exit(0);
					 }

				}

				time(&time_relative);
				zeit_angabe_sekunden = difftime(time_relative, time_start);

				printf(" For shot-nr. %d all traveltimes are calculated\n",data->shots[i]);
				printf("   (x=%f,y=%f,z=%f)\n",geo.x[(data->shots[i])-1], geo.y[(data->shots[i])-1], geo.z[(data->shots[i])-1]);
				printf("   Number of found receiver positions for the shot: %d\n",count);
				printf("   Computing time:%fh\n\n",zeit_angabe_sekunden/3600);

				/***************************************/
		   		/*Write out the traveltime cube*/
				/*out = fopen("tcalc.txt","wt");
					for(c=0; c<nz3;c++)
						for(b=0; b<ny3;b++)
							for(a=0; a<nx3; a++)
								fprintf(out,"%f\n",travel_t(a,b,c));
				fclose(out);*/
				/***************************************/


				/***************************************************************************************/
				/*Ray calculations (3-D version of the Aldridge&Oldenburg raypath-generation, 1993, Journal of seismic exploration,Vol.2,pages 257-274)*/

				/*Allocate temporary ray-structures for the specific shot*/
				raypath_tmp =(RP_STRUCT *)memory(NULL,data->lshots[i],sizeof(RP_STRUCT),"ForwardModRay");

				for(j=0;j<count;j++)
					raypath_tmp[j].n = nact_datapos[j];

				/*Calculate the rays*/
					RayCalc(tt,nx3,ny3,nz3, Xs,Ys,Zs,Xr,Yr,Zr,count,raypath_tmp);

				/*Copy the temporary raypath structures in structures that fit with the data structure*/
				for(j=0;j<count;j++)
				{
					raypath[nact_datapos[j]].n = raypath_tmp[j].n;
					raypath[nact_datapos[j]].nray = raypath_tmp[j].nray;

					if(raypath[nact_datapos[j]].nray != 0)
					{
						raypath[nact_datapos[j]].len = (double *) memory(NULL,raypath[nact_datapos[j]].nray,sizeof(double),"ForwardModRay");
						raypath[nact_datapos[j]].ele = (long *) memory(NULL,raypath[nact_datapos[j]].nray,sizeof(long),"ForwardModRay");
					}
					else
					{
						raypath[nact_datapos[j]].len = (double *) memory(NULL,1,sizeof(double),"ForwardModRay");
						raypath[nact_datapos[j]].ele = (long *) memory(NULL,1,sizeof(long),"ForwardModRay");
					}
					raypath[nact_datapos[j]].x = (double *) memory(NULL,raypath[nact_datapos[j]].nray + 1,sizeof(double),"ForwardModRay");
					raypath[nact_datapos[j]].y = (double *) memory(NULL,raypath[nact_datapos[j]].nray + 1,sizeof(double),"ForwardModRay");
					raypath[nact_datapos[j]].z = (double *) memory(NULL,raypath[nact_datapos[j]].nray + 1,sizeof(double),"ForwardModRay");

					for(k=0;k<raypath_tmp[j].nray;k++)
					{
						raypath[nact_datapos[j]].len[k] = raypath_tmp[j].len[k];
						raypath[nact_datapos[j]].ele[k] = raypath_tmp[j].ele[k];
						raypath[nact_datapos[j]].x[k] = raypath_tmp[j].x[k];
						raypath[nact_datapos[j]].y[k] = raypath_tmp[j].y[k];
						raypath[nact_datapos[j]].z[k] = raypath_tmp[j].z[k];
					}
					raypath[nact_datapos[j]].x[k] = raypath_tmp[j].x[k];
					raypath[nact_datapos[j]].y[k] = raypath_tmp[j].y[k];
					raypath[nact_datapos[j]].z[k] = raypath_tmp[j].z[k];
				}

				for(j=0;j<count;j++)
				{
					free(raypath_tmp[j].len);
					free(raypath_tmp[j].x);
					free(raypath_tmp[j].y);
					free(raypath_tmp[j].z);
					free(raypath_tmp[j].ele);
				}

				free(raypath_tmp);

				/****************************************************************************************/

				free(Xr);
				free(Yr);
				free(Zr);
				free(nact_rec);
				free(nact_datapos);

			}
		}

		/*End of the loop over all shots*/
		free(tmp_slow);
		free(tt);

	/*******************************************************************************************/
	/*******************************************************************************************/
	/*Check the modified structures*/
	for(i=0;i<data->ndata_seis;i++)
	{
		if(data->tcalc[i] == -1.0)
		{
			printf("For the shot-receiver combination %d no traveltime was calculated\n->Check the program\n",i+1);
			exit(0);
		}

		if(raypath[i].nray%1 != 0)
		{
			printf("For the shot-receiver combination %d no raypath was calculated\n->Check the program\n",i+1);
			exit(0);
		}
	}


	/***************************************/
	/*Write out the calculated traveltime cube*/
	 out = fopen("tcalc.txt","wt");
	 fprintf(out,"Shot-number    Tcalc(ms)\n");
		for(i=0; i<data->ndata_seis;i++)
					fprintf(out,"%d %d %d %f\n", i ,data->sno[i], data->rno[i], data->tcalc[i]);
	fclose(out);
	/***************************************/

	/****************************************************************************/
	/*Resort rays (If the ray runs along the boundary between to cells, the cell with the higher velocity will be considered)*/
	ResortRays(raypath,*data, grid);

	/*Write out the raypath*/
	WriteRayOut(raypath, *data, grid);

	/*Determine the number of ACTIVE rays and the RMS-value*/
	data->ndata_seis_act = 0;
	sum2 = 0.0;
	for(i=0;i<data->ndata_seis;i++)
	{
		if(raypath[i].nray != 0)
		{
			data->ndata_seis_act++;
			sum2 = sum2 + (data->tcalc[i] - data->tobs[i])*(data->tcalc[i] - data->tobs[i]);
		}
	}

	if(data->ndata_seis_act != 0)
		data->rms_seis = sqrt(sum2/data->ndata_seis_act);
	else
		data->rms_seis = -99999.9;

	printf("Seismic forward modeling is finished:\n");
	printf("RMS-values: %10.5f ms\n",data->rms_seis);
	printf("%d of %d shot-receiver combinations will be used for inversion\n", data->ndata_seis_act, data->ndata_seis);
	printf("----------------\n\n\n\n\n");


	return(1);
}

#undef travel_t

/*-------------------------------------------------------------*/
/*Bilinear interpolation in 3-D to determine the velocity from receiver/shot positions to next grid point:  */
/*Parameter:	x,y,z  := Position of the receiver/shots in "grid cells"  */
/*              *grid  := Pointer on the grid structure*/
/*              *data  := Data (usually velocities) from which the interpolated velocities are determined */
/*Output:  Determined velocity from receiver/shot positions to next grid point */


#define lo(val)    (int)floor((val))
#define hi(val)    (int)ceil((val))
#define dd(x,y,z)   data[nyz2*(x) + nz2*(y) + (z)]

float interpolate(float x,float y,float z,GRID_STRUCT *grid,float *data)
{
   float u,v,w;
   int ok,nx2,ny2,nz2,nyz2;
   float ival;

   nx2 = grid->nx +1 + 2*grid->nborder;
   ny2 = grid->ny +1 + 2*grid->nborder;
   nz2 = grid->nz +1 + 2*grid->nborder;
   nyz2 = ny2*nz2;

   /* Check, if point is in grid */
   ok = 1;
   if ( (x > nx2) || (x < 0) ) ok = 0;
   if ( (y > ny2) || (y < 0) ) ok = 0;
   if ( (z > nz2) || (z < 0) ) ok = 0;
   if (!ok)
   {
      printf("\nInterpolation point is out of the grid!\n");
      printf("x = %f y = %f z = %f\n",x,y,z);
      printf("nx = %d ny = %d nz = %d\n",grid->nx,grid->ny,grid->nz);
      printf("h = %f\n",grid->h);
      printf("orix = %f oriy = %f oriz = %f\n",grid->org[0],grid->org[1],grid->org[2]);
      exit(0);
   }

   /* Get interpolation distances */
   u = x - (float)floor(x);
   v = y - (float)floor(y);
   w = z - (float)floor(z);

   /* And now interpolate */
   ival = (1-u)*(1-v)*(1-w)*dd(lo(x),lo(y),lo(z)) +
          (u)*(1-v)*(1-w)*dd(hi(x),lo(y),lo(z)) +
          (u)*(v)*(1-w)*dd(hi(x),hi(y),lo(z)) +
          (1-u)*(v)*(1-w)*dd(lo(x),hi(y),lo(z)) +

          (1-u)*(1-v)*(w)*dd(lo(x),lo(y),hi(z)) +
          (u)*(1-v)*(w)*dd(hi(x),lo(y),hi(z)) +
          (u)*(v)*(w)*dd(hi(x),hi(y),hi(z)) +
          (1-u)*(v)*(w)*dd(lo(x),hi(y),hi(z));

   return(ival);
}

#undef lo
#undef hi
#undef dd

/*-------------------------------------------------------------*/
/*Calculation of the ray paths by the means of maximum traveltime descent*/
/*Parameter:	*tt       := Traveltimes at the grid cell edges  */
/*              nx,ny,nz  := Size of the grid*/
/*              Xs,Ys,Zs  := Normalized shot positions*/
/*              Xr,Yr,Zr  := Normalized receiver positions */
/*              nrec      := number of receivers per shot */
/*              *rp       := vector of the ray path structure*/

#define cell_index(x,y,z) ray_cell_index[nyz1*(x) + (nz1)*(y) + (z)]

int RayCalc(float *tt, int nx,int ny,int nz, float Xs, float Ys, float Zs,float *Xr,float *Yr, float *Zr,int nrec, RP_STRUCT *rp)
{
    int a,b,c;
	int i,count;
	int nx1,ny1,nz1;
	long nyz1;
	int *ray_cell_index; /*if 0=no ray in the cell; 1= ray path found in the cell*/
	long max_nr_of_ray_seg; /*max. number of ray-segments*/

	double *gradient; /*Components of the gradient*/
	CELL_STRUCT next_cell,cell;

	nx1=nx-1; /* nx:Number of nodes; nx1= number of cells in x-direction*/
	ny1=ny-1;
	nz1=nz-1;
	nyz1=ny1*nz1;
	max_nr_of_ray_seg = 2*(nx1+ny1+nz1);


	for(i=0;i<nrec;i++)
	{

		ray_cell_index = (int *)memory(NULL,(nx1)*(ny1)*(nz1),sizeof(int),"RayCalc");

		/*Set "boundary-index" to 1 for all grid cells:*/
		for(a=0; a<nx1; a++)
			for(b=0; b<ny1;b++)
				for(c=0; c<nz1;c++)
				{
					cell_index(a,b,c) = 0;
				}

		/****************************************************/
		/*Determine the rays starting at the receiver position:*/

		/*Check, if shot and receiver are in the same cell*/
		if(floor((double)Xs) == floor((double)Xr[i]) && floor((double) Ys) == floor((double) Yr[i]) && floor((double)Zs) == floor((double) Zr[i]))
		{
			rp[i].len = (double *)memory(NULL,1,sizeof(double),"RayCalc");
			rp[i].x = (double *)memory(NULL,2,sizeof(double),"RayCalc");
			rp[i].y = (double *)memory(NULL,2,sizeof(double),"RayCalc");
			rp[i].z = (double *)memory(NULL,2,sizeof(double),"RayCalc");
			rp[i].ele = (long *)memory(NULL,1,sizeof(long),"RayCalc");

			rp[i].len[0] = sqrt((double)(Xs-Xr[i])*(Xs-Xr[i]) + (Ys-Yr[i])*(Ys-Yr[i]) + (Zs-Zr[i])*(Zs-Zr[i])); /*Ray segment length*/
			rp[i].x[0] = (double) Xr[i];
			rp[i].y[0] = (double) Yr[i];
			rp[i].z[0] = (double) Zr[i];
			rp[i].x[1] = (double) Xs;
			rp[i].y[1] = (double) Ys;
			rp[i].z[1] = (double) Zs;

			cell.xno = (int)floor((double)Xr[i]);
			cell.yno = (int)floor((double)Yr[i]);
			cell.zno = (int)floor((double)Zr[i]);

			rp[i].ele[0] = nyz1*cell.xno + nz1*cell.yno + cell.zno; /*Determine the cell number, where the ray intersects*/
			rp[i].nray = 1;	 /*Number of segments of the ray*/

			goto fertig;
		}

		/*otherwise ...*/
		/*Set the first cell structure of the considered ray:*/
		cell.xno = (int)floor((double)Xr[i]);
		cell.yno = (int)floor((double)Yr[i]);
		cell.zno = (int)floor((double)Zr[i]);
		cell.xpos = Xr[i] - cell.xno;
		cell.ypos = Yr[i] - cell.yno;
		cell.zpos = Zr[i] - cell.zno;
		cell.dirx_i = 0;
		cell.diry_i = 0;
		cell.dirz_i = 0;

		/*Calculate the traveltime gradient*/
		gradient = TimeGrad(cell.xno,cell.yno,cell.zno,tt,ny,nz);

		/*Calculate the ray segment through the first grid cell*/
		 next_cell = RayBackTrace(gradient[0],gradient[1],gradient[2], cell, tt, ny, nz);
		 free(gradient);

				rp[i].len = (double *)memory(NULL,1,sizeof(double),"RayCalc");
				rp[i].x = (double *)memory(NULL,1,sizeof(double),"RayCalc");
				rp[i].y = (double *)memory(NULL,1,sizeof(double),"RayCalc");
				rp[i].z = (double *)memory(NULL,1,sizeof(double),"RayCalc");
				rp[i].ele = (long *)memory(NULL,1,sizeof(long),"RayCalc");

				rp[i].len[0] = sqrt((next_cell.xpos + next_cell.xno - cell.xpos -cell.xno)*(next_cell.xpos + next_cell.xno - cell.xpos - cell.xno) + (next_cell.ypos + next_cell.yno - cell.ypos -cell.yno)*(next_cell.ypos + next_cell.yno - cell.ypos - cell.yno) + (next_cell.zpos + next_cell.zno - cell.zpos - cell.zno)*(next_cell.zpos + next_cell.zno - cell.zpos -cell.zno)); /*Ray segment length*/
				rp[i].x[0] = (double) Xr[i];
				rp[i].y[0] = (double) Yr[i];
				rp[i].z[0] = (double) Zr[i];
				rp[i].ele[0] = nyz1*cell.xno + nz1*cell.yno + cell.zno; /*Determine the position number of the cell, which the ray intersects*/

				cell_index(cell.xno,cell.yno,cell.zno) = 1;

				/*Check, if the ray leave the cell*/
				if(next_cell.xno == 0 || next_cell.xno == nx1 || next_cell.yno == 0 || next_cell.yno == ny1 || next_cell.zno == 0 || next_cell.zno == nz1)
				{
					printf("The ray from the shot-receiver\ncombination %d leaves the model\n\n", rp[i].n +1);
					rp[i].nray = 0;
					goto fertig;
				}

				count =1;
				/******************************************************/
				/*Calculate iteratively the ray segments through the grid*/

				while(fabs(Xs - next_cell.xpos - next_cell.xno) > 1 || fabs(Ys - next_cell.ypos - next_cell.yno) > 1 || fabs(Zs -next_cell.zpos - next_cell.zno) > 1)
				{
						cell = next_cell;

						/*Calculate the traveltime gradient*/
						gradient = TimeGrad(cell.xno,cell.yno,cell.zno,tt,ny,nz);

						/*Calculate the ray segment through the corresponding grid cell*/
						next_cell = RayBackTrace(gradient[0],gradient[1],gradient[2], cell, tt, ny, nz);
						free(gradient);

						rp[i].len = (double *)memory((char *)rp[i].len,count+1,sizeof(double),"RayCalc");
						rp[i].x = (double *)memory((char *)rp[i].x,count+1,sizeof(double),"RayCalc");
						rp[i].y = (double *)memory((char *)rp[i].y,count+1,sizeof(double),"RayCalc");
						rp[i].z = (double *)memory((char *)rp[i].z,count+1,sizeof(double),"RayCalc");
						rp[i].ele = (long *)memory((char *)rp[i].ele,count+1,sizeof(long),"RayCalc");

						rp[i].len[count] = sqrt((next_cell.xpos + next_cell.xno - cell.xpos -cell.xno)*(next_cell.xpos + next_cell.xno - cell.xpos - cell.xno) + (next_cell.ypos + next_cell.yno - cell.ypos -cell.yno)*(next_cell.ypos + next_cell.yno - cell.ypos - cell.yno) + (next_cell.zpos + next_cell.zno - cell.zpos - cell.zno)*(next_cell.zpos + next_cell.zno - cell.zpos -cell.zno)); /*Ray segment length*/
						rp[i].x[count] = (double) cell.xpos + cell.xno;
						rp[i].y[count] = (double) cell.ypos + cell.yno;
						rp[i].z[count] = (double) cell.zpos + cell.zno;
						rp[i].ele[count] = nyz1*cell.xno + nz1*cell.yno + cell.zno; /*Determine the position number of the cell, which the ray intersects*/

						cell_index(cell.xno,cell.yno,cell.zno) = 1;

						/*Check, if the ray leave the cell*/
						if(next_cell.xno == 0 || next_cell.xno == nx1 || next_cell.yno == 0 || next_cell.yno == ny1 || next_cell.zno == 0 || next_cell.zno == nz1)
						{
							printf("The ray from the shot-receiver\ncombination %d leaves the model\n\n", rp[i].n+1);
							rp[i].nray = 0;
							goto fertig;
						}

						/*The ray will not be traced back, if the number of ray segments become too large; the ray is "fallen" probably in a local minima*/
						if(count >= max_nr_of_ray_seg)
						{
							printf("The discretized traveltime field of the ray from the shot-receiver\ncombination %d had a probably local minima\n\n",rp[i].n+1);
							rp[i].nray = 0;
							goto fertig;
						}


						count++;
				}

				/*Determine the last ray segment to the shot:*/
				rp[i].len = (double *)memory((char *)rp[i].len,count+1,sizeof(double),"RayCalc"); /*Normalized by the cell length*/
				rp[i].x = (double *)memory((char *)rp[i].x,count+2,sizeof(double),"RayCalc");
				rp[i].y = (double *)memory((char *)rp[i].y,count+2,sizeof(double),"RayCalc");
				rp[i].z = (double *)memory((char *)rp[i].z,count+2,sizeof(double),"RayCalc");
				rp[i].ele = (long *)memory((char *)rp[i].ele,count+1,sizeof(long),"RayCalc");

				rp[i].x[count] = (double) next_cell.xpos + next_cell.xno;
				rp[i].y[count] = (double) next_cell.ypos + next_cell.yno;
				rp[i].z[count] = (double) next_cell.zpos + next_cell.zno;

				rp[i].len[count] = sqrt((Xs - next_cell.xpos - next_cell.xno)*(Xs - next_cell.xpos - next_cell.xno) + (Ys - next_cell.ypos - next_cell.yno)*(Ys - next_cell.ypos - next_cell.yno) + (Zs - next_cell.zpos - next_cell.zno)*(Zs - next_cell.zpos - next_cell.zno)); /*Ray segment length*/
				rp[i].x[count+1] = (double) Xs;
				rp[i].y[count+1] = (double) Ys;
				rp[i].z[count+1] = (double) Zs;
				rp[i].ele[count] = nyz1*(int)floor((double)Xs) + nz1*(int)floor((double)Ys) + (int)floor((double)Zs); /*Determine the position number of the cell, which the ray intersects*/
				rp[i].nray = count+1; /*Number of the segments of the ray*/


		fertig:;

		free(ray_cell_index);

	}

	printf("All rays for the shot are calculated\n");
	printf("----------------\n\n");

	return(1);
}

#undef cell_index
/*-------------------------------------------------------------*/
/*Calculation the components of the time gradient*/
/*Parameter:	*tt			:= Traveltimes at the grid cell edges  */
/*				x,y,z	:= Position of the cell in the grid*/
/*				ny,nz	:= Size of the grid (number of grid cell EDGES)*/

/*	Output: Components of the gradient */

#define Traveltimes(a,b,c) tt[nyz*(a) + nz*(b) + (c)]

double *TimeGrad(int x, int y, int z, float *tt, int ny, int nz)
{
	int nyz;
	double *grad;

	nyz = ny*nz;

	grad = (double *)memory(NULL,3,sizeof(double),"TimeGrad");

	grad[0] = (- Traveltimes(x+1,y,z) - Traveltimes(x+1,y,z+1) - Traveltimes(x+1,y+1,z) - Traveltimes(x+1,y+1,z+1)
		      + Traveltimes(x,y,z) + Traveltimes(x,y,z+1) + Traveltimes(x,y+1,z) + Traveltimes(x,y+1,z+1))/4;		/*x-component*/

	grad[1] = (- Traveltimes(x,y+1,z) - Traveltimes(x,y+1,z+1) - Traveltimes(x+1,y+1,z) - Traveltimes(x+1,y+1,z+1)
		      + Traveltimes(x,y,z) + Traveltimes(x,y,z+1) + Traveltimes(x+1,y,z) + Traveltimes(x+1,y,z+1))/4;		/*y-component*/

	grad[2] = (- Traveltimes(x,y,z+1) - Traveltimes(x,y+1,z+1) - Traveltimes(x+1,y,z+1) - Traveltimes(x+1,y+1,z+1)
		      + Traveltimes(x,y,z) + Traveltimes(x,y+1,z) + Traveltimes(x+1,y,z) + Traveltimes(x+1,y+1,z))/4;		/*z-component*/

	return(grad);
}

#undef Traveltimes


/*-------------------------------------------------------------*/
/*Determine the ray segments through a cell in 3-D*/
/*Parameter:	gradx,grady,gradz	:= (negative) Traveltime-Gradient */
/*				cell				:= Cell structure of the actual cell*/
/*				*tt					:= traveltimes */
/*              ny,nz			:= grid size */

/*	Output: Cell structure of the next cell reached by the ray */

CELL_STRUCT RayBackTrace(double gradx, double grady, double gradz, CELL_STRUCT cell, float *tt,int ny,int nz)
{
	double eps = 0.01; /*Stabilize the program;*/
	double tmp_xpos, tmp_ypos, tmp_zpos, *gradient1, *gradient2;
	int diff;
	CELL_STRUCT next_cell;

	/*Safety settings to avoid instabilties (move the points sligthy away from the grid boundary)*/
	if(cell.dirx_i == 0 && cell.xpos == 0)
		cell.xpos = cell.xpos + eps;
	if(cell.dirx_i == 0 && cell.xpos == 1)
		cell.xpos = cell.xpos - eps;
	if(cell.diry_i == 0 && cell.ypos == 0)
		cell.ypos = cell.ypos + eps;
	if(cell.diry_i == 0 && cell.ypos == 1)
		cell.ypos = cell.ypos - eps;
	if(cell.dirz_i == 0 && cell.zpos == 0)
		cell.zpos = cell.zpos + eps;
	if(cell.dirz_i == 0 && cell.zpos == 1)
		cell.zpos = cell.zpos - eps;


	/*Ray enter the cell from the upper/lower x-plane*/
	/*CASE I: The ray run into the cell*/
	if(cell.dirx_i == 0 || (cell.dirx_i == 1 && gradx <= 0) || (cell.dirx_i == -1 && gradx >= 0))
	{
		if(gradx > 0)
		{

				tmp_ypos = cell.ypos - ((cell.xpos - 1)*(grady/gradx));
				tmp_zpos = cell.zpos - ((cell.xpos - 1)*(gradz/gradx));

			if(1 >= tmp_ypos && tmp_ypos >= 0 && 1 >= tmp_zpos && tmp_zpos >= 0)
			{
				next_cell.xpos = 0;
				next_cell.ypos = tmp_ypos;
				next_cell.zpos = tmp_zpos;
				next_cell.dirx_i = -1;
				next_cell.diry_i =  0;
				next_cell.dirz_i =  0;
				next_cell.xno = cell.xno + 1;
				next_cell.yno = cell.yno;
				next_cell.zno = cell.zno;

				goto bestimmt;

			}
		}

		if(gradx < 0)
		{
				tmp_ypos = cell.ypos - (cell.xpos *(grady/gradx));
				tmp_zpos = cell.zpos - (cell.xpos *(gradz/gradx));

			if( 1 >= tmp_ypos && tmp_ypos >= 0 && 1 >= tmp_zpos && tmp_zpos >= 0)
			{
				next_cell.xpos = 1;
				next_cell.ypos = tmp_ypos;
				next_cell.zpos = tmp_zpos;
				next_cell.dirx_i =  1;
				next_cell.diry_i =  0;
				next_cell.dirz_i =  0;
				next_cell.xno = cell.xno - 1;
				next_cell.yno = cell.yno;
				next_cell.zno = cell.zno;

				goto bestimmt;
			}

		}
	}


	/*Ray enter the cell from the positive or negative y-plane*/
	/*CASE I: The ray run into the cell*/
	if(cell.diry_i == 0 || (cell.diry_i == 1 && grady <= 0) || (cell.diry_i == -1 && grady >= 0))
	{
		if(grady > 0)
		{
				tmp_xpos = cell.xpos - ((cell.ypos - 1)*(gradx/grady));
				tmp_zpos = cell.zpos - ((cell.ypos - 1)*(gradz/grady));

			if(1 >= tmp_xpos && tmp_xpos >= 0 && 1 >= tmp_zpos && tmp_zpos >= 0)
			{
				next_cell.xpos = tmp_xpos;
				next_cell.ypos = 0;
				next_cell.zpos = tmp_zpos;
				next_cell.dirx_i =  0;
				next_cell.diry_i = -1;
				next_cell.dirz_i =  0;
				next_cell.xno = cell.xno;
				next_cell.yno = cell.yno +1;
				next_cell.zno = cell.zno;

				goto bestimmt;
			}
		}

		if(grady < 0)
		{
				tmp_xpos = cell.xpos - (cell.ypos *(gradx/grady));
				tmp_zpos = cell.zpos - (cell.ypos *(gradz/grady));

			if(1 >= tmp_xpos && tmp_xpos >= 0 && 1 >= tmp_zpos && tmp_zpos >= 0)
			{
				next_cell.xpos = tmp_xpos;
				next_cell.ypos = 1;
				next_cell.zpos = tmp_zpos;
				next_cell.dirx_i =  0;
				next_cell.diry_i =  1;
				next_cell.dirz_i =  0;
				next_cell.xno = cell.xno;
				next_cell.yno = cell.yno - 1;
				next_cell.zno = cell.zno;

				goto bestimmt;
			}

		}

	}


	/*Ray enter the cell from the positive or negative z-plane*/
	/*CASE I:The ray run into the cell*/
	if(cell.dirz_i == 0 || (cell.dirz_i == 1 && gradz <= 0) || (cell.dirz_i == -1 && gradz >= 0))
	{
		if(gradz > 0)
		{
				tmp_xpos = cell.xpos - ((cell.zpos - 1)*(gradx/gradz));
				tmp_ypos = cell.ypos - ((cell.zpos - 1)*(grady/gradz));

			if(1 >= tmp_xpos && tmp_xpos >= 0 && 1 >= tmp_ypos  && tmp_ypos >= 0)
			{
				next_cell.xpos = tmp_xpos;
				next_cell.ypos = tmp_ypos;
				next_cell.zpos = 0;
				next_cell.dirx_i =  0;
				next_cell.diry_i =  0;
				next_cell.dirz_i = -1;
				next_cell.xno = cell.xno;
				next_cell.yno = cell.yno;
				next_cell.zno = cell.zno +1;

				goto bestimmt;
			}
		}

		if(gradz < 0)
		{
				tmp_xpos = cell.xpos - (cell.zpos *(gradx/gradz));
				tmp_ypos = cell.ypos - (cell.zpos *(grady/gradz));

			if(1 >= tmp_xpos && tmp_xpos >= 0 && 1 >= tmp_ypos && tmp_ypos >= 0)
			{
				next_cell.xpos = tmp_xpos;
				next_cell.ypos = tmp_ypos;
				next_cell.zpos = 1;
				next_cell.dirx_i =  0;
				next_cell.diry_i =  0;
				next_cell.dirz_i =  1;
				next_cell.xno = cell.xno;
				next_cell.yno = cell.yno;
				next_cell.zno = cell.zno -1;

				goto bestimmt;
			}

		}
	}




	/*Ray enter the cell from the positive or negative x-plane*/
	/*CASE II: The ray run along the surface of the cell*/
	if(((cell.dirx_i == 1 && gradx > 0) || (cell.dirx_i == -1 && gradx < 0)))
	{
		if(grady < 0)
		{
			tmp_zpos = cell.zpos - cell.ypos*(gradz/grady);

			if(1 >= tmp_zpos && tmp_zpos >= 0)
			{

				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno + cell.dirx_i,cell.yno -1,cell.zno,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno,cell.yno -1,cell.zno,tt,ny,nz);
				if(gradient1[0] + gradient2[0] >= 0)
					diff = cell.dirx_i;
				else
					diff = - cell.dirx_i;

				if((gradient1[0] + gradient2[0]) >= 0)
					next_cell.xpos = 0;
				else
					next_cell.xpos = 1;
				next_cell.ypos = 1;
				next_cell.zpos = tmp_zpos;
				next_cell.dirx_i = 0;
				next_cell.diry_i = 1;
				next_cell.dirz_i = 0;
				if(diff>0)
					next_cell.xno = cell.xno + cell.dirx_i;
				else
					next_cell.xno = cell.xno;
				next_cell.yno = cell.yno -1;
				next_cell.zno = cell.zno;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

		if(grady > 0)
		{
			tmp_zpos = cell.zpos - (cell.ypos - 1)*(gradz/grady);

			if(1 >= tmp_zpos && tmp_zpos >= 0)
			{

				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno + cell.dirx_i,cell.yno +1,cell.zno,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno,cell.yno +1,cell.zno,tt,ny,nz);
				if(gradient1[0] + gradient2[0] >= 0)
					diff = cell.dirx_i;
				else
					diff = - cell.dirx_i;

				if(gradient1[0] + gradient2[0] >= 0)
					next_cell.xpos = 0;
				else
					next_cell.xpos = 1;
				next_cell.ypos = 0;
				next_cell.zpos = tmp_zpos;
				next_cell.dirx_i = 0;
				next_cell.diry_i = -1;
				next_cell.dirz_i = 0;
				if(diff>0)
					next_cell.xno = cell.xno + cell.dirx_i;
				else
					next_cell.xno = cell.xno;
				next_cell.yno = cell.yno +1;
				next_cell.zno = cell.zno;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

		if(gradz < 0)
		{
			tmp_ypos = cell.ypos - cell.zpos*(grady/gradz);

			if(1 >= tmp_ypos && tmp_ypos >= 0)
			{

				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno + cell.dirx_i,cell.yno ,cell.zno -1,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno,cell.yno, cell.zno -1,tt,ny,nz);
				if(gradient1[0] + gradient2[0] >= 0)
					diff = cell.dirx_i;
				else
					diff = -cell.dirx_i;

				if(gradient1[0] + gradient2[0] >= 0)
					next_cell.xpos = 0;
				else
					next_cell.xpos = 1;
				next_cell.ypos = tmp_ypos;
				next_cell.zpos = 1;
				next_cell.dirx_i = 0;
				next_cell.diry_i = 0;
				next_cell.dirz_i = 1;

				if(diff>0)
					next_cell.xno = cell.xno + cell.dirx_i;
				else
					next_cell.xno = cell.xno;
				next_cell.yno = cell.yno;
				next_cell.zno = cell.zno -1;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

		if(gradz > 0)
		{
			tmp_ypos = cell.ypos - (cell.zpos - 1)*(grady/gradz);

			if(1 >= tmp_ypos && tmp_ypos >= 0)
			{
				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno + cell.dirx_i,cell.yno ,cell.zno +1,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno,cell.yno, cell.zno +1,tt,ny,nz);
				if(gradient1[0] + gradient2[0] >= 0)
					diff = cell.dirx_i;
				else
					diff = -cell.dirx_i;

				if(gradient1[0] + gradient2[0] >= 0)
					next_cell.xpos = 0;
				else
					next_cell.xpos = 1;
				next_cell.ypos = tmp_ypos;
				next_cell.zpos = 0;
				next_cell.dirx_i = 0;
				next_cell.diry_i = 0;
				next_cell.dirz_i = -1;
				if(diff>0)
					next_cell.xno = cell.xno + cell.dirx_i;
				else
					next_cell.xno = cell.xno;
				next_cell.yno = cell.yno;
				next_cell.zno = cell.zno + 1;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

	}



	/*Ray enter the cell from the positive or negative y-plane*/
	/*CASE II: The ray run along the surface of the cell*/
	if(((cell.diry_i == 1 && grady > 0) || (cell.diry_i == -1 && grady < 0)))
	{
		if(gradx < 0)
		{
			tmp_zpos = cell.zpos - cell.xpos*(gradz/gradx);

			if(1 >= tmp_zpos && tmp_zpos >= 0)
			{
				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno -1,cell.yno + cell.diry_i,cell.zno,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno -1,cell.yno,cell.zno,tt,ny,nz);
				if(gradient1[1] + gradient2[1] >= 0)
					diff = cell.diry_i;
				else
					diff = -cell.diry_i;

				next_cell.xpos = 1;
				if(gradient1[1] + gradient2[1] >= 0)
					next_cell.ypos = 0;
				else
					next_cell.ypos = 1;
				next_cell.zpos = tmp_zpos;
				next_cell.dirx_i = 1;
				next_cell.diry_i = 0;
				next_cell.dirz_i = 0;
				next_cell.xno = cell.xno -1;
				if(diff>0)
					next_cell.yno = cell.yno + cell.diry_i;
				else
					next_cell.yno = cell.yno;
				next_cell.zno = cell.zno;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

		if(gradx > 0)
		{
			tmp_zpos = cell.zpos - (cell.xpos - 1)*(gradz/gradx);

			if(1 >= tmp_zpos && tmp_zpos >= 0)
			{
				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno +1,cell.yno + cell.diry_i,cell.zno,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno +1,cell.yno,cell.zno,tt,ny,nz);
				if(gradient1[1] + gradient2[1] >= 0)
					diff = cell.diry_i;
				else
					diff = -cell.diry_i;

				next_cell.xpos = 0;
				if(gradient1[1] + gradient2[1] >= 0)
					next_cell.ypos = 0;
				else
					next_cell.ypos = 1;
				next_cell.zpos = tmp_zpos;
				next_cell.dirx_i = -1;
				next_cell.diry_i = 0;
				next_cell.dirz_i = 0;
				next_cell.xno = cell.xno +1;
				if(diff>0)
					next_cell.yno = cell.yno + cell.diry_i;
				else
					next_cell.yno = cell.yno;
				next_cell.zno = cell.zno;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}
		if(gradz < 0)
		{
			tmp_xpos = cell.xpos - cell.zpos*(gradx/gradz);

			if(1 >= tmp_xpos && tmp_xpos >= 0)
			{
				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno,cell.yno + cell.diry_i,cell.zno -1,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno,cell.yno,cell.zno -1,tt,ny,nz);
				if(gradient1[1] + gradient2[1] >= 0)
					diff = cell.diry_i;
				else
					diff = -cell.diry_i;

				next_cell.xpos = tmp_xpos;
				if(gradient1[1] + gradient2[1] >= 0)
					next_cell.ypos = 0;
				else
					next_cell.ypos = 1;
				next_cell.zpos = 1;
				next_cell.dirx_i = 0;
				next_cell.diry_i = 0;
				next_cell.dirz_i = 1;

				next_cell.xno = cell.xno;
				if(diff>0)
					next_cell.yno = cell.yno + cell.diry_i;
				else
					next_cell.yno = cell.yno;
				next_cell.zno = cell.zno -1;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

		if(gradz > 0)
		{
			tmp_xpos = cell.xpos - (cell.zpos - 1)*(gradx/gradz);

			if(1 >= tmp_xpos && tmp_xpos >= 0)
			{


				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno,cell.yno + cell.diry_i,cell.zno +1,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno,cell.yno,cell.zno +1,tt,ny,nz);
				if(gradient1[1] + gradient2[1] >= 0)
					diff = cell.diry_i;
				else
					diff = -cell.diry_i;

				next_cell.xpos = tmp_xpos;
				if(gradient1[1] + gradient2[1] >= 0)
					next_cell.ypos = 0;
				else
					next_cell.ypos = 1;
				next_cell.zpos = 0;
				next_cell.dirx_i = 0;
				next_cell.diry_i = 0;
				next_cell.dirz_i = -1;
				next_cell.xno = cell.xno;
				if(diff>0)
					next_cell.yno = cell.yno + cell.diry_i;
				else
					next_cell.yno = cell.yno;
				next_cell.zno = cell.zno +1;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

	}




	/*Ray enter the cell from the positive or negative z-plane*/
	/*CASE II: The ray run along the surface of the cell*/
	if(((cell.dirz_i == 1 && gradz > 0) || (cell.dirz_i == -1 && gradz < 0)))
	{
		if(gradx < 0)
		{
			tmp_ypos = cell.ypos - cell.xpos*(grady/gradx);

			if(1 >= tmp_ypos && tmp_ypos >= 0)
			{
				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno -1,cell.yno,cell.zno +cell.dirz_i,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno -1,cell.yno,cell.zno,tt,ny,nz);
				if(gradient1[2] + gradient2[2] >= 0)
					diff = cell.dirz_i;
				else
					diff = -cell.dirz_i;

				next_cell.xpos = 1;
				next_cell.ypos = tmp_ypos;
				if(gradient1[2] + gradient2[2] >= 0)
					next_cell.zpos = 0;
				else
					next_cell.zpos = 1;
				next_cell.dirx_i = 1;
				next_cell.diry_i = 0;
				next_cell.dirz_i = 0;

				next_cell.xno = cell.xno -1;
				next_cell.yno = cell.yno;
				if(diff>0)
					next_cell.zno = cell.zno + cell.dirz_i;
				else
					next_cell.zno = cell.zno;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

		if(gradx > 0)
		{
			tmp_ypos = cell.ypos - (cell.xpos - 1)*(grady/gradx);

			if(1 >= tmp_ypos && tmp_ypos >= 0)
			{

				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno +1,cell.yno,cell.zno +cell.dirz_i,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno +1,cell.yno,cell.zno,tt,ny,nz);
				if(gradient1[2] + gradient2[2] >= 0)
					diff = cell.dirz_i;
				else
					diff = -cell.dirz_i;

				next_cell.xpos = 0;
				next_cell.ypos = tmp_ypos;
				if(gradient1[2] + gradient2[2] >= 0)
					next_cell.zpos = 0;
				else
					next_cell.zpos = 1;
				next_cell.dirx_i = -1;
				next_cell.diry_i = 0;
				next_cell.dirz_i = 0;
				next_cell.xno = cell.xno +1;
				next_cell.yno = cell.yno;
				if(diff>0)
					next_cell.zno = cell.zno + cell.dirz_i;
				else
					next_cell.zno = cell.zno;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}
		if(grady < 0)
		{
			tmp_xpos = cell.xpos - cell.ypos*(gradx/grady);

			if(1 >= tmp_xpos && tmp_xpos >= 0)
			{
				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno,cell.yno -1,cell.zno +cell.dirz_i,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno,cell.yno -1,cell.zno,tt,ny,nz);
				if(gradient1[2] + gradient2[2] >= 0)
					diff = cell.dirz_i;
				else
					diff = -cell.dirz_i;

				next_cell.xpos = tmp_xpos;
				next_cell.ypos = 1;
				if(gradient1[2] + gradient2[2] >= 0)
					next_cell.zpos = 0;
				else
					next_cell.zpos = 1;
				next_cell.dirx_i = 0;
				next_cell.diry_i = 1;
				next_cell.dirz_i = 0;
				next_cell.xno = cell.xno;
				next_cell.yno = cell.yno -1;
				if(diff>0)
					next_cell.zno = cell.zno + cell.dirz_i;
				else
					next_cell.zno = cell.zno;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

		if(grady > 0)
		{
			tmp_xpos = cell.xpos - (cell.ypos - 1)*(gradx/grady);

			if(1 >= tmp_xpos && tmp_xpos >= 0)
			{

				/*Comparison of the gradients(in the neighboring cells) to determine the next cell:*/
				/*Calculate the traveltime gradient*/
				gradient1 = TimeGrad(cell.xno,cell.yno +1,cell.zno +cell.dirz_i,tt,ny,nz);
				gradient2 = TimeGrad(cell.xno,cell.yno +1,cell.zno,tt,ny,nz);
				if(gradient1[2] + gradient2[2] >= 0)
				{
					diff = cell.dirz_i;
				}
				else
				{
					diff = -cell.dirz_i;
				}

				next_cell.xpos = tmp_xpos;
				next_cell.ypos = 0;
				if(gradient1[2] + gradient2[2] >= 0)
					next_cell.zpos = 0;
				else
					next_cell.zpos = 1;
				next_cell.dirx_i = 0;
				next_cell.diry_i = -1;
				next_cell.dirz_i = 0;

				next_cell.xno = cell.xno;
				next_cell.yno = cell.yno +1;
				if(diff>0)
					next_cell.zno = cell.zno + cell.dirz_i;
				else
					next_cell.zno = cell.zno;

				free(gradient1);
				free(gradient2);

				goto bestimmt;
			}
		}

	}




	/*Unexspected cases*/
	printf("Unexspected ray path in a grid cell\n");
	next_cell = cell;

	bestimmt:;

	return(next_cell);
}

/*-------------------------------------------------------------*/
/*Resort the ray segments; if a ray runs along the boundary of two cells this routine compares the velocity*/
/*in the two cells and assigns the ray segment to the cell with the higher velocity*/
/*Parameter:	*raypath			:= Raypath structures (Attention: x,y,z are normalized coordinates) */
/*				data				:= Data structure*/
/*				grid				:= Grid structure */

int ResortRays(RP_STRUCT *raypath,DATA_STRUCT data,GRID_STRUCT grid)
{
	long a,b,c,d,e;
	long nx,ny,nz,nyz,ny1,nz1,nyz1;

	double eps = 0.01;

	nx = grid.nx + 2*grid.nborder;
	ny = grid.ny + 2*grid.nborder;
	nz = grid.nz + 2*grid.nborder;
	nyz = ny*nz;

	ny1 = ny+1;
	nz1 = nz+1;
	nyz1 = ny1*nz1;


	/*Loop over all rays*/
	for(a=0;a<data.ndata_seis;a++)
		/*Loop over all segments of a ray*/
		for(b=0;b<raypath[a].nray;b++)
		{
					/*Find the right cell*/
					c = (long)floor((double)(raypath[a].ele[b]/nyz));
					d = (long)floor((double)((raypath[a].ele[b] - c*nyz)/ny));
					e = raypath[a].ele[b] - c*nyz - d*ny;
					if(raypath[a].ele[b] == (c*nyz + d*nz + e))
					{
						if(raypath[a].x[b]- (double)(c) <= eps && raypath[a].x[b+1] - (double)c <= eps)
						{
							if(grid.slow[c*nyz1 +d*nz1 +e] > grid.slow[(c-1)*nyz1 +d*nz1 +e])
							{
								raypath[a].ele[b] = ((c-1)*nyz + d*nz + e);
							}
						}

						if(raypath[a].x[b] - (double)c >= 1.0 - eps && raypath[a].x[b+1] - (double)c >= 1.0 - eps)
						{
							if (grid.slow[c*nyz1 +d*nz1 +e] > grid.slow[(c+1)*nyz1 +d*nz1 +e])
							{
								raypath[a].ele[b] = ((c+1)*nyz + d*nz + e);
							}
						}

						if(raypath[a].y[b]- (double)d <= eps && raypath[a].y[b+1] - (double)d <= eps)
						{
							if(grid.slow[c*nyz1 +d*nz1 +e] > grid.slow[c*nyz1 +(d-1)*nz1 +e])
							{
								raypath[a].ele[b] = (c*nyz + (d-1)*nz + e);
							}
						}

						if(raypath[a].y[b] - (double)d >= 1.0 - eps && raypath[a].y[b+1] - (double)d >= 1.0 - eps)
						{
							if (grid.slow[c*nyz1 +d*nz1 +e] > grid.slow[c*nyz1 +(d+1)*nz1 +e])
							{
								raypath[a].ele[b] = (c*nyz + (d+1)*nz + e);
							}
						}

						if(raypath[a].z[b]- (double)e <= eps && raypath[a].z[b+1] - (double)e <= eps)
						{
							if(grid.slow[c*nyz1 +d*nz1 +e] > grid.slow[c*nyz1 +d*nz1 +(e-1)])
							{
								raypath[a].ele[b] = (c*nyz + d*nz + (e-1));
							}
						}
						if(raypath[a].z[b] - (double)e >= 1.0 - eps && raypath[a].z[b+1] - (double)e >= 1.0 - eps)
						{
							if (grid.slow[c*nyz1 +d*nz1 +e] > grid.slow[c*nyz1 +d*nz1 +(e+1)])
							{
								raypath[a].ele[b] = (c*nyz + d*nz + (e+1));
							}
						}

						goto slutt;
					}
			slutt:;
		}


	return(1);
}
}
