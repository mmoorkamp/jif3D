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
/*Down-scaling of parallel rays following the idea of Hansruedi*/
/*		VERSUCH ZWEI		*/
/*Parameter: *raypath := Raypath structure*/
/*			 inv_cell  := Inversion cell structure*/
/*				ncell := Number of inversion cells*/
/*				nray  := Number of shot-receiver combinations*/
/*		angle_threshold:= maximun angle between rays for which the rays are considered to be "parallel"*/
/*     weighting_factor:= weighting factor of the rays (will be modified in this routine)*/
/*		kind_of_weight := Kind of weighting (absolut weighting = 1; relative weighting = 2)*/

/*Remark: The weightings are finally balanced such that the old and new weights have the same average value*/

int	ScalParallelRays( RP_STRUCT *raypath, BIDX_STRUCT *inv_cell, long ncell, long nray, double angle_threshold, double *weighting_factor, int kind_of_weight)
{
	long i,j,k;
	long *nr_of_rays_per_cell;		/*Number of rays per inversion cell*/
	long **ray_numbers_in_cell;		/*The number of the rays located in each cell*/
	long **ray_seg_numbers_in_cell; /*The number of the corresponding ray segments*/
	long *nr_parallel_rays;			/*Number of parallel rays for each ray*/
	long min_nr_parallel_rays;		/*Min.number of parallel rays for each ray*/
	double ***ray_vector;			/*vector of the ray segments in the cell*/
	double **ray_weight;			/*Weights of the rays in each cell*/
	double *summed_ray_weight;		/*Mean weights along the whole ray*/
	double summed_weighting_factor_old, summed_weighting_factor_new;
	double scalar_product;
	double mean_a, mean_b;
	double angle;					/*Angle between the vectors*/
	double length;					/*total length of the ray*/



	printf("----------------\n");
	printf("Start with re-adjusting the weight of parallel rays\n");
	printf("----------------\n\n");

	/*Allocate memory for the number of rays per cell*/
	nr_of_rays_per_cell = (long *)memory(NULL,ncell,sizeof(long),"ScalParallelRays");

	for(i=0;i<ncell;i++)
		nr_of_rays_per_cell[i] = 0;

	/*Allocate memory for location numbers of all rays within each cell*/
	ray_numbers_in_cell = (long **)memory(NULL,ncell,sizeof(long *),"ScalParallelRays");

	for(i=0;i<ncell;i++)
		ray_numbers_in_cell[i] = (long *)memory(NULL,1,sizeof(long),"ScalParallelRays");


	/*Allocate memory for location numbers of the corresponding ray segments*/
	ray_seg_numbers_in_cell = (long **)memory(NULL,ncell,sizeof(long *),"ScalParallelRays");

	for(i=0;i<ncell;i++)
		ray_seg_numbers_in_cell[i] = (long *)memory(NULL,1,sizeof(long),"ScalParallelRays");

	/*Allocate memory for the vector of the corresponding ray segments*/
	ray_vector = (double ***)memory(NULL,ncell,sizeof(double **),"ScalParallelRays");

	for(i=0;i<ncell;i++)
		ray_vector[i] = (double **)memory(NULL,3,sizeof(double *),"ScalParallelRays");

	for(i=0;i<ncell;i++)
		for(j=0;j<3;j++)
			ray_vector[i][j] = (double *)memory(NULL,1,sizeof(double),"ScalParallelRays");

	/*Allocate memory for the weight of the ray segments within each cell*/
	ray_weight = (double **)memory(NULL,ncell,sizeof(double *),"ScalParallelRays");

	for(i=0;i<ncell;i++)
		ray_weight[i] = (double *)memory(NULL,1,sizeof(double),"ScalParallelRays");
	

	/*Loop over all rays*/
	for(i=0;i<nray;i++)
	{
		/*Loop over all ray segments*/
		for(j=0;j<raypath[i].nray;j++)
		{

			if(raypath[i].ele[j] >= ncell)
			{
				printf("NOT enough memory is allocated: used %d, allocated %d\n",raypath[i].ele[j]+1, ncell);
				exit(0);
			}
			
			ray_numbers_in_cell[raypath[i].ele[j]] = (long *)memory((char *)ray_numbers_in_cell[raypath[i].ele[j]],nr_of_rays_per_cell[raypath[i].ele[j]]+1, sizeof(long),"ScalParallelRays");
			ray_seg_numbers_in_cell[raypath[i].ele[j]] = (long *)memory((char *) ray_seg_numbers_in_cell[raypath[i].ele[j]],nr_of_rays_per_cell[raypath[i].ele[j]]+1, sizeof(long),"ScalParallelRays");
			ray_vector[raypath[i].ele[j]][0] = (double *)memory((char *) ray_vector[raypath[i].ele[j]][0],nr_of_rays_per_cell[raypath[i].ele[j]]+1, sizeof(double),"ScalParallelRays");
			ray_vector[raypath[i].ele[j]][1] = (double *)memory((char *) ray_vector[raypath[i].ele[j]][1],nr_of_rays_per_cell[raypath[i].ele[j]]+1, sizeof(double),"ScalParallelRays");
			ray_vector[raypath[i].ele[j]][2] = (double *)memory((char *) ray_vector[raypath[i].ele[j]][2],nr_of_rays_per_cell[raypath[i].ele[j]]+1, sizeof(double),"ScalParallelRays");
			ray_weight[raypath[i].ele[j]] = (double *)memory((char *) ray_weight[raypath[i].ele[j]], nr_of_rays_per_cell[raypath[i].ele[j]]+1, sizeof(double),"ScalParallelRays");
			
			/*Make a list of rays crossing the corresponding cells*/
			ray_numbers_in_cell[raypath[i].ele[j]][nr_of_rays_per_cell[raypath[i].ele[j]]] = i;
			/* ... and make a list of the corresponding ray segments*/
			ray_seg_numbers_in_cell[raypath[i].ele[j]][nr_of_rays_per_cell[raypath[i].ele[j]]] = j;

			/*Determine the vectors of the ray segments*/
			ray_vector[raypath[i].ele[j]][0][nr_of_rays_per_cell[raypath[i].ele[j]]] = raypath[i].x[j];
			ray_vector[raypath[i].ele[j]][1][nr_of_rays_per_cell[raypath[i].ele[j]]] = raypath[i].y[j];
			ray_vector[raypath[i].ele[j]][2][nr_of_rays_per_cell[raypath[i].ele[j]]] = raypath[i].z[j];

			/*Set the weights to 1*/
			ray_weight[raypath[i].ele[j]][nr_of_rays_per_cell[raypath[i].ele[j]]] = 1.0;

			/*Determine the number of rays in the inversion cells*/
			nr_of_rays_per_cell[raypath[i].ele[j]]++;

		}
	}


	/***********************************************************/
	/*Calculating the angular dependent weighting factors for all rays in the different cells*/

	for(i=0;i<ncell;i++)
	{
		/*Only cells which are intersected by at least two (absolute weighting) or three (relative weighting) rays have to be considered*/
		if((nr_of_rays_per_cell[i] > 1 && kind_of_weight == 1) || (nr_of_rays_per_cell[i] > 2 && kind_of_weight == 2) )
		{

			/*Allocate memory for storing the number of parallel rays for each ray*/
			nr_parallel_rays = (long *)memory(NULL,nr_of_rays_per_cell[i],sizeof(long),"ScalParallelRays");

			for(j=0;j<nr_of_rays_per_cell[i];j++)
				nr_parallel_rays[j] = 1;

			/*Loop over all rays in the cell*/
			for(j=0;j<nr_of_rays_per_cell[i];j++)
			{
				/*Mean_value of the first vector*/
				mean_a = sqrt(ray_vector[i][0][j]*ray_vector[i][0][j] + ray_vector[i][1][j]*ray_vector[i][1][j] + ray_vector[i][2][j]*ray_vector[i][2][j]);

				/*Loop over all rays in the cell once again*/
				for(k=0;k<nr_of_rays_per_cell[i];k++)
				{

					/*Calculate the angles between the vectors:*/

					/*Scalar-product*/
					scalar_product = (ray_vector[i][0][j]*ray_vector[i][0][k])+(ray_vector[i][1][j]*ray_vector[i][1][k])+(ray_vector[i][2][j]*ray_vector[i][2][k]);

					/*Mean_value of the second vector*/
					mean_b = sqrt(ray_vector[i][0][k]*ray_vector[i][0][k] + ray_vector[i][1][k]*ray_vector[i][1][k] + ray_vector[i][2][k]*ray_vector[i][2][k]);

					if(mean_a != 0.0 && mean_b != 0.0 && j != k)
					{
//MODIFIED START!!
						/*Angle between the vectors:*/
						if((scalar_product/(mean_a * mean_b)) < 1.0 && (scalar_product/(mean_a * mean_b)) > -1.0)
							angle = (180/PI)*acos(scalar_product/(mean_a * mean_b));
						else
							angle = 0.0;
//MODIFIED END!!

						/*identify all "parallel" rays*/
						if(fabs(angle) <= angle_threshold)
							nr_parallel_rays[j]++;
					}

				}
			}

			/********************************/
			/*Absolute weighting*/
			if(kind_of_weight == 1)
			{
				/*Determine the weighting factor for the ray segments within each cell*/
				for(j=0;j<nr_of_rays_per_cell[i];j++)
					ray_weight[i][j] = ((1 + (double)nr_of_rays_per_cell[i] - (double)nr_parallel_rays[j])/(double)nr_of_rays_per_cell[i]);
			}

			/********************************/
			/*Relative weighting*/
			else
			{
				/*Find the factor for the re-weighting of the rays*/
				min_nr_parallel_rays = nr_of_rays_per_cell[i];

				for(j=0;j<nr_of_rays_per_cell[i];j++)
				{
					if(nr_parallel_rays[j] < min_nr_parallel_rays)
						min_nr_parallel_rays = nr_parallel_rays[j];
				}

				/*Apply the factor to the rays*/
				for(j=0;j<nr_of_rays_per_cell[i];j++)
					ray_weight[i][j] = min_nr_parallel_rays/nr_parallel_rays[j];
			}
			/********************************/


			free(nr_parallel_rays);

		}
	}


	/***********************************************************/

	/*Allocate memory for mean weights along the rays*/
	summed_ray_weight = (double *)memory(NULL,nray,sizeof(double),"ScalParallelRays");
	for(i=0;i<nray;i++)
		summed_ray_weight[i] = 0.0;


	/*Determine the mean weighting factors for each ray*/
	for(i=0;i<ncell;i++)
	{
		/*Only inversion cells will be considered in the calculation*/
		if(inv_cell[i].use == 1)
		{

			for(j=0;j<nr_of_rays_per_cell[i];j++)
			{
				/*Mean_value of the vector*/
				mean_a = sqrt(ray_vector[i][0][j]*ray_vector[i][0][j] + ray_vector[i][1][j]*ray_vector[i][1][j] + ray_vector[i][2][j]*ray_vector[i][2][j]);

				/*Sum all weights*/
				summed_ray_weight[ray_numbers_in_cell[i][j]] = summed_ray_weight[ray_numbers_in_cell[i][j]] + (mean_a * ray_weight[i][j]);
			}
		}
	}


	/*Summing of all "old" weighting factors*/
	summed_weighting_factor_old = 0.0;
	for(i=0;i<nray;i++)
		summed_weighting_factor_old = summed_weighting_factor_old + weighting_factor[i];



	/*Divide the weights by the ray lengths and modify the weighting factor:*/
	/*Loop over all rays*/
	for(i=0;i<nray;i++)
	{
		length = 0.0;

		/*Loop over all ray segments*/
		for(j=0;j<raypath[i].nray;j++)
		{
			/*Only inversion cells will be considered in the calculation*/
			if(inv_cell[raypath[i].ele[j]].use == 1)
				length = length + raypath[i].len[j];
		}

		/*Modify the weighting factors*/
		if(length != 0.0)
			weighting_factor[i] = (summed_ray_weight[i]/length)* weighting_factor[i];

	}


	/*Summing of all "new" weighting factors*/
	summed_weighting_factor_new = 0.0;
	for(i=0;i<nray;i++)
		summed_weighting_factor_new = summed_weighting_factor_new + weighting_factor[i];

	/*Adjust the individual weighting_factors*/
	for(i=0;i<nray;i++)
	{
		if(summed_weighting_factor_new != 0.0)
			weighting_factor[i] = (summed_weighting_factor_old/summed_weighting_factor_new)* weighting_factor[i];
	}
	

	printf("----------------\n");
	printf("Parallel rays are re-adjusted\n");
	printf("----------------\n\n");


	/*Free the memory*/
	for(i=0;i<ncell;i++)
	{
		for(j=0;j<3;j++)
		{
			free(ray_vector[i][j]);
		}
	}
	
	for(i=0;i<ncell;i++)
	{
		free(ray_numbers_in_cell[i]);
		free(ray_seg_numbers_in_cell[i]);
		free(ray_vector[i]);
		free(ray_weight[i]);
	}

	free(ray_vector);
	free(ray_numbers_in_cell);
	free(ray_seg_numbers_in_cell);
	free(nr_of_rays_per_cell);
	free(ray_weight);
	free(summed_ray_weight);

	return(1);
}

/*------------------------------------------------------------*/
/*Balance the entries of the different methods in the Jacobian matrix A*/
/*to avoid numerical problems during the inversion process.*/
/*Thereby, the scaling factors "inv->rel_scal_seis","inv->rel_scal_grav" and/or "inv->rel_scal_mt" are reset*/
/*Parameter:  flag := Structure with flags*/
/*			  data := Data structure*/
/*			    *A := sparse matrix*/
/*			   vec := Structure of the residual vector*/

int BalanceSparseMatrix(INV *inv, FLAG_STRUCT flag, SPMAT *A, double *res)
{
	long i,j,c;
	double mean_seis, mean_grav, mean_mt, val;

	/*Balancing is only required for JOINT inversions*/
	mean_seis = 0.0;
	mean_grav = 0.0;
	mean_mt = 0.0;

	printf("\n----------------\n");
	printf("Start balancing the different data sets\n");
	printf("----------------\n");


	/*Calculate the average size of the seismic derivatives per inversion cell*/
	if(flag.index_tseis != 0)
	{
		mean_seis = RegFactor(A,inv->seis_first_row,(inv->seis_first_row + inv->seis_nr_row),0,inv->nvel_used, 3, &c);
		printf("scaling factor (seismic): %lf\n\n", mean_seis);
	}
	/*Calculate the average size of the gravity derivatives per inversion cell*/
	if(flag.index_grav != 0)
	{
		mean_grav = RegFactor(A,inv->grav_first_row,(inv->grav_first_row + inv->grav_nr_row), 0,inv->nvel_used, 3, &c);
		printf("scaling factor (gravity): %lf\n\n", mean_grav);
	}
	/*Calculate the average size of the MT derivatives per inversion cell*/
	if(flag.index_mt != 0)
	{
		mean_mt = RegFactor(A,inv->mt_first_row,(inv->mt_first_row + inv->mt_nr_row), 0,inv->nvel_used, 3, &c);
		printf("scaling factor (resistivity): %lf\n\n", mean_mt);
	}

	/*Perform the balancing*/

	if(mean_seis >= mean_grav && mean_seis >= mean_mt)
	{
		if(mean_grav != 0)
		{
			inv->rel_scal_grav = mean_seis/mean_grav;

			for(i=inv->grav_first_row;i<(inv->grav_first_row + inv->grav_nr_row);i++)
			{
				for(j=0;j<inv->ncol;j++)
				{
					/*Balance the sparse matrix*/
					val= inv->rel_scal_grav * sp_get_val(A,i,j);
					if(val != 0)
						sp_set_val(A,i,j,val);
				}

				/*... and the residual data vector*/
				res[i] = inv->rel_scal_grav * res[i]; 
			}
		}

		if(mean_mt != 0)
		{
			inv->rel_scal_mt = mean_seis/mean_mt;

			for(i=inv->mt_first_row;i<(inv->mt_first_row + inv->mt_nr_row);i++)
			{
				for(j=0;j<inv->ncol;j++)
				{
					/*Balance the sparse matrix*/
					val= inv->rel_scal_mt * sp_get_val(A,i,j);
					if(val != 0)
						sp_set_val(A,i,j,val);
				}
				/*... and the residual data vector*/
				res[i] = inv->rel_scal_mt * res[i]; 
			}
		}

	}
	else if(mean_grav >= mean_seis && mean_grav >= mean_mt)
	{
		if(mean_seis != 0)
		{
			inv->rel_scal_seis = mean_grav/mean_seis;

			for(i=inv->seis_first_row;i<(inv->seis_first_row + inv->seis_nr_row);i++)
			{
				for(j=0;j<inv->ncol;j++)
				{
					/*Balance the sparse matrix*/
					val= inv->rel_scal_seis * sp_get_val(A,i,j);
					if(val != 0)
						sp_set_val(A,i,j,val);
				}
				/*... and the residual data vector*/
				res[i] = inv->rel_scal_seis * res[i]; 
			}
		}

		if(mean_mt != 0)
		{
			inv->rel_scal_mt = mean_grav/mean_mt;

			for(i=inv->mt_first_row;i<(inv->mt_first_row + inv->mt_nr_row);i++)
			{
				for(j=0;j<inv->ncol;j++)
				{
					/*Balance the sparse matrix*/
					val= inv->rel_scal_mt * sp_get_val(A,i,j);
					if(val != 0)
						sp_set_val(A,i,j,val);
				}
				/*... and the residual data vector*/
				res[i] = inv->rel_scal_mt * res[i]; 
			}

		}
		
	}
	else
	{
		if(mean_seis != 0)
		{
			inv->rel_scal_seis = mean_mt/mean_seis;

			for(i=inv->seis_first_row;i<(inv->seis_first_row + inv->seis_nr_row);i++)
			{
				for(j=0;j<inv->ncol;j++)
				{
					/*Balance the sparse matrix*/
					val = inv->rel_scal_seis * sp_get_val(A,i,j);
					if(val != 0)
						sp_set_val(A,i,j,val);
				}
				/*... and the residual data vector*/
				res[i] = inv->rel_scal_seis * res[i]; 
			}
		}

		if(mean_grav != 0)
		{
			inv->rel_scal_grav = mean_mt/mean_grav;

			for(i=inv->grav_first_row;i<(inv->grav_first_row + inv->grav_nr_row);i++)
			{
				for(j=0;j<inv->ncol;j++)
				{
					/*Balance the sparse matrix*/
					val= inv->rel_scal_grav * sp_get_val(A,i,j);
					if(val != 0)
						sp_set_val(A,i,j,val);
				}
				/*... and the residual data vector*/
				res[i] = inv->rel_scal_grav * res[i]; 
			}
		}
	}

	return(1);
}

/*------------------------------------------------------------*/
/*Balance the MT entries in the Jacobian matrix A for the different frequencies*/
/*(for each frequency the entries are equally weighted)*/
/*Parameter:   data := Data structure*/
/*			    *A := sparse matrix*/
/*			   vec := Structure of the residual vector*/

int BalanceMTfrequencies(INV *inv, SPMAT *A, double *res)
{
	double *index;
	long i,j,k,nfreq;
	long c;
	double *mean_mt_freq,val;

	printf("\n----------------\n");
	printf("Start balancing the MT data for the different frequencies\n");
	printf("----------------\n");

	nfreq = (inv->mt_nr_row/2);


	if(nfreq != 0)
	{
		mean_mt_freq = (double *)memory(NULL,nfreq+1,sizeof(double),"BalanceMTfrequencies");
		index = (double *)memory(NULL,nfreq+1,sizeof(double),"BalanceMTfrequencies");
	}
	else
	{
		mean_mt_freq = (double *)memory(NULL,1,sizeof(double),"BalanceMTfrequencies");
		index = (double *)memory(NULL,nfreq+1,sizeof(double),"BalanceMTfrequencies");
	}

	for(i=0;i<(nfreq+1);i++)
	{
		mean_mt_freq[i] = 0.0;
		index[i] = i;
	}

	if(nfreq != 0)
	{
		for(i=0;i<nfreq;i++)
		{
			/*Determing the size of the summed sensitivities*/
			mean_mt_freq[i+1] = RegFactor(A,(inv->mt_first_row + (2*i)), (inv->mt_first_row + (2*i) + 2), 0,inv->nvel_used, 3, &c);
		}

		/*sort them with size*/
		sort2(nfreq,mean_mt_freq,index);

		/*Determing the weighting factors*/
		for(i=1;i<nfreq;i++)
		{
			if(mean_mt_freq[i] != 0)
				mean_mt_freq[i] = mean_mt_freq[nfreq]/mean_mt_freq[i];
		}

		mean_mt_freq[nfreq] = 1;

		/*Assign the weighting factors to the matrix*/
		k = 0;

		for(i=inv->mt_first_row;i<(inv->mt_first_row + inv->mt_nr_row);i++)
		{
			for(j=0;j<inv->ncol;j++)
				{
					/*Balance the sparse matrix*/
					val= mean_mt_freq[(long)index[k+1]] * sp_get_val(A,i,j);
					if(val != 0)
						sp_set_val(A,i,j,val);
				}
				/*... and the residual data vector*/
				res[i] = mean_mt_freq[(long)index[k+1]] * res[i];

				if((i - inv->mt_first_row)%2 == 1)
					k=k+1;
		}


	}

	free(mean_mt_freq);
	free(index);

	return(1);
}


/*------------------------------------------------------------*/
/*Multiply the external weighting to the sparse matrix entries*/
/*Parameter:    *fatray			:= Fat-ray structure*/
/*				*raypath		:= Ray structure*/
/*				*grav			:= Gravity structure*/
/*				*mt				:= MT structure*/
/*				*inv_cell		:= Structure of the inversion cells*/
/*				 nray			:= Number of shot-receiver combinations*/
/*				 *A				:= Sparse matrix*/
/*				flag			:= General flags*/

int ReweightSparseMatrix(FLAG_STRUCT flag, F_RP_STRUCT *fatray,RP_STRUCT *raypath, GRAV_STRUCT *grav, MT_STRUCT *mt, BIDX_STRUCT *inv_cell, DATA_STRUCT data, SPMAT *A)
{
	long i,k,j,m;
	double val;

	/*Seismic*/
	if(flag.index_tseis != 0)
	{
		if(flag.kind_of_rays == 1)
		{
			/*for conventional rays*/
			/*Loop over all rays*/
			for(i=0;i<data.ndata_seis;i++)
				/*Loop over all ray segments*/
				for(j=0;j<raypath[i].nray;j++)
				{
					/*Add weights for the sparse matrix entries*/
					val= sp_get_val(A,raypath[i].n,inv_cell[raypath[i].ele[j]].used_nr);
					sp_set_val(A,raypath[i].n,inv_cell[raypath[i].ele[j]].used_nr,val*data.weigth_seis[i]);
				}
		}
		else
		{
			/*for fat-rays*/
			/*Loop over all fatrays*/
			for(i=0;i<data.ndata_seis;i++)
				/*Loop over all fat-ray segments*/
				for(j=0;j<fatray[i].ncell;j++)
				{
					/*Add weights for the sparse matrix entries*/
					val= sp_get_val(A,fatray[i].n,inv_cell[fatray[i].ele[j]].used_nr); 
					sp_set_val(A,fatray[i].n,inv_cell[fatray[i].ele[j]].used_nr,val*data.weigth_seis[i]);
				}	
		}
	}


	/*Gravity*/
	if(flag.index_grav != 0)
	{
		/*Loop over all stations*/
			for(i=0;i<data.ndata_grav;i++)
			{
				/*Loop over all inversion cells*/
				for(j=0;j<grav[i].ncell;j++)
				{
					/*Add weights for the sparse matrix entries*/
					val= sp_get_val(A,grav[i].n,inv_cell[grav[i].ele[j]].used_nr); 
					sp_set_val(A,grav[i].n,inv_cell[grav[i].ele[j]].used_nr,val*data.weigth_grav[i]);
				}
			}
	}


	/*MT*/
	if(flag.index_mt != 0)
	{
		/*Loop over all stations*/
		for(i=0;i<data.ndata_mt;i++)
		{
			/*Loop over all cells that affects the MT measurements*/
			for(j=0;j<mt[i].ncell;j++)
			{
				/*Only one mode is used*/
				if(flag.nr_of_modes == 1)
				{
					/*Loop over all frequencies*/
					for(k=0;k<(2*mt[i].nfreq);k++)
					{
						 m = k/2;

						/*Add weights for the sparse matrix entries*/
						val= sp_get_val(A,mt[i].n[k],inv_cell[mt[i].ele[j]].used_nr); 
						sp_set_val(A,mt[i].n[k],inv_cell[mt[i].ele[j]].used_nr,val*data.weigth_mt[i][m]);
					}
				}
				/*Both modes are used*/
				else
				{
					/*Loop over all frequencies*/
					for(k=0;k<(4*mt[i].nfreq);k++)
					{
						 m = k/4;

						/*Add weights for the sparse matrix entries*/
						val= sp_get_val(A,mt[i].n[k],inv_cell[mt[i].ele[j]].used_nr); 
						sp_set_val(A,mt[i].n[k],inv_cell[mt[i].ele[j]].used_nr,val*data.weigth_mt[i][m]);
					}
				}
			}
		}
	}


	return(1);
}

/*------------------------------------------------------------*/
/*Multiply the external weighting to the residual data vector*/
/*Parameter:    *fatray			:= Fat-ray structure*/
/*				*raypath		:= Ray structure*/
/*				*grav			:= Gravity structure*/
/*				*mt				:= MT structure*/
/*				data			:= Data structure*/
/*				flag			:= General flags*/
/*				vec				:= Structure for the residual vector*/

int ReweightResdDataVector(FLAG_STRUCT flag,VEC *vec, DATA_STRUCT data, RP_STRUCT *raypath, F_RP_STRUCT *fatray, GRAV_STRUCT *grav, MT_STRUCT *mt)
{
	long i,j;

	/*Seismic: traveltimes in s*/
	if(flag.index_tseis != 0)
	{
		/*Conventional rays*/
		if(flag.kind_of_rays == 1)
		{
			for(i=0;i<data.ndata_seis;i++) 
				vec->ve[raypath[i].n] = vec->ve[raypath[i].n]*data.weigth_seis[i];
		}
		else
		{
			for(i=0;i<data.ndata_seis;i++) 
				vec->ve[fatray[i].n] = vec->ve[fatray[i].n]*data.weigth_seis[i];
		}
	}

	/*Gravity: acceleration in mgal*/
	if(flag.index_grav != 0)
	{
		for(i=0;i<data.ndata_grav;i++)
			vec->ve[grav[i].n] = vec->ve[grav[i].n]*data.weigth_grav[i];
	}

	/*MT:Real- and imaginary part of the impedance in ohm*/
	if(flag.index_mt != 0)
	{
		for(i=0;i<data.ndata_mt;i++)
		{
			/*Only one mode is used*/
			if(flag.nr_of_modes == 1)
			{
				/*Loop over the frequencies*/
				for(j=0;j<(data.nfreq_mt[i]);j++)
				{
					/*Real part*/
					vec->ve[mt[i].n[2*j]] = vec->ve[mt[i].n[2*j]]*data.weigth_mt[i][j];
					/*Imaginary part*/
					vec->ve[mt[i].n[2*j+1]] = vec->ve[mt[i].n[2*j+1]]*data.weigth_mt[i][j];
				}
			}
			/*Both modes are used*/
			else
			{
				/*Loop over the frequencies*/
				for(j=0;j<(data.nfreq_mt[i]);j++)
				{
					/*TE-mode*/
					/*Real part*/
					vec->ve[mt[i].n[4*j]] = vec->ve[mt[i].n[4*j]]*data.weigth_mt[i][j];
					/*Imaginary part*/
					vec->ve[mt[i].n[4*j+1]] = vec->ve[mt[i].n[4*j+1]]*data.weigth_mt[i][j];
					/*TM-mode*/
					/*Real part*/
					vec->ve[mt[i].n[4*j+2]] = vec->ve[mt[i].n[4*j+2]]*data.weigth_mt[i][j];
					/*Imaginary part*/
					vec->ve[mt[i].n[4*j+3]] = vec->ve[mt[i].n[4*j+3]]*data.weigth_mt[i][j];					
				}
			}
		}
	}

	return(1);
}