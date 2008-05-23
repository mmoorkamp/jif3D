#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "dcomplex.h"
#include "inv3d.h"
#include "fkt.h"

/*-------------------------------------------------------------*/
/*Performing the 1D MT modeling:*/
/*Parameter:	geo  := Geometry structure  */
/*              grid  := Grid structure */
/*              *data  := Pointer on data structure (the MT impedances will be (re-)determined in this routine)*/
/*       kind_of_data  := Kind of data 1= TE-mode data; 2= TM-mode data; 3= Berdichewsky average*/

/*Remark: the magnetic permeability is considered to correspond to the one of the free air*/

#define res(x,y,z) grid.res[((grid.ny)*(grid.nz))*(x) + (grid.nz)*(y) + (z)] 
#define M_KONST 4*PI*1E-7 /*magnetic permeability of the air( considered for the earth) in [H/m]*/

int ForwardModMT_1D(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, int kind_of_data)
{
	int nr_xcells, nr_ycells;
	long i,j,k,m,n,r,r_index, count, nr_layer, nr_layer1;
	long int_tmp, nx_int[2], ny_int[2], nz_int;
	float xs,ys,zs;
	double imped_r, imped_i; /*real- and imaginary part of the impedance*/
	double nx_float, ny_float, nz_float, dz;
	double *thickness; /*Thicknesses of the "layers"*/
	double *resist; /*Resistivities in the considered column*/
	double ang_freq; /*angular frequency*/
 	double app_res; /*apparent resistivity*/
	double berd_real, berd_imag; /*real and imaginary part of the Berdichewsky average of the observed data*/

	FILE *out;
	out = fopen("MT_response_calc.txt","w");

	fprintf(out,"station frequency[Hz] Re(Z) Im(Z)  app.res[Ohmm]\n");
	
	printf("1-D\n\n");

	for(i=0;i<data->ndata_mt;i++)
	{
		for(j=0;j<data->nfreq_mt[i];j++)
		{
			data->calc_real_mt_TE[i][j] = 0.0; 
			data->calc_imag_mt_TE[i][j] = 0.0;
		}
	}

	if(data->ndata_mt == 0)
		printf("!!WARNING!! NO MT data exists but the resistivity forward modelling is activated\n\n");


	count = 0;
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

			/********************/
			/*Determine the thicknesses of the layers*/
			thickness = (double *)memory(NULL,1, sizeof(double),"ForwardModMT_1D");
			thickness[0] = dz;
			nr_layer1 = 1;

			for(r=1;r<grid.nz - nz_int -1;r++)
			{
				thickness = (double *)memory((char *)thickness, r+1, sizeof(double),"ForwardModMT_1D");
				thickness[r] = grid.h;

				nr_layer1++;
			}

			/********************************************/
			/*Calculate the MT response*/

			/*Relate the stations to the corresponding MT measurements*/
			for(m=0;m<data->ndata_mt;m++)
			{
				if(data->mts[i] == data->mno[m])
				{
					count++;
//					fprintf(out,"%d %d %d %d %d\n",count,data->mno[m],-9999,3,data->nfreq_mt[m]);


					/*Loop over all frequencies*/
					for(n=0;n<data->nfreq_mt[m];n++)
					{

						/*Use angular frequencies 2*PI*f*/
						ang_freq = 2*PI*(data->freq_mt[m][n]);

						r_index = 0;

						imped_r = 0.0;
						imped_i = 0.0;

						/*Loop over all used cells in xy direction*/
						for(j=0;j<nr_xcells;j++)
							for(k=0;k<nr_ycells;k++)
							{


								/*Determine the resistivities in the considered column*/
								resist = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_1D");
								nr_layer = 0;

								for(r=0;r<grid.nz - nz_int;r++)
								{
									resist = (double *)memory((char *)resist, r+1, sizeof(double),"ForwardModMT_1D");
									resist[r] = res(nx_int[j],ny_int[k],r + nz_int);
									nr_layer++;
								}

								if(nr_layer != (nr_layer1 + 1))
								{
									printf("The number of determined thicknesses do not corresponds to the number\nof the layer in the 1-D MT modelling\n");
									exit(0);
								}

								/*Calculate the 1-D MT response (TM-mode)*/
								MT_1D_CALC(ang_freq, &imped_r, &imped_i, thickness, nr_layer, resist, 0.0);

								/*If the station are located exactly at the border between to cells, the results are summed*/
								data->calc_real_mt_TE[m][n] = data->calc_real_mt_TE[m][n] + imped_r;
								data->calc_imag_mt_TE[m][n] = data->calc_imag_mt_TE[m][n] + imped_i;

								r_index++;

								free(resist);
							}

							data->calc_real_mt_TE[m][n] = data->calc_real_mt_TE[m][n]/r_index;
							data->calc_imag_mt_TE[m][n] = data->calc_imag_mt_TE[m][n]/r_index;

							/*Determine the apparent resistivity*/
							app_res = ((data->calc_real_mt_TE[m][n])*(data->calc_real_mt_TE[m][n]) + (data->calc_imag_mt_TE[m][n])*(data->calc_imag_mt_TE[m][n]))*(1/(ang_freq*(M_KONST)));

							/*Write out the MT response */
							  fprintf(out,"%d %f %f %f %f\n",m,data->freq_mt[m][n],data->calc_real_mt_TE[m][n],data->calc_imag_mt_TE[m][n], app_res);
//							  fprintf(out,"%f %f %f\n",data->freq_mt[m][n],data->calc_real_mt_TE[m][n],data->calc_imag_mt_TE[m][n]);

					}
					
				}
			}

		if(i%10 == 0)
			printf("Modeling for %d of %d stations is completed\n", i+1, geo.nstat_mt);

		free(thickness);
		}

		/***************************************************************/
		/*Calculate the RMS values of the real and imaginary part of Z*/
		data->rms_mt_real = 0.0;
		data->rms_mt_imag = 0.0;

		k=0;

		/*Loop over all data*/
		for(m=0;m<data->ndata_mt;m++)
		{
			/*TE-mode data*/
			if(kind_of_data == 1)
			{
				/*Loop over all frequencies*/
				for(n=0;n<data->nfreq_mt[m];n++)
				{
					data->rms_mt_real = data->rms_mt_real + (((data->calc_real_mt_TE[m][n] - data->real_mt_TE[m][n]) * (data->calc_real_mt_TE[m][n] - data->real_mt_TE[m][n])));
					data->rms_mt_imag = data->rms_mt_imag + (((data->calc_imag_mt_TE[m][n] - data->imag_mt_TE[m][n]) * (data->calc_imag_mt_TE[m][n] - data->imag_mt_TE[m][n])));
							
					k++;
				}
			}
			/*TM-mode data*/
			else if(kind_of_data == 2)
			{
				for(n=0;n<data->nfreq_mt[m];n++)
				{
					data->rms_mt_real = data->rms_mt_real + (((data->calc_real_mt_TE[m][n] + data->real_mt_TM[m][n]) * (data->calc_real_mt_TE[m][n] + data->real_mt_TM[m][n])));
					data->rms_mt_imag = data->rms_mt_imag + (((data->calc_imag_mt_TE[m][n] + data->imag_mt_TM[m][n]) * (data->calc_imag_mt_TE[m][n] + data->imag_mt_TM[m][n])));
							
					k++;
				}
			}
			/*Berdichewsky average*/
			else
			{
				for(n=0;n<data->nfreq_mt[m];n++)
				{
					berd_real = 0.5*(data->real_mt_TE[m][n] - data->real_mt_TM[m][n]);
					berd_imag = 0.5*(data->imag_mt_TE[m][n] - data->imag_mt_TM[m][n]);
					data->rms_mt_real = data->rms_mt_real + (((data->calc_real_mt_TE[m][n] - berd_real) * (data->calc_real_mt_TE[m][n] - berd_real)));
					data->rms_mt_imag = data->rms_mt_imag + (((data->calc_imag_mt_TE[m][n] - berd_imag) * (data->calc_imag_mt_TE[m][n] - berd_imag)));
							
					k++;
				}
			}
		}

		if(data->ndata_mt != 0)
		{
			data->rms_mt_real = sqrt(data->rms_mt_real/k);
			data->rms_mt_imag = sqrt(data->rms_mt_imag/k);
		}
		else
		{
			data->rms_mt_real = -99999.9;
			data->rms_mt_imag = -99999.9;
		}

	printf("MT modeling is finished:\n");
	printf("RMS-values: Re-part: %10.5f ohm; Im-part: %10.5f ohm\n",data->rms_mt_real, data->rms_mt_imag);
	printf("----------------\n\n\n");

	
	fclose(out);

	return(1);
}

/*-------------------------------------------------------------*/
/*Calculate analytically the MT 1D response (TE-Mode) for a layered media*/
/*Parameters:		  ang_freq = angular frequency*/
/*	*z_real_part, *z_imag_part = real and imaginary part of the impedance (E/H)*/
/*				            *dz = thickness of the layers*/
/*					  nr_layer = Nr of layer*/
/*					   *resist = resistivities within the layers*/
/*                           p = horizontal wavernumber (if the "waves" are running vertical it is 0.0)*/

int	MT_1D_CALC(double ang_freq, double *z_real_part, double *z_imag_part, double *dz, long nr_layer, double *resist, double p)
{
	long r;
	double alpha_square_imag; /*imaginary part of the squared spatial wavenumber*/
	dcomplex alpha; /*spatial wavenumber in z-direction; because only the damped term of the telegraphic equation is considered, it is strongly complex*/
	dcomplex Q_value; /*Ratio: -nu*E_y/i*w*B_x = nu/(sqrt(i*w*nu*sigma))*/
	dcomplex cmplx_nu, imag_one, nominator, denominator;
	dcomplex h_times_alpha, nu_div_alpha, alpha_div_nu, tanh_div_alpha, tanh_times_Qvalue;

	cmplx_nu = Complex((M_KONST),0);
	imag_one = Complex(0,1);

	/*calculate the response from the deepest layer (homogeneous half-space)*/
	alpha_square_imag = (M_KONST)*ang_freq*(1/resist[nr_layer - 1]); /*Imag.part of squared spatial wavenumber a*/
	alpha = Complex(p*p,alpha_square_imag); /*Make a complex number*/
	alpha = Csqrt(alpha);				  /*Calculate the complex spatial wavenumber a*/

	Q_value = Cdiv(cmplx_nu,alpha);		/*Calculate the ratio: nu*E_y/i*w*B_x*/

	/**********************/
	/*calculate the recursively the response from the other layers*/
	/*(Backward)loop over the used cells in z-direction (ONLY cells benneath the station are considered)*/
	for(r=(nr_layer -2);r >= 0; r--)
	{
		/*Recursive calculation of the Q-ratio*/

		alpha_square_imag = (M_KONST)*ang_freq*(1/resist[r]); /*Imag.part of squared spatial wavenumber a*/
		alpha = Complex(0,alpha_square_imag); /*Make a complex number*/
		alpha = Csqrt(alpha);				  /*Calculate the complex spatial wavenumber a*/

		/*Determine the tanh*/
			h_times_alpha = RCmul(dz[r],alpha);	/*Determine the product of slowness and thickness of the layer/cells*/
			h_times_alpha = CBtanh(h_times_alpha);		/*Calculate the tangens hyperbolicus*/

		/*Calculate the nominator*/
			nu_div_alpha = Cdiv(cmplx_nu,alpha);				/*Determine the quotient of nu_0 and the spatial wavenumber a*/
			tanh_div_alpha = Cmul(h_times_alpha,nu_div_alpha);	/*Determine the product of tanh(...)*nu_0/a */
			nominator = Cadd(tanh_div_alpha,Q_value);

		/*Calculate the denominator*/
			alpha_div_nu = Cdiv(alpha,cmplx_nu);				/*Determine the quotient of the spatial wavenumber a and nu_0*/
			tanh_times_Qvalue = Cmul(h_times_alpha,Q_value);	/*Determine the product of tanh(...)*Q_value (of the former iteration)*/
			denominator = Cmul(tanh_times_Qvalue,alpha_div_nu);
			denominator.r = denominator.r + 1.0;

		/*Determine the new Q-ratio*/
			Q_value = Cdiv(nominator,denominator);

	}

	/*Determine the impedance from the Q-ratio*/
	Q_value = Cmul(Q_value, imag_one);
	*z_real_part = ((-1)*ang_freq*Q_value.r);
	*z_imag_part = ((-1)*ang_freq*Q_value.i);

	return(1);
}


/*-------------------------------------------------------------*/
/*Performing the 2D MT modeling:*/
/*Parameter:	geo  := Geometry structure  */
/*              grid  := Grid structure */
/*              *data  := Pointer on data structure (the MT impedances will be (re-)determined in this routine)*/
/*       2D_mt_struct  := Structure including 2D MT specific parameters (output data)*/
/*                flag := Structure of important flags */
/*				  topo := Topography structure*/

/*Output: Number of 2D slices*/

/*Remark: the magnetic permeability is considered to correspond to the one of free air*/

/****************************************************/
/****************************************************/
#define NR_OF_1D_CALCULATION 5		/*Nr of equally spaced positions, where the 1-D response was calculated to determine the number of cells in the air by means*/
									/*of the largest apparent resistivity*/

#define GLOBAL_FACTOR_REFINEMENT 0.3		/*This Factor will be multiplied by the minimum skin depth found in the grid. The resulting measure
											is an estimate for the chosen grid cell size (see Weaver, 1994 or, Simpson and Bahr, 2005)*/
#define LOCAL_FACTOR_REFINEMENT 0.03       /*Local refinement of the grid used only in the vicinity of the stations*/
#define LOCAL_AREA_FACTOR 0.25               /*This factor multiplied with the minimum skin depth obtained from the 1-D response specify the thickness of the highly refined grid around the stations*/

/****************** IONOSPHERE **********************/
#define RES_IONOS	1				/*Resistivity of the ionosphere in ohmm*/

#define DZ_IONOS 1000				/*Thickness of the ionosphere in m (following the code from Tartis)*/

/****************** ATHMOSPHERE **********************/
#define RES_AIR	1E10				/*considered resistivity of air in OHMm; (used for the calculation of the max. skin depths) <-EVENTUELL NOCH MAL NEU FESTLEGEN*/

#define FACTOR_AIR 1.5                /*This factor will multiplied with the maximum skin depth obtained from the app.resistivities to determine the extension of the grid
									in the air*/
#define FACTOR_SCAL_AIR 2			/*Factor specifying the "grow" of the cells in the air (size is used in the code from Tartis)*/

#define MIN_SIZE_AIR 1E5			/*Minimum thickness of the air layer in m*/
/************* BASEMENT **********************/
#define FACTOR_BASEMENT 30			/*This factor will multiplied with the maximum skin depth at the basement to determine the extension of the grid
									toward larger depths*/
#define FACTOR_SCAL_BAS 1.25		/*Factor specifying the "grow" of the cells in the basement (size is used in the code from Tartis)*/

/************* LATERAL BORDER ****************/
#define FACTOR_BORDER_LEFT_RIGHT 4	/*This factor will be multiplied with the maximum skin depth at the left and right border to determine the extension
									of the grid towards the left and right*/
#define FACTOR_SCAL_BORDER 1.5 	   /*Factor specifying the "grow" of the cells within the left and right border (size is used in the code from Tartis)*/
/*********************************************/
/*********************************************/

int ForwardModMT_2D(GEOMETRY geo, GRID_STRUCT grid, DATA_STRUCT *data, CALC_2D_MT *mt_2D_struct_out, FLAG_STRUCT flag, TOPO_STRUCT topo)
{
	int kind_of_data;
	long nr_of_slices;       /*Number of slices*/
	long a,b,c,i,j,k,m,n,nsamples;
	long first_sample;
	long ix,iz;
	long nx,nz;
	long ny2,nz2,nyz2;
	long nz_air,nz_ionos;	 /*Nr of additional cells in the air*/
	long nx_left_border, nx_right_border; /*Nr of cells at the left and right border*/
	long nz_basement;	/*Nr of cells at the basement*/
	double berd_real, berd_imag, berd_real_calc, berd_imag_calc;
	double min_res;				/*Minimum resistivity in the complete slice*/
	double max_res_basement;	/*Maximum resistivity at the basement*/
	double max_res_left, max_res_right; /*Maximum resistivity at the left and right border of the grid*/
	double  *org_res_slice, *org_index_slice, *tmp_3D_b_index;
	double  T, res_ionos;
	double *grid_points; /*z-coordinate of the intersections of the middle of the slices with the topography*/
	double *Ex_r,*Ex_i,*Hy_r,*Hy_i; /*Field components at the upper interface (NOT THE CENTER !!!!) of the cells*/ 
	double *Hx_r,*Hx_i,*Ey_r,*Ey_i; /*Field components at the upper interface (NOT THE CENTER !!!!) of the cells*/ 
	double *Ez_r,*Ez_i;
	double E_abs, H_abs,E_phase, H_phase, app_res1;

	CALC_2D_MT mt_2D_struct[1]; /*MT structure that will be determined in this routine*/

//char fname1[32],fname2[32];
//FILE *out,*out1;
	
	printf("2-D\n\n");

	for(i=0;i<data->ndata_mt;i++)
	{
		for(j=0;j<data->nfreq_mt[i];j++)
		{
			data->calc_real_mt_TE[i][j] = 0.0; 
			data->calc_imag_mt_TE[i][j] = 0.0;
			data->calc_real_mt_TM[i][j] = 0.0; 
			data->calc_imag_mt_TM[i][j] = 0.0;
		}
	}

	if(data->ndata_mt == 0)
		printf("!!WARNING!! NO MT data exists but the resistivity forward modelling is activated\n\n");


	kind_of_data = flag.kind_of_data_mt;
	first_sample = 0;

	/************************************************************/		
	/*Specify the number of 2-D slices*/
		/*Slice along x-direction*/
	if(flag.direc_2D_mt != 2) 
		nr_of_slices = (long) ceil((double)grid.ny/(double)(flag.nr_cells_2d_mt)); /*Number of slices*/
		/*Slice along y-direction*/
	else
		nr_of_slices = (long) ceil((double)grid.nx/(double)(flag.nr_cells_2d_mt)); /*Number of slices*/

	/*******************************************************************************************/
	/*******************************************************************************************/
								/*START LOOP OVER ALL SLICES*/
	/*******************************************************************************************/
	/*******************************************************************************************/

	/*Start loop over all slices*/
	for(i=0;i<nr_of_slices;i++)
	{
		
		/*Specify the number of cells in a slice*/
		mt_2D_struct->nr_cells_2d_mt = flag.nr_cells_2d_mt;

		/*Specify the 2-D parameters of the 2D slices*/

			/*Slice along x-direction*/
		if(flag.direc_2D_mt != 2) 
			nx = grid.nx;
			/*Slice along y-direction*/
		else
			nx = grid.ny;

		nz = grid.nz;
		
		/************************************************************/
		/*Determine the "thickness" of the slice*/
		DetThickness2DSlice(mt_2D_struct,grid,i,flag.direc_2D_mt, flag.nr_cells_2d_mt);
		/*******************************************************************************************/
		/*Check, if any MT-stations are located within the area/tube of the slice*/
		if(data->ndata_mt == 0)
			printf("!!WARNING!! NO MT data exists but the resistivity forward modelling is activated\n\n");

		/*Allocate the memory for the list of station numbers, and the x,z coordinates*/
		mt_2D_struct->mts = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");
		mt_2D_struct->x = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
		mt_2D_struct->z = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");

		DetStation2DSlice(data,geo,mt_2D_struct, flag.direc_2D_mt);

		/*******************************************************************************************/
		/*Allocate the memory for the resistivities/indeces along the slice*/
		org_res_slice = (double *)memory(NULL,nx*nz,sizeof(double),"ForwardModMT_2D");
		org_index_slice = (double *)memory(NULL,nx*nz,sizeof(double),"ForwardModMT_2D");

		tmp_3D_b_index = (double *)memory(NULL,(grid.nx*grid.ny*grid.nz),sizeof(double),"ForwardModMT_2D");

		ny2 = 2*grid.nborder + grid.ny;
		nz2 = 2*grid.nborder + grid.nz;
		nyz2 = ny2*nz2;

		/*Adjusting the index that specify border cells*/
		for(a=0;a<grid.nx;a++)
			for(b=0;b<grid.ny;b++)
				for(c=0;c<grid.nz;c++)
				{
					tmp_3D_b_index[a*(grid.ny*grid.nz) + b*grid.nz + c] = (double)grid.border_index[(a+grid.nborder)*nyz2 + (b+grid.nborder)*nz2 +(c+grid.nborder)];
				}
		
		/*Calculate the resistivity values/indeces within each cell by averaging the cells in the considered rod vertical two the considered 2D-plane*/
		CalcAveValues2DSlice(grid, first_sample, flag.direc_2D_mt, flag.nr_cells_2d_mt, nx, nz, org_res_slice, grid.res);
		CalcAveValues2DSlice(grid, first_sample, flag.direc_2D_mt, flag.nr_cells_2d_mt, nx, nz, org_index_slice, tmp_3D_b_index);

		free(tmp_3D_b_index);

		grid_points = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
		/*******************************************************************************************/

		/*Allocate the memory for the frequencies used for the stations in the slice*/
		mt_2D_struct->freq = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
		mt_2D_struct->freq_ndata = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");
		mt_2D_struct->freq_data = (long **)memory(NULL,1,sizeof(long *),"ForwardModMT_2D");
		mt_2D_struct->freq_stat = (long **)memory(NULL,1,sizeof(long *),"ForwardModMT_2D");

		mt_2D_struct->Ex_r = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
		mt_2D_struct->Ex_i = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
		mt_2D_struct->Ey_r = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
		mt_2D_struct->Ey_i = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
		mt_2D_struct->Hx_r = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
		mt_2D_struct->Hx_i = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
		mt_2D_struct->Hy_r = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
		mt_2D_struct->Hy_i = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");

		/*Determine the frequencies (+ indeces of the measurements and the station numbers)*/
		AdjustFreq2Dslice(mt_2D_struct, data);

		/*******************************************************************************************/
		/*Find the minimum resistivity in the complete slice*/
		min_res = org_res_slice[0]; /*in OHMm*/

			/*Loop over all cells*/
		for(k=0;k<(nx*nz);k++)
		{
			if(min_res > org_res_slice[k])
				min_res = org_res_slice[k];
		}

		if(min_res < 0.1)
		{
			printf("WARNING!!! There are very low resistivities < 0.1 OHMm found in the grid!!!\n");
			printf("Are these values okay?\n\n");
		}

		/*A zero resistivity will lead to an infinite small min. skin depth*/
		if(min_res < 0.0)
			min_res = 0.0000001;

		/********************************/
		/*Find the maximum resistivity at the basement*/
		max_res_basement = org_res_slice[nz -1]; /*in OHMm*/

			/*Loop over all cells at the basement*/
		for(ix=0;ix<nx;ix++)
		{
			if(max_res_basement < org_res_slice[(ix+1)*nz -1])
				max_res_basement = org_res_slice[(ix+1)*nz -1];
		}

		if(max_res_basement > RES_AIR) /*Check if the resistivity is larger than the one of air*/
			max_res_basement = RES_AIR;

		/********************************/
		/*Find the maximum resistivity at the left and right border*/
		max_res_left = org_res_slice[0];
		max_res_right = org_res_slice[(nx -1)*(nz)];

			/*Loop over all cells at the left and right border*/
		for(iz=0;iz<nz;iz++)
		{
			if(max_res_left < org_res_slice[iz])
				max_res_left = org_res_slice[iz];

			if(max_res_right < org_res_slice[(nx -1)*nz + iz])
				max_res_right = org_res_slice[(nx -1)*nz + iz];
		}

		/*The resistivity of "air" is the maximum resistivity*/

		if(max_res_left > RES_AIR) /*Check if the resistivity is larger than the one of air*/
			max_res_left = RES_AIR;
		if(max_res_right > RES_AIR) /*Check if the resistivity is larger than the one of air*/
			max_res_right = RES_AIR;

		/********************************/
			/*Allocate some memory for parameters of the MT structure*/

		if(mt_2D_struct->nfreq != 0)
		{
			mt_2D_struct->org[0] = (double *)memory(NULL,mt_2D_struct->nfreq,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->org[1] = (double *)memory(NULL,mt_2D_struct->nfreq,sizeof(double),"ForwardModMT_2D");

			mt_2D_struct->nx = (long *)memory(NULL,mt_2D_struct->nfreq,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->nz = (long *)memory(NULL,mt_2D_struct->nfreq,sizeof(long),"ForwardModMT_2D");

			mt_2D_struct->nz_atm = (long *)memory(NULL,mt_2D_struct->nfreq,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->nz_ionos = (long *)memory(NULL,mt_2D_struct->nfreq,sizeof(long),"ForwardModMT_2D");

			mt_2D_struct->nz_basement = (long *)memory(NULL,mt_2D_struct->nfreq,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->nx_right_border = (long *)memory(NULL,mt_2D_struct->nfreq,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->nx_left_border = (long *)memory(NULL,mt_2D_struct->nfreq,sizeof(long),"ForwardModMT_2D");

			mt_2D_struct->hx = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"ForwardModMT_2D");
			mt_2D_struct->hz = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"ForwardModMT_2D");

			mt_2D_struct->res_slice = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"ForwardModMT_2D");
			mt_2D_struct->index_slice = (int **)memory(NULL,mt_2D_struct->nfreq,sizeof(int *),"ForwardModMT_2D");

			mt_2D_struct->cell_scal_factor = (int *)memory(NULL,mt_2D_struct->nfreq,sizeof(int),"ForwardModMT_2D");

			for(j=0;j<mt_2D_struct->nfreq;j++)
			{
				mt_2D_struct->res_slice[j] = (double *)memory(NULL,nx*nz,sizeof(double),"ForwardModMT_2D");
				mt_2D_struct->index_slice[j] = (int *)memory(NULL,nx*nz,sizeof(int),"ForwardModMT_2D");
				mt_2D_struct->hx[j] = (double *)memory(NULL,nx,sizeof(double),"ForwardModMT_2D");
				mt_2D_struct->hz[j] = (double *)memory(NULL,nz,sizeof(double),"ForwardModMT_2D");
			}
		}
		else
		{
			mt_2D_struct->org[0] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->org[1] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");

			mt_2D_struct->nx = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->nz = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");
			
			mt_2D_struct->nz_atm = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->nz_ionos = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");

			mt_2D_struct->nz_basement = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->nx_right_border = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->nx_left_border = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");

			mt_2D_struct->hx = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
			mt_2D_struct->hz = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");

			mt_2D_struct->res_slice = (double **)memory(NULL,1,sizeof(double *),"ForwardModMT_2D");
			mt_2D_struct->index_slice = (int **)memory(NULL,1,sizeof(int *),"ForwardModMT_2D");

			mt_2D_struct->cell_scal_factor = (int *)memory(NULL,1,sizeof(int),"ForwardModMT_2D");
			
		}

		printf("Start the calculation for the layer %d\n",i+1);


		/*******************************************************************************************/
		/*******************************************************************************************/
									/*START LOOP OVER ALL FREQUENCIES*/
		/*******************************************************************************************/
		/*******************************************************************************************/

		for(j=0;j<mt_2D_struct->nfreq;j++)
		{

			nz_air = 0;
			nz_ionos = 0;
			nx_left_border = 0;
			nx_right_border = 0;
			nz_basement = 0;

//sprintf(fname1,"test_grid%7f.txt",mt_2D_struct->freq[j]);
//sprintf(fname2,"res_value%7f.txt",mt_2D_struct->freq[j]);
//out = fopen(fname1,"w");
//out1 = fopen(fname2,"w");

			/*Re-set all parameters for the iteration*/

			/*Reset the number of cells*/
				/*Slice along x-direction*/
			if(flag.direc_2D_mt != 2) 
			{
				mt_2D_struct->nx[j] = grid.nx;
				mt_2D_struct->org[0][j] = grid.org[0] - 0.5*grid.h;
			}
				/*Slice along y-direction*/
			else
			{
				mt_2D_struct->nx[j] = grid.ny;
				mt_2D_struct->org[0][j] = grid.org[1] - 0.5*grid.h;
			}

			mt_2D_struct->nz[j] = grid.nz;
			mt_2D_struct->org[1][j] = grid.org[2] - 0.5*grid.h;

			for(k=0;k<nx;k++)
				mt_2D_struct->hx[j][k] = grid.h;
			for(k=0;k<nz;k++)
				mt_2D_struct->hz[j][k] = grid.h;

			/*******************************************************************************************/
			/* Make the global refined MT-Grid governed by the minimum resistivity in the grid*/
			MakeRefinedGridMT(grid.h, mt_2D_struct, min_res, org_res_slice, org_index_slice, j);
			
			/*****************************************************************/
			/*Make an additional refinement close to the stations governed by the minimum resistivity
			  in the grid and the skin depth of the 1D MT calculation at the stations*/
			MakeLocalRefinedGridMT(grid.h, mt_2D_struct, min_res, i, j);

			/*******************************************************************************************/
			/*Add the topography information:*/
			/*Determine the z-coordinates of intersections of the slice (in the middle) with the topography*/
			grid_points = (double *)memory((char *)grid_points,(mt_2D_struct->nx[j])+1,sizeof(double),"ForwardModMT_2D");
			DetTopography2DSlice(grid, first_sample, flag.direc_2D_mt, flag.nr_cells_2d_mt, mt_2D_struct->nx[j], mt_2D_struct->hx[j] ,topo, grid_points);

			/*Add the topography information, which comes from the positions of the MT stations*/
			DetTopographyMTStat2DSlice(geo, grid_points, flag.direc_2D_mt, mt_2D_struct->nx[j], mt_2D_struct->hx[j], mt_2D_struct->org[0][j], mt_2D_struct->nstat_mt, mt_2D_struct->mts);

			/*Assign the topography information to the cells*/
			AddTopography2DSlice(mt_2D_struct->nx[j], mt_2D_struct->nz[j], mt_2D_struct->hx[j], grid_points, mt_2D_struct->org[1][j], mt_2D_struct->res_slice[j], grid.res_air_water ,mt_2D_struct->index_slice[j]);

			/*******************************************************************************************/
			/*Add boundary cells at the sides and at the basement*/
			/*Obtained from the maximum skin depths at the boundaries*/
			/*and the skin depths from the resistivity at the bottom of the grid*/
			AddCellsAtBoundaryMT(mt_2D_struct, max_res_basement, max_res_left, max_res_right,j, &nx_left_border, &nx_right_border, &nz_basement);

			/*******************************************************************************************/
			/*******************************************************************************************/
				/*TM-MODE*/
			if(kind_of_data != 1)
			{

				/*Allocate memory for the Ex and Hy components of the fields*/
				Ex_r = (double *)memory(NULL,mt_2D_struct->nx[j]*mt_2D_struct->nz[j], sizeof(double),"ForwardModMT_2D");
				Ex_i = (double *)memory(NULL,mt_2D_struct->nx[j]*mt_2D_struct->nz[j], sizeof(double),"ForwardModMT_2D");
				Ez_r = (double *)memory(NULL,mt_2D_struct->nx[j]*mt_2D_struct->nz[j], sizeof(double),"ForwardModMT_2D");
				Ez_i = (double *)memory(NULL,mt_2D_struct->nx[j]*mt_2D_struct->nz[j], sizeof(double),"ForwardModMT_2D");
				Hy_r = (double *)memory(NULL,mt_2D_struct->nx[j]*mt_2D_struct->nz[j], sizeof(double),"ForwardModMT_2D");
				Hy_i = (double *)memory(NULL,mt_2D_struct->nx[j]*mt_2D_struct->nz[j], sizeof(double),"ForwardModMT_2D");

				for(k=0;k<(mt_2D_struct->nx[j]*mt_2D_struct->nz[j]);k++)
				{
					Ex_r[k]= 0.0;
					Ex_i[k]= 0.0;
					Ez_r[k]= 0.0;
					Ez_i[k]= 0.0;
					Hy_r[k]= 0.0;
					Hy_i[k]= 0.0;
				}

				T = 1/mt_2D_struct->freq[j]; /*Periode in s*/

				/*Perform the calculation for the TM-Mode*/
				/*Using routine from 2-D finite difference code (frequency domain) of Tarits,1987 */
				/*ATTENTION !!! The response will be determined at the upper interface of the cells*/
				int temp_nx = mt_2D_struct->nx[j];
				int temp_nz = mt_2D_struct->nz[j];
				hpol_(&T, &temp_nx, &temp_nz, mt_2D_struct->hx[j], mt_2D_struct->hz[j], mt_2D_struct->res_slice[j], Hy_r, Hy_i, Ex_r, Ex_i, Ez_r, Ez_i);

				mt_2D_struct->nx[j] = temp_nx;
				mt_2D_struct->nz[j] = temp_nz;
				printf("The TM-Mode is calculated for the %d frequency: %fHz\n",j+1,mt_2D_struct->freq[j]);
                 
				/***************************************************/
				/*Calculate the apparent resistivities*/
				for(k=0;k<mt_2D_struct->nx[j]*mt_2D_struct->nz[j];k++)
				{
                   E_abs = sqrt(Ex_r[k]*Ex_r[k] + Ex_i[k]*Ex_i[k]); /*absolute values of the fields*/
				   H_abs = sqrt(Hy_r[k]*Hy_r[k] + Hy_i[k]*Hy_i[k]);
				   E_phase = (180/PI) * atan2(Ex_i[k],Ex_r[k]); /*phases of the fields*/
				   H_phase = (180/PI) * atan2(Hy_i[k],Hy_r[k]);

				   app_res1 = ((E_abs*E_abs)/(H_abs*H_abs))*(1/(2*PI*mt_2D_struct->freq[j]*(M_KONST))); /*app.resistivities*/
				}

				/*Determine the impedances at the station positions*/
				DetImpAtStations(mt_2D_struct,data,Hy_r,Hy_i,Ex_r,Ex_i,2,0,j);

				/*Write the resistivity models in a temporary file (only if 2D inversion will be performed lateron)*/
				nsamples = (mt_2D_struct->nx[j]*(mt_2D_struct->nz[j] - nz_air - nz_ionos));

				if(flag.dimension_mt == 2) /*Using the analytical 2-D formula for the derivatives (see Jegen)*/
				{
					WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, Ex_r, 1);
					WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, Ex_i, 2);

					WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, Ez_r, 9);
					WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, Ez_i, 10);
				}
				else if(flag.dimension_mt == 4) /*Using the generalized RRI (see Yamane et al.1996)*/
				{
					WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, Hy_r, 7);
					WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, Hy_i, 8);
				}

				free(Ex_r);
				free(Ex_i);
				free(Ez_r);
				free(Ez_i);
				free(Hy_r);
				free(Hy_i);
			}

			/*******************************************************************************************/
			/*******************************************************************************************/
				/*TE-MODE*/
			if(kind_of_data != 2)
			{

				/*******************************************************************************************/
				/*Add boundary cells in the air*/
				/*Obtained from the maximum skin depths from 1D calculations*/
				nz_air = 0;
				nz_ionos = 0;
				AddCellsAtBoundaryAirMT(mt_2D_struct, j, &nz_air, &nz_ionos);

				/******************************************************************************/

				/*Allocate memory for the Ex and Hy components of the fields*/
				Hx_r = (double *)memory(NULL,mt_2D_struct->nx[j]*(mt_2D_struct->nz[j] - nz_air - nz_ionos), sizeof(double),"ForwardModMT_2D");
				Hx_i = (double *)memory(NULL,mt_2D_struct->nx[j]*(mt_2D_struct->nz[j] - nz_air - nz_ionos), sizeof(double),"ForwardModMT_2D");
				Ey_r = (double *)memory(NULL,mt_2D_struct->nx[j]*(mt_2D_struct->nz[j] - nz_air - nz_ionos), sizeof(double),"ForwardModMT_2D");
				Ey_i = (double *)memory(NULL,mt_2D_struct->nx[j]*(mt_2D_struct->nz[j] - nz_air - nz_ionos), sizeof(double),"ForwardModMT_2D");

				for(k=0;k<(mt_2D_struct->nx[j]*(mt_2D_struct->nz[j] - nz_air - nz_ionos));k++)
				{
					Ey_r[k]= 0.0;
					Ey_i[k]= 0.0;
					Hx_r[k]= 0.0;
					Hx_i[k]= 0.0;
				}

				T = 1/mt_2D_struct->freq[j]; /*Periode in s*/
				res_ionos = RES_IONOS; /*Assumed resistivity of the ionosphere in ohmm*/

				/*Perform the calculation for the TE-Mode*/
				/*Using routine from 2-D finite difference code (frequency domain) of Tarits,1987 */
				/*ATTENTION !!! The response will be determined at the upper interface of the cells*/
				epol_(&T, &mt_2D_struct->nx[j], &mt_2D_struct->nz[j], &nz_ionos, &nz_air, mt_2D_struct->hx[j], mt_2D_struct->hz[j], mt_2D_struct->res_slice[j], &res_ionos, Hx_r, Hx_i, Ey_r, Ey_i);

				printf("The TE-Mode is calculated for the %d frequency: %fHz\n",j+1,mt_2D_struct->freq[j]);

				/***************************************************/
				/*Calculate the apparent resistivities*/
				for(k=0;k<mt_2D_struct->nx[j]*(mt_2D_struct->nz[j]- nz_air - nz_ionos);k++)
				{
                   E_abs = sqrt(Ey_r[k]*Ey_r[k] + Ey_i[k]*Ey_i[k]); /*absolute values of the fields*/
				   H_abs = sqrt(Hx_r[k]*Hx_r[k] + Hx_i[k]*Hx_i[k]);
				   E_phase = (180/PI) * atan2(Ey_i[k],Ey_r[k]); /*phases of the fields*/
				   H_phase = (180/PI) * atan2(Hx_i[k],Hx_r[k]);

				   app_res1 = ((E_abs*E_abs)/(H_abs*H_abs))*(1/(2*PI*mt_2D_struct->freq[j]*(M_KONST))); /*app.resistivities*/

/*************/
 //fprintf(out1,"%d %12.10f %12.10f %f %f %f %f\n",k,E_abs,H_abs,E_phase,H_phase,app_res1,mt_2D_struct->res_slice[j][k]);
/*************/

				}

/*************/
/*fprintf(out,"%f %d %d %d+1\n\n",mt_2D_struct->freq[j], mt_2D_struct->nx[j], mt_2D_struct->nz[j], nz_air);

for(k=0;k<(mt_2D_struct->nx[j]);k++)
	fprintf(out,"%d %f\n",k,mt_2D_struct->hx[j][k]);

fprintf(out,"\n");

for(k=0;k<(mt_2D_struct->nz[j]);k++)
	fprintf(out,"%d %f\n",k,mt_2D_struct->hz[j][k]);
/*************/

				/*Determine the impedances at the station positions*/
				DetImpAtStations(mt_2D_struct,data,Hx_r,Hx_i,Ey_r,Ey_i,1,(nz_air + nz_ionos),j);

				/*Write the resistivity models in a temporary file (only if afterwards a 2-D inversion is performed)*/
				nsamples = (mt_2D_struct->nx[j]*(mt_2D_struct->nz[j] - nz_air - nz_ionos));

				if(flag.dimension_mt == 2 || flag.dimension_mt == 4)/*Using the analytical 2-D formula (see Jegen)*/
																	/*or the RRI (see Smith and Booker,1991) for the derivatives*/
				{
					WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, Ey_r, 3);
					WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, Ey_i, 4);
				}

				free(Hx_r);
				free(Hx_i);
				free(Ey_r);
				free(Ey_i);
			}

//fclose(out);
//fclose(out1);

			mt_2D_struct->nz_atm[j] = nz_air;
			mt_2D_struct->nz_ionos[j] = nz_ionos;
			mt_2D_struct->nx_left_border[j] = nx_left_border;
			mt_2D_struct->nx_right_border[j] = nx_right_border;
			mt_2D_struct->nz_basement[j] = nz_basement;

			/*Write the resistivity models in a temporary file*/
			nsamples = (mt_2D_struct->nx[j]*(mt_2D_struct->nz[j] - nz_air - nz_ionos));

			if(flag.dimension_mt == 2 || flag.dimension_mt == 4)
			{
				WriteMTResistivityOut(i, mt_2D_struct->freq[j], nsamples, mt_2D_struct->res_slice[j], 0);
				WriteMTResistivityIndexOut(i, mt_2D_struct->freq[j], nsamples, mt_2D_struct->index_slice[j]);
			}

			/*Transfer the local MT structure to the global MT structure*/
			MakeGlobMTStruct(i,j,mt_2D_struct,mt_2D_struct_out);
		
		}
		/*******************************************************************************************/
		/*******************************************************************************************/
								/*END OF THE LOOP OVER ALL FREQUENCIES*/
		/*******************************************************************************************/
		/*******************************************************************************************/

		first_sample = first_sample + mt_2D_struct->nr_cells_2d_mt;

		if(mt_2D_struct->nfreq == 0)
			/*Transfer the local MT structure to the global MT structure*/
			MakeGlobMTStruct(i,0,mt_2D_struct,mt_2D_struct_out);

		for(j=0;j<mt_2D_struct->nfreq;j++)
		{
			free(mt_2D_struct->freq_data[j]);
			free(mt_2D_struct->freq_stat[j]);
			free(mt_2D_struct->hx[j]);
			free(mt_2D_struct->hz[j]);
			free(mt_2D_struct->res_slice[j]);
			free(mt_2D_struct->index_slice[j]);

			free(mt_2D_struct->Ex_r[j]);
			free(mt_2D_struct->Ex_i[j]);
			free(mt_2D_struct->Ey_r[j]);
			free(mt_2D_struct->Ey_i[j]);
			free(mt_2D_struct->Hx_r[j]);
			free(mt_2D_struct->Hx_i[j]);
			free(mt_2D_struct->Hy_r[j]);
			free(mt_2D_struct->Hy_i[j]);

		}
		
		free(mt_2D_struct->freq);
		free(mt_2D_struct->freq_ndata);
		free(mt_2D_struct->freq_data);
		free(mt_2D_struct->freq_stat);
		free(mt_2D_struct->res_slice);
		free(mt_2D_struct->index_slice);
		free(mt_2D_struct->mts);
		free(mt_2D_struct->x);
		free(mt_2D_struct->z);
		free(mt_2D_struct->nx);
		free(mt_2D_struct->nz);
		free(mt_2D_struct->nz_atm);
		free(mt_2D_struct->nz_ionos);
		free(mt_2D_struct->nz_basement);
		free(mt_2D_struct->nx_left_border);
		free(mt_2D_struct->nx_right_border);
		free(mt_2D_struct->hx);
		free(mt_2D_struct->hz);
		free(mt_2D_struct->org[0]);
		free(mt_2D_struct->org[1]);
		free(mt_2D_struct->cell_scal_factor);

		free(mt_2D_struct->Ex_r);
		free(mt_2D_struct->Ex_i);
		free(mt_2D_struct->Ey_r);
		free(mt_2D_struct->Ey_i);
		free(mt_2D_struct->Hx_r);
		free(mt_2D_struct->Hx_i);
		free(mt_2D_struct->Hy_r);
		free(mt_2D_struct->Hy_i);

		free(org_res_slice);
		free(org_index_slice);
		free(grid_points);
	}
	/*******************************************************************************************/
	/*******************************************************************************************/
							/*END OF THE LOOP OVER ALL SLICES*/
	/*******************************************************************************************/
	/*******************************************************************************************/

	
	/*Calculate the RMS values of the real and imaginary part of Z*/
	data->rms_mt_real = 0.0;
	data->rms_mt_imag = 0.0;

	k=0;

	/*Loop over all data*/
	for(m=0;m<data->ndata_mt;m++)
	{
		/*TE-mode data*/
		if(kind_of_data == 1)
		{
			/*Loop over all frequencies*/
			for(n=0;n<data->nfreq_mt[m];n++)
			{
				data->rms_mt_real = data->rms_mt_real + (((data->calc_real_mt_TE[m][n] - data->real_mt_TE[m][n]) * (data->calc_real_mt_TE[m][n] - data->real_mt_TE[m][n])));
				data->rms_mt_imag = data->rms_mt_imag + (((data->calc_imag_mt_TE[m][n] - data->imag_mt_TE[m][n]) * (data->calc_imag_mt_TE[m][n] - data->imag_mt_TE[m][n])));
					
				k++;
			}
		}
		/*TM-mode data*/
		else if(kind_of_data == 2)
		{
			for(n=0;n<data->nfreq_mt[m];n++)
			{
				data->rms_mt_real = data->rms_mt_real + (((data->calc_real_mt_TM[m][n] - data->real_mt_TM[m][n]) * (data->calc_real_mt_TM[m][n] - data->real_mt_TM[m][n])));
				data->rms_mt_imag = data->rms_mt_imag + (((data->calc_imag_mt_TM[m][n] - data->imag_mt_TM[m][n]) * (data->calc_imag_mt_TM[m][n] - data->imag_mt_TM[m][n])));
						
				k++;
			}
		}
		/*Berdichewsky average*/
		else
		{
			for(n=0;n<data->nfreq_mt[m];n++)
			{
				berd_real = 0.5*(data->real_mt_TE[m][n] - data->real_mt_TM[m][n]);
				berd_imag = 0.5*(data->imag_mt_TE[m][n] - data->imag_mt_TM[m][n]);
				berd_real_calc = 0.5*(data->calc_real_mt_TE[m][n] - data->calc_real_mt_TM[m][n]);
				berd_imag_calc = 0.5*(data->calc_imag_mt_TE[m][n] - data->calc_imag_mt_TM[m][n]);

				data->rms_mt_real = data->rms_mt_real + ((berd_real_calc - berd_real) * (berd_real_calc - berd_real));
				data->rms_mt_imag = data->rms_mt_imag + ((berd_imag_calc - berd_imag) * (berd_imag_calc - berd_imag));
							
				k++;
			}
		}
	}

	if(data->ndata_mt != 0)
	{
		data->rms_mt_real = sqrt(data->rms_mt_real/k);
		data->rms_mt_imag = sqrt(data->rms_mt_imag/k);
	}
	else
	{
		data->rms_mt_real = -99999.9;
		data->rms_mt_imag = -99999.9;
	}

	printf("MT modeling is finished:\n");
	printf("RMS-values: Re-part: %10.5f ohm; Im-part: %10.5f ohm\n",data->rms_mt_real, data->rms_mt_imag);
	printf("----------------\n\n\n");

	/****************************************/
	/*Write out the impedances at the stations*/
/*	out2 = fopen("MT_response_calc.txt","w");

	if(kind_of_data != 2)
	{
		for(m=0;m<data->ndata_mt;m++)
			for(k=0;k<data->nfreq_mt[m];k++)
			{
				fprintf(out2,"%d %f %10.8f %10.8f\n",m,data->freq_mt[m][k],data->calc_real_mt_TE[m][k],data->calc_imag_mt_TE[m][k]);
			}
	}

	fprintf(out2,"\n");

	if(kind_of_data != 1)
	{
		for(m=0;m<data->ndata_mt;m++)
			for(k=0;k<data->nfreq_mt[m];k++)
			{
				fprintf(out2,"%d %f %10.8f %10.8f\n",m,data->freq_mt[m][k],data->calc_real_mt_TM[m][k],data->calc_imag_mt_TM[m][k]);
			}
	}

	fclose(out2);
	/****************************************/


	return(nr_of_slices);
}


/*-------------------------------------------------------------*/
/*Determine the thickness of the considered 2-D slice*/
/*Parameter:     grid  := Grid structure */
/*       2D_mt_struct  := Structure including 2D MT specific parameters*/
/*             layer   := nr of layer*/
/*       direc_2D_mt   := direction of the slices (1=x-axis, 2=y-axis)*/
/*     nr_cells_2D_mt  := Number of cells per layer (specified in the input file)*/

int DetThickness2DSlice(CALC_2D_MT *mt_2D_struct, GRID_STRUCT grid, long layer, int direc_2D_mt, long nr_cells_2D_mt)
{
	/*Slice along x-direction*/
	if(direc_2D_mt != 2)
	{
		mt_2D_struct->min_y = (grid.org[1] - 0.5*grid.h) + (layer*nr_cells_2D_mt*grid.h);
		if((layer+1)*nr_cells_2D_mt <= grid.ny)
			mt_2D_struct->max_y = (grid.org[1] - 0.5*grid.h) + ((layer+1)*nr_cells_2D_mt*grid.h);
		else
		{
		    mt_2D_struct->max_y = (grid.org[1] - 0.5*grid.h) + grid.ny*grid.h;
				/*Modifying the width of the last layer*/
			mt_2D_struct->nr_cells_2d_mt = grid.ny - (layer*nr_cells_2D_mt);
		}
	}
		/*Slice along y-direction*/
	else
	{
		mt_2D_struct->min_y = (grid.org[0] - 0.5*grid.h) + (layer*nr_cells_2D_mt*grid.h);
		if((layer+1)*nr_cells_2D_mt <= grid.nx)
			mt_2D_struct->max_y = (grid.org[0] - 0.5*grid.h) + ((layer+1)*nr_cells_2D_mt*grid.h);
		else
		{
		    mt_2D_struct->max_y = (grid.org[0] - 0.5*grid.h) + grid.nx*grid.h;
				/*Modifying the width of the last layer*/
			mt_2D_struct->nr_cells_2d_mt = grid.nx - ((layer+1)*nr_cells_2D_mt);
		}
	}

	return(1);
}

/*-------------------------------------------------------------*/
/*Determine the station + coordinates in the considered 2D slice*/
/*Parameter:   2D_mt_struct  := Structure including 2D MT specific parameters*/
/*				geo  := Geometry structure  */
/*              *data  := Pointer on data structure (the MT impedances will be (re-)determined in this routine)*/
/*       direc_2D_mt   := direction of the slices (1=x-axis, 2=y-axis)*/

int DetStation2DSlice(DATA_STRUCT *data,GEOMETRY geo, CALC_2D_MT *mt_2D_struct, int direc_2D_mt)
{
	long j;
	float xs,ys,zs;

	mt_2D_struct->nstat_mt = 0;
		
	/*Loop over all stations*/
	for(j=0;j<geo.nstat_mt;j++)
	{
		/*Positions of the MT stations:*/
		/*Slice along x-direction*/
		if(direc_2D_mt != 2)
		{
				xs = (float)geo.x[(data->mts[j]-1)];
				ys = (float)geo.y[(data->mts[j]-1)];
		}
		/*Slice along y-direction*/
		else
		{
				xs = (float)geo.y[(data->mts[j]-1)];
				ys = (float)geo.x[(data->mts[j]-1)];
		}
				zs = (float)geo.z[(data->mts[j]-1)]; 
			
		/*Find the stations in the considered slice*/
		if(ys <= mt_2D_struct->max_y && ys >= mt_2D_struct->min_y)
		{

			/*Allocate the memory for the list of station numbers, and the x,z coordinates*/
			mt_2D_struct->mts = (long *)memory((char *)mt_2D_struct->mts,(mt_2D_struct->nstat_mt)+1,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->x = (double *)memory((char *)mt_2D_struct->x,(mt_2D_struct->nstat_mt)+1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->z = (double *)memory((char *)mt_2D_struct->z,(mt_2D_struct->nstat_mt)+1,sizeof(double),"ForwardModMT_2D");

			/*List of station numbers*/
			mt_2D_struct->mts[mt_2D_struct->nstat_mt] = data->mts[j];
			/*x,z-coordinates of the stations*/
			mt_2D_struct->x[mt_2D_struct->nstat_mt] = xs;
			mt_2D_struct->z[mt_2D_struct->nstat_mt] = zs;

			mt_2D_struct->nstat_mt++;
		}
	}		

	return(1);
}

/*-------------------------------------------------------------*/
/*Fill a slice with values by averaging the values in the rod vertical to the considered plane:*/
/*Parameter:     grid  := Grid structure */
/*       first_sample   := first sample of the layer*/
/*       direc_2D_mt   := direction of the slices (1=x-axis, 2=y-axis)*/
/*     nr_cells_2D_mt  := Number of cells per layer (specified in the input file)*/
/*              nx,nz  : Number of cells in x and z direction*/
/*               cube  := Cube that includes the parameters that should be averaged*/
/*				slice := Resistivities of the slices and indeces of the slices*/

#define cube(x,y,z) cube[((grid.ny)*(grid.nz))*(x) + (grid.nz)*(y) + (z)]

int CalcAveValues2DSlice(GRID_STRUCT grid, long first_sample, int direc_2D_mt, long nr_cells_2D_mt, long nx, long nz, double *slice, double *cube)
{
	long j,ix,iz;
	long index, index_air_water;
	long nyz2,ny2,nz2;
	double air_water;

	ny2 = grid.ny + 2*grid.nborder;
	nz2 = grid.nz + 2*grid.nborder;
	nyz2 = ny2*nz2;

	for(j=0;j<nx*nz;j++)
		slice[j]= 0.0;

	/*Loop over all cells*/
	for(ix=0;ix<nx;ix++)
		for(iz=0;iz<nz;iz++)
		{
			index = 0;
			index_air_water = 0;

			air_water = 0.0;

			/*Slice along x-direction*/
			if(direc_2D_mt != 2)
			{
				/*Calculate the average resistivity values in a cell of the slice*/
				for(j=first_sample;j< (first_sample + nr_cells_2D_mt) && j<grid.ny ;j++)
				{	
					/*Cells in the air/water will not be not considered in the averaging*/
					if(grid.border_index[nyz2*(ix+grid.nborder) + nz2*(j+grid.nborder) + (iz+grid.nborder)] != 0)
					{
						slice[ix*nz + iz] = cube(ix,j,iz) + slice[ix*nz + iz];
						index++;
					}
					else
					{
						air_water = cube(ix,j,iz)+ air_water; 
						index_air_water++;
					}
				}
			}
			/*Slice along y-direction*/
			else
			{
				/*Calculate the average resistivity values in a cell of the slice*/
				for(j=first_sample;j<(first_sample + nr_cells_2D_mt) && j<grid.nx ;j++)
				{
					/*Cells in the air/water will not be not considered in the averaging*/
					if(grid.border_index[nyz2*(j+grid.nborder) + nz2*(ix+grid.nborder) + (iz+grid.nborder)] != 0)
					{
						slice[ix*nz + iz] = cube(j,ix,iz) + slice[ix*nz + iz];

						index++;
					}
					else
					{
						air_water = cube(j,ix,iz) + air_water;
						index_air_water++;
					}
				}
			}

			/*Only if at least half part of the cells are NOT border cells, the cells will be averaged*/
			if(index >= index_air_water)
				slice[ix*nz + iz] = slice[ix*nz + iz]/index;
			/*Otherwise they will be filled by cells with air/water entries*/
			else
				slice[ix*nz + iz] = air_water/index_air_water;

		}

	return(1);
}

#undef cube


/*-------------------------------------------------------------*/
/*Determine the cutting points of the topography at the middel of the 2-D slice :*/
/*Parameter:     grid  := Grid structure */
/*       first_sample   := first sample of the layer*/
/*       direc_2D_mt   := direction of the slices (1=x-axis, 2=y-axis)*/
/*     nr_cells_2D_mt  := Number of cells per layer (specified in the input file)*/
/*                 nx   := Number of cells in x direction*/
/*				  *hx   := Thickness of the cells in x-direction*/
/*				topo	:= Topography structure*/
/*      *grid_points    := z-coordinate of the intersection of the topography with the border of the cells*/

int DetTopography2DSlice(GRID_STRUCT grid, long first_sample, int direc_2D_mt, long nr_cells_2d_mt, long nx, double *hx,TOPO_STRUCT topo, double *grid_points)
{
	long ix;
	double *xgrid_edges, *ygrid_edges, *zgrid_edges;
	double org_xy;

	/*Allocate the memory for the edges*/
	xgrid_edges = (double *)memory(NULL,(nx+1),sizeof(double),"DetTopography2DSlice");
	ygrid_edges = (double *)memory(NULL,(nx+1),sizeof(double),"DetTopography2DSlice");
	zgrid_edges = (double *)memory(NULL,(nx+1),sizeof(double),"DetTopography2DSlice");
	
	/*Slice along x-direction*/
	if(direc_2D_mt != 2)
	{
		/*Re-specify the thickness of the slice close to "upper" border*/
		if(grid.ny < (first_sample + nr_cells_2d_mt))
			nr_cells_2d_mt = grid.ny - first_sample;

		org_xy = grid.org[0] - (grid.h/2.0);

		/*Loop over all cells in x-direction*/
		for(ix=0;ix<(nx+1);ix++)
		{
			/*Determine the x-coordinate*/
			xgrid_edges[ix] = org_xy;
			/*Determine the y-coordinate*/
			ygrid_edges[ix] = grid.org[1] + grid.h*((double)first_sample + ((double)nr_cells_2d_mt - 1.0)/2.0);

			if(ix != nx)
				org_xy = org_xy + hx[ix];
		}
	}

	/*Slice along y-direction*/
	else
	{
		/*Re-specify the thickness of the slice close to "upper" border*/
		if(grid.nx < (first_sample + nr_cells_2d_mt))
			nr_cells_2d_mt = grid.nx - first_sample;

		org_xy = grid.org[1] - (grid.h/2.0);

		/*Loop over all cells in y-direction*/
		for(ix=0;ix<(nx+1);ix++)
		{
			/*Determine the y-coordinate*/
			ygrid_edges[ix] = org_xy;
			/*Determine the x-coordinate*/
			xgrid_edges[ix] = grid.org[0] + grid.h*((double)first_sample + ((double)nr_cells_2d_mt - 1.0)/2.0);

			if(ix != nx)
				org_xy = org_xy + hx[ix];
		}
	}

	/*Calculate the intersections*/
	TriangleInterpol(xgrid_edges,ygrid_edges,zgrid_edges,(nx+1),topo.x,topo.y,topo.z,topo.nr_of_topo_points,topo.index,topo.nr_of_triangles);

	/*Assign the intersections with the topography*/
	for(ix=0;ix<(nx+1);ix++)
		grid_points[ix] = zgrid_edges[ix];

	free(xgrid_edges);
	free(ygrid_edges);
	free(zgrid_edges);

	return(1);
}

/*-------------------------------------------------------------*/
/*Modify the cutting points close to the MT stations within the slice*/
/*		geometry	:= Geometry structure*/
/*	*grid_points    := z-coordinate of the intersection of the topography with the border of the cells*/
/*   direc_2d_mt    := direction of the slices (1=x-axis, 2=y-axis)*/
/*      nr_stat_mt   := Nr of MT stations*/
/*            *mts  : station indeces of the MT stations*/
/*           nx		:= Nr of samples in x/y-direction*/
/*           hx		:= Sample intervals in the x/y direction*/
/*          orgx    := Origin of the grid in x-direction*/

int DetTopographyMTStat2DSlice(GEOMETRY geo, double *grid_points, int direc_2D_mt, long nx, double *hx, double orgx, long nr_stat_mt, long *mts)
{

	long *index_cell;
	long j,i;
	double pos_x1, pos_x2, delta_x;

	/*Allocate the memory for the edges*/
	index_cell = (long *)memory(NULL,(nx+1),sizeof(long),"DetTopographyMTStat2DSlice");

	for(i=0;i<(nx+1);i++)
		index_cell[i] = -1;

	/*Loop over all MT stations in the 2D slice*/
	for(j=0;j<nr_stat_mt;j++)
	{
		/*If the station location has topography information at the same time*/
		if(geo.coor_info[(mts[j] - 1)] != 0)
		{
			/*Determine the direction of the 2-D slice*/
				/*Slice along the x-direction*/
			if(direc_2D_mt != 2)
			{
				pos_x1 = orgx;
				pos_x2 = orgx + hx[0];

				/*Loop over all samples in x-direction*/
				for(i=0;i<nx;i++)
				{
					/*Find the position of the station in the grid*/
					if(geo.x[(mts[j] - 1)] >= pos_x1 && geo.x[(mts[j] - 1)] <= pos_x2)
					{

						/*Check, if there is NOT an other station in the grid cell*/
						if(index_cell[i] == -1)
						{
							index_cell[i] = mts[j];
							/*Reset the intersection points*/
							grid_points[i] = geo.z[(mts[j] - 1)];
						}
						else
						{
							if(geo.x[index_cell[i] -1] > pos_x1 )
								delta_x = geo.x[index_cell[i] -1] - pos_x1;
							else
								delta_x = pos_x1 - geo.x[index_cell[i] -1];

							if(delta_x > (geo.x[(mts[j] - 1)] - pos_x1))
							{
								index_cell[i] = mts[j];
								/*Reset the intersection points*/
								grid_points[i] = geo.z[(mts[j] - 1)];
							}

						}

						if(index_cell[i+1] == -1)
						{
							index_cell[i+1] = mts[j];
							/*Reset the intersection points*/
							grid_points[i+1] = geo.z[(mts[j] - 1)];
						}
						else
						{
							if(geo.x[index_cell[i+1] -1] > pos_x2)
								delta_x = geo.x[index_cell[i+1] -1] - pos_x2;
							else
								delta_x = pos_x2 - geo.x[index_cell[i+1] -1];

							if(delta_x > (pos_x2 - geo.x[(mts[j] - 1)]))
							{
								index_cell[i+1] = mts[j];
								/*Reset the intersection points*/
								grid_points[i+1] = geo.z[(mts[j] - 1)];
							}
						}

						break;
					}

					if(i!=(nx-1))
					{
						pos_x1 = pos_x1 + hx[i];
						pos_x2 = pos_x2 + hx[i+1];
					}
				}
			}
				/*Slice along the y-direction*/
			else
			{
				pos_x1 = orgx;
				pos_x2 = orgx + hx[0];

				/*Loop over all samples in y-direction*/
				for(i=0;i<nx;i++)
				{
					/*Find the position of the station in the grid*/
					if(geo.y[(mts[j] - 1)] >= pos_x1 && geo.y[(mts[j] - 1)] <= pos_x2)
					{

						/*Check, if there is NOT an other station in the grid cell*/
						if(index_cell[i] == -1)
						{
							index_cell[i] = mts[j];
							/*Reset the intersection points*/
							grid_points[i] = geo.z[(mts[j] - 1)];
						}
						else
						{
							if(geo.y[index_cell[i] -1] > pos_x1)
								delta_x = geo.y[index_cell[i] -1] - pos_x1;
							else
								delta_x = pos_x1 - geo.y[index_cell[i] -1];

							if(delta_x > (geo.y[(mts[j] - 1)] - pos_x1))
							{
								index_cell[i] = mts[j];
								/*Reset the intersection points*/
								grid_points[i] = geo.z[(mts[j] - 1)];
							}
						}

						if(index_cell[i+1] == -1)
						{
							index_cell[i+1] = mts[j];
							/*Reset the intersection points*/
							grid_points[i+1] = geo.z[(mts[j] - 1)];
						}
						else
						{
							if(geo.y[index_cell[i+1] -1] > pos_x2)
								delta_x = geo.y[index_cell[i+1] -1] - pos_x2;
							else
								delta_x = pos_x2 - geo.y[index_cell[i+1] -1];

							if(delta_x > (pos_x2 - geo.y[(mts[j] - 1)]))
							{
								index_cell[i+1] = mts[j];
								/*Reset the intersection points*/
								grid_points[i+1] = geo.z[(mts[j] - 1)];
							}
						}
	
						break;
					}

					if(i!=(nx-1))
					{
						pos_x1 = pos_x1 + hx[i];
						pos_x2 = pos_x2 + hx[i+1];
					}
				}
			}
		}
	}

	free(index_cell);

	return(1);
}

/*-------------------------------------------------------------*/
/*Add the topography information to the refined grid; reset the */
/*resistivity values immediately above and below the topography*/
/*Parameter:       nx,nz   := Number of cells in the x/z-direction*/
/*				  *hz  := Thickness of the cells in z-direction*/
/*		*grid_points    := z-coordinate of the intersection of the topography with the border of the cells*/
/*                orgz  := Origin of the data cube (z-component; upper edge and NOT center of the cell)*/
/*		   *res_value   := resistivity values in the slice*/
/*   res_value_water_air := Resistivity value used from the input parameter file (only used, if the resistivities can not be related fro the neighboring cells)*/
/*         *index_topo  := Index specifying if the point belongs to the "air/water" (==0) or the inverted part (==1)*/
/*REMARK: DetTopography2DSlice have to be applied before to determine the *grid_points */

int AddTopography2DSlice(long nx, long nz, double *hz, double *grid_points, double orgz, double *res_value , float res_water_air, int *index_topo)
{
	long i,j,k,j_down, j_up ,j_new, r_index;
	double average_z, hz_l, hz_u, res_below_topo, res_above_topo;

	/*Loop over all cells in slice-direction*/
	for(i=0; i<nx; i++)
	{
		/*Identify, if topography exists for these points*/
		if(grid_points[i]!= -99999.9 || grid_points[i+1]!= -99999.9)
		{

			/*Average coordinate of intersection*/
			average_z = (grid_points[i] + grid_points[i+1])/2.0;

			/*Specify the lower and upper border of the cells*/
			hz_l = orgz;
			hz_u = orgz + hz[0];

			/*Loop over all depths*/
			for(j=0; j<nz; j++)
			{
				/*Find the cell that will be intersected by the topography*/
				if(hz_l <= average_z && hz_u >= average_z)
				{

					if((average_z - hz_l) < (hz_u - average_z))
						j_new = j;
					else
						j_new = j+1;

					/*Point belong to the air before*/
					if(index_topo[i*nz + j_new] == 0)
					{
						j_down = j_new; 

						r_index = 0;

						/*Make the cells ACTIVE below the topography*/
						while(index_topo[i*nz + j_down] == 0 && (j_down < nz))
						{
							index_topo[i*nz + j_down] = 1;
							j_down++;
							r_index++;
						}

						/*Fill the modified cells with the resistivties that are found before directly below the topography*/
						if(j_down != nz)
						{
							res_below_topo = res_value[i*nz + j_down];

							for(k=j_new;k<j_down;k++)
								res_value[i*nz + k] = res_below_topo;
						}
						else if(r_index != 0)
							printf("\nWARNING!!! There are NO active cells below the topography\n");

					}
					/*Point belong to the active region before*/
					else
					{
						j_up = j_new - 1;
						
						r_index = 0;

						/*Make the cells DE-ACTIVE above the topography*/
						while(index_topo[i*nz + j_up] != 0 && (j_up >= 0))
						{
							index_topo[i*nz + j_up] = 0;
							j_up--;
							r_index++;
						}

						/*Fill the DE-ACTIVATED cells with the water/air resistivity above the cell*/
						if(j_up != -1)
						{
							res_above_topo = res_value[i*nz + j_up];

							for(k=j_up+1;k<j_new;k++)
								res_value[i*nz + k] = res_above_topo;
						}
						/*Fill the DE-ACTIVATED cells with the water/air resistivity from the parameterfile*/
						/*Required because NO cells with water/air resistivity exists above the now DE-ACTIVATED CELLS*/
						else if(r_index != 0)
						{
							printf("\nWARNING!!! There are NO de-active cells above the topography!!\n");
							printf("The water/air resistivity values from the parameter-file are used:\n");
							printf("%f ohmm\n\n",res_water_air);

							for(k=j_up+1;k<j_new;k++)
								res_value[i*nz + k] = (double)res_water_air;
						}

					}

					goto next_nx;
				}
			
				if(j != nz)
				{
					hz_l = hz_l + hz[j];
					hz_u = hz_u + hz[j+1];
				}

			}

		}

		next_nx:;

	}

	return(1);
}

/*-------------------------------------------------------------*/
/*Determine the frequencies used for the 2D slice and adjust the 2D_mt structure such that it can be used for the 2D forward calculation*/
/*Parameter:     *data  := Data structure */
/*       2D_mt_struct  := Structure including 2D MT specific parameters*/

int AdjustFreq2Dslice(CALC_2D_MT *mt_2D_struct, DATA_STRUCT *data)
{

	long j,k,m,n;
	int freq_index, counter;

	mt_2D_struct->freq_ndata[0] = 0;

	/*Make a list of the used frequencies (+ the corresponding measurement indices) in the slice*/
	mt_2D_struct->nfreq = 0;

	for(j=0;j<mt_2D_struct->nstat_mt;j++)
		for(m=0;m<data->ndata_mt;m++)
		{
			if(mt_2D_struct->mts[j] == data->mno[m])
			{
				/*Loop over all frequencies*/
				for(k=0;k<data->nfreq_mt[m];k++)
				{
					freq_index = 0;

					/*Check if the frequency has been already found*/
					for(n=0;n<mt_2D_struct->nfreq;n++)
					{
						if(mt_2D_struct->freq[n] == data->freq_mt[m][k])
						{
							mt_2D_struct->freq_ndata[n]++;

							freq_index = 1;
							goto supercool;
						}
					}
						
					supercool:;
						
					/*Make a list of all frequencies in the slice*/
					if(freq_index == 0)
					{
						/*Reallocate the memory for the frequencies used for the stations in the slice*/
						mt_2D_struct->freq = (double *)memory((char *)mt_2D_struct->freq,(mt_2D_struct->nfreq)+1,sizeof(double),"ForwardModMT_2D");
						mt_2D_struct->freq_ndata = (long *)memory((char *)mt_2D_struct->freq_ndata,(mt_2D_struct->nfreq)+1,sizeof(long),"ForwardModMT_2D");
						mt_2D_struct->freq_data = (long **)memory((char *)mt_2D_struct->freq_data,(mt_2D_struct->nfreq)+1,sizeof(long *),"ForwardModMT_2D");
						mt_2D_struct->freq_stat = (long **)memory((char *)mt_2D_struct->freq_stat,(mt_2D_struct->nfreq)+1,sizeof(long *),"ForwardModMT_2D");

						/*Reallocate memory for the field components at the stations*/
						mt_2D_struct->Ex_r = (double **)memory((char *)mt_2D_struct->Ex_r,(mt_2D_struct->nfreq)+1,sizeof(double *),"ForwardModMT_2D");
						mt_2D_struct->Ex_i = (double **)memory((char *)mt_2D_struct->Ex_i,(mt_2D_struct->nfreq)+1,sizeof(double *),"ForwardModMT_2D");
						mt_2D_struct->Ey_r = (double **)memory((char *)mt_2D_struct->Ey_r,(mt_2D_struct->nfreq)+1,sizeof(double *),"ForwardModMT_2D");
						mt_2D_struct->Ey_i = (double **)memory((char *)mt_2D_struct->Ey_i,(mt_2D_struct->nfreq)+1,sizeof(double *),"ForwardModMT_2D");
						mt_2D_struct->Hx_r = (double **)memory((char *)mt_2D_struct->Hx_r,(mt_2D_struct->nfreq)+1,sizeof(double *),"ForwardModMT_2D");
						mt_2D_struct->Hx_i = (double **)memory((char *)mt_2D_struct->Hx_i,(mt_2D_struct->nfreq)+1,sizeof(double *),"ForwardModMT_2D");
						mt_2D_struct->Hy_r = (double **)memory((char *)mt_2D_struct->Hy_r,(mt_2D_struct->nfreq)+1,sizeof(double *),"ForwardModMT_2D");
						mt_2D_struct->Hy_i = (double **)memory((char *)mt_2D_struct->Hy_i,(mt_2D_struct->nfreq)+1,sizeof(double *),"ForwardModMT_2D");

						mt_2D_struct->freq[mt_2D_struct->nfreq] = data->freq_mt[m][k];
						mt_2D_struct->freq_ndata[mt_2D_struct->nfreq] = 1;
						mt_2D_struct->nfreq++;
					}
				}
			}
		}


	for(n=0;n<mt_2D_struct->nfreq;n++)
	{
		counter = 0;

		/*Allocate memory for the list of measurements for the different frequencies*/
		if(mt_2D_struct->freq_ndata[n] != 0)
		{
			mt_2D_struct->freq_data[n] = (long *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->freq_stat[n] = (long *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(long),"ForwardModMT_2D");

			mt_2D_struct->Ex_r[n] = (double *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Ex_i[n] = (double *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Ey_r[n] = (double *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Ey_i[n] = (double *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Hx_r[n] = (double *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Hx_i[n] = (double *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Hy_r[n] = (double *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Hy_i[n] = (double *)memory(NULL,mt_2D_struct->freq_ndata[n],sizeof(double),"ForwardModMT_2D");

		}
		else
		{
			mt_2D_struct->freq_data[n] = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");
			mt_2D_struct->freq_stat[n] = (long *)memory(NULL,1,sizeof(long),"ForwardModMT_2D");

			mt_2D_struct->Ex_r[n] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Ex_i[n] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Ey_r[n] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Ey_i[n] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Hx_r[n] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Hx_i[n] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Hy_r[n] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->Hy_i[n] = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");

		}

		for(j=0;j<mt_2D_struct->nstat_mt;j++)
			for(m=0;m<data->ndata_mt;m++)
			{
				if(mt_2D_struct->mts[j] == data->mno[m])
				{
					/*Loop over all frequencies*/
					for(k=0;k<data->nfreq_mt[m];k++)
					{
						if(mt_2D_struct->freq[n] == data->freq_mt[m][k])
						{
							mt_2D_struct->freq_data[n][counter] = m;
							mt_2D_struct->freq_stat[n][counter] = mt_2D_struct->mts[j];
							counter++;

							if(counter > mt_2D_struct->freq_ndata[n])
							{
								printf("An error occurs during making the 2D_MT structure!!\n ");
								exit(0);
							}
						}
					}
				}
			}
	}

	return(1);
}

/*-------------------------------------------------------------*/
/*Setup the global refined MT grid*/
/*- Global refinement (governed by the smallest resistivity in the grid and corresponding skin depth)*/
/*Parameter: grid_h = cell size of the original grid*/
/*           *mt_2D_struct = structure for the 2D MT grid*/
/*            min_res = Minimum resistivity in the grid*/
/*      *org_res_slice = Original resistivity slice*/
/*      *org_index_slice = Original indeces in the slice*/
/*            freq_index  = index specifying the used frequency */

int MakeRefinedGridMT(double grid_h, CALC_2D_MT *mt_2D_struct, double min_res, double *org_res_slice, double *org_index_slice, long freq_index)
{
	long k;
	double min_skin_depth; /*skin depth for the smallest resistivity in the grid*/
	double tmp_cell_size, *tmp_res_slice, *tmp_index_slice;


	/*Find factor for refinement of the grid*/
	min_skin_depth = (1/(2*PI))*sqrt((1.0E7*min_res)/mt_2D_struct->freq[freq_index]); /*minimum skin depth*/
	tmp_cell_size = min_skin_depth * GLOBAL_FACTOR_REFINEMENT;
		/*Scaling factor*/
	mt_2D_struct->cell_scal_factor[freq_index] = (int)ceil(grid_h/tmp_cell_size);
	/********************************/
	/*Refine the total grid*/

	/*Allocate the memory for the grid cell sizes*/
	mt_2D_struct->hx[freq_index] = (double *)memory((char *) mt_2D_struct->hx[freq_index],(mt_2D_struct->cell_scal_factor[freq_index])* (mt_2D_struct->nx[freq_index]),sizeof(double),"ForwardModMT_2D");
	mt_2D_struct->hz[freq_index] = (double *)memory((char *) mt_2D_struct->hz[freq_index],(mt_2D_struct->cell_scal_factor[freq_index])* (mt_2D_struct->nz[freq_index]),sizeof(double),"ForwardModMT_2D");

	/*Allocate the memory for a temporary grid*/
	tmp_res_slice = (double *)memory(NULL,(mt_2D_struct->cell_scal_factor[freq_index])*(mt_2D_struct->cell_scal_factor[freq_index])*(mt_2D_struct->nx[freq_index])*(mt_2D_struct->nz[freq_index]),sizeof(double),"ForwardModMT_2D");
	tmp_index_slice = (double *)memory(NULL,(mt_2D_struct->cell_scal_factor[freq_index])*(mt_2D_struct->cell_scal_factor[freq_index])*(mt_2D_struct->nx[freq_index])*(mt_2D_struct->nz[freq_index]),sizeof(double),"ForwardModMT_2D");

	/*Grid refinement*/
		/*resistivities+border indeces*/
	grid_refinement(&mt_2D_struct->nx[freq_index], &mt_2D_struct->nz[freq_index], mt_2D_struct->hx[freq_index], mt_2D_struct->hz[freq_index], grid_h, mt_2D_struct->cell_scal_factor[freq_index], org_res_slice, tmp_res_slice, org_index_slice, tmp_index_slice);

	mt_2D_struct->res_slice[freq_index] = (double *)memory((char *)mt_2D_struct->res_slice[freq_index],(mt_2D_struct->nx[freq_index])*(mt_2D_struct->nz[freq_index]),sizeof(double),"ForwardModMT_2D");
	mt_2D_struct->index_slice[freq_index] = (int *)memory((char *)mt_2D_struct->index_slice[freq_index],(mt_2D_struct->nx[freq_index])*(mt_2D_struct->nz[freq_index]),sizeof(int),"ForwardModMT_2D");
				
	for(k=0;k<((mt_2D_struct->nx[freq_index])*(mt_2D_struct->nz[freq_index]));k++)
	{
		mt_2D_struct->res_slice[freq_index][k] = tmp_res_slice[k];
		mt_2D_struct->index_slice[freq_index][k] = (int)tmp_index_slice[k];
	}

	free(tmp_res_slice);
	free(tmp_index_slice);

	return(0);
}


/*-------------------------------------------------------------*/
/*Refine the grid in a specified region*/
/*Parameter: *nx= number of cells in x-direction (will be modified)*/
/*           *nz= number of cells in z-direction (will be modified)*/
/*           *hx= List of cell intervals in x-direction (will be modified)*/
/*           *hz= List of cell intervals in z-direction (will be modified)*/
/*           org_cell_size = original cell size*/
/*           cell_scal_factor = scaling factor*/
/*           *org_slice1/2 = values in the original slice*/
/*           *new_slice1/2 = values in the new slice*/

int grid_refinement(long *nx, long *nz, double *hx, double *hz, double org_cell_size, long cell_scal_factor, double *org_slice1, double *new_slice1, double *org_slice2, double *new_slice2)
{
	long k,m,ix,iz;
	double tmp_grid_h;

	tmp_grid_h = org_cell_size/cell_scal_factor;

			/*Modify the cell sizes*/
			for(k=0;k<(cell_scal_factor)*(*nx);k++)
				hx[k] = tmp_grid_h;
			for(k=0;k<(cell_scal_factor)*(*nz);k++)
				hz[k] = tmp_grid_h;

						
			for(ix=0;ix<(*nx);ix++)
				for(iz=0;iz<(*nz);iz++)
					for(k=0;k<cell_scal_factor;k++)
						for(m=0;m<cell_scal_factor;m++)
						{
							/*Fill the refined grid with parameters*/
							new_slice1[((ix*cell_scal_factor)+k)*((*nz)*cell_scal_factor) + ((iz*cell_scal_factor)+m)] = org_slice1[ix*(*nz) + iz];
							new_slice2[((ix*cell_scal_factor)+k)*((*nz)*cell_scal_factor) + ((iz*cell_scal_factor)+m)] = org_slice2[ix*(*nz) + iz];
						}

			/*Reset the number of cells in x and z diretion*/
			(*nx) = (cell_scal_factor)*(*nx);
			(*nz) = (cell_scal_factor)*(*nz);


	return(0);
}

/*-------------------------------------------------------------*/
/*Additional local refinement close to the stations; governed by the smallest resistivity
/*  in the grid and the 1D skin depths at a specified number of stations)*/
/*Parameter:        grid_h = cell size of the original grid*/
/*           *mt_2D_struct = structure for the 2D MT grid*/
/*                min_res  = Minimum resistivity in the grid*/
/*               nr_layer  = layer of the 2-D calculation*/
/*            freq_index  = index specifying the used frequency */

int MakeLocalRefinedGridMT(double grid_h, CALC_2D_MT *mt_2D_struct, double min_res, long nr_layer, long freq_index)
{
	int *tmp_index_slice;
	long k,n,ix,iz;
	long local_cell_scal_factor; /*Scaling factor of the grid refinement*/
	long min_iz, max_iz;
	long station_nz, station_ix, station_iz;
	long nr_station,*used_station;  /*List of the used stations, where the skin depth as determined*/
	long factor;
	long tmp_nr_hz;
	double min_z, max_z; /*Minimum and maximum depth of the stations*/
	double ang_freq; /*angular frequency*/
	double tmp_grid_h, *tmp_hz, *tmp_res_slice;
	double app_res_max; /*maximum apparent resistivity at the considered stations*/
	double app_res; /*apparent resistivity at the stations*/
	double pos_x, pos_z;
	double station_dz;
	double *station_hz, *station_res;
	double imped_r, imped_i; /*Impedances obtained from the 1D calculations*/
	double skin_depth_station; /*skin depth at the stations*/
	double min_skin_depth;     /*minimum skin depths in the grid*/
	double tmp_cell_size;

	/*****************************************************************/
	/*Determine the z-positions of the stations*/
	min_z = mt_2D_struct->z[0];
	max_z = mt_2D_struct->z[0];

	for(k=0;k<mt_2D_struct->nstat_mt;k++)
	{
		/*Find the minimum and maximum z*/
		if(mt_2D_struct->z[k] < min_z);
			min_z = mt_2D_struct->z[k];
		if(mt_2D_struct->z[k] > max_z);
			max_z = mt_2D_struct->z[k];
	}

	/*****************************************************************/
	/*Determine the maximum skin depth for a number of mt-station locations:*/

		/*Determine the position of the considered stations*/
	if(mt_2D_struct->nstat_mt == 0)
	{
		printf("NO MT-station exists in the %d layer\n",nr_layer);
		printf("But calculations are performed!!!\n");
		exit(0);
	}
	else if(NR_OF_1D_CALCULATION >= mt_2D_struct->nstat_mt)
	{
		used_station = (long *)memory(NULL, mt_2D_struct->nstat_mt, sizeof(long), "ForwardModMT_2D");

		for(k=0;k<mt_2D_struct->nstat_mt;k++)
			used_station[k] = k;

			nr_station = mt_2D_struct->nstat_mt;
	}
	else
	{
		used_station = (long *)memory(NULL, NR_OF_1D_CALCULATION, sizeof(long), "ForwardModMT_2D");

		for(k=0;k<NR_OF_1D_CALCULATION;k++)
			used_station[k] = (long) ceil((double)k*(mt_2D_struct->nstat_mt/(NR_OF_1D_CALCULATION - 1)));

		nr_station = NR_OF_1D_CALCULATION;
	}

	  /*Calculate the skin depths*/
	ang_freq = 2*PI*mt_2D_struct->freq[freq_index];

	app_res_max = 0.0;

	for(k=0;k<nr_station;k++)
	{
			
		/*Find the z-location of the station*/
		pos_z = mt_2D_struct->org[1][freq_index];
		pos_x = mt_2D_struct->org[0][freq_index];
		station_dz = 0.0;
		station_iz = 0;
		station_ix = 0;

		for(n=0;n<(mt_2D_struct->nz[freq_index] - 1);n++)
		{
			if(pos_z <= mt_2D_struct->z[k] && pos_z + mt_2D_struct->hz[freq_index][n] >= mt_2D_struct->z[k])
			{
				station_dz = (pos_z + mt_2D_struct->hz[freq_index][n]) - mt_2D_struct->z[k];
				station_iz = n;
				break;
			}
			pos_z = pos_z + mt_2D_struct->hz[freq_index][n];	
		}

		for(n=0;n<(mt_2D_struct->nx[freq_index] - 1);n++)
		{
			if(pos_x <= mt_2D_struct->x[k] && pos_x + mt_2D_struct->hx[freq_index][n] >= mt_2D_struct->x[k])
			{
				station_ix = n;
				break;
			}
		   pos_x = pos_x + mt_2D_struct->hx[freq_index][n];	
		}


		/*Determine the cells below the station*/
		station_nz = mt_2D_struct->nz[freq_index] - station_iz;

		station_hz = (double *)memory(NULL,station_nz,sizeof(double),"ForwardModMT_2D");
		station_res = (double *)memory(NULL,station_nz,sizeof(double),"ForwardModMT_2D");

		station_hz[0] = station_dz;
		station_res[0] = mt_2D_struct->res_slice[freq_index][station_ix*mt_2D_struct->nz[freq_index] + station_iz];

		for(n=1;n<station_nz;n++)
		{
			station_hz[n] = mt_2D_struct->hz[freq_index][n+station_iz];
			station_res[n] = mt_2D_struct->res_slice[freq_index][station_ix*mt_2D_struct->nz[freq_index] +n+station_iz];
		}

		imped_r = 0;
		imped_i = 0;

		/*1D calculation*/
		MT_1D_CALC(ang_freq, &imped_r, &imped_i, station_hz, station_nz, station_res, 0.0);

		/*Determine the apparent resistivity*/
		app_res = ((imped_r)*(imped_r) + (imped_i)*(imped_i))*(1/(ang_freq*(M_KONST)));

		if(app_res > app_res_max)
			app_res_max = app_res;

		free(station_hz);
		free(station_res);

	}

	free(used_station);

		/*calculate the maximum skin depth at these stations*/
	skin_depth_station = (1/(2*PI))*sqrt((1.0E7*app_res_max)/mt_2D_struct->freq[freq_index]);
	skin_depth_station = LOCAL_AREA_FACTOR * skin_depth_station;

	/*****************************************************************/
	/*Determine the cells that should be refined*/
	min_z = min_z - skin_depth_station;
	max_z = max_z + skin_depth_station;


	pos_z = mt_2D_struct->org[1][freq_index];
	min_iz = 0;

		/*Loop over all depths*/
	for(n=0;n< mt_2D_struct->nz[freq_index];n++)
	{
		if(pos_z >= min_z && n%mt_2D_struct->cell_scal_factor[freq_index] == 0)
		{
			/*First cell that will be refined*/
			min_iz = n;
			break;
		}

		pos_z = pos_z + mt_2D_struct->hz[freq_index][n];	
	}

	max_iz = min_iz;

	for(n=min_iz;n<mt_2D_struct->nz[freq_index];n++)
	{
		if(pos_z >= max_z && n%mt_2D_struct->cell_scal_factor[freq_index] == 0)
		{
			/*Last cell that will be refined*/
			max_iz = n;
			goto hier_weiter;
		}

		pos_z = pos_z + mt_2D_struct->hz[freq_index][n];
	}

	max_iz = mt_2D_struct->nz[freq_index];

	hier_weiter:;

	/*****************************************************************/
	/*Make the local refinement of the cells*/
	if(LOCAL_FACTOR_REFINEMENT < GLOBAL_FACTOR_REFINEMENT)
	{
		min_skin_depth = (1/(2*PI))*sqrt((1.0E7*min_res)/mt_2D_struct->freq[freq_index]); /*minimum skin depth*/
		tmp_cell_size = min_skin_depth * LOCAL_FACTOR_REFINEMENT;
			/*Scaling factor*/
		local_cell_scal_factor = (int)ceil(grid_h/tmp_cell_size);
			/*cell sizes*/
		tmp_grid_h = grid_h/local_cell_scal_factor;

		/*Remark:After the first iteration the loop will be used to increase gradually the cell size from the local to the global system*/
		while(local_cell_scal_factor > mt_2D_struct->cell_scal_factor[freq_index])
		{

			/*Number of cells in z-direction*/					
			factor = (max_iz - min_iz)/mt_2D_struct->cell_scal_factor[freq_index];
			factor = local_cell_scal_factor*factor;

			tmp_nr_hz = min_iz + factor + (mt_2D_struct->nz[freq_index] - max_iz);

			/*Modify the cell sizes*/
			tmp_hz = (double *)memory(NULL,tmp_nr_hz,sizeof(double),"ForwardModMT_2D");

			for(n=0;n<min_iz;n++)
				tmp_hz[n] = mt_2D_struct->hz[freq_index][n];
			for(n=0;n<factor;n++)
				tmp_hz[n + min_iz] = tmp_grid_h;
			for(n=0;n<(mt_2D_struct->nz[freq_index] - max_iz);n++)
				tmp_hz[n + min_iz + factor] = mt_2D_struct->hz[freq_index][max_iz + n];

			mt_2D_struct->hz[freq_index] = (double *)memory((char *) mt_2D_struct->hz[freq_index], tmp_nr_hz, sizeof(double),"ForwardModMT_2D");

			for(n=0;n<tmp_nr_hz;n++)
				mt_2D_struct->hz[freq_index][n] = tmp_hz[n];

			/*Refill the cells*/
			if(tmp_nr_hz != 0)
			{
				tmp_res_slice = (double *)memory(NULL,(mt_2D_struct->nx[freq_index]*tmp_nr_hz),sizeof(double),"ForwardModMT_2D"); 
				tmp_index_slice = (int *)memory(NULL,(mt_2D_struct->nx[freq_index]*tmp_nr_hz),sizeof(int),"ForwardModMT_2D"); 
			}
			else
			{
				tmp_res_slice = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D"); 
				tmp_index_slice = (int *)memory(NULL,1,sizeof(int),"ForwardModMT_2D"); 
			}

			for(ix=0;ix<mt_2D_struct->nx[freq_index];ix++)
				for(iz=0;iz<min_iz;iz++)
				{
					tmp_res_slice[ix*tmp_nr_hz + iz] = mt_2D_struct->res_slice[freq_index][ix*(mt_2D_struct->nz[freq_index]) + iz];
					tmp_index_slice[ix*tmp_nr_hz + iz] = mt_2D_struct->index_slice[freq_index][ix*(mt_2D_struct->nz[freq_index]) + iz];
				}

			for(ix=0;ix<mt_2D_struct->nx[freq_index];ix++)
				for(iz=0;iz<(max_iz- min_iz)/mt_2D_struct->cell_scal_factor[freq_index];iz++)
					for(n=0;n<local_cell_scal_factor;n++)
					{
						/*Fill the refined grid with parameters*/
						tmp_res_slice[ix*tmp_nr_hz + min_iz + ((iz*local_cell_scal_factor)+n)] = mt_2D_struct->res_slice[freq_index][ix*(mt_2D_struct->nz[freq_index]) + (mt_2D_struct->cell_scal_factor[freq_index])*iz + min_iz];
						tmp_index_slice[ix*tmp_nr_hz + min_iz + ((iz*local_cell_scal_factor)+n)] = mt_2D_struct->index_slice[freq_index][ix*(mt_2D_struct->nz[freq_index]) + (mt_2D_struct->cell_scal_factor[freq_index])*iz + min_iz];
					}

			for(ix=0;ix<mt_2D_struct->nx[freq_index];ix++)
				for(iz=0;iz<(mt_2D_struct->nz[freq_index] - max_iz);iz++)
				{
					tmp_res_slice[ix*tmp_nr_hz + min_iz + factor + iz] = mt_2D_struct->res_slice[freq_index][ix*(mt_2D_struct->nz[freq_index]) + iz + max_iz];
					tmp_index_slice[ix*tmp_nr_hz + min_iz + factor + iz] = mt_2D_struct->index_slice[freq_index][ix*(mt_2D_struct->nz[freq_index]) + iz + max_iz];
				}	

			mt_2D_struct->res_slice[freq_index] = (double *)memory((char *) mt_2D_struct->res_slice[freq_index],(mt_2D_struct->nx[freq_index]*tmp_nr_hz),sizeof(double),"ForwardModMT_2D");
			mt_2D_struct->index_slice[freq_index] = (int *)memory((char *) mt_2D_struct->index_slice[freq_index],(mt_2D_struct->nx[freq_index]*tmp_nr_hz),sizeof(int),"ForwardModMT_2D");

			for(n=0;n<(mt_2D_struct->nx[freq_index]*tmp_nr_hz);n++)
			{
				mt_2D_struct->res_slice[freq_index][n] = tmp_res_slice[n];
				mt_2D_struct->index_slice[freq_index][n] = tmp_index_slice[n];
			}

			/*modify the number of cells in z-direction*/
			mt_2D_struct->nz[freq_index] = tmp_nr_hz;

			free(tmp_res_slice);
			free(tmp_index_slice);
			free(tmp_hz);

			/************/
			/*Set the parameters for the next iteration:*/
			/************/
			min_iz = min_iz + factor;
			max_iz = min_iz + mt_2D_struct->cell_scal_factor[freq_index];

			local_cell_scal_factor--;
			tmp_grid_h = grid_h/local_cell_scal_factor;

		}

	}

	return(0);
}

/*-------------------------------------------------------------*/
/* Adding cells at the left and right boundary (governed by the largest 1D skin depth at the boundaries)*/
/* and at the basement (governed by the resistivity and corresponding skin depth for the deepest cell of the original grid)*/
/*Parameter: *mt_2D_struct = structure for the 2D MT grid*/
/*            max_res_basement = max. resistivity at the basement*/
/*             freq_index  = index specifying the used frequency */
/*            max_res_left, max_res_right = Maximum resistivity at the left and right side*/
/*    *nx_right_border, *nz_left_border = Nr of cells at the left and right border*/
/*                          *nz_basement= Nr of cells at the basement*/

int AddCellsAtBoundaryMT(CALC_2D_MT *mt_2D_struct, double max_res_basement, double max_res_left, double max_res_right, long freq_index, long *nx_left_border, long *nx_right_border, long *nz_basement)
{
	int *tmp_index_slice;
	long ix,iz,k;
	double skin_depth_basement; /*Ski depth at the basement*/
	double skin_depth_left, skin_depth_right; /*skin depth at the left and the right side*/
	double sum_z;
	double sum_x_left, sum_x_right;
	double *tmp_hx, *tmp_left_hx, *tmp_right_hx, *tmp_res_slice;

	/*******************************************************************************************/
	/*Specify the grid parameters at the basement*/
	skin_depth_basement = (1/(2*PI))*sqrt((1.0E7*max_res_basement)/mt_2D_struct->freq[freq_index]);	/*Max skin depth at the basement*/
	skin_depth_basement = skin_depth_basement * FACTOR_BASEMENT;							/*Size of the extension of the basement*/

	sum_z = 0.0;

	while(skin_depth_basement > sum_z)
	{
		/*Determine the thickness of the cells in the basement*/
		mt_2D_struct->hz[freq_index] = (double *)memory((char *) mt_2D_struct->hz[freq_index],(mt_2D_struct->nz[freq_index] + (*nz_basement) + 2),sizeof(double),"ForwardModMT_2D");
		mt_2D_struct->hz[freq_index][mt_2D_struct->nz[freq_index] + (*nz_basement)] = FACTOR_SCAL_BAS * mt_2D_struct->hz[freq_index][mt_2D_struct->nz[freq_index] -1 + (*nz_basement)];

		sum_z = mt_2D_struct->hz[freq_index][mt_2D_struct->nz[freq_index] + (*nz_basement)] + sum_z;

		(*nz_basement)++; /*Nr of cells in the basement*/

	}

	mt_2D_struct->hz[freq_index][mt_2D_struct->nz[freq_index] + (*nz_basement)] = 0.0; /*<- For some unknown reason the TM-Mode calculation requires this*/

	/*******************************************************************************************/
	/*Specify the grid parameters at the left and right border*/
	skin_depth_left = (1/(2*PI))*sqrt((1.0E7*max_res_left)/mt_2D_struct->freq[freq_index]);	/*Max. skin depth at the left side*/
	skin_depth_right= (1/(2*PI))*sqrt((1.0E7*max_res_right)/mt_2D_struct->freq[freq_index]);	/*Max. skin depth at the right side*/
	skin_depth_left = skin_depth_left * FACTOR_BORDER_LEFT_RIGHT;				/*Size of the extension at the left side*/
	skin_depth_right = skin_depth_right * FACTOR_BORDER_LEFT_RIGHT;				/*Size of the extension at the right side*/

	/**************/
	/*left*/
		/*Allocate temporary memory for the cells at the left border*/
	tmp_left_hx = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");

		/*border cell closest to the grid*/
	tmp_left_hx[0] = FACTOR_SCAL_BORDER * mt_2D_struct->hx[freq_index][0];
	sum_x_left = tmp_left_hx[0];
	(*nx_left_border) = 1;
			
	while(skin_depth_left > sum_x_left)
	{
		/*Determine the thickness of the cells at the left border*/
		tmp_left_hx = (double *)memory((char *) tmp_left_hx, (*nx_left_border)+1 ,sizeof(double),"ForwardModMT_2D");
		tmp_left_hx[(*nx_left_border)] = FACTOR_SCAL_BORDER * tmp_left_hx[(*nx_left_border) - 1];

		sum_x_left = tmp_left_hx[(*nx_left_border)] + sum_x_left;

		(*nx_left_border)++; /*Nr of cells at the left border*/
	}

	/*Move the origin*/
	mt_2D_struct->org[0][freq_index] =  mt_2D_struct->org[0][freq_index] - sum_x_left; 

	/**************/
	/*right*/
		/*Allocate temporary memory for the cells at the right border*/
	tmp_right_hx = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");

		/*border cell closest to the grid*/
	tmp_right_hx[0] = FACTOR_SCAL_BORDER * mt_2D_struct->hx[freq_index][mt_2D_struct->nx[freq_index] - 1];
	sum_x_right = tmp_right_hx[0];
	(*nx_right_border) = 1;
			
	while(skin_depth_right > sum_x_right)
	{
		/*Determine the thickness of the cells at the right border*/
		tmp_right_hx = (double *)memory((char *) tmp_right_hx, (*nx_right_border)+1 ,sizeof(double),"ForwardModMT_2D");
		tmp_right_hx[(*nx_right_border)] = FACTOR_SCAL_BORDER * tmp_right_hx[(*nx_right_border) - 1];

		sum_x_right = tmp_right_hx[(*nx_right_border)] + sum_x_right;

		(*nx_right_border)++; /*Nr of cells at the right border*/
	}

	/**************/
	/*Adjust the array including the cell sizes in x-direction*/
	tmp_hx = (double *)memory(NULL,mt_2D_struct->nx[freq_index], sizeof(double),"ForwardModMT_2D");
	for(k=0;k<mt_2D_struct->nx[freq_index];k++)
		tmp_hx[k] = mt_2D_struct->hx[freq_index][k];

	/*Re-allocate the memory for the grid cell sizes*/
	mt_2D_struct->hx[freq_index] = (double *)memory((char *) mt_2D_struct->hx[freq_index], ((*nx_left_border) + mt_2D_struct->nx[freq_index] + (*nx_right_border) + 1), sizeof(double),"ForwardModMT_2D");
	/*and re-fill*/
	for(k=0;k<(*nx_left_border); k++)				/*left border*/
		mt_2D_struct->hx[freq_index][(*nx_left_border) - 1 - k] = tmp_left_hx[k];
	for(k=0;k<mt_2D_struct->nx[freq_index]; k++)	/*central area*/
		mt_2D_struct->hx[freq_index][(*nx_left_border) + k] = tmp_hx[k];
	for(k=0;k<(*nx_right_border); k++)			/*right border*/
		mt_2D_struct->hx[freq_index][(*nx_left_border) + mt_2D_struct->nx[freq_index] + k] = tmp_right_hx[k];
			
	mt_2D_struct->hx[freq_index][(*nx_left_border) + mt_2D_struct->nx[freq_index] + (*nx_right_border)] = 0.0; /*<- For some unknown reason the TM-Mode calculation requires this*/

	free(tmp_left_hx);
	free(tmp_right_hx);
	free(tmp_hx);
			
	/*******************************************************************************************/
	/*Make the final grid file for the TM-mode*/
	/*Allocate the memory for a temporary grid*/
	tmp_res_slice = (double *)memory(NULL,(mt_2D_struct->nx[freq_index])*(mt_2D_struct->nz[freq_index]),sizeof(double),"ForwardModMT_2D");
	tmp_index_slice = (int *)memory(NULL,(mt_2D_struct->nx[freq_index])*(mt_2D_struct->nz[freq_index]),sizeof(int),"ForwardModMT_2D");

		
	/*Loop over all cells*/
	for(k=0;k<(mt_2D_struct->nx[freq_index])*(mt_2D_struct->nz[freq_index]);k++)
	{
		tmp_res_slice[k] = mt_2D_struct->res_slice[freq_index][k];	
		tmp_index_slice[k] = mt_2D_struct->index_slice[freq_index][k];	
	}

	mt_2D_struct->res_slice[freq_index] = (double *)memory((char *)mt_2D_struct->res_slice[freq_index],((*nx_left_border) + mt_2D_struct->nx[freq_index] + (*nx_right_border))*((*nz_basement) + mt_2D_struct->nz[freq_index]),sizeof(double),"ForwardModMT_2D");
	mt_2D_struct->index_slice[freq_index] = (int *)memory((char *)mt_2D_struct->index_slice[freq_index],((*nx_left_border) + mt_2D_struct->nx[freq_index] + (*nx_right_border))*((*nz_basement) + mt_2D_struct->nz[freq_index]),sizeof(int),"ForwardModMT_2D");


	/*Refill the structure*/
		/*central area*/
	for(ix=0;ix<mt_2D_struct->nx[freq_index];ix++)
		for(iz=0;iz<mt_2D_struct->nz[freq_index];iz++)
		{
			mt_2D_struct->res_slice[freq_index][(ix + (*nx_left_border))*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz] = tmp_res_slice[ix*mt_2D_struct->nz[freq_index] + iz];
			mt_2D_struct->index_slice[freq_index][(ix + (*nx_left_border))*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz] = tmp_index_slice[ix*mt_2D_struct->nz[freq_index] + iz];
		}

		/*basement*/
	for(ix=0;ix<mt_2D_struct->nx[freq_index];ix++)
		for(iz=0;iz<(*nz_basement);iz++)
		{
			mt_2D_struct->res_slice[freq_index][(ix + (*nx_left_border))*((*nz_basement) + mt_2D_struct->nz[freq_index]) + (iz + mt_2D_struct->nz[freq_index])] = tmp_res_slice[ix*mt_2D_struct->nz[freq_index] + mt_2D_struct->nz[freq_index] -1];
			mt_2D_struct->index_slice[freq_index][(ix + (*nx_left_border))*((*nz_basement) + mt_2D_struct->nz[freq_index]) + (iz + mt_2D_struct->nz[freq_index])] = tmp_index_slice[ix*mt_2D_struct->nz[freq_index] + mt_2D_struct->nz[freq_index] -1];
		}

		/*left border*/
	for(ix=0;ix<(*nx_left_border);ix++)
		for(iz=0;iz<(mt_2D_struct->nz[freq_index] + (*nz_basement));iz++)
		{
			mt_2D_struct->res_slice[freq_index][ix*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz] = mt_2D_struct->res_slice[freq_index][((*nx_left_border))*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz]; 
			mt_2D_struct->index_slice[freq_index][ix*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz] = mt_2D_struct->index_slice[freq_index][((*nx_left_border))*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz]; 
		}

		/*right border*/
	for(ix=0;ix<(*nx_right_border);ix++)
		for(iz=0;iz<(mt_2D_struct->nz[freq_index] + (*nz_basement));iz++)
		{
			mt_2D_struct->res_slice[freq_index][(ix + (*nx_left_border) + mt_2D_struct->nx[freq_index])*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz] = mt_2D_struct->res_slice[freq_index][((*nx_left_border) + mt_2D_struct->nx[freq_index] -1)*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz]; 
			mt_2D_struct->index_slice[freq_index][(ix + (*nx_left_border) + mt_2D_struct->nx[freq_index])*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz] = mt_2D_struct->index_slice[freq_index][((*nx_left_border) + mt_2D_struct->nx[freq_index] -1)*((*nz_basement) + mt_2D_struct->nz[freq_index]) + iz]; 
		}


	mt_2D_struct->nz[freq_index] = (*nz_basement) + mt_2D_struct->nz[freq_index];			/*Extent grid with the cells at the basement*/
	mt_2D_struct->nx[freq_index] = (*nx_left_border) + mt_2D_struct->nx[freq_index] + (*nx_right_border);	/*Extent grid with the cells at the left and right border*/		

	free(tmp_res_slice);
	free(tmp_index_slice);

	return(0);
}

/*-------------------------------------------------------------*/
/*Add cells in the air for the TE-Mode (governed by the skin depth*/
/*of the 1D calculation)*/
/*Parameter:	mt_2D_struct = structure for the 2D MT slices*/
/*             freq_index  = index specifying the used frequency */
/*             nz_air,nz_ionos = Number of cells in the air and in the ionosphere*/

int AddCellsAtBoundaryAirMT(CALC_2D_MT *mt_2D_struct, long freq_index, long *nz_air, long *nz_ionos)
{
	long m,k;
	long *pos_of_1d_calculations;
	double ang_freq; /*angular frequency*/
	double app_res_max, app_res; /*Apparent resistivity*/ 
	double imped_r, imped_i, *resist;
	double skin_depth_air; /*Skin depth required for the air*/
	double *tmp_air_hz, *tmp_hz, sum_z_air;

	/*Specify the grid parameters in the air (required for the TE-Mode)*/

	/*Determine the max. apparent resistivity by means of 1-D MT calculations*/
		/*Number of positions, for which the apparent resistivity will be calculated*/
	if(NR_OF_1D_CALCULATION < 2)
		pos_of_1d_calculations = (long *)memory(NULL,2, sizeof(long),"ForwardModMT_2D");
	else
		pos_of_1d_calculations = (long *)memory(NULL,NR_OF_1D_CALCULATION, sizeof(long),"ForwardModMT_2D");

	/*Determine the positions (cell number) in the grid where the 1-D calculation will be performed*/
	pos_of_1d_calculations[0] = 0;												/*left border*/
	pos_of_1d_calculations[NR_OF_1D_CALCULATION -1] = mt_2D_struct->nx[freq_index] - 1;		/*right border*/

	for(k=1;k<NR_OF_1D_CALCULATION-1;k++)
		pos_of_1d_calculations[k] = (long)(((double)mt_2D_struct->nx[freq_index]/(NR_OF_1D_CALCULATION - 1))*k); /*other positions*/


	/******************************************************************************/
	/*Perform the 1-D calculations*/
	ang_freq = 2*PI*mt_2D_struct->freq[freq_index];

	app_res_max = 0.0;

	for(k=0;k<NR_OF_1D_CALCULATION;k++)
	{
		imped_r = 0;
		imped_i = 0;

		resist = (double *)memory(NULL,mt_2D_struct->nz[freq_index],sizeof(double),"ForwardModMT_2D");

		for(m=0;m<mt_2D_struct->nz[freq_index];m++)
			resist[m] = mt_2D_struct->res_slice[freq_index][pos_of_1d_calculations[k]*mt_2D_struct->nz[freq_index] + m];

		/*MT 1D calculation*/
		MT_1D_CALC(ang_freq, &imped_r, &imped_i, mt_2D_struct->hz[freq_index], mt_2D_struct->nz[freq_index], resist, 0.0);

		/*Determine the apparent resistivity*/
		app_res = ((imped_r)*(imped_r) + (imped_i)*(imped_i))*(1/(ang_freq*(M_KONST)));

		if(app_res > app_res_max)
			app_res_max = app_res;

		free(resist);
	}

	free(pos_of_1d_calculations);

	/******************************************************************************/
	/*Specyfing the grid parameters in the atmosphere*/
	skin_depth_air = (1/(2*PI))*sqrt((1.0E7*app_res_max)/mt_2D_struct->freq[freq_index]);	/*Skin depth related from the apparent resistivity*/
	skin_depth_air = skin_depth_air * FACTOR_AIR;									/*Size of the extension in the air*/

	/*Allocate temporary memory for the cells in the air*/
	tmp_air_hz = (double *)memory(NULL,1,sizeof(double),"ForwardModMT_2D");

	/*cell in the air closest to the usual grid*/
	tmp_air_hz[0] = FACTOR_SCAL_AIR * mt_2D_struct->hz[freq_index][0];
	sum_z_air = tmp_air_hz[0];
	(*nz_air) = 1;

	while(skin_depth_air > sum_z_air || MIN_SIZE_AIR > sum_z_air)
	{
	/*Determine the thickness of the cells in the air*/
		tmp_air_hz = (double *)memory((char *) tmp_air_hz, (*nz_air)+1 ,sizeof(double),"ForwardModMT_2D");
		tmp_air_hz[(*nz_air)] = FACTOR_SCAL_AIR * tmp_air_hz[(*nz_air) - 1];

		sum_z_air = tmp_air_hz[(*nz_air)] + sum_z_air;

		(*nz_air)++;
	}

	/*Add also ONE layer Ionosphere (following the 2-D code from Tarits)*/
	tmp_air_hz = (double *)memory((char *) tmp_air_hz, (*nz_air)+1 ,sizeof(double),"ForwardModMT_2D");
	tmp_air_hz[(*nz_air)] = DZ_IONOS;
	sum_z_air = tmp_air_hz[(*nz_air)] + sum_z_air;

	(*nz_ionos) = 1; /*Nr of cells in the ionosphere*/
			
	/******************************************************************************/

	/*Adjust the array including the cell sizes of the air in z-direction*/
	tmp_hz = (double *)memory(NULL,mt_2D_struct->nz[freq_index], sizeof(double),"ForwardModMT_2D");
	for(k=0;k<mt_2D_struct->nz[freq_index];k++)
		tmp_hz[k] = mt_2D_struct->hz[freq_index][k];

	/*Re-allocate the memory for the grid cell sizes*/
	mt_2D_struct->hz[freq_index] = (double *)memory((char *) mt_2D_struct->hz[freq_index], ((*nz_air) + (*nz_ionos) + mt_2D_struct->nz[freq_index]), sizeof(double),"ForwardModMT_2D");
	/*and re-fill*/
	for(k=0;k<((*nz_air) + (*nz_ionos)); k++)				/*air + ionosphere*/
		mt_2D_struct->hz[freq_index][((*nz_air) + (*nz_ionos)) - 1 - k] = tmp_air_hz[k];
	for(k=0;k<mt_2D_struct->nz[freq_index]; k++)	/*central area*/
		mt_2D_struct->hz[freq_index][((*nz_air) + (*nz_ionos) )+ k] = tmp_hz[k];

	mt_2D_struct->nz[freq_index] = (*nz_air) + (*nz_ionos) + mt_2D_struct->nz[freq_index];			/*Extent grid with the cells in the air*/

	free(tmp_hz);
	free(tmp_air_hz);

	return(0);
}


/*-------------------------------------------------------------*/
/*Determine the impedances at the positions of the stations*/
/*Parameter:     mt_structure  := MT 2D structure */
/*                       data  := data structure   */
/*             H_r,H_i,E_r,E_i := real and imaginary part of the H- and E-field*/
/*                index_mode   := specifying the considered mode (1=TE, 2=TM)*/
/*                   nz_atm    := Number of cells in the atmosphere*/
/*             frequency_nr    := Frequency number*/

int DetImpAtStations(CALC_2D_MT *mt_structure, DATA_STRUCT *data, double *H_r, double *H_i, double *E_r, double *E_i, int index_mode, long nz_atm, long frequency_nr)
{
	int index_x, index_z;
	long i,j,k,index_freq, index_stat;
	double pos_x, pos_z, int_z1, int_z2;
	dcomplex H,E,Imped;


	/*Loop over all stations in the 2-D layer*/
    for(i=0;i<mt_structure->freq_ndata[frequency_nr];i++)
	{

		pos_x = mt_structure->org[0][frequency_nr];
		index_x = 0;

		pos_z = mt_structure->org[1][frequency_nr];
		int_z1 = 0.0;
		int_z2 = 0.5*mt_structure->hz[frequency_nr][nz_atm];
		index_z = 0;

		for(k=0;k<mt_structure->nstat_mt;k++)
		{
			if(mt_structure->freq_stat[frequency_nr][i] == mt_structure->mts[k])
			{
				index_stat = k;
				break;
			}
		}

		/*Find the cell, in which the station is located*/
		/************************************************/
		/*Loop over all cells in the x-direction*/
		for(j=0;j<mt_structure->nx[frequency_nr];j++)
		{
			/*Find the x-position (relative to the grid cell centers)*/
			if((mt_structure->x[k] >= pos_x) && (mt_structure->x[k] <= (pos_x + mt_structure->hx[frequency_nr][j])))
			{
				goto funnet_x;
			}

			pos_x = mt_structure->hx[frequency_nr][j] + pos_x;
			index_x++;
		}

		/*Station is not found within the grid*/
		printf("The MT-station %d is NOT located within the grid!!\n", mt_structure->mts[k]);
		exit(0);

		funnet_x:;

		/************************************************/
		/*Loop over all cells in the z-direction*/
		for(j=nz_atm;j<mt_structure->nz[frequency_nr];j++)
		{
			/*Find the z-position (relative to the grid cell border !!!)*/
			if((mt_structure->z[k] >= (pos_z - int_z1)) && (mt_structure->z[k] <= (pos_z + int_z2)))
			{
				goto funnet_z;
			}

			pos_z = mt_structure->hz[frequency_nr][j] + pos_z;
			int_z1 = int_z2;
			int_z2 = 0.5*mt_structure->hz[frequency_nr][j+1];

			if(j == (mt_structure->nz[frequency_nr] - 2))
				int_z2 = int_z2 + 0.5*mt_structure->hz[frequency_nr][j+1];

			index_z++;
		}

		/*Station is not found within the grid*/
		printf("The MT-station %d is NOT located within the grid!!\n", mt_structure->mts[k]);
		exit(0);

		funnet_z:;

		/************************************************/
		/*Determine the fields at the station locations:*/
		H = Complex(H_r[index_x*(mt_structure->nz[frequency_nr] - nz_atm) + index_z], H_i[index_x*(mt_structure->nz[frequency_nr] - nz_atm) + index_z]);
    	E = Complex(E_r[index_x*(mt_structure->nz[frequency_nr] - nz_atm) + index_z], E_i[index_x*(mt_structure->nz[frequency_nr] - nz_atm) + index_z]);

		/*Determine Impedances*/
		Imped = Cdiv(E,H);

		/*Find the indeces of the frequency for the data measurement*/
		for(k=0;k<data->nfreq_mt[mt_structure->freq_data[frequency_nr][i]];k++)
		{
			if(data->freq_mt[mt_structure->freq_data[frequency_nr][i]][k] == mt_structure->freq[frequency_nr])
			{
				index_freq = k;
				goto immer_weiter;
			}
		}

		printf("The frequency %f will not be used from station %d\n!!", mt_structure->freq[frequency_nr],data->mno[mt_structure->freq_data[frequency_nr][i]]);
		exit(0);

		immer_weiter:;


		/*TE-Mode*/
		if(index_mode == 1)
		{
			data->calc_real_mt_TE[mt_structure->freq_data[frequency_nr][i]][index_freq] = (double)Imped.r;
			data->calc_imag_mt_TE[mt_structure->freq_data[frequency_nr][i]][index_freq] = (double)Imped.i;

			/*Fields at the stations*/
			mt_structure->Ey_r[frequency_nr][i] = (double)E.r;
			mt_structure->Ey_i[frequency_nr][i] = (double)E.i;
			mt_structure->Hx_r[frequency_nr][i] = (double)H.r;
			mt_structure->Hx_i[frequency_nr][i] = (double)H.i;

		}
		/*TM-Mode*/
		else if(index_mode == 2)
		{
			data->calc_real_mt_TM[mt_structure->freq_data[frequency_nr][i]][index_freq] = (double)Imped.r;
			data->calc_imag_mt_TM[mt_structure->freq_data[frequency_nr][i]][index_freq] = (double)Imped.i;

			/*Fields at the stations*/
			mt_structure->Ex_r[frequency_nr][i] = (double)E.r;
			mt_structure->Ex_i[frequency_nr][i] = (double)E.i;
			mt_structure->Hy_r[frequency_nr][i] = (double)H.r;
			mt_structure->Hy_i[frequency_nr][i] = (double)H.i;

		}
		else
		{
			printf("Whether TE- nor TM-mode is specified in the DetImpAtStations-Routine\n");
			exit(0);
		}

	}

	return(1);
}


/*-------------------------------------------------------------*/
/*Transfer the data from the local to the global MT structure*/
/*(write the resistivity values in a temporary file)*/
/*Parameter:	layer_i = specify the layer*/
/*               freq_i = specify the frequency */
/*             mt_2D_struct = local MT structure*/
/*             mt_2D_struct_out = global MZ structures*/

int MakeGlobMTStruct(long layer_i,long freq_i, CALC_2D_MT *mt_2D_struct, CALC_2D_MT *mt_2D_struct_out)
{
	long i;

	/*Allocate memory and assign parameters to the new structures*/
	if(freq_i == 0)
	{
		mt_2D_struct_out[layer_i].nr_cells_2d_mt = mt_2D_struct->nr_cells_2d_mt;
		mt_2D_struct_out[layer_i].max_y = mt_2D_struct->max_y;
		mt_2D_struct_out[layer_i].min_y = mt_2D_struct->min_y;

		if(mt_2D_struct->nfreq == 0)
		{
			mt_2D_struct_out[layer_i].org[0] = (double *)memory(NULL, 1, sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].org[1] = (double *)memory(NULL, 1, sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].nx = (long *)memory(NULL, 1, sizeof(long),"MakeGlobMTStruct");	
			mt_2D_struct_out[layer_i].nz = (long *)memory(NULL, 1, sizeof(long),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].hx = (double **)memory(NULL, 1, sizeof(double *),"MakeGlobMTStruct");	
			mt_2D_struct_out[layer_i].hz = (double **)memory(NULL, 1, sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].nz_atm = (long *)memory(NULL, 1, sizeof(long),"MakeGlobMTStruct");	
			mt_2D_struct_out[layer_i].nz_ionos = (long *)memory(NULL, 1, sizeof(long),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].nz_basement = (long *)memory(NULL, 1, sizeof(long),"MakeGlobMTStruct");	
			mt_2D_struct_out[layer_i].nx_left_border = (long *)memory(NULL, 1, sizeof(long),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].nx_right_border = (long *)memory(NULL, 1, sizeof(long),"MakeGlobMTStruct");
		}
		else
		{
			mt_2D_struct_out[layer_i].org[0] = (double *)memory(NULL, mt_2D_struct->nfreq, sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].org[1] = (double *)memory(NULL, mt_2D_struct->nfreq, sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].nx = (long *)memory(NULL, mt_2D_struct->nfreq, sizeof(long),"MakeGlobMTStruct");	
			mt_2D_struct_out[layer_i].nz = (long *)memory(NULL, mt_2D_struct->nfreq, sizeof(long),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].hx = (double **)memory(NULL, mt_2D_struct->nfreq, sizeof(double *),"MakeGlobMTStruct");	
			mt_2D_struct_out[layer_i].hz = (double **)memory(NULL, mt_2D_struct->nfreq, sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].nz_atm = (long *)memory(NULL, mt_2D_struct->nfreq, sizeof(long),"MakeGlobMTStruct");	
			mt_2D_struct_out[layer_i].nz_ionos = (long *)memory(NULL, mt_2D_struct->nfreq, sizeof(long),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].nz_basement = (long *)memory(NULL,mt_2D_struct->nfreq, sizeof(long),"MakeGlobMTStruct");	
			mt_2D_struct_out[layer_i].nx_left_border = (long *)memory(NULL, mt_2D_struct->nfreq, sizeof(long),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].nx_right_border = (long *)memory(NULL, mt_2D_struct->nfreq, sizeof(long),"MakeGlobMTStruct");
		}
	}

	if(mt_2D_struct->nfreq != 0)
	{
		mt_2D_struct_out[layer_i].hx[freq_i] = (double *)memory(NULL, mt_2D_struct->nx[freq_i],sizeof(double),"MakeGlobMTStruct");	
		mt_2D_struct_out[layer_i].hz[freq_i] = (double *)memory(NULL, mt_2D_struct->nz[freq_i],sizeof(double),"MakeGlobMTStruct");

		mt_2D_struct_out[layer_i].org[0][freq_i] = mt_2D_struct->org[0][freq_i];
		mt_2D_struct_out[layer_i].org[1][freq_i] = mt_2D_struct->org[1][freq_i];
		mt_2D_struct_out[layer_i].nx[freq_i] = mt_2D_struct->nx[freq_i];
		mt_2D_struct_out[layer_i].nz[freq_i] = mt_2D_struct->nz[freq_i];

		for(i=0;i<mt_2D_struct->nx[freq_i];i++)
			mt_2D_struct_out[layer_i].hx[freq_i][i] = mt_2D_struct->hx[freq_i][i];
		for(i=0;i<mt_2D_struct->nz[freq_i];i++)
			mt_2D_struct_out[layer_i].hz[freq_i][i] = mt_2D_struct->hz[freq_i][i];

		mt_2D_struct_out[layer_i].nz_atm[freq_i] = mt_2D_struct->nz_atm[freq_i];
		mt_2D_struct_out[layer_i].nz_ionos[freq_i] = mt_2D_struct->nz_ionos[freq_i];
		mt_2D_struct_out[layer_i].nz_basement[freq_i] = mt_2D_struct->nz_basement[freq_i];
		mt_2D_struct_out[layer_i].nx_left_border[freq_i] = mt_2D_struct->nx_left_border[freq_i];
		mt_2D_struct_out[layer_i].nx_right_border[freq_i] = mt_2D_struct->nx_right_border[freq_i];
	}

	if(freq_i == 0)
	{
		mt_2D_struct_out[layer_i].nstat_mt = mt_2D_struct->nstat_mt;

		if(mt_2D_struct->nstat_mt != 0)
		{
			mt_2D_struct_out[layer_i].x = (double *)memory(NULL,mt_2D_struct->nstat_mt,sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].z = (double *)memory(NULL,mt_2D_struct->nstat_mt,sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].mts = (long *)memory(NULL,mt_2D_struct->nstat_mt,sizeof(long),"MakeGlobMTStruct");
		}
		else
		{
			mt_2D_struct_out[layer_i].x = (double *)memory(NULL,1,sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].z = (double *)memory(NULL,1,sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].mts = (long *)memory(NULL,1,sizeof(long),"MakeGlobMTStruct");
		}

		for(i=0;i<mt_2D_struct->nstat_mt;i++)
		{
			mt_2D_struct_out[layer_i].x[i] = mt_2D_struct->x[i];
			mt_2D_struct_out[layer_i].z[i] = mt_2D_struct->z[i];
			mt_2D_struct_out[layer_i].mts[i] = mt_2D_struct->mts[i];
		}

		mt_2D_struct_out[layer_i].nfreq = mt_2D_struct->nfreq;
	
		if(mt_2D_struct->nfreq == 0)
		{
			mt_2D_struct_out[layer_i].freq = (double *)memory(NULL,1,sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].freq_ndata = (long *)memory(NULL,1,sizeof(long),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].cell_scal_factor = (int *)memory(NULL,1,sizeof(int),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].freq_data = (long **)memory(NULL,1,sizeof(long *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].freq_stat = (long **)memory(NULL,1,sizeof(long *),"MakeGlobMTStruct");
	
			mt_2D_struct_out[layer_i].Ex_r = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Ex_i = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Ey_r = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Ey_i = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Hx_r = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Hx_i = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Hy_r = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Hy_i = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");


		}
		else
		{
			mt_2D_struct_out[layer_i].freq = (double *)memory(NULL,mt_2D_struct->nfreq,sizeof(double),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].freq_ndata = (long *)memory(NULL,mt_2D_struct->nfreq,sizeof(long),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].cell_scal_factor = (int *)memory(NULL,mt_2D_struct->nfreq,sizeof(int),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].freq_data = (long **)memory(NULL,mt_2D_struct->nfreq,sizeof(long *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].freq_stat = (long **)memory(NULL,mt_2D_struct->nfreq,sizeof(long *),"MakeGlobMTStruct");

			mt_2D_struct_out[layer_i].Ex_r = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Ex_i = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Ey_r = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Ey_i = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Hx_r = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Hx_i = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Hy_r = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].Hy_i = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
		}
	}

	if(mt_2D_struct->nfreq != 0)
	{
		mt_2D_struct_out[layer_i].freq_data[freq_i] = (long *)memory(NULL, mt_2D_struct->freq_ndata[freq_i],sizeof(long),"MakeGlobMTStruct");	
		mt_2D_struct_out[layer_i].freq_stat[freq_i] = (long *)memory(NULL, mt_2D_struct->freq_ndata[freq_i],sizeof(long),"MakeGlobMTStruct");

		mt_2D_struct_out[layer_i].Ex_r[freq_i] = (double *)memory(NULL,mt_2D_struct->freq_ndata[freq_i],sizeof(double),"MakeGlobMTStruct");
		mt_2D_struct_out[layer_i].Ex_i[freq_i] = (double *)memory(NULL,mt_2D_struct->freq_ndata[freq_i],sizeof(double),"MakeGlobMTStruct");
		mt_2D_struct_out[layer_i].Ey_r[freq_i] = (double *)memory(NULL,mt_2D_struct->freq_ndata[freq_i],sizeof(double),"MakeGlobMTStruct");
		mt_2D_struct_out[layer_i].Ey_i[freq_i] = (double *)memory(NULL,mt_2D_struct->freq_ndata[freq_i],sizeof(double),"MakeGlobMTStruct");
		mt_2D_struct_out[layer_i].Hx_r[freq_i] = (double *)memory(NULL,mt_2D_struct->freq_ndata[freq_i],sizeof(double),"MakeGlobMTStruct");
		mt_2D_struct_out[layer_i].Hx_i[freq_i] = (double *)memory(NULL,mt_2D_struct->freq_ndata[freq_i],sizeof(double),"MakeGlobMTStruct");
		mt_2D_struct_out[layer_i].Hy_r[freq_i] = (double *)memory(NULL,mt_2D_struct->freq_ndata[freq_i],sizeof(double),"MakeGlobMTStruct");
		mt_2D_struct_out[layer_i].Hy_i[freq_i] = (double *)memory(NULL,mt_2D_struct->freq_ndata[freq_i],sizeof(double),"MakeGlobMTStruct");

		mt_2D_struct_out[layer_i].freq[freq_i] = mt_2D_struct->freq[freq_i];
		mt_2D_struct_out[layer_i].freq_ndata[freq_i] = mt_2D_struct->freq_ndata[freq_i];
		mt_2D_struct_out[layer_i].cell_scal_factor[freq_i] = mt_2D_struct->cell_scal_factor[freq_i];

		for(i=0;i<mt_2D_struct->freq_ndata[freq_i];i++)
		{
			mt_2D_struct_out[layer_i].freq_data[freq_i][i] = mt_2D_struct->freq_data[freq_i][i];
			mt_2D_struct_out[layer_i].freq_stat[freq_i][i] = mt_2D_struct->freq_stat[freq_i][i]; 
			mt_2D_struct_out[layer_i].Ex_r[freq_i][i] = mt_2D_struct->Ex_r[freq_i][i];
			mt_2D_struct_out[layer_i].Ex_i[freq_i][i] = mt_2D_struct->Ex_i[freq_i][i];
			mt_2D_struct_out[layer_i].Ey_r[freq_i][i] = mt_2D_struct->Ey_r[freq_i][i];
			mt_2D_struct_out[layer_i].Ey_i[freq_i][i] = mt_2D_struct->Ey_i[freq_i][i];
			mt_2D_struct_out[layer_i].Hx_r[freq_i][i] = mt_2D_struct->Hx_r[freq_i][i];
			mt_2D_struct_out[layer_i].Hx_i[freq_i][i] = mt_2D_struct->Hx_i[freq_i][i];
			mt_2D_struct_out[layer_i].Hy_r[freq_i][i] = mt_2D_struct->Hy_r[freq_i][i];
			mt_2D_struct_out[layer_i].Hy_i[freq_i][i] = mt_2D_struct->Hy_i[freq_i][i];
		}
	}

		
	/**************************************************************************/
	/*The resistivity calues will be NOT stored in the parameter, because they maybe are too huge for the memory*/
	/*Instead they will be written out in external tenporary files*/

	if(freq_i == 0)
	{
		if(mt_2D_struct->nfreq == 0)
		{
			mt_2D_struct_out[layer_i].res_slice = (double **)memory(NULL,1,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].index_slice = (int **)memory(NULL,1,sizeof(int *),"MakeGlobMTStruct");
		}
		else
		{
			mt_2D_struct_out[layer_i].res_slice = (double **)memory(NULL,mt_2D_struct->nfreq,sizeof(double *),"MakeGlobMTStruct");
			mt_2D_struct_out[layer_i].index_slice = (int **)memory(NULL,mt_2D_struct->nfreq,sizeof(int *),"MakeGlobMTStruct");
		}
	}

	if(mt_2D_struct->nfreq != 0)
	{
		mt_2D_struct_out[layer_i].res_slice[freq_i] = (double *)memory(NULL,1,sizeof(double),"MakeGlobMTStruct");
		mt_2D_struct_out[layer_i].index_slice[freq_i] = (int *)memory(NULL,1,sizeof(int),"MakeGlobMTStruct");
	}

	
return(0);
}



#undef M_KONST

#undef RES_AIR
#undef NR_OF_1D_CALCULATION
#undef GLOBAL_FACTOR_REFINEMENT
#undef LOCAL_FACTOR_REFINEMENT
#undef LOCAL_AREA_FACTOR
#undef RES_IONOS
#undef DZ_IONOS
#undef FACTOR_AIR
#undef FACTOR_SCAL_AIR
#undef MIN_SIZE_AIR
#undef FACTOR_BASEMENT
#undef FACTOR_SCAL_BAS
#undef FACTOR_BORDER_LEFT_RIGHT
#undef FACTOR_SCAL_BORDER

#undef res