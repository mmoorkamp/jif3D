#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "inv3d.h"
#include "fkt.h"

/*---------------------------------------------------------------------------*/
/* Build the Data and the Geometry structure from the the DAT- and GEO-files*/
/* parameters: *geo       := Pointer on the unmerged geometry structures */
/*             *files_geo := Pointer on the filenames of the GEO files */
/*             *geo_all   := Pointer on the merged geometry structure*/
/*	    	   *data       := Pointer on the unmerged data structures */
/*             *files_data := Pointer on the filenames of the DAT files */
/*             *data_all   := Pointer on the merged data structure*/
/*				ngeo_files := Number of GEO resp. DAT files*/
/*				eps        := Distance between to shot positions, so that they are considered as the same position*/
/*        *kind_of_data_mt := Kind of data used for MT (TE-Mode =1, TM-Mode=2, 1D:Berdichewsky-average=3, 2D:Both TE- and TM-Mode=3)*/

/*OUTPUT: Kind of MT-data format: 0 = NO MT data, 1 = Only one impedance per frequency, 2 = Two impedances per frequency*/

int MakeGeoandDataStructure(GEOMETRY *geo,FILENAME *files_geo,GEOMETRY *geo_all,DATA_STRUCT *data,FILENAME *files_dat, DATA_STRUCT *data_all, int ngeo_file, float eps, int *kind_of_data_mt)
{
	int i,j;
	long *count, *count_tmp, mt_data_format, mt_data_format_tmp; /*specify the kind of found MT data format (0= no mt-data, 1= only one impedance per frequneny, 2 = two impedances per frequency)*/

	/* Read geometry files*/
		
	printf("Read in GEO-inputfiles:\n");	
	printf("----------------\n\n");
		
	for(i=0;i<ngeo_file;i++)
	{
		ReadGeoFiles(files_geo[i].name, &(geo[i]));
		printf("%s is read in\n",files_geo[i].name);
	}
	printf("----------------\n\n");

	free(files_geo);


	/*-----------------------------------------------------------*/
	/* Read data files*/
	printf("Read in DAT-inputfiles:\n");
	printf("----------------\n\n");

    count = (long *)memory(NULL,ngeo_file,sizeof(long),"MakeGeoandDataStructure");
	count_tmp = (long *)memory(NULL,ngeo_file,sizeof(long),"MakeGeoandDataStructure");
	
	mt_data_format_tmp = 0;
	mt_data_format = 0;

	for(i=0;i<ngeo_file;i++)
	{
		count[i] =0;

		/*For the OLD Hanruedi-Format activate "ReadDataFilesOLD" for the Bjoern-Format activate "ReadDataFiles"*/
		mt_data_format = ReadDataFiles(files_dat[i].name, &(data[i]), &(geo[i]), &(count_tmp[i]), kind_of_data_mt);
		/*ReadDataFilesOLD(files_dat[i].name, &(data[i]), &(geo[i]));*/
		
		/*Check if the MT data format is not changing in the routine*/
		if(mt_data_format != 0)
		{
			if(i>0)
			{
				if(mt_data_format_tmp != mt_data_format && mt_data_format_tmp != 0)
				{
					printf("The kind of MT-data formats vary in the different dat-files\n");
					printf("Adjust the different formats!!\n");
					exit(0);
				}
			}
			mt_data_format_tmp = mt_data_format;
		}

		count[i] = count_tmp[i] + count[i]; 
		printf("%s is read in\n",files_dat[i].name);
	}
	printf("----------------\n\n");

	free(files_dat);

	/*-----------------------------------------------------------*/
	/* Merge geometry structures(LoadData HAS TO be performed before) */
	*geo_all = MergeGeometry(geo,ngeo_file);

	/*-----------------------------------------------------------*/
	/* Merge data structures(ReadDataFiles and ReadGeoFiles HAVE TO be performed before))*/
	*data_all = MergeData(data,ngeo_file, geo);
 
	/*Optimize data structure (Shots at the same positions are merged)*/
	OptiData(data_all,geo_all, eps);

	/*Calculate the shot-receiver distances for all shot-receiver combinations*/
	CalcShotRecDist(data_all,*geo_all);

	/*Determine the number of active receivers/shots for a shot/receiver position:*/
	StatData(data_all, *geo_all);

	for(i=0;i<ngeo_file;i++)
	{
		free(geo[i].x);
		free(geo[i].y);
		free(geo[i].z); 
		free(geo[i].coor_info);
	}
	free(geo);

	for(i=0;i<ngeo_file;i++)
	{
		for(j=0;j<count[i];j++)
		{
			free(data[i].freq_mt[j]);
			free(data[i].real_mt_TE[j]);
			free(data[i].imag_mt_TE[j]);
			free(data[i].real_mt_TM[j]);
			free(data[i].imag_mt_TM[j]);
			free(data[i].weigth_mt[j]);
		}

		free(data[i].sno);
		free(data[i].rno);
		free(data[i].tobs);
		free(data[i].weigth_seis);
		free(data[i].shots);
		free(data[i].recs);

		free(data[i].obs_grav);
		free(data[i].gno);
		free(data[i].gravs);
		free(data[i].weigth_grav);

		free(data[i].nfreq_mt);
		free(data[i].mno);
		free(data[i].mts);
		free(data[i].freq_mt);
		free(data[i].real_mt_TE);
		free(data[i].imag_mt_TE);
		free(data[i].real_mt_TM);
		free(data[i].imag_mt_TM);
		free(data[i].weigth_mt);
	}
	free(data);
	
	free(count);
	free(count_tmp);

	return(mt_data_format);
}


/*---------------------------------------------------------------------------*/
/* read in geometry files*/
/* parameters: *fname     := Pointer on the geo filenames */
/*             *geo       := Pointer on the structure of the receiver/shot positions and locations of gravity and MT measurement */

int ReadGeoFiles(char *fname, GEOMETRY *geo)
{
	int zero_line;
	long a,no;
	FILE *inf;
	char line[128], *last_trace;

	geo->x = (float *)memory(NULL,1,sizeof(float),"ReadGeoFiles");
	geo->y = (float *)memory(NULL,1,sizeof(float),"ReadGeoFiles");
	geo->z = (float *)memory(NULL,1,sizeof(float),"ReadGeoFiles");
	geo->coor_info = (int *)memory(NULL,1,sizeof(int),"ReadGeoFiles");

	a=0;
		inf = fopen(fname,"rt");
		if (inf == NULL)
		{
			fprintf(stderr,"Unable to open %s\n", fname);
			exit(0);
		}

		while(!feof(inf))
		{
			last_trace = fgets(line,128,inf);

			if(last_trace == NULL)
				break;

			zero_line = (int)strlen(line);

			if(zero_line > 1)
			{
				geo->x = (float *)memory((char *)geo->x,a+1,sizeof(float),"ReadGeoFiles");
				geo->y = (float *)memory((char *)geo->y,a+1,sizeof(float),"ReadGeoFiles");
				geo->z = (float *)memory((char *)geo->z,a+1,sizeof(float),"ReadGeoFiles");
				geo->coor_info = (int *)memory((char *)geo->coor_info,a+1,sizeof(int),"ReadGeoFiles");
				sscanf(line,"%d %f %f %f %d\n",&no, &(geo->x[a]),&(geo->y[a]),&(geo->z[a]),&(geo->coor_info[a]));
			
				if(no != a+1)
				{
					printf("Error while reading geometry file no = %d a = %d\n in file %s\n", no, a+1, *fname);
					exit(0);
				}
			a++;
			}
		}

		fclose(inf);

	geo->nstat = a;
	/*Seismic*/
	geo->nshot = 0; 	/*REMARK:The number of seismic shots and receivers locations, gravimetric and MT stations will be determined while including the Data-Files (in ReadDataFiles)*/
	geo->nrec = 0;
	/*Gravimety*/
	geo->nstat_grav =0; 
	/*MT*/
	geo->nstat_mt = 0;

	geo->eps = 0.0;		/*REMARK: eps will be set in the routine OptiData*/

	if(geo->nstat == 0)
	{
		printf("!!! The geometry file %s is probably empty !!!\n", fname);
		exit(0);
	}

	return(1);
}

/*---------------------------------------------------------------------------*/
/* merge several GEOMETRY structures to one GEOMETRY structure*/
/* parameters:   *geo       := Pointer on the structures of the receiver/shot locations and gravity/MT stations that should be merged */
/*               num_geo    := number of structures that should be merged*/
/* Output: The merged output GEOMETRY structure*/
/* LoadData HAVE TO be performed before!!*/

 GEOMETRY MergeGeometry(GEOMETRY *geo, int num_geo)
{
  int i,tmp_nstat;
  GEOMETRY geo_all;	
  FILE *out;

  geo_all.nstat=0;

  geo_all.nshot=0;
  geo_all.nrec=0;
  geo_all.nstat_grav=0;
  geo_all.nstat_mt=0;

  /*Determine number of rows in the final GEOMETRY structure*/
  for(i=0;i<num_geo;i++)
  {
	geo_all.nstat = geo_all.nstat + geo[i].nstat;

	geo_all.nshot = geo_all.nshot + geo[i].nshot;
	geo_all.nrec  = geo_all.nrec  + geo[i].nrec;
	geo_all.nstat_grav = geo_all.nstat_grav + geo[i].nstat_grav;
	geo_all.nstat_mt = geo_all.nstat_mt + geo[i].nstat_mt;
  }


  /*Allocate memory*/
				geo_all.x = (float *)memory(NULL,geo_all.nstat,sizeof(float),"MergeGeometry");
				geo_all.y = (float *)memory(NULL,geo_all.nstat,sizeof(float),"MergeGeometry");
				geo_all.z = (float *)memory(NULL,geo_all.nstat,sizeof(float),"MergeGeometry");
				geo_all.coor_info = (int *)memory(NULL,geo_all.nstat,sizeof(int),"MergeGeometry");

  tmp_nstat=0;
  
  /* merge structures */
  for(i=0;i<num_geo;i++)
  {
		memcpy(&(geo_all.x[tmp_nstat]),&(geo[i].x[0]),geo[i].nstat*sizeof(float));
		memcpy(&(geo_all.y[tmp_nstat]),&(geo[i].y[0]),geo[i].nstat*sizeof(float));
		memcpy(&(geo_all.z[tmp_nstat]),&(geo[i].z[0]),geo[i].nstat*sizeof(float));
		memcpy(&(geo_all.coor_info[tmp_nstat]),&(geo[i].coor_info[0]),geo[i].nstat*sizeof(int));
		tmp_nstat = geo[i].nstat + tmp_nstat;
  }

    
   printf("All structures from the GEO-files are merged together\n");
   printf(" Number of shot and receiver locations: %d, %d\n", geo_all.nshot, geo_all.nrec);
   printf(" Number of the gravity and MT stations: %d, %d\n", geo_all.nstat_grav, geo_all.nstat_mt);
   printf(" Number of ALL station locations: %d\n", geo_all.nstat);
   printf("----------------\n\n");

   			/***************************************/
			 out = fopen("merged.geo","wt");
			 for(i=0;i<geo_all.nstat;i++)
			   	fprintf(out,"%d %f %f %f %d\n", i+1, geo_all.x[i], geo_all.y[i], geo_all.z[i], geo_all.coor_info[i]);
			 fclose(out);
			/***************************************/

	geo_all.eps =0.0; /*REMARK: eps will be set in the routine OptiData*/
 

  return(geo_all);
}

/*---------------------------------------------------------------------------*/
/*Read in the "Hansruedi" dat-files and build up the data structure ttdata*/
/* Parameter: *ttdata := Pointer on the data structure*/
/*            *geo    := Pointer on the geometry structure*/
/* REMARK: Also the number of shots and receivers will be determined by this routine (and will be written 
/*		   in the corresponding Geometry structure)*/

 int ReadDataFilesOLD(char *infile,DATA_STRUCT *ttdata,GEOMETRY *geo)
{

   FILE *inf;
   int a,b,c,d;
   char line[128];
   int shnr,frec,lrec;
   int found;
   char hcr;
   char mess_str[128];
/* double minx,minz;*/
   int nitems;
   
   inf = fopen(infile,"rt");
   ttdata->ndata_seis = 0;
   geo->nshot = geo->nrec = 0;
   ttdata->tobs = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD");
   ttdata->weigth_seis = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD");
   ttdata->rno = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD");
   ttdata->sno = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD");
   if (ttdata->shots == NULL) ttdata->shots = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD");
   if (ttdata->recs == NULL) ttdata->recs = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD");

   ttdata->obs_grav = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD"); /*<- NOT used*/
   ttdata->gravs = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD"); /*<- NOT used*/
   ttdata->gno = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD"); /*<-NOT used */
   ttdata->weigth_grav = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD"); /*<- NOT used*/

   ttdata->nfreq_mt = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD"); /*<-NOT used */
   ttdata->mno = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD"); /*<-NOT used */
   ttdata->mts = (int *)memory(NULL,1,sizeof(int),"LoadDataOLD"); /*<- NOT used*/
   ttdata->freq_mt = (double **)memory(NULL,1,sizeof(double *),"LoadDataOLD"); /*<-NOT used */
   ttdata->freq_mt[0] = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD");
   ttdata->weigth_mt = (double **)memory(NULL,1,sizeof(double *),"LoadDataOLD");/*<-NOT used */
   ttdata->weigth_mt[0] = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD");
   ttdata->real_mt_TE = (double **)memory(NULL,1,sizeof(double *),"LoadDataOLD");/*<-NOT used */
   ttdata->real_mt_TE[0] = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD");
   ttdata->imag_mt_TE = (double **)memory(NULL,1,sizeof(double *),"LoadDataOLD");/*<-NOT used */
   ttdata->imag_mt_TE[0] = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD");
   ttdata->real_mt_TM = (double **)memory(NULL,1,sizeof(double *),"LoadDataOLD");/*<-NOT used */
   ttdata->real_mt_TM[0] = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD");
   ttdata->imag_mt_TM = (double **)memory(NULL,1,sizeof(double *),"LoadDataOLD");/*<-NOT used */
   ttdata->imag_mt_TM[0] = (double *)memory(NULL,1,sizeof(double),"LoadDataOLD");
	/*ATTENTION: The "dummy" gravity and MT allocation is NOT tested until now*/

   while(!feof(inf))
   {
      memset(line,0,80);
      fgets(line,82,inf);
      if (feof(inf))
      {
	 fclose(inf);
	 break;
      }
      
      if (!strncmp(line,"       99999       99999",24))
      {
		printf("End of file marker detected! in %s\n",infile);
		fclose(inf);
		break;
      }
      if (!strncmp(line+7,"99999",5))
      {
	 /* Read header line */
	 sscanf(line+13,"%d %d %d",&shnr,&frec,&lrec);
	 
	 /* Update shot list */
	 found = 0;
	 for (a=0;a<geo->nshot;a++)
	 {
	    if (ttdata->shots[a] == shnr)
	    {
	       found = 1;
	       break;
	    }
	 }
	 if (!found)
	 {
        ttdata->shots = (int *)memory((char *)ttdata->shots,geo->nshot+1,sizeof(int),"LoadDataOLD");
	    ttdata->shots[geo->nshot] = shnr;
	    geo->nshot++;
	 }
	 
	 /* Update receiver list */
	 for (b=frec;b<=lrec;b++)
	 {
	    found = 0;
	    for (a=0;a<geo->nrec;a++)
	    {
	       if (ttdata->recs[a] == b)
	       {
	          found = 1;
	          break;
	       }
	    }
	    if (!found)
	    {
           ttdata->recs = (int *)memory((char *)ttdata->recs,geo->nrec+1,sizeof(int),"LoadDataOLD");
	       ttdata->recs[geo->nrec] = b;
	       (geo->nrec)++;
	    }
	 }
	 
	 /* Read traveltime data */
	 b = (ttdata->ndata_seis); 
	 c = 0;
	 d = 0;
	 for (a=(ttdata->ndata_seis);a<((ttdata->ndata_seis)+lrec-frec+1);a++) 
	 {
            ttdata->tobs = (double *)memory((char *)ttdata->tobs,b+1,sizeof(double),"LoadDataOLD");
			ttdata->weigth_seis = (double *)memory((char *)ttdata->weigth_seis,b+1,sizeof(double),"LoadDataOLD");
            ttdata->rno = (int *)memory((char *)ttdata->rno,b+1,sizeof(int),"loadDataOLD");
            ttdata->sno = (int *)memory((char *)ttdata->sno,b+1,sizeof(int),"LoadDataOLD");
	    nitems = fscanf(inf,"%lf",&ttdata->tobs[b]);
	    d++;
	    if (nitems != 1)
	    {
	       printf("Error wihile reading data file!\n");
	       printf("Shot no: %d\n",shnr);
	       printf("First rec %d   last rec %d\n",frec,lrec);
	       printf("Erroneous element %d\n",a-(ttdata->ndata_seis)+1);
               exit(0);
	    }
	    ttdata->sno[b] = shnr;
	    ttdata->rno[b] = a - (ttdata->ndata_seis) + frec;
		ttdata->weigth_seis[b] = 1.0; /*Remark: In the Hansruedi Format is no weighting implemented*/
	    if (ttdata->tobs[b] > 500000.0) 
	    {
			printf(mess_str,"Shot No. %d:  Station No. %d Suspicious value (%lf)\n", shnr,a-(ttdata->ndata_seis)+frec,ttdata->tobs[b]);
		}
	    if (ttdata->tobs[b] > 0.0) 
	    {
	       b++; 
	    }
	    else 
	    {
	       c++;
	    }
	 }
	 if ((d-1) % 10)
	 {
	    hcr = 0;
	    while (hcr != '\n') hcr = getc(inf);
	 }
	 ttdata->ndata_seis += lrec-frec+1 - c;
      }
   }     
   fclose(inf);

   /* Apply timeshift */
   if (ttdata->timeshift) for (a=0;a<ttdata->ndata_seis;a++) {ttdata->tobs[a] -= ttdata->timeshift; if (ttdata->tobs[a] < 0.0) ttdata->tobs[a] = 0.0;}

   /*Set gravity and MT stations and data to zero:*/
   geo->nstat_grav = 0;		
   ttdata->ndata_grav = 0;

   geo->nstat_mt = 0;
   ttdata->ndata_mt = 0;

   return(1);

}

/*---------------------------------------------------------------------------*/
/*Read in the "Bjoerns" dat-files and build up the data structure ttdata*/
/* Parameter: *ttdata := Pointer on the data structure*/
/*            *geo    := Pointer on the geometry structure*/
/*			  *count  := Minimum (!!!) number of gravity measurements (Required to free the memory lateron) */
/*        *kind_of_data_mt := Kind of data used for MT (TE-Mode =1, TM-Mode=2, 1D:Berdichewsky-average=3, 2D:Both TE- and TM-Mode=3)*/

/* REMARK: Also the number of shots/receivers, gravity and MT locations will be determined by this routine (and will be written */
/*		   in the corresponding Geometry structure)*/

/*OUTPUT: The kind of identied MT data format: 0= no MT data, 1= Only one impedance per frequency, 2= Two impedances per frequency*/

int ReadDataFiles(char *infile,DATA_STRUCT *ttdata,GEOMETRY *geo, long *count, int *kind_of_data_mt)
{
	 FILE *inf;
	 int zero_line;
     int b,c,d,e,f,g,h,j,fr_i,no,nr_col, nr_col2;
	 int found_sou, found_rec, found_grav, found_mt;
	 int sou_index, rec_index, index_input_mt, tmp_index_input_mt, counter_mt;
	 int kind_of_measure; /*Specify, if the measurement is seismic (=1), gravity (=2), or MT (=3)*/
	 float ttime,measured_data,measured_data_weight, measured_data_weight2; /*weighting of the data*/
	 float freq_mt, real, imag, real2, imag2; /*Real and imaginary part in the MT measurements*/
     char line[128], *last_trace;
  
	/*Seismic parameters*/
	ttdata->ndata_seis = 0;
    geo->nshot = geo->nrec = 0;
	/*Gravitry parameters*/
	ttdata->ndata_grav = 0;
	/*MT parameters*/
	ttdata->ndata_mt = 0;
	ttdata->ndata_freq_mt = 0;

	b=0;
	d=0;
	e=0;
	f=0;
	g=0;
	h=0;
	j=0;

	 /* Allocate memory*/
	 /*Seismic parameters*/
     ttdata->tobs = (double *)memory(NULL,1,sizeof(double),"LoadData");
     ttdata->rno = (int *)memory(NULL,1,sizeof(int),"LoadData");
     ttdata->sno = (int *)memory(NULL,1,sizeof(int),"LoadData");

	 ttdata->shots = (int *)memory(NULL,1,sizeof(int),"LoadData");
     ttdata->recs = (int *)memory(NULL,1,sizeof(int),"LoadData");
	 ttdata->weigth_seis = (double *)memory(NULL,1,sizeof(double),"LoadData");

	 /*Gravimetry parameters*/
	 ttdata->obs_grav = (double *)memory(NULL,1,sizeof(double),"LoadData");
	 ttdata->gno = (int *)memory(NULL,1,sizeof(int),"LoadData");
	 ttdata->gravs = (int *)memory(NULL,1,sizeof(int),"LoadData");
	 ttdata->weigth_grav = (double *)memory(NULL,1,sizeof(double),"LoadData");

	 /*MT parameters*/
	 ttdata->nfreq_mt = (int *)memory(NULL,1,sizeof(int),"LoadData");
	 ttdata->mno = (int *)memory(NULL,1,sizeof(int),"LoadData");
	 ttdata->mts = (int *)memory(NULL,1,sizeof(int),"LoadData");


	 /*--------------------------------------------------*/
	 /*Loop to estimate the minimum size of the MT arrays in the data files*/
	 /*(Required because the first index cannot be reallocated in C)*/
	 inf = fopen(infile,"rt");
	 *count = 1;
	 
	 if (inf == NULL)
	 {
 		fprintf(stderr,"Unable to open %s\n", infile);
 		exit(0);
	 }

	
	 while(!feof(inf))
	 {
		last_trace = fgets(line,128,inf);

		if(last_trace == NULL)
			break;

		zero_line = (int)strlen(line); /*Looking for empty traces*/

		if(zero_line > 1)
		{
			sscanf(line,"%d %d %d %d %f %f\n",&no,&sou_index,&rec_index,&kind_of_measure, &measured_data, &measured_data_weight);
			
			if(kind_of_measure == 3)
				(*count)++;
		}
	 }

	 fclose(inf);
	/*--------------------------------------------------*/

	 ttdata->freq_mt = (double **)memory(NULL,(*count),sizeof(double *),"LoadData");
     ttdata->real_mt_TE = (double **)memory(NULL,(*count),sizeof(double *),"LoadData");
     ttdata->imag_mt_TE = (double **)memory(NULL,(*count),sizeof(double *),"LoadData");
     ttdata->real_mt_TM = (double **)memory(NULL,(*count),sizeof(double *),"LoadData");
     ttdata->imag_mt_TM = (double **)memory(NULL,(*count),sizeof(double *),"LoadData");
	 ttdata->weigth_mt = (double **)memory(NULL,(*count),sizeof(double *),"LoadData");

	 for(c=0;c<(*count);c++)
	 {
		 	ttdata->freq_mt[c] = (double *)memory(NULL,1,sizeof(double),"LoadData");
			ttdata->real_mt_TE[c] = (double *)memory(NULL,1,sizeof(double),"LoadData");
			ttdata->imag_mt_TE[c] = (double *)memory(NULL,1,sizeof(double),"LoadData");
			ttdata->real_mt_TM[c] = (double *)memory(NULL,1,sizeof(double),"LoadData");
			ttdata->imag_mt_TM[c] = (double *)memory(NULL,1,sizeof(double),"LoadData");
			ttdata->weigth_mt[c] = (double *)memory(NULL,1,sizeof(double),"LoadData");
	 }
	 /*--------------------------------------------------*/
	 
	 
	 inf = fopen(infile,"rt");

		if (inf == NULL)
		{
			fprintf(stderr,"Unable to open %s\n", infile);
			exit(0);
		}

		index_input_mt = 0;
		tmp_index_input_mt = 0;
		counter_mt = 0;

		while(!feof(inf))
			{
				last_trace = fgets(line,128,inf);

				if(last_trace == NULL)
					break;
				
				zero_line = (int)strlen(line); 

				if(zero_line > 1)
				{
					nr_col = sscanf(line,"%d %d %d %d %f %f\n",&no,&sou_index,&rec_index,&kind_of_measure, &measured_data, &measured_data_weight);


					/*---------------------------MT parameters-----------------------*/
					/*Only if the sou_index is specified in the geofile, the line in the dat were considered*/
					if(kind_of_measure == 3 && sou_index <= (geo->nstat)) 
					{
						ttdata->nfreq_mt = (int *)memory((char *)ttdata->nfreq_mt,g+1,sizeof(int),"LoadData");
						ttdata->mno = (int *)memory((char *)ttdata->mno,g+1,sizeof(int),"LoadData");
						ttdata->nfreq_mt[g] = (int) measured_data;
						ttdata->mno[g] = sou_index;

						for(fr_i=0;fr_i<ttdata->nfreq_mt[g];fr_i++)
						{

							last_trace = fgets(line,128,inf);

							zero_line = (int)strlen(line); /*Looking for empty traces*/

							if(zero_line <= 1 || last_trace == NULL)
								goto jump;
							
							nr_col2 = sscanf(line,"%f %f %f %f %f %f\n",&freq_mt, &real, &imag, &real2, &imag2, &measured_data_weight2);
							
							/*If input data consists ONLY of one impedance*/
							if(nr_col2 == 3 || nr_col2 ==4)
								index_input_mt = 1;
							/*If input data consists of two impedances*/
							else
								index_input_mt = 2;


							if(counter_mt > 0)
							{
								if(index_input_mt != tmp_index_input_mt)
								{
									printf("The kind of MT-data format varies within a dat-file\n");
									printf("Adjust the different formats!!\n");
									exit(0);
								}
							}

							counter_mt++;
							tmp_index_input_mt = index_input_mt;

							if(nr_col2 == 3 || nr_col2 == 5)
								measured_data_weight2 = 1.0; /*If the weighting is not specified, it is set to 1.0*/
							if(nr_col2 == 4)
								measured_data_weight2 = real2;

							if(nr_col2 < 3 || nr_col2 > 6)
							{
								printf("Input format is not correct while reading the MT-data\nNumber of columns is %d\n",nr_col);
								exit(0);
							}


							/*Read in the frequency, real/imaginary part and the weighting factor for the MT data*/
							ttdata->freq_mt[g] = (double *)memory((char *)ttdata->freq_mt[g],fr_i+1,sizeof(double),"LoadData");
							ttdata->real_mt_TE[g] = (double *)memory((char *)ttdata->real_mt_TE[g],fr_i+1,sizeof(double),"LoadData");
							ttdata->imag_mt_TE[g] = (double *)memory((char *)ttdata->imag_mt_TE[g],fr_i+1,sizeof(double),"LoadData");
							ttdata->real_mt_TM[g] = (double *)memory((char *)ttdata->real_mt_TM[g],fr_i+1,sizeof(double),"LoadData");
							ttdata->imag_mt_TM[g] = (double *)memory((char *)ttdata->imag_mt_TM[g],fr_i+1,sizeof(double),"LoadData");
							ttdata->weigth_mt[g] = (double *)memory((char *)ttdata->weigth_mt[g],fr_i+1,sizeof(double),"LoadData");


							ttdata->freq_mt[g][fr_i] = freq_mt;

							if(index_input_mt != 1)
							{
								ttdata->real_mt_TE[g][fr_i] = real;
								ttdata->imag_mt_TE[g][fr_i] = imag;

								ttdata->real_mt_TM[g][fr_i] = real2;
								ttdata->imag_mt_TM[g][fr_i] = imag2;
							}
							else if(index_input_mt == 1 && *kind_of_data_mt != 2)
							{
								ttdata->real_mt_TE[g][fr_i] = real;
								ttdata->imag_mt_TE[g][fr_i] = imag;

								ttdata->real_mt_TM[g][fr_i] = 0.0;
								ttdata->imag_mt_TM[g][fr_i] = 0.0;
							}
							else
							{
								ttdata->real_mt_TM[g][fr_i] = real;
								ttdata->imag_mt_TM[g][fr_i] = imag;

								ttdata->real_mt_TE[g][fr_i] = 0.0;
								ttdata->imag_mt_TE[g][fr_i] = 0.0;
							}


							ttdata->weigth_mt[g][fr_i] = measured_data_weight2;
						}

						g++;
					
						/*Determine the list of MT stations:*/	
						found_mt = 0;
						
						for(c=0;c<(g-1);c++)
						{
							if(ttdata->mno[c] == sou_index)
							{
								printf(" WARNING !! The MT position %d in %s is used more than one time !!!\n", sou_index, infile);
								found_mt = 1;
							}
						}
							
						if(found_mt != 1)
						{
							/*Number of MT stations*/
						    (geo->nstat_mt)++;

							ttdata->mts = (int *)memory((char *)ttdata->mts,j+1,sizeof(int),"LoadData");
							(ttdata->mts[j]) = sou_index;
							j++;
						}

					}

					
					/*If the last column (including the weights of the measured data) is empty*/
					if(nr_col == 5)
						measured_data_weight = 1.0; /*Data are NOT specific weighted*/



					/*---------------------------Seismic parameters------------------*/
					/*If the ttime < 0, the corresponding first break was not picked!*/
					/*Only if the sou_index and the rec_index are specified in the geofile, the corresponding line in the data were considered*/
					if(kind_of_measure == 1 && measured_data > 0.0 && sou_index <= (geo->nstat) && rec_index <= (geo->nstat))
					{   
						ttime = measured_data;

						ttdata->tobs = (double *)memory((char *)ttdata->tobs,b+1,sizeof(double),"LoadData");
						ttdata->rno = (int *)memory((char *)ttdata->rno,b+1,sizeof(int),"LoadData");
						ttdata->sno = (int *)memory((char *)ttdata->sno,b+1,sizeof(int),"LoadData");
						ttdata->weigth_seis = (double *)memory((char *)ttdata->weigth_seis,b+1,sizeof(double),"LoadData");


						ttdata->sno[b] = sou_index;
						ttdata->rno[b] = rec_index;
						ttdata->tobs[b] = ttime; /*Scaled to seconds*/
						ttdata->weigth_seis[b] = measured_data_weight;

						 if (ttdata->tobs[b] > 500000.0) 
						{
							printf("Shot No. %d:  Receiver No. %d Suspicious value %f s in file %s\n\n",ttdata->sno[b],ttdata->rno[b],(ttdata->tobs[b]),infile);
						}

						b++;

					
					/*Determine the number of shots and receiver:*/			
						found_sou = 0;
						found_rec = 0;

						for(c=0;c<(b-1);c++)
						{
							if((ttdata->sno[c]) == sou_index) 
								found_sou = 1;
							if((ttdata->rno[c]) == rec_index) 
								found_rec = 1;
						}

						if(found_sou != 1)
						{
							(geo->nshot)++;

							ttdata->shots = (int *)memory((char *)ttdata->shots,d+1,sizeof(int),"LoadData");
							(ttdata->shots[d]) = sou_index;
							d++;
						}
						if(found_rec != 1)
						{
							(geo->nrec)++;

							ttdata->recs = (int *)memory((char *)ttdata->recs,e+1,sizeof(int),"LoadData");
							(ttdata->recs[e]) = rec_index;
							e++;
						}
					}


					/*---------------------------Gravity parameters------------------*/
					/*Only if the sou_index is specified in the geofile, the corresponding line in the data were considered*/
					else if(kind_of_measure == 2 && sou_index <= (geo->nstat))
					{
						ttdata->obs_grav = (double *)memory((char *)ttdata->obs_grav,f+1,sizeof(double),"LoadData");
						ttdata->gno = (int *)memory((char *)ttdata->gno,f+1,sizeof(int),"LoadData");
						ttdata->weigth_grav = (double *)memory((char *)ttdata->weigth_grav,f+1,sizeof(double),"LoadData");

						ttdata->obs_grav[f] = measured_data;
						ttdata->gno[f] = sou_index;
						ttdata->weigth_grav[f] = measured_data_weight;

						f++;

						/*Determine the list of gravity stations:*/	
						found_grav = 0;
						
						for(c=0;c<(f-1);c++)
						{
							if(ttdata->gno[c] == sou_index)
							{
								printf(" WARNING !! The gravity position %d in %s is used more than one time !!!\n", sou_index, infile);
								found_grav = 1;
							}
						}
							
						if(found_grav != 1)
						{
							/*Number of gravity stations*/
						    (geo->nstat_grav)++;

							ttdata->gravs = (int *)memory((char *)ttdata->gravs,h+1,sizeof(int),"LoadData");
							(ttdata->gravs[h]) = sou_index;
							h++;
						}

					}
				}
			}

	  jump:;

	 /*If Only one impedance is available from the data it will be considered as TE-Mode*/
     if(index_input_mt == 1 && *kind_of_data_mt !=1 && *kind_of_data_mt != 2)
	 {
		*kind_of_data_mt = 1;
		printf("\nATTENTION: Because ONLY one impedance is found, the MT data will be\nconsidered as TE-Mode data!!!\n");
		printf("The former specification of the MT-data input is NOT relevant any more\n\n");
	 }

	  /*Seismic parameters*/
      ttdata->ndata_seis = b;
	  ttdata->timeshift = 0.0;	/*TIME SHIFTS CAN BE IMPLEMENTED HER !!!*/
	  /*Gravimetry parameters*/
	  ttdata->ndata_grav = f;
	  /*MT parameters*/
	  ttdata->ndata_mt = g;

	  for(c=0;c<ttdata->ndata_mt;c++)
		ttdata->ndata_freq_mt = ttdata->ndata_freq_mt + ttdata->nfreq_mt[c]; 

	 fclose(inf);

	 return(index_input_mt);
}

/*---------------------------------------------------------------------------*/
/* merge several DATA_STRUCT structures to one DATA_STRUCT structure*/
/* parameters:   *data       := Pointer on the data structures that should be merged*/
/*               num_dat     := number of structures that should be merged*/
/*				 *geo        := the UNMERGED geometry structure */
/*        kind_of_data_mt := Kind of data used for MT (TE-Mode =1, TM-Mode=2, 1D:Berdichewsky-average=3, 2D:Both TE- and TM-Mode=3)*/

/* Output: The merged output DATA_STRUCT structure*/

DATA_STRUCT MergeData(DATA_STRUCT *data, int num_dat, GEOMETRY *geo)
{
	int i,j,k; /*number of all stations*/
	int nstat_all,  nrecs_all, nshots_all, nstat_grav_all, nstat_mt_all;
	int tmp_nstat, tmp_nshots, tmp_nrecs, tmp_nstat_grav, tmp_nstat_mt;
	int tmp_ndata_seis, tmp_ndata_grav, tmp_ndata_mt;
	DATA_STRUCT data_all;


	data_all.ndata_seis = 0;
	data_all.timeshift = 0; /* NOT IMPLEMENTED UNTIL NOW !!!*/
	data_all.ndata_grav = 0;
	data_all.ndata_mt = 0;
	data_all.ndata_freq_mt = 0;

	nstat_all = 0;
	nshots_all = 0;
	nrecs_all = 0;
	nstat_grav_all = 0;
	nstat_mt_all = 0;

    /*Determine the number of the different stations*/
	for(i=0;i<num_dat;i++)
	{
		nstat_all = nstat_all + geo[i].nstat;
		nshots_all = nshots_all + geo[i].nshot;
		nrecs_all = nrecs_all + geo[i].nrec;
		nstat_grav_all = nstat_grav_all + geo[i].nstat_grav;
		nstat_mt_all = nstat_mt_all + geo[i].nstat_mt;
	
		data_all.ndata_seis = data_all.ndata_seis + data[i].ndata_seis;
		data_all.ndata_grav = data_all.ndata_grav + data[i].ndata_grav;
		data_all.ndata_mt = data_all.ndata_mt + data[i].ndata_mt;
		data_all.ndata_freq_mt = data_all.ndata_freq_mt + data[i].ndata_freq_mt;
	}

	    /*Allocate memory for the merged data structure*/
		/*Seismic*/
		if(data_all.ndata_seis != 0)
		{
			data_all.tobs = (double *)memory( NULL,data_all.ndata_seis,sizeof(double),"MergeData");
			data_all.tcalc = (double *)memory( NULL,data_all.ndata_seis,sizeof(double),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
			data_all.rno = (int *)memory( NULL,data_all.ndata_seis,sizeof(int),"MergeData");
			data_all.sno = (int *)memory( NULL,data_all.ndata_seis,sizeof(int),"MergeData");
			data_all.weigth_seis = (double *)memory( NULL,data_all.ndata_seis,sizeof(double),"MergeData");
			data_all.shots = (int *)memory( NULL, nshots_all, sizeof(int),"MergeData");
			data_all.recs= (int *)memory( NULL, nrecs_all, sizeof(int),"MergeData");
		}
		else /*Case that no seismic data exist*/
		{
			data_all.tobs = (double *)memory( NULL,1,sizeof(double),"MergeData");
			data_all.tcalc = (double *)memory( NULL,1,sizeof(double),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
			data_all.rno = (int *)memory( NULL,1,sizeof(int),"MergeData");
			data_all.sno = (int *)memory( NULL,1,sizeof(int),"MergeData");
			data_all.weigth_seis = (double *)memory( NULL,1,sizeof(double),"MergeData");
			data_all.shots = (int *)memory( NULL,1, sizeof(int),"MergeData");
			data_all.recs= (int *)memory( NULL,1, sizeof(int),"MergeData");
		}
	
		/*Gravimetry*/
		if(data_all.ndata_grav != 0)
		{
			data_all.obs_grav = (double *)memory(NULL,data_all.ndata_grav,sizeof(double),"MergeData");
			data_all.calc_grav = (double *)memory(NULL,data_all.ndata_grav,sizeof(double),"MergeData");/*Just memory allocated but not used until the forward calculation was performed*/
			data_all.gno = (int *)memory(NULL,data_all.ndata_grav,sizeof(int),"MergeData");
			data_all.gravs = (int *)memory(NULL,nstat_grav_all,sizeof(int),"MergeData");
			data_all.weigth_grav = (double *)memory(NULL,data_all.ndata_grav,sizeof(double),"MergeData");
		}
		else /*Case that no gravity data exist*/
		{
			data_all.obs_grav = (double *)memory(NULL,1,sizeof(double),"MergeData");
			data_all.calc_grav = (double *)memory(NULL,1,sizeof(double),"MergeData");/*Just memory allocated but not used until the forward calculation was performed*/
			data_all.gno = (int *)memory(NULL,1,sizeof(int),"MergeData");
			data_all.gravs = (int *)memory(NULL,1,sizeof(int),"MergeData");
			data_all.weigth_grav = (double *)memory(NULL,1,sizeof(double),"MergeData");
		}
		/*MT*/
		if(data_all.ndata_mt != 0)
		{
			data_all.nfreq_mt = (int *)memory(NULL,data_all.ndata_mt,sizeof(int),"MergeData");
			data_all.mno = (int *)memory(NULL,data_all.ndata_mt,sizeof(int),"MergeData");
			data_all.mts = (int *)memory(NULL,nstat_mt_all,sizeof(int),"MergeData");

			data_all.freq_mt = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData");
			data_all.real_mt_TE = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData");
			data_all.imag_mt_TE = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData");
			data_all.calc_real_mt_TE = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
			data_all.calc_imag_mt_TE = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
			data_all.real_mt_TM = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData");
			data_all.imag_mt_TM = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData");
			data_all.calc_real_mt_TM = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
			data_all.calc_imag_mt_TM = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
	
			data_all.weigth_mt = (double **)memory(NULL,data_all.ndata_mt,sizeof(double *),"MergeData");

		}
		else /*Case that no MT data exist*/
		{
			data_all.nfreq_mt = (int *)memory(NULL,1,sizeof(int),"MergeData");
			data_all.mno = (int *)memory(NULL,1,sizeof(int),"MergeData");
			data_all.mts = (int *)memory(NULL,1,sizeof(int),"MergeData");

			data_all.freq_mt = (double **)memory(NULL,1,sizeof(double *),"MergeData");
			data_all.real_mt_TE = (double **)memory(NULL,1,sizeof(double *),"MergeData");
			data_all.imag_mt_TE = (double **)memory(NULL,1,sizeof(double *),"MergeData");
			data_all.calc_real_mt_TE = (double **)memory(NULL,1,sizeof(double *),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
			data_all.calc_imag_mt_TE = (double **)memory(NULL,1,sizeof(double *),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
			data_all.real_mt_TM = (double **)memory(NULL,1,sizeof(double *),"MergeData");
			data_all.imag_mt_TM = (double **)memory(NULL,1,sizeof(double *),"MergeData");
			data_all.calc_real_mt_TM = (double **)memory(NULL,1,sizeof(double *),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
			data_all.calc_imag_mt_TM = (double **)memory(NULL,1,sizeof(double *),"MergeData"); /*Just memory allocated but not used until the forward calculation was performed*/
	
			data_all.weigth_mt = (double **)memory(NULL,1,sizeof(double *),"MergeData");
		}



		tmp_nstat = 0;
		tmp_nshots = 0;
		tmp_nrecs = 0;
		tmp_nstat_grav = 0;
		tmp_nstat_mt = 0;

		tmp_ndata_seis = 0;
		tmp_ndata_grav = 0;
		tmp_ndata_mt = 0;

	for(i=0;i<num_dat;i++)
	{

		/*Merge all single data structures*/
		/*Seismic*/
		if(data_all.ndata_seis != 0)
		{
			memcpy(&(data_all.tobs[tmp_ndata_seis]),&(data[i].tobs[0]),data[i].ndata_seis*sizeof(double));
			memcpy(&(data_all.rno[tmp_ndata_seis]),&(data[i].rno[0]),data[i].ndata_seis*sizeof(int));
			memcpy(&(data_all.sno[tmp_ndata_seis]),&(data[i].sno[0]),data[i].ndata_seis*sizeof(int));
			memcpy(&(data_all.weigth_seis[tmp_ndata_seis]),&(data[i].weigth_seis[0]),data[i].ndata_seis*sizeof(double));
			memcpy(&(data_all.shots[tmp_nshots]),&(data[i].shots[0]),geo[i].nshot*sizeof(int));
			memcpy(&(data_all.recs[tmp_nrecs]),&(data[i].recs[0]),geo[i].nrec*sizeof(int));
		}
		/*Gravimetry*/
		if(data_all.ndata_grav != 0)
		{
			memcpy(&(data_all.obs_grav[tmp_ndata_grav]),&(data[i].obs_grav[0]),data[i].ndata_grav*sizeof(double));
			memcpy(&(data_all.gno[tmp_ndata_grav]),&(data[i].gno[0]),data[i].ndata_grav*sizeof(int));
			memcpy(&(data_all.weigth_grav[tmp_ndata_grav]),&(data[i].weigth_grav[0]),data[i].ndata_grav*sizeof(double));
			memcpy(&(data_all.gravs[tmp_nstat_grav]),&(data[i].gravs[0]),geo[i].nstat_grav*sizeof(int));
		}
		/*MT*/
		if(data_all.ndata_mt != 0)
		{
			memcpy(&(data_all.nfreq_mt[tmp_ndata_mt]),&(data[i].nfreq_mt[0]), data[i].ndata_mt*sizeof(int));
			memcpy(&(data_all.mno[tmp_ndata_mt]),&(data[i].mno[0]), data[i].ndata_mt*sizeof(int));
			memcpy(&(data_all.mts[tmp_nstat_mt]), &(data[i].mts[0]), geo[i].nstat_mt*sizeof(int));
		}

		/*Give the stations new numbers*/
		/*Seismic*/
		for(j=0; j<data[i].ndata_seis;j++)
		{
			data_all.sno[tmp_ndata_seis + j] = tmp_nstat + data_all.sno[tmp_ndata_seis + j];
			data_all.rno[tmp_ndata_seis + j] = tmp_nstat + data_all.rno[tmp_ndata_seis + j];
		}
		for(j=0; j<geo[i].nshot;j++)
			data_all.shots[tmp_nshots + j] = tmp_nstat + data_all.shots[tmp_nshots + j];
		for(j=0; j<geo[i].nrec;j++)
			data_all.recs[tmp_nrecs + j]   = tmp_nstat + data_all.recs[tmp_nrecs + j];

		/*Gravity*/
		for(j=0; j<data[i].ndata_grav;j++)
			data_all.gno[tmp_ndata_grav + j] = tmp_nstat + data_all.gno[tmp_ndata_grav + j];
		for(j=0; j<geo[i].nstat_grav;j++)
			data_all.gravs[tmp_nstat_grav + j] = tmp_nstat + data_all.gravs[tmp_nstat_grav +j];

		/*MT*/
		for(j=0; j<data[i].ndata_mt;j++)
			data_all.mno[tmp_ndata_mt + j] = tmp_nstat + data_all.mno[tmp_ndata_mt +j];
		for(j=0; j<geo[i].nstat_mt;j++)
			data_all.mts[tmp_nstat_mt + j] = tmp_nstat + data_all.mts[tmp_nstat_mt + j];


		tmp_nstat = geo[i].nstat + tmp_nstat;
		tmp_nshots = geo[i].nshot + tmp_nshots;
		tmp_nrecs = geo[i].nrec + tmp_nrecs;
		tmp_nstat_grav = geo[i].nstat_grav + tmp_nstat_grav;
		tmp_nstat_mt = geo[i].nstat_mt + tmp_nstat_mt;

		tmp_ndata_seis = data[i].ndata_seis + tmp_ndata_seis;
		tmp_ndata_grav = data[i].ndata_grav + tmp_ndata_grav;
		tmp_ndata_mt = data[i].ndata_mt + tmp_ndata_mt;


	}

	/*Allocate the memory for the arrays of the MT data including the frequencies, real and imaginary parts*/
	for(i=0;i<data_all.ndata_mt;i++)
	{
		if(data_all.nfreq_mt != 0)
		{
			data_all.freq_mt[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");
			data_all.real_mt_TE[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");
			data_all.imag_mt_TE[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");
			data_all.calc_real_mt_TE[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");/*Just memory allocated but not used until the forward calculation was performed*/
			data_all.calc_imag_mt_TE[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");/*Just memory allocated but not used until the forward calculation was performed*/
			data_all.real_mt_TM[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");
			data_all.imag_mt_TM[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");
			data_all.calc_real_mt_TM[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");/*Just memory allocated but not used until the forward calculation was performed*/
			data_all.calc_imag_mt_TM[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");/*Just memory allocated but not used until the forward calculation was performed*/

			data_all.weigth_mt[i] = (double *)memory(NULL,data_all.nfreq_mt[i],sizeof(double),"MergeData");
		}
		else
		{
			data_all.freq_mt[i] = (double *)memory(NULL,1,sizeof(double),"MergeData");
			data_all.real_mt_TE[i] = (double *)memory(NULL,1,sizeof(double),"MergeData");
			data_all.imag_mt_TE[i] = (double *)memory(NULL,1,sizeof(double),"MergeData");
			data_all.calc_real_mt_TE[i] = (double *)memory(NULL,1,sizeof(double),"MergeData");
			data_all.calc_imag_mt_TE[i] = (double *)memory(NULL,1,sizeof(double),"MergeData");
			data_all.weigth_mt[i] = (double *)memory(NULL,1,sizeof(double),"MergeData");
		}
	}
	

	/*Merge the arrays of the MT data including the frequencies, real and imaginary parts*/
	tmp_ndata_mt = 0;

	for(i=0;i<num_dat;i++)
	{
		for(j=0;j<data[i].ndata_mt;j++)
			for(k=0;k<data[i].nfreq_mt[j];k++)
			{
				data_all.freq_mt[j + tmp_ndata_mt][k] = data[i].freq_mt[j][k];
				data_all.real_mt_TE[j + tmp_ndata_mt][k] = data[i].real_mt_TE[j][k];
				data_all.imag_mt_TE[j + tmp_ndata_mt][k] = data[i].imag_mt_TE[j][k];
				data_all.real_mt_TM[j + tmp_ndata_mt][k] = data[i].real_mt_TM[j][k];
				data_all.imag_mt_TM[j + tmp_ndata_mt][k] = data[i].imag_mt_TM[j][k];

				data_all.weigth_mt[j + tmp_ndata_mt][k] = data[i].weigth_mt[j][k];
			}

				tmp_ndata_mt = data[i].ndata_mt + tmp_ndata_mt;
	}

	printf("All structures from the DAT-files are merged together\n");
	printf(" Number of source-receiver combinations: %d\n", data_all.ndata_seis);
	printf(" Number of gravity measurements: %d\n", data_all.ndata_grav);
	printf(" Number of MT measurements: %d\n", data_all.ndata_mt);
	printf(" Number of ALL station locations: %d\n", nstat_all);
    printf("----------------\n\n");


	return(data_all);
}

/*---------------------------------------------------------------------------*/
/* Optimize DATA_STRUCT: Looking for shots/receivers/gravity and MT stations with the same coordinates and */
/* merge them together afterwards*/
/* parameters:   *data       := Pointer on the data structures that should be merged*/
/*               eps		 := Maximum distance between two coordinates location, */
/*								so that they are assigned to the same shot resp. receiver location */ 
/*				 *geo        := the geometry structure */

#define xs(a) geo->x[(data->shots[a])-1] /*x-coordinate of the shot locations*/
#define ys(a) geo->y[(data->shots[a])-1] /*y-coordinate of the shot locations*/
#define zs(a) geo->z[(data->shots[a])-1] /*z-coordinate of the shot locations*/
#define xr(a) geo->x[(data->recs[a])-1] /*x-coordinate of the receiver locations*/
#define yr(a) geo->y[(data->recs[a])-1] /*y-coordinate of the receiver locations*/
#define zr(a) geo->z[(data->recs[a])-1] /*z-coordinate of the receiver locations*/
#define xg(a) geo->x[(data->gravs[a])-1] /*x-coordinate of the gravity station location*/
#define yg(a) geo->y[(data->gravs[a])-1] /*y-coordinate of the gravity station location*/
#define zg(a) geo->z[(data->gravs[a])-1] /*z-coordinate of the gravity station location*/
#define xm(a) geo->x[(data->mts[a])-1] /*x-coordinate of the MT station location*/
#define ym(a) geo->y[(data->mts[a])-1] /*y-coordinate of the MT station location*/
#define zm(a) geo->z[(data->mts[a])-1] /*z-coordinate of the MT station location*/

int OptiData(DATA_STRUCT *data,GEOMETRY *geo, float eps)
{

	int i,j,k;
	float distS; /*Distances between the source locations*/
	float distR; /*Distances between the receiver locations*/
	float distG; /*Distances between the locations of the gravity stations*/
	float distM; /*Distances between the locations of the MT stations*/
	FILE *out;

	geo->eps = (float)fabs(eps);

    printf("Optimize data structure:\n");
	printf("Start comparing the shot locations:\n");
	printf("----------------\n\n");

	/*Compare shot locations:*/
	for(i=0;i<(geo->nshot);i++)
		for(j=i+1;j<(geo->nshot);j++)
		{
		  distS = (float)sqrt((xs(i)-xs(j))*(xs(i)-xs(j)) + (ys(i)-ys(j))*(ys(i)-ys(j)) + (zs(i)-zs(j))*(zs(i)-zs(j)));
		  if(distS <= eps)
		  { 
			  /*Exchange the shot-index*/
			  for(k=0;k<(data->ndata_seis);k++)
				  if(data->sno[k] == data->shots[j])
						data->sno[k] = data->shots[i];

			    printf("Shot location %d has the coordinates: x=%f, y=%f, z=%f\n",(data->shots[i]),xs(i),ys(i),zs(i));
				printf("Shot location %d has the coordinates: x=%f, y=%f, z=%f\n",(data->shots[j]),xs(j),ys(j),zs(j));
				printf("Distance: %f; Shot location: %d -> %d\n\n",distS,(data->shots[j]),(data->shots[i]));
				
				memmove(&(data->shots[j]),&(data->shots[j+1]),((geo->nshot)-j)*sizeof(int));
				(geo->nshot) = (geo->nshot) - 1;


		  }
		}



	/*Compare receiver locations:*/
   	for(i=0;i<(geo->nrec);i++)
		for(j=i+1;j<(geo->nrec);j++)
		{
		  distR = (float)sqrt((xr(i)-xr(j))*(xr(i)-xr(j)) + (yr(i)-yr(j))*(yr(i)-yr(j)) + (zr(i)-zr(j))*(zr(i)-zr(j)));
		  if(distR <= eps)
		  {
				/*Exchange the receiver-index*/
				for(k=0;k<(data->ndata_seis);k++)
					if(data->rno[k] == data->recs[j])
						data->rno[k] = data->recs[i];

					printf("Rec. location %d has the coordinates: x=%f, y=%f, z=%f\n",(data->recs[i]),xr(i),yr(i),zr(i));
					printf("Rec. location %d has the coordinates: x=%f, y=%f, z=%f\n",(data->recs[j]),xr(j),yr(j),zr(j));
					printf("Distance: %f; Receiver location: %d -> %d\n\n",distR,(data->recs[j]),(data->recs[i]));

					memmove(&(data->recs[j]),&(data->recs[j+1]),((geo->nrec)-j)*sizeof(int));
					(geo->nrec) = (geo->nrec) - 1;

		  }
		}
   

	/*Compare the locations of the gravity stations:*/
	for(i=0;i<(geo->nstat_grav);i++)
		for(j=i+1;j<(geo->nstat_grav);j++)
		{
			distG = (float)sqrt((xg(i)-xg(j))*(xg(i)-xg(j)) + (yg(i)-yg(j))*(yg(i)-yg(j)) + (zg(i)-zg(j))*(zg(i)-zg(j)));
			if(distG <= eps)
			{
				for(k=0;k<(data->ndata_grav);k++)
					if(data->gno[k] == data->gravs[j])
						data->gno[k] = data->gravs[i];

					printf("Gravity station %d has the coordinates: x=%f, y=%f, z=%f\n",(data->gravs[i]),xg(i),yg(i),zg(i));
					printf("Gravity station %d has the coordinates: x=%f, y=%f, z=%f\n",(data->gravs[j]),xg(j),yg(j),zg(j));
					printf("Distance: %f; Location of gravity station: %d -> %d\n\n",distG,(data->gravs[j]),(data->gravs[i]));

					memmove(&(data->gravs[j]),&(data->gravs[j+1]),((geo->nstat_grav)-j)*sizeof(int));
					(geo->nstat_grav) = (geo->nstat_grav) -1;
			}

		}

	/*Compare the locations of the MT stations:*/
	 for(i=0;i<(geo->nstat_mt);i++)
		 for(j=i+1;j<(geo->nstat_mt);j++)
		 {
			 distM = (float)sqrt((xm(i)-xm(j))*(xm(i)-xm(j)) + (ym(i)-ym(j))*(ym(i)-ym(j)) + (zm(i)-zm(j))*(zm(i)-zm(j)));
			 if(distM <= eps)
			 {
				 for(k=0;k<(data->ndata_mt);k++)
					 if(data->mno[k] == data->mts[j])
						 data->mno[k] = data->mts[i];

					printf("MT station %d has the coordinates: x=%f, y=%f, z=%f\n",(data->mts[i]),xm(i),ym(i),zm(i));
					printf("MT station %d has the coordinates: x=%f, y=%f, z=%f\n",(data->mts[j]),xm(j),ym(j),zm(j));
					printf("Distance: %f; Location of gravity station: %d -> %d\n\n",distM,(data->mts[j]),(data->mts[i]));

					memmove(&(data->mts[j]),&(data->mts[j+1]),((geo->nstat_mt)-j)*sizeof(int));
					(geo->nstat_mt) = (geo->nstat_mt) -1;					 
			 }
		 }

		printf("\nThe data structure is optimized\n");
		printf(" Number of shot and receiver locations: %d, %d\n", geo->nshot, geo->nrec);
		printf(" Number of gravity and MT station locations: %d %d\n", geo->nstat_grav, geo->nstat_mt);

	
			/***************************************/
			/* Write out used receiver/shot positions and locations of the gravity and MT data*/
			 out = fopen("station-statistic-opt.txt","wt");
			 fprintf(out,"# shot-positions: no: shots %d recs %d\n",geo->nshot, geo->nrec);
			 for(i=0;i<geo->nshot;i++)
			 {
			   	fprintf(out,"%d %d\n", i+1,data->shots[i]);
			 }

			 fprintf(out,"\n# receiver-positions:\n");
			 for(i=0;i<geo->nrec;i++)
			 {
			   	fprintf(out,"%d %d\n", i+1,data->recs[i]);
			 }

			 fprintf(out,"\n# gravity stations: no: %d\n", geo->nstat_grav);
			 for(i=0;i<geo->nstat_grav;i++)
			 {
				 fprintf(out,"%d %d\n", i+1, data->gravs[i]);
			 }

			 fprintf(out,"\n# MT stations: no: %d\n", geo->nstat_mt);
			 for(i=0;i<geo->nstat_mt;i++)
			 {
				 fprintf(out,"%d %d\n", i+1, data->mts[i]);
			 }

			 fclose(out);
			 /***************************************/



	return(1);
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

/*---------------------------------------------------------------------------*/
/* Calculate the shot-receiver distance for all shot receiver-combinations*/
/* parameters:   *data       := Pointer on the data structure*/ 
/*				 geo        := the geometry structure */

int CalcShotRecDist(DATA_STRUCT *data,GEOMETRY geo)
{
	long i;

		if(data->ndata_seis != 0)
		{
			data->xdist = (float *)memory( NULL,data->ndata_seis,sizeof(float),"MergeData");
			data->ydist = (float *)memory( NULL,data->ndata_seis,sizeof(float),"MergeData");
			data->zdist = (float *)memory( NULL,data->ndata_seis,sizeof(float),"MergeData");
		}
		else
		{
			data->xdist = (float *)memory( NULL,1,sizeof(float),"MergeData");
			data->ydist = (float *)memory( NULL,1,sizeof(float),"MergeData");
			data->zdist = (float *)memory( NULL,1,sizeof(float),"MergeData");
		}

		for(i=0;i<data->ndata_seis;i++)
		{
			data->xdist[i] = geo.x[(data->rno[i])-1] - geo.x[(data->sno[i])-1];
			data->ydist[i] = geo.y[(data->rno[i])-1] - geo.y[(data->sno[i])-1];
			data->zdist[i] = geo.z[(data->rno[i])-1] - geo.z[(data->sno[i])-1];
		}

		return(1);
}

/*---------------------------------------------------------------------------*/
/* Calculate for every shot(receiver) position in DATA_STRUCT the number of active receivers (used shots)*/
/* parameters:   *data       := Pointer on the data structure*/ 
/*				 geo        := the geometry structure */

int StatData(DATA_STRUCT *data,GEOMETRY geo)
{

	int i,j;

	if(geo.nshot != 0)
		data->lshots = (int *)memory(NULL,geo.nshot,sizeof(int),"StatData");
	else
		data->lshots = (int *)memory(NULL,1,sizeof(int),"StatData");

		/*Determine the number of shot positions*/
	  for(i=0; i<geo.nshot; i++)
	  {
		data->lshots[i] = 0;
		for(j=0; j<data->ndata_seis; j++)
		{
			if(data->shots[i] == data->sno[j])
			{
				data->lshots[i]++;	
			}
		}
	  }

	  if(geo.nrec != 0)
		data->lrecs = (int *)memory(NULL,geo.nrec,sizeof(int),"StatData");
	  else
		data->lrecs = (int *)memory(NULL,1,sizeof(int),"StatData");

	/*Determine the number of receiver positions*/
	  for(i=0; i<geo.nrec; i++)
	  {
		data->lrecs[i] = 0;
		for(j=0; j<data->ndata_seis; j++)
		{
			if(data->recs[i] == data->rno[j])
			{
				data->lrecs[i]++;	
			}
		}
	  }

	return(1);
}
