#ifndef MODELING_SEISMIC_H
#define MODELING_SEISMIC_H
#include "../Global/convert.h"
#include "inv3d.h"
#include <string>
#include "Podvin.h"
/*! \file
 * The file containing all the function declarations */
/*--------------------------------------------------------------*/
/* Define function:*/

namespace jiba
  {
/*in file toolbox2.c:*/
/*! Allocate new memory with additional safeguards and error information in case of failure */
char *memory (char *prev_addr, size_t n, size_t size, std::string progname)
  {
    //if there is nothing to allocate
    if (n == 0)
      {
        std::cerr << "\n\n!!!WARNING!!!: In " << progname << ", no memory was allocated, n = " << n << std::endl << std::endl;
        return((char *) NULL); /* Take care of n = 0 */
      }

    char *tmp;
    if (prev_addr) // if we have to reallocate
      {
        if ((tmp = (char *)realloc ( prev_addr, n * size)) == NULL)
          {
            std::string ErrorString = "Fatal Error: " + progname + " could not reallocate more memory, n = " + stringify(n);
            throw std::runtime_error(ErrorString);

          }
      }
    else //if we allocate new memory
      {
        if ((tmp = (char *)calloc ( n, size)) == NULL)
          {
            std::string ErrorString = "Fatal Error: " + progname + " could not allocate more memory, n = " + stringify(n);
            throw std::runtime_error(ErrorString);
          }
      }
    return (tmp);
  }
//int GetNext(char *line,FILE *inf);


//int WriteModSeisOut(int ninv, DATA_STRUCT data, GEOMETRY geo, GRID_STRUCT grid);
//int WriteRayOut(RP_STRUCT *rp, DATA_STRUCT data, GRID_STRUCT grid);
//int WriteFatRayOut(F_RP_STRUCT *fr,DATA_STRUCT data,GRID_STRUCT grid);
//int WriteTravelTimeModOut(float *tt, int loc,int nx,int ny,int nz);
//float *ReadTravelTimeModIn(char *fname);
//int WriteRayDenseOut(GRID_STRUCT grid, int ninv, long ninv_cell, long *nrays, double *length, BIDX_STRUCT *inv_cell);

/** \addtogroup seismod Seismic tomography modelling */
/* @{ */
int ForwardModRay(GEOMETRY geo,GRID_STRUCT grid,DATA_STRUCT *data, RP_STRUCT *raypath, time_t time_start);
/*in podvin-lecomte-3D.c*/

/* @} */
  }
#endif
