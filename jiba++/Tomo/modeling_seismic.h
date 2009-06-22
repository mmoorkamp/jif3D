#ifndef MODELING_SEISMIC_H
#define MODELING_SEISMIC_H


#include "inv3d.h"


/*! \file
 * The file containing all the function declarations */
/*--------------------------------------------------------------*/
/* Define function:*/

namespace jiba
  {
    /** \addtogroup seismod Seismic tomography modelling */
    /* @{ */
    int ForwardModRay(const GEOMETRY &geo, const GRID_STRUCT &grid,
        DATA_STRUCT *data, RP_STRUCT *raypath);
    /*in podvin-lecomte-3D.c*/
    float interpolate(float x, float y, float z, const GRID_STRUCT &grid,
        float *data);
  /* @} */
  }
#endif
