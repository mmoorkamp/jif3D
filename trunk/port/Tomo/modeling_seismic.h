#ifndef MODELING_SEISMIC_H
#define MODELING_SEISMIC_H


#include <vector>

#include "../Global/Jif3DGlobal.h"
#include "tomo_types.h"

/*! \file
 * This file contains some data structures and function definitions for tomography modeling.
 * These are based on the code by B. Heincke with some modifications to avoid c-style allocations
 * and screen output.
 *  */
/*--------------------------------------------------------------*/
/* Define function:*/

namespace jif3D
  {
    /** \addtogroup tomo Seismic tomography classes and functions */
    /* @{ */

    //! The basic forward modeling function, uses the Podvin-Lecomte algorithm to calculate traveltimes and then constructs rays from these times
    J3DEXPORT int ForwardModRay(const GEOMETRY &geo, const GRID_STRUCT &grid, DATA_STRUCT &data,
        std::vector<RP_STRUCT> &raypath);
    /*in podvin-lecomte-3D.c*/
    J3DEXPORT float interpolate(float x, float y, float z, const GRID_STRUCT &grid, float *data);

  /* @} */
  }



J3DEXPORT jif3D::RayResult ForwardModShot(const jif3D::GEOMETRY &geo,
    const jif3D::GRID_STRUCT &grid);





#endif
