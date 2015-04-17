#ifndef MODELING_SEISMIC_H
#define MODELING_SEISMIC_H

#ifdef HAVEHPX
#include <hpx/config.hpp>
#include <hpx/include/actions.hpp>
#endif
#include <vector>
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
    int ForwardModRay(const GEOMETRY &geo, const GRID_STRUCT &grid, DATA_STRUCT &data,
        std::vector<RP_STRUCT> &raypath);
    /*in podvin-lecomte-3D.c*/
    float interpolate(float x, float y, float z, const GRID_STRUCT &grid, float *data);

  /* @} */
  }



jif3D::RayResult ForwardModShot(const jif3D::GEOMETRY &geo,
    const jif3D::GRID_STRUCT &grid);


inline jif3D::RayResult ForwardModShotArray(const std::vector<jif3D::GEOMETRY> &geo,
    const jif3D::GRID_STRUCT &grid)
  {
    jif3D::RayResult Result;
    for (size_t i = 0; i < geo.size(); ++i)
      {
        jif3D::RayResult CurrResult = ForwardModShot(geo[i], grid);
        std::copy(CurrResult.raypath.begin(),CurrResult.raypath.end(),back_inserter(Result.raypath));
        std::copy(CurrResult.tcalc.begin(),CurrResult.tcalc.end(),back_inserter(Result.tcalc));
      }
    return Result;
  }

#ifdef HAVEHPX
HPX_DEFINE_PLAIN_ACTION(ForwardModShotArray, ForwardModShotArray_action);
#endif

#endif
