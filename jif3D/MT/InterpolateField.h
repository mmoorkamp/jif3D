//============================================================================
// Name        : InterpolateFields.h
// Author      : 21 Jun 2013
// Version     :
// Copyright   : 2013, mm489
//============================================================================

#ifndef INTERPOLATEFIELDS_H_
#define INTERPOLATEFIELDS_H_

#include "../Global/Jif3DGlobal.h"
#include "X3DModel.h"
#include <complex>
#include <vector>

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! Interpolate magnetic and electric fields horizontally to allow for arbitrary horizontal site positions
    J3DEXPORT std::complex<double> InterpolateField(std::vector<std::complex<double> > &Field,
        const X3DModel &Model, size_t MeasIndex,
        const std::vector<size_t> &MeasDepthIndices);
  /* @} */
  }

#endif /* INTERPOLATEFIELDS_H_ */
