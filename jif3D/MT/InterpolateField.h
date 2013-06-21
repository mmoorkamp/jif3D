//============================================================================
// Name        : InterpolateFields.h
// Author      : 21 Jun 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================


#ifndef INTERPOLATEFIELDS_H_
#define INTERPOLATEFIELDS_H_

#include <complex>
#include <vector>
#include "X3DModel.h"

namespace jif3D {

  std::complex<double> InterpolateField(std::vector<std::complex<double> > &Field,
      const X3DModel &Model, size_t MeasIndex, const std::vector<size_t> &MeasDepthIndices);
}


#endif /* INTERPOLATEFIELDS_H_ */
