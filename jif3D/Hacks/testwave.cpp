//============================================================================
// Name        : testwave.cpp
// Author      : 7 Mar 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include "../Global/Wavelet.h"
#include "../Global/VecMat.h"
#include "../ModelTransforms/WaveletModelTransform.h"
#include <gsl/gsl_wavelet2d.h>
#include <boost/multi_array.hpp>
#include <cstdlib>
#include <algorithm>
#include <iostream>

int main()
  {
    const size_t n = 16;
    boost::multi_array<double, 2> Data(boost::extents[n][ n]);
    std::generate_n(Data.origin(), n * n, drand48);

    gsl_wavelet *w;
    gsl_wavelet_workspace *work;

    w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 4);
    work = gsl_wavelet_workspace_alloc(n * n);

    double *gsldata = (double *)malloc(n * n * sizeof(double));
    std::copy(Data.origin(), Data.origin() + Data.num_elements(), gsldata);
    //gsl_wavelet_transform_forward(w, gsldata, 1, n*n, work);
    gsl_wavelet2d_transform_forward(w, gsldata, n, n, n, work);

    jif3D::WaveletTransform(Data);

    for (size_t i = 0; i < n * n; ++i)
      {
        std::cout << i << "  " << Data.data()[i] << " " << gsldata[i] << std::endl;
      }

  }
