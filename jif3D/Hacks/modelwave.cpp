//============================================================================
// Name        : modelwave.cpp
// Author      : 5 Mar 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================

#include <iostream>
#include <algorithm>
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Global/FileUtil.h"
#include "../Global/Wavelet.h"
#include "../Global/VecMat.h"
#include "../ModelTransforms/WaveletModelTransform.h"

int main()
  {
    jif3D::ThreeDSeismicModel Model;
    std::string InfileName = jif3D::AskFilename("Input Model: ");
    Model.ReadNetCDF(InfileName);

    const size_t nx = Model.GetXCellSizes().size();
    const size_t ny = Model.GetYCellSizes().size();
    const size_t nz = Model.GetZCellSizes().size();

    jif3D::WaveletTransform(Model.SetSlownesses());

    jif3D::ThreeDSeismicModel WaveHigh(Model);
    jif3D::ThreeDSeismicModel WaveLow(Model);

    for (size_t i = 3 * nx / 4; i < nx; ++i)
      {
        for (size_t j = 3 * ny / 4; j < ny; ++j)
          {
            for (size_t k = 3 * nz / 4; k < nz; ++k)
              {
                WaveHigh.SetData()[i][j][k] = 0;
              }
          }
      }

    for (size_t i = 0; i < nx / 4; ++i)
      {
        for (size_t j = 0; j < ny / 4; ++j)
          {
            for (size_t k = 0; k < nz / 4; ++k)
              {
                WaveLow.SetData()[i][j][k] = 0;
              }
          }
      }

    jif3D::InvWaveletTransform(WaveHigh.SetSlownesses());
    jif3D::InvWaveletTransform(WaveLow.SetSlownesses());
    WaveHigh.WriteNetCDF(InfileName + ".waveh.nc");
    WaveLow.WriteNetCDF(InfileName + ".wavel.nc");
    Model.WriteNetCDF(InfileName + ".wcoeff.nc");

    jif3D::rvec Elements(Model.GetSlownesses().num_elements());
    std::copy(Model.GetSlownesses().origin(),
        Model.GetSlownesses().origin() + Model.GetSlownesses().num_elements(),
        Elements.begin());
    const size_t nElem = Model.GetSlownesses().num_elements() / 1000;
    std::nth_element(Elements.begin(), Elements.begin() + nElem, Elements.end(),
        std::greater<double>());
    std::cout << "The highest elements are: ";
    std::copy(Elements.begin(), Elements.begin() + nElem,
        std::ostream_iterator<double>(std::cout, " "));
    std::cout << std::endl;
    for (size_t i = 0; i < nx; ++i)
      {
        for (size_t j = 0; j < ny; ++j)
          {
            for (size_t k = 0; k < nz; ++k)
              {
                if (Model.GetData()[i][j][k] >= Elements(nElem - 1))
                  {
                    std::cout << i << " " << j << " " << k << " "
                        << Model.GetData()[i][j][k] << std::endl;
                  }

              }
          }
      }
  }

