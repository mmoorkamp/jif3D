/*
 * CudaMagneticSusceptibilityImp.cpp
 *
 *  Created on: Jul 19, 2022
 *      Author: max
 */

#include "CudaMagneticSusceptibilityImp.h"
#include "../Gravity/BasicGravElements.h"
#include <numeric>
#include <cuda_runtime.h>
#include <boost/math/constants/constants.hpp>

// Calculate the effect of a single measurement, it returns the gridded part of the sensitivities
extern "C" void SingleFTGMeas(const double x_meas, const double y_meas,
    const double z_meas, double *d_xcoord, double *d_ycoord, double *d_zcoord,
    double *d_xsize, double *d_ysize, double *d_zsize, double *d_result,
    const unsigned int nx, const unsigned int ny, const unsigned int nz,
    double *returnvalue, const unsigned int BLOCK_SIZE);
//Perform the allocation of arrays on the GPU and copy the values
extern "C" void PrepareData(double **d_xcoord, double **d_ycoord, double **d_zcoord,
    double **d_xsize, double **d_ysize, double **d_zsize, double **d_result,
    const double *xcoord, const double *ycoord, const double *zcoord, const double *xsize,
    const double *ysize, const double *zsize, unsigned int nx, unsigned int ny,
    unsigned int nz);
// Free the allocated data on the GPU
extern "C" void FreeData(double **d_xcoord, double **d_ycoord, double **d_zcoord,
    double **d_xsize, double **d_ysize, double **d_zsize, double **d_result);

//the default block size for a CUDA execution block
const size_t defaultblocksize = 128;

namespace jif3D
  {

  rvec CudaMagneticSusceptibilityImp::Calculate(const ThreeDSusceptibilityModel &Model, const TotalFieldMagneticData &Data,
      ThreeDGravMagCalculator<TotalFieldMagneticData> &Calculator)
    {

      const unsigned int nx = Model.GetSusceptibilities().shape()[0];
      const unsigned int ny = Model.GetSusceptibilities().shape()[1];
      const unsigned int nz = Model.GetSusceptibilities().shape()[2];
      //allocate memory
      PrepareData(&d_xcoord, &d_ycoord, &d_zcoord, &d_xsize, &d_ysize,
          &d_zsize, &d_result, Model.GetXCoordinates().data(),
          Model.GetYCoordinates().data(), Model.GetZCoordinates().data(),
          Model.GetXCellSizes().data(), Model.GetYCellSizes().data(),
          Model.GetZCellSizes().data(), nx, ny, nz);
      // call the base class that coordinates the calculation of gridded and background parts
      rvec result(ThreeDGravMagImplementation<TotalFieldMagneticData>::Calculate(Model, Data, Calculator));
      // free memory
      FreeData(&d_xcoord, &d_ycoord, &d_zcoord, &d_xsize, &d_ysize, &d_zsize,
          &d_result);
      return result;
    }


    rvec CudaMagneticSusceptibilityImp::CalcGridded(const size_t measindex,
        const ThreeDSusceptibilityModel &Model, const TotalFieldMagneticData &Data,
        rmat &Sensitivities)
      {
        //first we define some constants and abbreviations
        const double x_meas = Data.GetMeasPosX()[measindex];
        const double y_meas = Data.GetMeasPosY()[measindex];
        const double z_meas = Data.GetMeasPosZ()[measindex];

        //calculate the directional angles
        const double BxComp = cos(GetInclination()) * cos(GetDeclination());
        const double ByComp = cos(GetInclination()) * sin(GetDeclination());
        const double BzComp = sin(GetInclination());
        //the conversion factor to go from tensor gravity sensitivities to magnetics
        const double factor = 1.0
            / (4 * boost::math::constants::pi<double>() * jif3D::Grav_const);

        const size_t nbglayers = Model.GetBackgroundThicknesses().size();
        const size_t ngrid = Model.GetNModelElements();
        //the CUDA routines return the 6 independent elements of the FTG tensor (we consider Uzz independent)
        const size_t nreturnelements = 6;
        // we determine whether there are enough elements in the sensitivity matrix
        const bool storesens = (Sensitivities.size1() >= RawDataPerMeasurement())
            && (Sensitivities.size2() >= ngrid + nbglayers);
        //if currsens does not have sufficient size
        if (currsenssize != ngrid)
          {
            //check whether it has been allocated before
            if (currsens != NULL)
              {
                cudaFreeHost(currsens);
              }
            //allocate memory and store how much
            cudaMallocHost((void**) &currsens, ngrid * sizeof(double) * nreturnelements);
            currsenssize = ngrid;
          }


        std::fill_n(currsens, ngrid * nreturnelements, 0.0);
        //This call goes into the GPU, implementation in gravcuda.cu
        SingleFTGMeas(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord, d_xsize,
            d_ysize, d_zsize, d_result, Model.GetSusceptibilities().shape()[0],
            Model.GetSusceptibilities().shape()[1], Model.GetSusceptibilities().shape()[2], currsens,
            blocksize);

        double Bx = 0.0, By = 0.0, Bz = 0.0;
        const double fs = GetFieldStrength();
#pragma omp parallel default(shared)  reduction(+:Bx,By,Bz)
          {
#pragma omp for
            for (int offset = 0; offset < ngrid; ++offset)
              {

                //project the tensor sensitivities to magnetic field components
                const double BxSens = (currsens[offset] * BxComp
                    + currsens[offset + ngrid] * ByComp
                    + currsens[offset + 2 * ngrid] * BzComp) * fs * factor;
                const double BySens = (currsens[offset + ngrid] * BxComp
                    + currsens[offset + 3 * ngrid] * ByComp
                    + currsens[offset + 4 * ngrid] * BzComp) * fs * factor;
                const double BzSens =
                    (currsens[offset + 2 * ngrid] * BxComp
                        + currsens[offset + 4 * ngrid] * ByComp
                        + currsens[offset + 5 * ngrid]) * fs * factor;
                //times the Susceptibility
                const double Susceptibility = 1000.0;
                    //Model.GetSusceptibilities().origin()[offset];
                Bx += BxSens * Susceptibility;
                By += BySens * Susceptibility;
                Bz += BzSens * Susceptibility;
                //if we need the sensitivity values for each cell in the inversion
                //we store it in the appropriate component
                if (storesens)
                  {
                    Sensitivities(0, offset) = BxSens;
                    Sensitivities(1, offset) = BySens;
                    Sensitivities(2, offset) = BzSens;
                  }
              }
          }

        rvec returnvalue(RawDataPerMeasurement());
        returnvalue(0) = Bx;
        returnvalue(1) = By;
        returnvalue(2) = Bz;

        return returnvalue;
      }

    CudaMagneticSusceptibilityImp::CudaMagneticSusceptibilityImp(double Inc, double Dec,
        double Fs) :
        OMPMagneticSusceptibilityImp(Inc, Dec, Fs), d_xcoord(NULL), d_ycoord(NULL), d_zcoord(
        NULL), d_xsize(NULL), d_ysize(
        NULL), d_zsize(NULL), d_result(NULL), currsens(NULL), currsenssize(0), blocksize(
            defaultblocksize)
      {
        // TODO Auto-generated constructor stub

      }

    CudaMagneticSusceptibilityImp::~CudaMagneticSusceptibilityImp()
      {
        //if we allocated memory we have to free it
        if (currsens != NULL)
          {
            cudaFreeHost(currsens);
          }
        cudaDeviceReset();
        currsens = NULL;
      }

  } /* namespace jif3D */
