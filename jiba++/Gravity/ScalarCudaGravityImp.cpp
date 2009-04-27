//============================================================================
// Name        : ScalarCudaGravityImp.cpp
// Author      : Dec 10, 2008
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#include <numeric>
#include <cuda_runtime.h>
#include "ScalarCudaGravityImp.h"
#include "ThreeDGravityCalculator.h"
#include "BasicGravElements.h"

namespace jiba
  {
    // These are forward calculations for the functions declared in gravcuda.cu that
    // need to be called from here, they provide the interface to CUDA and purely
    // operate with basic C types

    // Calculate the effect of a single measurement, it returns the gridded part of the sensitivities
    extern "C"
void      SingleScalarMeas(const double x_meas, const double y_meas,
          const double z_meas, double *d_xcoord, double *d_ycoord,
          double *d_zcoord, double *d_xsize, double *d_ysize, double *d_zsize,
          double *d_result, const unsigned int nx, const unsigned int ny, const unsigned int nz,
          double *returnvalue);
      //Perform the allocation of arrays on the GPU and copy the values
      extern "C" void PrepareData(double **d_xcoord, double **d_ycoord, double **d_zcoord,
          double **d_xsize, double **d_ysize, double **d_zsize, double **d_result,
          const double *xcoord,const double *ycoord,const double *zcoord,
          const double *xsize,const double *ysize,const double *zsize, unsigned int nx,unsigned int ny,unsigned int nz);
      // Free the allocated data on the GPU
      extern "C" void FreeData(double **d_xcoord, double **d_ycoord, double **d_zcoord,
          double **d_xsize, double **d_ysize, double **d_zsize, double **d_result);

      ScalarCudaGravityImp::ScalarCudaGravityImp()
        {
          // we have to do some raw pointer operations for handling sensitivities with CUDA
          currsens = NULL;
          currsenssize = 0;
        }

      ScalarCudaGravityImp::~ScalarCudaGravityImp()
        {
          //if we allocated memory we have to free it
          if (currsens != NULL)
            {
              cudaFreeHost( currsens );
              //delete []currsens;
            }
          cudaThreadExit();
          currsens = NULL;
        }
      /*! Calculate the effect of the background layers for a single measurement
       * @param measindex The index of the measurement in the Model
       * @param xwidth The total width of the gridded domain in x-direction in m
       * @param ywidth The total width of the gridded domain in y-direction in m
       * @param zwidth The total width of the gridded domain in z-direction in m
       * @param Model The Gravity model
       * @param Sensitivities If the matrix passed here holds \f$ 1 \times nbg+ngrid \f$ or more elements, store sensitivity information in the right fields
       * @return A vector with a single component that contains the gravitational effect of the background
       */
      rvec ScalarCudaGravityImp::CalcBackground(const size_t measindex, const double xwidth, const double ywidth,
          const double zwidth, const ThreeDGravityModel &Model,
          rmat &Sensitivities)
        {
          const double x_meas = Model.GetMeasPosX()[measindex];
          const double y_meas = Model.GetMeasPosY()[measindex];
          const double z_meas = Model.GetMeasPosZ()[measindex];
          const size_t nbglayers = Model.GetBackgroundThicknesses().size();
          double result = 0.0;
          double currtop = 0.0;
          double currvalue = 0.0;
          double currbottom = 0.0;
          const size_t modelsize = Model.GetDensities().num_elements();
          const bool storesens = (Sensitivities.size1() >= ndatapermeas)
          && (Sensitivities.size2() >= modelsize + nbglayers);
          // for all layers of the background
          for (size_t j = 0; j < nbglayers; ++j)
            {
              const double currthick = Model.GetBackgroundThicknesses()[j];
              currbottom = currtop + currthick;
              // first assume an infinite sheet for the current layer
              currvalue = CalcInfSheetTerm(z_meas,currtop,currbottom);
              // and then subtract the value for the modelling domain, as this is already calculated in the discretized routine
              // if the background layer complete coincides with the discretized area
              if (currtop < zwidth && (currbottom <= zwidth))

                {
                  currvalue -= CalcGravBoxTerm(x_meas, y_meas, z_meas, 0.0, 0.0,
                      currtop, xwidth, ywidth, currthick);
                }
              //if some of the background coincides and some is below
              if (currtop < zwidth && currbottom> zwidth)

                {
                  currvalue -= CalcGravBoxTerm(x_meas, y_meas, z_meas, 0.0, 0.0,
                      currtop, xwidth, ywidth, (zwidth - currtop));
                }
              if (storesens)
                {
                  Sensitivities(0, modelsize + j) = currvalue;
                }
              result += currvalue * Model.GetBackgroundDensities()[j];
              currtop += currthick;
            }
          rvec returnvector(1);
          returnvector(0) = result;
          return returnvector;
        }

      /*! Calculate the effect of the gridded domain on the GPU for a single measurement
       * @param measindex The index of the measurement
       * @param Model The gravity model
       * @param Sensitivities If this parameters has \f$ 1 \times nbg+ngrid \f$ or more elements, store sensitivity information
       * @return A single component vector with the accelerational effect of the gridded domain
       */
      rvec ScalarCudaGravityImp::CalcGridded(const size_t measindex, const ThreeDGravityModel &Model,
          rmat &Sensitivities)
        {
          //first we define some constants and abbreviations
          const double x_meas = Model.GetMeasPosX()[measindex];
          const double y_meas = Model.GetMeasPosY()[measindex];
          const double z_meas = Model.GetMeasPosZ()[measindex];
          const size_t nbglayers = Model.GetBackgroundThicknesses().size();
          const size_t ngrid = Model.GetDensities().num_elements();
          // we determine whether there are enough elements in the sensitivity matrix
          const bool storesens = (Sensitivities.size1() >= ndatapermeas)
          && (Sensitivities.size2() >= ngrid + nbglayers);
          //if currsens does not have sufficient size
          if (currsenssize != ngrid)
            {
              //check whether it has been allocated before
              if (currsens != NULL)
                {
                  cudaFreeHost( currsens );
                }
              //allocate memory and store how much
              cudaMallocHost( (void**) &currsens, ngrid *sizeof(double) );
              currsenssize = ngrid;
            }
          std::fill_n(currsens,currsenssize,0.0);
          //This call goes into the GPU, implementation in gravcuda.cu
          SingleScalarMeas(x_meas,y_meas,z_meas,d_xcoord,d_ycoord,d_zcoord,d_xsize,d_ysize,d_zsize,d_result,
              Model.GetDensities().shape()[0],Model.GetDensities().shape()[1],Model.GetDensities().shape()[2],currsens);
          rvec result(ndatapermeas);
          //the GPU only calculates the sensitivities, we calculate the acceleration with the densities
          result(0) = std::inner_product(currsens,
              currsens+ngrid,Model.GetDensities().origin(),0.0);
          if (storesens)
            {
              std::copy(currsens,currsens+ngrid,Sensitivities.data().begin());
            }

          return result;
        }
      /*! The reimplementation for CUDA mainly manages the allocation of memory on the GPU
       * @param Model The gravity model
       * @param Calculator The calculator object
       * @return The vector holding the gravitational acceleration at each measurement site
       */
      rvec ScalarCudaGravityImp::Calculate(const ThreeDGravityModel &Model,ThreeDGravityCalculator &Calculator)
        {

          const unsigned int nx = Model.GetDensities().shape()[0];
          const unsigned int ny = Model.GetDensities().shape()[1];
          const unsigned int nz = Model.GetDensities().shape()[2];
          CheckTransform();
          //allocate memory
          PrepareData(&d_xcoord,&d_ycoord,&d_zcoord,&d_xsize,&d_ysize,&d_zsize,&d_result,Model.GetXCoordinates().data(),Model.GetYCoordinates().data(),
              Model.GetZCoordinates().data(),Model.GetXCellSizes().data(),Model.GetYCellSizes().data(),Model.GetZCellSizes().data(),nx,ny,nz);
          // call the base class that coordinates the calculation of gridded and background parts
          rvec result(ThreeDGravityImplementation::Calculate(Model,Calculator));
          // free memory
          FreeData(&d_xcoord,&d_ycoord,&d_zcoord,&d_xsize,&d_ysize,&d_zsize,&d_result);
          return result;
        }

    }
