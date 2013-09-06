
#include <stdio.h>
#include "gravcuda_kernel.cu"
#include "ftgcuda_kernel.cu"
#include "../Global/CudaError.h"

// Free all the memory allocated on the GPU
extern "C"
void  FreeData(double **d_xcoord, double **d_ycoord, double **d_zcoord,
      double **d_xsize, double **d_ysize, double **d_zsize, double **d_result)
    {
      CUDA_SAFE_CALL(cudaFree(*d_xcoord));
      CUDA_SAFE_CALL(cudaFree(*d_ycoord));
      CUDA_SAFE_CALL(cudaFree(*d_zcoord));
      CUDA_SAFE_CALL(cudaFree(*d_xsize));
      CUDA_SAFE_CALL(cudaFree(*d_ysize));
      CUDA_SAFE_CALL(cudaFree(*d_zsize));
      CUDA_SAFE_CALL(cudaFree(*d_result));
      *d_xcoord = NULL;
      *d_ycoord = NULL;
      *d_zcoord = NULL;
      *d_xsize = NULL;
      *d_ysize = NULL;
      *d_zsize = NULL;
      *d_result = NULL;
    }

/* Initialize the GPU, allocate memory on the GPU and transfer coordinate
 * information from CPU to GPU. The variables starting with d_ point to memory
 * on the GPU while the corresponding variables without the prefix point to
 * main memory.
 * 
 * xcoord, ycoord, zcoord: The coordinates of the left, upper front corner for each cell in m
 * xsize, ysize, zsize: The size of each cell in the respective direction in m
 * nx, ny, nz the number of cells in x-direction, y-direction, z-direction, respectively
 * 
 */
extern "C"
void  PrepareData(double **d_xcoord, double **d_ycoord, double **d_zcoord,
      double **d_xsize, double **d_ysize, double **d_zsize, double **d_result,
      const double *xcoord,const double *ycoord,const double *zcoord,
      const double *xsize,const double *ysize,const double *zsize, unsigned int nx,unsigned int ny,unsigned int nz)
    {
      //initialize the device and perform basic checks
      // we always use the first device
      int deviceCount;
      CUDA_SAFE_CALL(cudaGetDeviceCount(&deviceCount));
      if (deviceCount == 0)
        {
          fprintf(stderr, "CUTIL CUDA error: no devices supporting CUDA.\n");
          exit(-1);
        }
      int dev = 0;
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, dev);
      if (deviceProp.major < 1 || deviceProp.major == 9999 )
        {
          fprintf(stderr, "Cutil error: Device does not support CUDA.\n");
          exit(-1);
        }
      int version = deviceProp.major * 10 + deviceProp.minor;
      if (version < 13)
        {
          fprintf(stderr, "Cutil error: Device does not support double precision.\n");
          exit(-1);
        }

      cudaSetDevice(dev);
      cudaGetLastError();
      //determine the memory requirements in bytes
      size_t xmem_size = sizeof(double) * nx;
      size_t ymem_size = sizeof(double) * ny;
      size_t zmem_size = sizeof(double) * nz;
      size_t result_size = sizeof(double) * nx * ny * nz;
      //allocate memory
      CUDA_SAFE_CALL(cudaMalloc((void**) d_xcoord, xmem_size));
      CUDA_SAFE_CALL(cudaMalloc((void**) d_ycoord, ymem_size));
      CUDA_SAFE_CALL(cudaMalloc((void**) d_zcoord, zmem_size));
      CUDA_SAFE_CALL(cudaMalloc((void**) d_xsize, xmem_size));
      CUDA_SAFE_CALL(cudaMalloc((void**) d_ysize, ymem_size));
      CUDA_SAFE_CALL(cudaMalloc((void**) d_zsize, zmem_size));
      CUDA_SAFE_CALL(cudaMalloc((void**) d_result, result_size));
      cudaMemset((void *)*d_result,0.0,result_size);

      // copy host memory to device
      CUDA_SAFE_CALL(cudaMemcpy(*d_xcoord, xcoord,
              xmem_size, cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(*d_ycoord, ycoord,
              ymem_size, cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(*d_zcoord, zcoord,
              zmem_size, cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(*d_xsize, xsize, xmem_size,
              cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(*d_ysize, ysize, ymem_size,
              cudaMemcpyHostToDevice));
      CUDA_SAFE_CALL(cudaMemcpy(*d_zsize, zsize, zmem_size,
              cudaMemcpyHostToDevice));
    }

/* Calculate the sensitivities for scalar gravity data for a single measurement position
 *
 * x_meas, y_meas, z_meas the measurement coordinates in m
 * d_xcoord, d_ycoord, d_zcoord, ..., d_result, the pointers to GPU memory as initialized by
 * PrepareData.
 * returnvalue will contain an array of size nx*ny*nz with the sensitivities in c-storage order
 * BLOCK_SIZE the number of threads per CUDA execution block
 */
extern "C"
void  SingleScalarMeas(const double x_meas, const double y_meas,
      const double z_meas, double *d_xcoord, double *d_ycoord,
      double *d_zcoord, double *d_xsize, double *d_ysize, double *d_zsize,
      double *d_result, const unsigned int nx, const unsigned int ny, const unsigned int nz,
      double *returnvalue, const unsigned int BLOCK_SIZE)
    {
      dim3 threads( BLOCK_SIZE);
      const unsigned int nelements = nx * ny * nz;
      const unsigned int nblocks = (nelements/threads.x)+1;
      dim3 grid(nblocks);
      CalcScalarMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
          d_xsize, d_ysize, d_zsize, nx, ny,
          nz, d_result);
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue, d_result, nelements * sizeof(double),
              cudaMemcpyDeviceToHost));
    }

/* Calculate the sensitivities for scalar gravity data for a single measurement position
 *
 * x_meas, y_meas, z_meas the measurement coordinates in m
 * d_xcoord, d_ycoord, d_zcoord, ..., d_result, the pointers to GPU memory as initialized by
 * PrepareData.
 * returnvalue will contain an array of size 6*nx*ny*nz with the sensitivities for ftg, with
 * the sensitivities for Uxx first in c-storage order, then Uxy, Uxz, Uyy, Uyz, Uzz
 * BLOCK_SIZE the number of threads per CUDA execution block
 */

extern "C"
void  SingleFTGMeas(const double x_meas, const double y_meas,
      const double z_meas, double *d_xcoord, double *d_ycoord,
      double *d_zcoord, double *d_xsize, double *d_ysize, double *d_zsize,
      double *d_result, const unsigned int nx, const unsigned int ny, const unsigned int nz,
      double *returnvalue, const unsigned int BLOCK_SIZE)
    {
      dim3 threads( BLOCK_SIZE);
      const unsigned int nelements = nx * ny * nz;
      const unsigned int nblocks = (nelements/threads.x)+1;
      const unsigned int mem_elements = nelements * sizeof(double);
      dim3 grid(nblocks);
      CalcUxxMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
          d_xsize, d_ysize, d_zsize, nx, ny,
          nz, d_result);
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      CalcUxyMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
          d_xsize, d_ysize, d_zsize, nx, ny,
          nz, d_result);
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      //now  Uxz
      CalcUxzMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
          d_xsize, d_ysize, d_zsize, nx, ny,
          nz, d_result);
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + 2*nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      //now Uyy
      CalcUyyMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
                d_xsize, d_ysize, d_zsize, nx, ny,
                nz, d_result);
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + 3* nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      //now  Uyz
      CalcUyzMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
                d_xsize, d_ysize, d_zsize, nx, ny,
                nz, d_result);
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + 4*nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      //and finally Uzz
      CalcUzzMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
                d_xsize, d_ysize, d_zsize, nx, ny,
                nz, d_result);
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + 5* nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

    }
