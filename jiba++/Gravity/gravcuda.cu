#include <cutil.h>
#include <cutil_inline.h>
#include <stdio.h>
#include "gravcuda_kernel.cu"
#include "ftgcuda_kernel.cu"
extern "C"
void  FreeData(double **d_xcoord, double **d_ycoord, double **d_zcoord,
      double **d_xsize, double **d_ysize, double **d_zsize, double **d_result)
    {
      CUT_CHECK_ERROR("Kernel execution failed");
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

extern "C"
void  PrepareData(double **d_xcoord, double **d_ycoord, double **d_zcoord,
      double **d_xsize, double **d_ysize, double **d_zsize, double **d_result,
      const double *xcoord,const double *ycoord,const double *zcoord,
      const double *xsize,const double *ysize,const double *zsize, unsigned int nx,unsigned int ny,unsigned int nz)
    {

      int deviceCount;
      cutilSafeCallNoSync(cudaGetDeviceCount(&deviceCount));
      if (deviceCount == 0)
        {
          fprintf(stderr, "CUTIL CUDA error: no devices supporting CUDA.\n");
          exit(-1);
        }
      int dev = 0;
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, dev);
      //printf("\nThe Properties of the Device with ID %d are\n",dev);
      //printf("\tDevice Name : %s",deviceProp.name);
      //printf("\n\tDevice Memory Size (in bytes) : %u",deviceProp.bytes);
      //printf("\n\tDevice Major Revision Numbers : %d",deviceProp.major);
      //printf("\n\tDevice Minor Revision Numbers : %d",deviceProp.minor);
      if (deviceProp.major < 1)
        {
          fprintf(stderr, "cutil error: device does not support CUDA.\n");
          exit(-1);
        }

      cudaSetDevice(dev);
      cudaGetLastError();

      size_t xmem_size = sizeof(double) * nx;
      size_t ymem_size = sizeof(double) * ny;
      size_t zmem_size = sizeof(double) * nz;
      size_t result_size = sizeof(double) * nx * ny * nz;

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
      cutilCheckMsg("Allocation and memory copy failed");
    }

extern "C"
void  SingleScalarMeas(const double x_meas, const double y_meas,
      const double z_meas, double *d_xcoord, double *d_ycoord,
      double *d_zcoord, double *d_xsize, double *d_ysize, double *d_zsize,
      double *d_result, const unsigned int nx, const unsigned int ny, const unsigned int nz,
      double *returnvalue)
    {
      const int BLOCK_SIZE = 128;
      dim3 threads( BLOCK_SIZE);
      const unsigned int nelements = nx * ny * nz;
      const unsigned int nblocks = (nelements/threads.x)+1;
      dim3 grid(nblocks);
      CalcScalarMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
          d_xsize, d_ysize, d_zsize, nx, ny,
          nz, d_result);
      cutilCheckMsg("Kernel execution failed");
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue, d_result, nelements * sizeof(double),
              cudaMemcpyDeviceToHost));

      CUT_CHECK_ERROR("Result copy failed");
    }

extern "C"
void  SingleFTGMeas(const double x_meas, const double y_meas,
      const double z_meas, double *d_xcoord, double *d_ycoord,
      double *d_zcoord, double *d_xsize, double *d_ysize, double *d_zsize,
      double *d_result, const unsigned int nx, const unsigned int ny, const unsigned int nz,
      double *returnvalue)
    {
      const int BLOCK_SIZE = 128;
      dim3 threads( BLOCK_SIZE);
      const unsigned int nelements = nx * ny * nz;
      const unsigned int nblocks = (nelements/threads.x)+1;
      const unsigned int mem_elements = nelements * sizeof(double);
      dim3 grid(nblocks);
      CalcUxxMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
          d_xsize, d_ysize, d_zsize, nx, ny,
          nz, d_result);
      cutilCheckMsg("Kernel execution failed for Uxx");
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      CalcUxyMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
          d_xsize, d_ysize, d_zsize, nx, ny,
          nz, d_result);
      cutilCheckMsg("Kernel execution failed for Uxy");
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      //now  Uxz
      CalcUxzMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
          d_xsize, d_ysize, d_zsize, nx, ny,
          nz, d_result);
      cutilCheckMsg("Kernel execution failed for Uxz");
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + 2*nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      //now Uyy
      CalcUyyMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
                d_xsize, d_ysize, d_zsize, nx, ny,
                nz, d_result);
      cutilCheckMsg("Kernel execution failed for Uyy");
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + 3* nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      //now  Uyz
      CalcUyzMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
                d_xsize, d_ysize, d_zsize, nx, ny,
                nz, d_result);
      cutilCheckMsg("Kernel execution failed for Uyz");
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + 4*nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      //and finally Uzz
      CalcUzzMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
                d_xsize, d_ysize, d_zsize, nx, ny,
                nz, d_result);
      cutilCheckMsg("Kernel execution failed for Uyy");
      CUDA_SAFE_CALL(cudaMemcpy(returnvalue + 5* nelements, d_result, mem_elements,
              cudaMemcpyDeviceToHost));

      CUT_CHECK_ERROR("Result copy failed");
    }
