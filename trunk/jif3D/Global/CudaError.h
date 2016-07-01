//============================================================================
// Name        : CudaError.h
// Author      : 6 Sep 2013
// Version     : 
// Copyright   : 2013, mm489
//============================================================================


#ifndef CUDAERROR_H_
#define CUDAERROR_H_


#define CUDA_SAFE_CALL(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#endif /* CUDAERROR_H_ */
