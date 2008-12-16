#include <cutil.h>
#include "gravcuda_kernel.cu"


extern "C" void FreeData(double *d_xcoord, double *d_ycoord, double *d_zcoord,
		double *d_xsize, double *d_ysize, double *d_zsize, double *d_result)
{
	CUT_CHECK_ERROR("Kernel execution failed");
	CUDA_SAFE_CALL(cudaFree(d_xcoord));
	CUDA_SAFE_CALL(cudaFree(d_ycoord));
	CUDA_SAFE_CALL(cudaFree(d_zcoord));
	CUDA_SAFE_CALL(cudaFree(d_xsize));
	CUDA_SAFE_CALL(cudaFree(d_ysize));
	CUDA_SAFE_CALL(cudaFree(d_zsize));
	CUDA_SAFE_CALL(cudaFree(d_result));
}

extern "C" void PrepareData(double *d_xcoord, double *d_ycoord, double *d_zcoord,
		double *d_xsize, double *d_ysize, double *d_zsize, double *d_result,
		const double *xcoord,const  double *ycoord,const  double *zcoord,
		const double *xsize,const  double *ysize,const  double *zsize, unsigned int nx,unsigned int ny,unsigned int nz) {
	CUT_DEVICE_INIT(argc, argv);

	unsigned int xmem_size = sizeof(double) * nx;
	unsigned int ymem_size = sizeof(double) * ny;
	unsigned int zmem_size = sizeof(double) * nz;

	CUDA_SAFE_CALL(cudaMalloc((void**) &d_xcoord, xmem_size));
	CUDA_SAFE_CALL(cudaMalloc((void**) &d_ycoord, ymem_size));
	CUDA_SAFE_CALL(cudaMalloc((void**) &d_zcoord, zmem_size));
	CUDA_SAFE_CALL(cudaMalloc((void**) &d_xsize, xmem_size));
	CUDA_SAFE_CALL(cudaMalloc((void**) &d_ysize, ymem_size));
	CUDA_SAFE_CALL(cudaMalloc((void**) &d_zsize, zmem_size));
	CUDA_SAFE_CALL(cudaMalloc((void**) &d_result, xmem_size * ymem_size
					* zmem_size));

	// copy host memory to device
	CUDA_SAFE_CALL(cudaMemcpy(d_xcoord, xcoord,
					xmem_size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_ycoord, ycoord,
					ymem_size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_zcoord, zcoord,
					zmem_size, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_xsize, xsize, xmem_size,
					cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_ysize, ysize, ymem_size,
					cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_zsize, zsize, zmem_size,
					cudaMemcpyHostToDevice));

}

extern "C" void SingleScalarMeas(const double x_meas, const double y_meas,
		const double z_meas, double *d_xcoord, double *d_ycoord,
		double *d_zcoord, double *d_xsize, double *d_ysize, double *d_zsize,
		double *d_result, const int nx, const int ny, const int nz,
		double *returnvalue) {
	const int BLOCK_SIZE = 16;
	dim3 threads( BLOCK_SIZE);
	int nelements = nx * ny * nz;
	dim3 grid(nelements / threads.x+1);

	CalcScalarMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
			d_xsize, d_ysize, d_zsize, nx, ny,
			nz, d_result);
	cudaMemcpy(returnvalue, d_result, nx * sizeof(double)*ny * sizeof(double)*nz * sizeof(double),
			cudaMemcpyDeviceToHost);
}
