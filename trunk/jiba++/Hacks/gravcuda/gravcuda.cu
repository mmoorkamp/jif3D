/*
 * Copyright 1993-2007 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and
 * international Copyright laws.  Users and possessors of this source code
 * are hereby granted a nonexclusive, royalty-free license to use this code
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL,
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOURCE CODE.
 *
 * U.S. Government End Users.   This source code is a "commercial item" as
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of
 * "commercial computer  software"  and "commercial computer software
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995)
 * and is provided to the U.S. Government only as a commercial end item.
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the
 * source code with only those rights set forth herein.
 *
 * Any use of this source code in individual and commercial software must
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

/* Template project which demonstrates the basics on how to setup a project
 * example application.
 * Host code.
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <numeric>

// includes, project
#include <cutil.h>

extern "C" void CalcScalarMeas(const double x_meas, const double y_meas,
		const double z_meas, const double *XCoord,
		const double *YCoord, const double *ZCoord, const double *XSizes,
		const double *YSizes, const double *ZSizes, const double *Densities,
		const int nx, const int ny, const int nz, double *returnvalue);
// includes, kernels
#include <gravcuda_kernel.cu>

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void runTest(int argc, char** argv);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) {
	runTest(argc, argv);
}

////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void runTest(int argc, char** argv) {

	CUT_DEVICE_INIT(argc, argv);

	unsigned int nx = 100;
	unsigned int ny = 60;
	unsigned int nz = 60;
	unsigned int xmem_size = sizeof(double) * nx;
	unsigned int ymem_size = sizeof(double) * ny;
	unsigned int zmem_size = sizeof(double) * nz;
	unsigned int grid_size = sizeof(double) * nx * ny * nz;
	// allocate host memory
	double* xcoord = (double*) malloc(xmem_size);
	double* ycoord = (double*) malloc(ymem_size);
	double* zcoord = (double*) malloc(zmem_size);
	double* xsize = (double*) malloc(xmem_size);
	double* ysize = (double*) malloc(ymem_size);
	double* zsize = (double*) malloc(zmem_size);
	double* density = (double*) malloc(grid_size);
	double* result = (double*) malloc(grid_size);
	double returnvalue;
	// initalize the memory
	for (unsigned int i = 0; i < nx; ++i) {
		xcoord[i] = (double) i;
		xsize[i] = 1.0;
	}
	for (unsigned int i = 0; i < ny; ++i) {
		ycoord[i] = (double) i;
		ysize[i] = 1.0;
	}
	for (unsigned int i = 0; i < nz; ++i) {
		zcoord[i] = (double) i;
		zsize[i] = 1.0;
	}
	for (unsigned int i = 0; i < nx * ny * nz; ++i)
		density[i] = 1.0;
	double x_meas = -1.5;
	double y_meas = 1.5;
	double z_meas = -1.0;

	returnvalue = 0.0;
	// allocate device memory
	double *d_xcoord, *d_ycoord, *d_zcoord;
	double *d_xsize, *d_ysize, *d_zsize;
	double *d_density, *d_result;

	CUDA_SAFE_CALL(cudaMalloc( (void**) &d_xcoord, xmem_size));
	CUDA_SAFE_CALL(cudaMalloc( (void**) &d_ycoord, ymem_size));
	CUDA_SAFE_CALL(cudaMalloc( (void**) &d_zcoord, zmem_size));
	CUDA_SAFE_CALL(cudaMalloc( (void**) &d_xsize, xmem_size));
	CUDA_SAFE_CALL(cudaMalloc( (void**) &d_ysize, ymem_size));
	CUDA_SAFE_CALL(cudaMalloc( (void**) &d_zsize, zmem_size));
	CUDA_SAFE_CALL(cudaMalloc( (void**) &d_density, grid_size));
	CUDA_SAFE_CALL(cudaMalloc( (void**) &d_result, grid_size));

	// copy host memory to device
	CUDA_SAFE_CALL(cudaMemcpy(d_xcoord, xcoord, xmem_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL(cudaMemcpy(d_ycoord, ycoord, ymem_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL(cudaMemcpy(d_zcoord, zcoord, zmem_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL(cudaMemcpy(d_xsize, xsize, xmem_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL(cudaMemcpy(d_ysize, ysize, ymem_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL(cudaMemcpy(d_zsize, zsize, zmem_size, cudaMemcpyHostToDevice) );
	CUDA_SAFE_CALL(cudaMemcpy(d_density, density, grid_size, cudaMemcpyHostToDevice) );

	// setup execution parameters
	//dim3 grid( 1, 1, 1);
	//dim3 threads( num_threads,
	//		num_threads, num_threads);

	const int BLOCK_SIZE = 16;
	dim3 threads( BLOCK_SIZE);
	int nelements = nx * ny * nz;
	dim3 grid(nelements / threads.x+1);

	CalcScalarMeas<<< grid, threads >>>(x_meas, y_meas, z_meas, d_xcoord, d_ycoord, d_zcoord,
			d_xsize, d_ysize, d_zsize, d_density, nx, ny,
			nz, d_result);
	cudaMemcpy(result, d_result, grid_size,
			cudaMemcpyDeviceToHost);
	for (size_t i = 0; i < nelements; ++i) {
		//printf("Device: %e \n", result[i]);
		returnvalue += result[i];
	}

	CUT_CHECK_ERROR("Kernel execution failed");

	printf("Result: %e \n", returnvalue);
	//double reduceTime = cutGetAverageTimerValue(timer);
	//printf("Average time: %f ms\n", reduceTime);
	// cleanup memory
	free(xcoord);
	free(ycoord);
	free(zcoord);
	free(xsize);
	free(ysize);
	free(zsize);
	free(density);
	CUDA_SAFE_CALL(cudaFree(d_xcoord));
	CUDA_SAFE_CALL(cudaFree(d_ycoord));
	CUDA_SAFE_CALL(cudaFree(d_zcoord));
	CUDA_SAFE_CALL(cudaFree(d_xsize));
	CUDA_SAFE_CALL(cudaFree(d_ysize));
	CUDA_SAFE_CALL(cudaFree(d_zsize));
	CUDA_SAFE_CALL(cudaFree(density));
}
