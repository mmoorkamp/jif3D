#ifndef _TEMPLATE_KERNEL_H_
#define _TEMPLATE_KERNEL_H_

#define Grav_const 6.67428e-8

__device__ void OffsetToIndex(const int offset,const int ny,const int nz, int &xi, int &yi, int &zi)
{
	zi = offset % nz;
	xi = (offset - zi) / nz;
	yi = xi % ny;
	xi = (xi-yi) / ny;
}

__device__ float CalcGravTerm(const float x, const float y, const float z) {
	//if the distance r between the measurement point and one of the border points is very small
	//the log(...) expressions become undefined, so we shift the measurement point by a tiny amount
	const float r = sqrt(x * x + y * y + z * z);
	float rvalue = x * log(y + r) + y * log(x + r);
	// atan2 takes care of small denominators
	rvalue += z * atan2(z * r, x * y);

	return rvalue;
}

__device__ float CalcGravBoxTerm(const float meas_x, const float meas_y,
		const float meas_z, const float ul_corner_x, const float ul_corner_y,
		const float ul_corner_z, const float x_size, const float y_size,
		const float z_size) {
	//we need 8 terms, one for each corner of the prism, they have varying signs
	float returnvalue = -CalcGravTerm(meas_x - ul_corner_x, meas_y
			- ul_corner_y, meas_z - ul_corner_z);
	// one corner shifted -> positive sign
	returnvalue += CalcGravTerm(meas_x - (ul_corner_x + x_size), meas_y
			- ul_corner_y, meas_z - ul_corner_z);
	returnvalue += CalcGravTerm(meas_x - ul_corner_x, meas_y - (ul_corner_y
			+ y_size), meas_z - ul_corner_z);
	returnvalue += CalcGravTerm(meas_x - ul_corner_x, meas_y - ul_corner_y,
			meas_z - (ul_corner_z + z_size));
	// two corners shifted -> negative sign
	returnvalue -= CalcGravTerm(meas_x - (ul_corner_x + x_size), meas_y
			- (ul_corner_y + y_size), meas_z - ul_corner_z);
	returnvalue -= CalcGravTerm(meas_x - (ul_corner_x + x_size), meas_y
			- ul_corner_y, meas_z - (ul_corner_z + z_size));
	returnvalue -= CalcGravTerm(meas_x - ul_corner_x, meas_y - (ul_corner_y
			+ y_size), meas_z - (ul_corner_z + z_size));
	//three corners shifted -> positive sign
	returnvalue += CalcGravTerm(meas_x - (ul_corner_x + x_size), meas_y
			- (ul_corner_y + y_size), meas_z - (ul_corner_z + z_size));
	// we multiply the geometric term by the gravitational constant and the density
	return -returnvalue * Grav_const;
}

__global__ void CalcScalarMeas(const float x_meas, const float y_meas,
		const float z_meas, const float *XCoord, const float *YCoord,
		const float *ZCoord, const float *XSizes, const float *YSizes,
		const float *ZSizes, const float *Densities, const int nx,
		const int ny, const int nz, float *returnvalue) {
	const unsigned int offset = blockIdx.x * blockDim.x + threadIdx.x;
	int xindex, yindex, zindex;

	if (offset < nx * ny * nz) {
		OffsetToIndex(offset, ny, nz, xindex, yindex, zindex);
		//we store the current value for possible sensitivity calculations
		//currvalue contains the geometric term, i.e. the sensitivity

		float currvalue = CalcGravBoxTerm(x_meas, y_meas, z_meas,
				XCoord[xindex], YCoord[yindex], ZCoord[zindex], XSizes[xindex],
				YSizes[yindex], ZSizes[zindex]);
		returnvalue[offset] = currvalue * Densities[offset];
	}

}

#endif // #ifndef _TEMPLATE_KERNEL_H_
