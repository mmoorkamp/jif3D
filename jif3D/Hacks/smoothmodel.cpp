/*
 * smoothmodel.cpp
 *
 *  Created on: Jul 17, 2014
 *      Author: mmoorkamp
 */

#include <iostream>
#include <algorithm>
#include "../MT/X3DModel.h"
#include "../Global/FileUtil.h"
#include "../Global/VecMat.h"
#include "../ModelTransforms/WaveletModelTransform.h"

int main()
{
	jif3D::X3DModel Model;
	std::string InfileName = jif3D::AskFilename("Input Model: ");
	Model.ReadNetCDF(InfileName);

	jif3D::X3DModel NewModel(Model);

	const size_t nx = Model.GetXCellSizes().size();
	const size_t ny = Model.GetYCellSizes().size();
	const size_t nz = Model.GetZCellSizes().size();

	const int padding = 2;

	for (size_t i = padding; i < nx - padding; ++i)
	{
		for (size_t j = padding; j < ny - padding; ++j)
		{
			for (size_t k = 0; k < nz; ++k)
			{
/*				typedef boost::multi_array_types::index_range range;
				jif3D::ThreeDModelBase::t3DModelData::index_gen indices;
				jif3D::ThreeDModelBase::t3DModelData::const_array_view<2>::type myview =
						Model.GetConductivities()[indices[range(i - padding,
								i + padding)][range(j - padding, j + padding)][k]];
				const size_t nv = myview.num_elements();

				std::vector<double> values(nv);
				std::copy(myview.origin(), myview.origin() + nv,
						values.begin());	*/
				std::vector<double> values;
				for (int l = -padding; l < padding; ++l)
				{
					for (int m = -padding; m < padding; ++m)
					{
						values.push_back(Model.GetConductivities()[i+l][j+m][k]);
					}
				}
				std::sort(values.begin(), values.end());
				NewModel.SetConductivities()[i][j][k] = values[values.size() / 2];
			}
		}
	}

	NewModel.WriteNetCDF(InfileName + "_smooth.nc");
	NewModel.WriteVTK(InfileName + "_smooth.vtk");
}

