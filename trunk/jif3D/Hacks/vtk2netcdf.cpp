//============================================================================
// Name        : vtk2netcdf.cpp
// Author      : 15 Jun 2013
// Version     :
// Copyright   : 2013, mm489
//============================================================================

#include <iostream>
#include <algorithm>
#include "../MT/X3DModel.h"
#include "../Global/FileUtil.h"
#include "../ModelBase/VTKTools.h"

int main()
{
	jif3D::X3DModel Model;
	std::string Infilename = jif3D::AskFilename("Filename: ");
	jif3D::ThreeDModelBase::t3DModelDim XCellSizes, YCellSizes;
	jif3D::Read3DModelFromVTK(Infilename, XCellSizes, YCellSizes,
			Model.SetZCellSizes(), Model.SetConductivities());
	Model.SetHorizontalCellSize(XCellSizes[1], YCellSizes[1],
			XCellSizes.num_elements(), YCellSizes.num_elements());
	const size_t nz = Model.GetZCellSizes().num_elements();
	std::vector<double> bg_cond(nz), bg_thick(nz);
	//std::copy(Model.GetZCellSizes().begin(), Model.GetZCellSizes().end(),
	//		bg_thick.end());

	for (size_t i = 0; i < nz; ++i)
	{
		bg_thick.at(i) = Model.GetZCellSizes()[i];
		bg_cond.at(i) = Model.GetConductivities()[0][0][i];
	}

	Model.SetBackgroundConductivities(bg_cond);
	Model.SetBackgroundThicknesses(bg_thick);
	Model.WriteNetCDF(Infilename + ".nc");
}

