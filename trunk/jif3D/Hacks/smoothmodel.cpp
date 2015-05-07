/*
 * smoothmodel.cpp
 *
 *  Created on: Jul 17, 2014
 *      Author: mmoorkamp
 */

#include <iostream>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include "../MT/X3DModel.h"
#include "../Global/FileUtil.h"
#include "../Global/VecMat.h"
#include "../ModelTransforms/WaveletModelTransform.h"

namespace po = boost::program_options;
using namespace boost::accumulators;

int main(int argc, char *argv[])
{

	int xpadding = 2;
	int ypadding = 2;
	int zpadding = 0;

	int startlayer = 0;
	int endlayer = std::numeric_limits<int>::max();
	po::options_description desc("General options");
	desc.add_options()("help", "produce help message")("xsize",
			po::value<int>(&xpadding)->default_value(2),
			"Size of smoothing area in x-direction")("ysize",
			po::value<int>(&ypadding)->default_value(2),
			"Size of smoothing area in y-direction")("zsize",
			po::value<int>(&zpadding)->default_value(0),
			"Size of smoothing area in z-direction")("startlayer",
			po::value<int>(&startlayer)->default_value(0),
			"First layer to apply the smoothing")(
			"endlayer",
			po::value<int>(&endlayer)->default_value(
					std::numeric_limits<int>::max()),
			"Last layer to apply the smoothing");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	if (vm.count("help"))
	{
		std::cout << desc << "\n";
		return 1;
	}

	jif3D::X3DModel Model;
	std::string InfileName = jif3D::AskFilename("Input Model: ");
	Model.ReadNetCDF(InfileName);

	jif3D::X3DModel NewModel(Model);

	const int nx = Model.GetXCellSizes().size();
	const int ny = Model.GetYCellSizes().size();
	const int nz = Model.GetZCellSizes().size();

	startlayer = std::max(0, startlayer);
	endlayer = std::min(endlayer, nz);
	for (int i = xpadding; i < nx - xpadding; ++i)
	{
		for (int j = ypadding; j < ny - ypadding; ++j)
		{
			for (int k = startlayer; k < endlayer - zpadding; ++k)
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
				//std::vector<double> values;
				accumulator_set<double, stats<tag::median> > acc;
				for (int l = -xpadding; l <= xpadding; ++l)
				{
					for (int m = -ypadding; m <= ypadding; ++m)
					{
						for (int n = -zpadding; n <= zpadding; ++n)
						{
							acc(Model.GetConductivities()[i + l][j + m][k + n]);
						}
					}
				}
				//std::sort(values.begin(), values.end());
				NewModel.SetConductivities()[i][j][k] = median(acc);
			}
		}
	}

	NewModel.WriteNetCDF(InfileName + "_smooth.nc");
	NewModel.WriteVTK(InfileName + "_smooth.vtk");
}

