//============================================================================
// Name        : reldiff.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2009, MM
//============================================================================

#include <iostream>
#include <string>
#include <boost/bind.hpp>
#include "ThreeDSeismicModel.h"
#include "../Global/FileUtil.h"

/*! \file reldiff.cpp
 * Calculate the relative difference between two models and output it as a modelfile
 * and as a single column.
 */

int main()
  {
    std::string firstname = jif3D::AskFilename("First model file: ", true);
    std::string secondname = jif3D::AskFilename("Second model file: ", true);
    jif3D::ThreeDSeismicModel FirstModel;
    jif3D::ThreeDSeismicModel SecondModel;

    FirstModel.ReadNetCDF(firstname);
    SecondModel.ReadNetCDF(secondname);
    size_t ngrid = FirstModel.GetSlownesses().num_elements();
    if (SecondModel.GetSlownesses().num_elements() != ngrid)
      {
        std::cerr << "Number of model parameter does not match." << std::endl;
        return 100;
      }
    std::transform(FirstModel.GetSlownesses().origin(),
        FirstModel.GetSlownesses().origin() + ngrid,
        SecondModel.GetSlownesses().origin(),
        SecondModel.SetSlownesses().origin(), boost::bind(
            std::divides<double>(), boost::bind(std::minus<double>(), _1, _2),
            _1));
    SecondModel.WriteNetCDF("diff.nc");
    SecondModel.WriteVTK("diff.vtk");
    std::ofstream outfile("diff.hist");
    std::copy(SecondModel.GetSlownesses().origin(),
        SecondModel.GetSlownesses().origin() + ngrid, std::ostream_iterator<
            double>(outfile, "\n"));
  }
