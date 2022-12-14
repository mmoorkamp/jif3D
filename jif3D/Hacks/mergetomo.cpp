//============================================================================
// Name        : mergetomo.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2009, MM
//============================================================================


#include <iostream>
#include <string>
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Tomo/ReadWriteTomographyData.h"
#include "../ModelBase/VTKTools.h"
#include "../Global/FileUtil.h"
#include <boost/cast.hpp>

using namespace std;

int main()
  {

    jif3D::ThreeDSeismicModel BaseModel;
    std::string ModelFilename = jif3D::AskFilename("Base Filename: ");
    BaseModel.ReadNetCDF(ModelFilename);

    jif3D::ThreeDSeismicModel AnomalyModel;
    ModelFilename = jif3D::AskFilename("Anomaly Filename: ");
    AnomalyModel.ReadNetCDF(ModelFilename);

    const size_t ngrid = BaseModel.GetSlownesses().num_elements();

    const double threshold = 0.00027;
    for (size_t i = 0; i < ngrid; ++i)
      {
        const double value = *(AnomalyModel.GetSlownesses().origin() + i);
        if (value > threshold)
          *(BaseModel.SetSlownesses().origin() + i) = value;
      }

    BaseModel.WriteNetCDF("ano.nc");
    BaseModel.WriteVTK("ano.vtk");
  }
