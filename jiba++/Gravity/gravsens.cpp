//============================================================================
// Name        : gravsens.cpp
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

/*! \file Calculate the sensitivity matrix for a given model and data geometry and write out the model eigenvectors.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include "../Global/convert.h"
#include "ThreeDGravityModel.h"
#include "ReadWriteGravityData.h"
#include "../Inversion/MatrixTools.h"
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include "../ModelBase/VTKTools.h"
namespace atlas = boost::numeric::bindings::atlas;

int main(int argc, char *argv[])
  {
    jiba::ThreeDGravityModel Model(true);

    std::string modelfilename, datafilename;
    //get the name of the file containing the mesh information
    std::cout << "Mesh Filename: ";
    std::cin >> modelfilename;
    //we read in a complete modelfile, but we only use the mesh information
    Model.ReadNetCDF(modelfilename);

    //we also read in data, but we only use the information about the measurements
    jiba::ThreeDGravityModel::tScalarMeasVec Data;
    jiba::ThreeDGravityModel::tMeasPosVec PosX, PosY, PosZ;
    std::cout << "Data Filename: ";
    std::cin >> datafilename;
    jiba::ReadScalarGravityMeasurements(datafilename, Data, PosX, PosY, PosZ);
    const size_t nmeas = PosX.size();
    //set the measurement points in the model to those of the data
    Model.ClearMeasurementPoints();
    for (size_t i = 0; i < nmeas; ++i)
      {
        Model.AddMeasurementPoint(PosX.at(i), PosY.at(i), PosZ.at(i));
      }

    //calculate the response, we don't actually care about the densities
    //and the model response, but we need the sensitivity matrix
    Model.CalcGravity();
    jiba::rmat Sensitivities(Model.GetScalarSensitivities());
    //do the SVD to get the eigenvalues and eigenvectors
    jiba::rvec s;
    jiba::rmat u, vt;
    jiba::SVD(Sensitivities, s, u, vt);
    //write out the eigenvalues
    std::ofstream evalfile((modelfilename + ".eval").c_str());
    std::copy(s.begin(), s.end(), std::ostream_iterator<double>(evalfile, "\n"));

    //print some information and ask which eigenvectors to write out
    std::cout << "There are " << s.size() << " eigenvalues.\n";
    size_t startindex = 0;
    std::cout << "Enter the first index for the eigenvalues to plot: ";
    std::cin >> startindex;
    size_t endindex = 0;
    std::cout << "Enter the last index for the eigenvalues to plot: ";
    std::cin >> endindex;
    //make sure we stay in the right range
    startindex = std::max(size_t(0),startindex);
    endindex = std::min(s.size(),endindex);
    const size_t xsize = Model.GetDensities().shape()[0];
    const size_t ysize = Model.GetDensities().shape()[1];
    const size_t zsize = Model.GetDensities().shape()[2];
    //make a structure with the shape of the model to plot the values
    jiba::ThreeDModelBase::t3DModelData sens(
        boost::extents[xsize][ysize][zsize]);
    //output the eigenvectors for the chosen indices
    //each eigenvector gets its own file
    for (size_t currindex = startindex; currindex < endindex; ++currindex)
      {
        const double maxvalue = *std::max_element(row(vt, currindex).begin(), row(vt, currindex).end(),boost::bind(fabs,_1) < boost::bind(fabs,_2));
        std::transform(row(vt, currindex).begin(), row(vt, currindex).end(), sens.data(),boost::bind(std::multiplies<double>(),_1,1./maxvalue));
        jiba::Write3DModelToVTK(modelfilename + ".sens."+stringify(currindex)+".vtk", "scalar_sens",
            Model.GetXCellSizes(), Model.GetYCellSizes(),
            Model.GetZCellSizes(), sens);
      }
  }
