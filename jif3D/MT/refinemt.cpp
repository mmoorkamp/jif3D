//============================================================================
// Name        : refinemt.cpp
// Author      : Apr 27, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../ModelBase/ModelRefiner.h"
#include "../Global/FileUtil.h"
#include "X3DModel.h"

/*! \file refinemt.cpp
 * Refine an existing MT modeling mesh. The program asks for a file with a conductivity model
 * (original model) and another file with a refined model (refinement mesh).
 * From this second file only the mesh
 * information is used. The output file has the conductivities of the original model projected
 * onto the refinement mesh. Note that the program does not perform any interpolation. Also, it
 * does not check whether the new model fulfills the grid requirements of the forward code.
 */

int main()
  {
    jif3D::X3DModel CoarseMod, FineMesh;
    //read the file with the original conductivity model
    std::string coarsemodname = jif3D::AskFilename("Original model: ");
    CoarseMod.ReadNetCDF(coarsemodname);
    //read the file with the mesh for refinement
    std::string finemeshname = jif3D::AskFilename("Refinement mesh: ");
    FineMesh.ReadNetCDF(finemeshname);
    //create a ModelRefiner object and pass the information about
    //the coordinates of the fine model
    jif3D::ModelRefiner Refiner;
    Refiner.SetXCoordinates(FineMesh.GetXCoordinates());
    Refiner.SetYCoordinates(FineMesh.GetYCoordinates());
    Refiner.SetZCoordinates(FineMesh.GetZCoordinates());
    //project the values from the coarse model onto the fine model
    Refiner.RefineModel(CoarseMod, FineMesh);
    //for the background we simply copy the configuration from the original model
    //as the background layers do not have to have the same thickness as the mesh
    FineMesh.SetBackgroundConductivities(
        CoarseMod.GetBackgroundConductivities());
    FineMesh.SetBackgroundThicknesses(CoarseMod.GetBackgroundThicknesses());
    //ask for the name of the output file and write a new netcdf file
    std::string outname = jif3D::AskFilename("Outputfile: ", false);
    FineMesh.WriteNetCDF(outname);
    FineMesh.WriteVTK(outname+".vtk");
  }
