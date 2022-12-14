//============================================================================
// Name        : refinemt.cpp
// Author      : Apr 27, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#include "../ModelBase/ModelRefiner.h"
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "X3DModel.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Magnetics/ThreeDMagneticModel.h"


#include <boost/variant.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>

#include <netcdf>
#include <typeinfo>
typedef boost::mpl::vector<jif3D::ThreeDMagneticModel, jif3D::X3DModel,
    jif3D::ThreeDSeismicModel, jif3D::ThreeDGravityModel> vecModelType;
//boost::mpl::at_c<vecType, 3>::type hi = 3;
typedef boost::make_variant_over<vecModelType>::type data_type;
data_type Model;

struct model_test
  {
  std::string filename;

  template<typename U>
  void operator()(U x)
    {
      U test;

      try
        {
          test.ReadNetCDF(filename);
          Model = test;
        } catch (jif3D::FatalException &e)
        {
          // ignore
        } catch(const netCDF::exceptions::NcException &ex)
        {
          // ignore
        }

    }
  model_test(std::string f) :
      filename(f)
    {
    }
  };

class write_model: public boost::static_visitor<>
  {
public:
  template<typename T>
  void operator()(T & Model) const
    {
      try {
        const std::string finemeshname = jif3D::AskFilename("Refinement mesh: ");

        //read the file with the mesh for refinement
        auto FineMesh(Model);
        FineMesh.ReadNetCDF(finemeshname);
        //create a ModelRefiner object and pass the information about
        //the coordinates of the fine model
        jif3D::ModelRefiner Refiner;
        Refiner.SetXCoordinates(FineMesh.GetXCoordinates());
        Refiner.SetYCoordinates(FineMesh.GetYCoordinates());
        Refiner.SetZCoordinates(FineMesh.GetZCoordinates());
        //project the values from the coarse model onto the fine model
        Refiner.RefineModel(Model, FineMesh);
        //for the background we simply copy the configuration from the original model
        //as the background layers do not have to have the same thickness as the mesh
        /*FineMesh.SetBackgroundConductivities(
         CoarseMod.GetBackgroundConductivities());
         FineMesh.SetBackgroundThicknesses(CoarseMod.GetBackgroundThicknesses());*/
        //ask for the name of the output file and write a new netcdf file
        const std::string outname = jif3D::AskFilename("Outputfile: ", false);
        FineMesh.WriteNetCDF(outname);
        FineMesh.WriteVTK(outname + ".vtk");
      } catch(netCDF::exceptions::NcException &ex) {
        // ignore
      }
    }

  };

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
    using boost::mpl::for_each;
    using boost::mpl::range_c;

//read the file with the original conductivity model
    const std::string coarsemodname = jif3D::AskFilename("Original model: ");
    model_test CoarseMod(coarsemodname);
    for_each<vecModelType>(CoarseMod);

    boost::apply_visitor(write_model(), Model);
  }
