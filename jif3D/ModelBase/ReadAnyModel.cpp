/*
 * ReadAnyModel.cpp
 *
 *  Created on: 7 Apr 2015
 *      Author: mm489
 */

#include "ReadAnyModel.h"
#include "../ModelBase/ModelRefiner.h"
#include "../Global/FileUtil.h"
#include "../Global/FatalException.h"
#include "../MT/X3DModel.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Magnetics/ThreeDMagneticModel.h"
#include "../SurfaceWaves/SurfaceWaveModel.h"
#include "../DCResistivity/ThreeDDCResistivityModel.h"


#include <boost/variant.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/for_each.hpp>
#include <netcdf>
#include <typeinfo>

namespace jif3D
  {

    typedef boost::mpl::vector<jif3D::SurfaceWaveModel, jif3D::ThreeDMagneticModel, jif3D::X3DModel,
        jif3D::ThreeDSeismicModel, jif3D::ThreeDGravityModel, jif3D::ThreeDDCResistivityModel> vecModelType;
//boost::mpl::at_c<vecType, 3>::type hi = 3;
    typedef boost::make_variant_over<vecModelType>::type data_type;
    data_type Model;
    boost::shared_ptr<jif3D::ThreeDModelBase> ReturnModel;

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
            } catch (const jif3D::FatalException &e)
            {
              // ignore
            } catch(const netCDF::exceptions::NcException &ex) {
              // ignore
            }

        }

      model_test(std::string f) :
          filename(f)
        {
        }
      };

    class assign_model: public boost::static_visitor<>
      {
    public:
      template<typename T>
      void operator()(T & MyModel) const
        {
          ReturnModel = boost::shared_ptr<jif3D::ThreeDModelBase>(
              new jif3D::ThreeDModelBase(MyModel));
        }

      };

    boost::shared_ptr<jif3D::ThreeDModelBase> ReadAnyModel(const std::string &Filename)
      {
        using boost::mpl::for_each;
        using boost::mpl::range_c;

//read the file with the original conductivity model

        model_test CoarseMod(Filename);
        for_each<vecModelType>(CoarseMod);

        boost::apply_visitor(assign_model(), Model);
        return ReturnModel;
      }

  }

