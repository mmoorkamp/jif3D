//============================================================================
// Name        : SetupMT.h
// Author      : Mar 1, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPMT_H_
#define SETUPMT_H_

#include "../MT/X3DModel.h"
#include "../MT/X3DMTCalculator.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/JointObjective.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the MT part of the joint inversion
    /*!This class reads the information about grid sizes, measurement configuration, data etc.
     * from the command line and through interactive input. It configures the objective
     * function for MT data and adds it to the joint objective.
     */
    class SetupMT
      {
    private:
      //! The objective function object for magnetotelluric data using X3D as a forward engine
      boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> > MTObjective;
      //! The relative data error to assume for construction of the data variance
      double relerr;
      //! The file name for a model with cell refinements
      std::string FineModelName;
      //! The object containing the geometry of the MT inversion model, has to match the geometry of the starting model
      jif3D::X3DModel MTModel;
      //! The name of the netcdf file containing information about the inverse covariance matrix
      std::string MTInvCovarName;
    public:
      //! Get read-only access to the objective function object, for example to output misfit information
      const jif3D::ThreeDModelObjective<jif3D::X3DMTCalculator> &GetMTObjective()
        {
          return *MTObjective;
        }
      //! Setup the program options for the MT part of the inversion
      po::options_description SetupOptions();
      //! Setup the MT objective function
      /*! Setup the objective function for inverting magnetotelluric data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created MT objective is added to
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param xorigin The origin for the inversion grid in x-direction
       * @param yorigin The origin for the inversion grid in y-direction
       * @param TempDir Set the directory to which all temporary files are written, this directory must exist
       * @return True if the weight for  the MT objective is greater zero, i.e. we added an objective function to JointObjective, false otherwise
       */
      bool
      SetupObjective(const po::variables_map &vm, jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::GeneralModelTransform> Transform, double xorigin = 0.0,
          double yorigin = 0.0, boost::filesystem::path TempDir =
              boost::filesystem::current_path());
          //! Return the MT model that has been set for the inversion
          const
      jif3D::X3DModel &GetModel() const
        {
          return MTModel;
        }
      SetupMT();
      virtual ~SetupMT();
      };
  /* @} */
  }

#endif /* SETUPMT_H_ */
