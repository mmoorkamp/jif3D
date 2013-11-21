//============================================================================
// Name        : SetupMagnetics.h
// Author      : Mar 1, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPMAGNETICS_H_
#define SETUPMAGNETICS_H_

#include "../Inversion/ThreeDModelObjective.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "../Magnetics/ThreeDMagneticModel.h"
#include "../Inversion/JointObjective.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the Magnetics objective function for joint inversion
    /*! This class reads the information about grid sizes, measurement configuration etc.
     * from the command line and through interactive input.
     */
    class SetupMagnetics
      {
    public:
      typedef typename jif3D::DiskGravMagCalculator<jif3D::ThreeDMagneticModel> CalculatorType;
    private:
      double inclination;
      double declination;
      double fieldstrength;
      //! The relative error for the scalar data to assume for construction of the data variance
      double relerr;
      //! The minimum error for the scalar data to assume for construction of the data variance
      double minerr;
      //! Stores the grid for the scalar Magnetics model and the starting model
      jif3D::ThreeDMagneticModel Model;
      //! Possible pointer to the scalar Magnetics objective function, gets assigned below depending on user input
      boost::shared_ptr<jif3D::ThreeDModelObjective<CalculatorType> > MagObjective;

    public:

      //! read-only access to the objective function for scalar Magnetics data
      const jif3D::ThreeDModelObjective<CalculatorType> &GetObjective()
        {
          return *MagObjective;
        }

      //! Return the scalar Magnetics starting model and measurement positions  that were read in
      /*! The model object also contains the measurement positions. As
       * we might have different stations for scalar and FTG Magnetics
       * data, we return two model objects. The densities will be identical for both of them.
       * @return The model object containing the starting model and measurement positions.
       */
      const jif3D::ThreeDMagneticModel &GetModel() const
        {
          return Model;
        }

      //! Return an options descriptions object for boost::program_options that contains information about Magnetics options
      po::options_description SetupOptions();
      //! Setup the objective function and add to the joint objective
      /*! Setup the objective function for inverting scalar and tensorial data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created Magnetics objective(s) are added to
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param xorigin The origin for the inversion grid in x-direction
       * @param yorigin The origin for the inversion grid in y-direction
       * @param TempDir A directory to store temporary files with sensitivity information
       * @return True if the weight for one of the Magnetics objectives is greater zero, i.e. we added an objective function to JointObjective, false otherwise
       */
      bool
      SetupObjective(const po::variables_map &vm, jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::GeneralModelTransform> &Transform,
          double xorigin = 0.0, double yorigin = 0.0, boost::filesystem::path TempDir =
              boost::filesystem::current_path());
      SetupMagnetics();
      virtual ~SetupMagnetics();
      };
  /* @} */
  }

#endif /* SETUPMAGNETICS_H_ */
