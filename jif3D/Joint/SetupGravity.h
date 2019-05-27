//============================================================================
// Name        : SetupGravity.h
// Author      : Mar 1, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPGRAVITY_H_
#define SETUPGRAVITY_H_

#include "../Global/Jif3DGlobal.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/JointObjective.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ScalarGravityData.h"
#include "../Gravity/TensorGravityData.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the gravity objective function for joint inversion
    /*! This class reads the information about grid sizes, measurement configuration etc.
     * from the command line and through interactive input.
     */
    class J3DEXPORT SetupGravity
      {
    public:
      typedef typename jif3D::DiskGravMagCalculator<jif3D::ScalarGravityData> ScalarCalculatorType;
      typedef typename jif3D::DiskGravMagCalculator<jif3D::TensorGravityData> TensorCalculatorType;
    private:
      //! The relative error for the scalar data to assume for construction of the data variance
      double scalrelerr;
      //! The relative error for the ftg data to assume for construction of the data variance
      double ftgrelerr;
      //! The minimum error for the scalar data to assume for construction of the data variance
      double scalminerr;
      //! The minimum error for the ftg data to assume for construction of the data variance
      double ftgminerr;
      //! Stores the grid for the scalar gravity model and the starting model
      jif3D::ThreeDGravityModel ScalGravModel;
      //! Stores the grid for the FTG gravity model and the starting model
      jif3D::ThreeDGravityModel FTGGravModel;
      //! Possible pointer to the scalar gravity objective function, gets assigned below depending on user input
      boost::shared_ptr<jif3D::ThreeDModelObjective<ScalarCalculatorType> > ScalGravObjective;
      //! Possible pointer to the tensor gravity objective function, gets assigned below depending on user input
      boost::shared_ptr<jif3D::ThreeDModelObjective<TensorCalculatorType> > FTGObjective;
      //! Does the user want scalar gravity calculations and have we set up everything?
      bool HaveScal;
      //! Does the user want tensor gravity calculations and have we set up everything?
      bool HaveFTG;
    public:
      //! Does the user want scalar gravity calculations and have we set up everything?
      bool GetHaveScal() const
        {
          return HaveScal;
        }
      //! Does the user want tensor gravity calculations and have we set up everything?
      bool GetHaveFTG() const
        {
          return HaveFTG;
        }
      //! read-only access to the objective function for scalar gravity data
      const jif3D::ThreeDModelObjective<ScalarCalculatorType> &GetScalGravObjective()
        {
          return *ScalGravObjective;
        }
      //! read-only access to the objective function for tensor gravity data
      const jif3D::ThreeDModelObjective<TensorCalculatorType> &GetFTGObjective()
        {
          return *FTGObjective;
        }
      //! Return the scalar gravity starting model and measurement positions  that were read in
      /*! The model object also contains the measurement positions. As
       * we might have different stations for scalar and FTG gravity
       * data, we return two model objects. The densities will be identical for both of them.
       * @return The model object containing the starting model and measurement positions.
       */
      const jif3D::ThreeDGravityModel &GetScalModel() const
        {
          return ScalGravModel;
        }
      //! Return the FTG gravity starting model and measurement positions  that were read in
      /*! The model object also contains the measurement positions. As
       * we might have different stations for scalar and FTG gravity
       * data, we return two model objects. The densities will be identical for both of them.
       * @return The model object containing the starting model and measurement positions.
       */
      const jif3D::ThreeDGravityModel &GetFTGModel() const
        {
          return FTGGravModel;
        }
      //! Return an options descriptions object for boost::program_options that contains information about gravity options
      po::options_description SetupOptions();
      //! Setup the objective function and add to the joint objective
      /*! Setup the objective function for inverting scalar and tensorial data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created gravity objective(s) are added to
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param xorigin The origin for the inversion grid in x-direction
       * @param yorigin The origin for the inversion grid in y-direction
       * @param TempDir A directory to store temporary files with sensitivity information
       * @return True if the weight for one of the gravity objectives is greater zero, i.e. we added an objective function to JointObjective, false otherwise
       */
      bool
      SetupObjective(const po::variables_map &vm, jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::GeneralModelTransform> &Transform, double xorigin = 0.0,
          double yorigin = 0.0, boost::filesystem::path TempDir =
              boost::filesystem::current_path());
      SetupGravity();
      virtual ~SetupGravity();
      };
  /* @} */
  }

#endif /* SETUPGRAVITY_H_ */
