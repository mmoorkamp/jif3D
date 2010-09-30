//============================================================================
// Name        : SetupGravity.h
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPGRAVITY_H_
#define SETUPGRAVITY_H_

#include "../Gravity/ThreeDGravityModel.h"
#include "../Inversion/JointObjective.h"
#include <boost/program_options.hpp>

namespace jiba
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the gravity objective function for joint inversion
    /*! This class reads the information about grid sizes, measurement configuration etc.
     * from the command line and through interactive input.
     */
    class SetupGravity
      {
    private:
      // The relative error for the scalar data to assume for construction of the data variance
      double scalrelerr;
      // The relative error for the ftg data to assume for construction of the data variance
      double ftgrelerr;
      // The minimum error for the scalar data to assume for construction of the data variance
      double scalminerr;
      // The minimum error for the ftg data to assume for construction of the data variance
      double ftgminerr;
      jiba::ThreeDGravityModel GravModel;
    public:
      //! Return an options descriptions object for boost::program_options that contains information about gravity options
      po::options_description SetupOptions();
      //! Setup the objective function and add to the joint objective
      /*! Setup the objective function for inverting scalar and tensorial data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created gravity objective(s) are added to
       * @param StartModel The starting model geometry
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param NeedStartModel Do we need to ask the user for a starting model (true) or do we generate the geometry from the parameter StartModel
       */
      void
      SetupObjective(const po::variables_map &vm,
          jiba::JointObjective &Objective, const ThreeDModelBase &StartModel,
          boost::shared_ptr<jiba::GeneralModelTransform> Transform);
      //! Return the gravity model that was read in
      const jiba::ThreeDGravityModel &GetModel()
        {
          return GravModel;
        }
      SetupGravity();
      virtual ~SetupGravity();
      };
  /* @} */
  }

#endif /* SETUPGRAVITY_H_ */
