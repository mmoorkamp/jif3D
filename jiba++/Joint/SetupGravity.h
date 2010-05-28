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
#include "../Tomo/ThreeDSeismicModel.h"
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
      jiba::ThreeDGravityModel GravModel;
    public:
      //! Return an options descriptions object for boost::program_options that contains information about gravity options
      po::options_description SetupOptions();
      //! Setup the objective function and add to the joint objective
      void
      SetupObjective(const po::variables_map &vm,
          jiba::JointObjective &Objective,
          const ThreeDSeismicModel &StartModel, boost::shared_ptr<
              jiba::GeneralModelTransform> Transform);
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
