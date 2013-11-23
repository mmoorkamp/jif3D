//============================================================================
// Name        : SetupInversion.h
// Author      : Mar 2, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPINVERSION_H_
#define SETUPINVERSION_H_

#include "../Inversion/ObjectiveFunction.h"
#include "../Inversion/GradientBasedOptimization.h"
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the core optimization
    /*! Depending on command line options this class sets up different optimization algorithm
     * that can be used for joint inversion.
     */
    class SetupInversion
      {
    private:
    	//! Scale the initial gradient for the L-BFGS optimization algorithm in the hope to improve convergence
    	bool scalegrad;
	//! The number of correction pairs for the L-BFGS optimization algorithm
      int corrpairs;
    public:
      //! Setup the possible command line options
      /*! This function creates an options_description object for boost::program
       * options that describes all possible options for the inversion. This can
       * then be used to parse the command line or a configuration file etc.
       * @return The options_description object for the inversion options
       */
      po::options_description SetupOptions();
      //! Setup the optimizer object so that it is ready to use for the inversion
      /*! We setup the optimizer object based on the program options and pass all
       * information about objective functions, starting model and covariance to it
       * so that it can be directly used for inversion.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param ObjFunction The objective function to use for the optimization
       * @param InvModel The vector containing the current model parameters
       * @param CovModVec The vector of covariances
       * @return A shared pointer to a new optimizer object
       */
      boost::shared_ptr<jif3D::GradientBasedOptimization> ConfigureInversion(
          const po::variables_map &vm, boost::shared_ptr<
              jif3D::ObjectiveFunction> ObjFunction, const jif3D::rvec &InvModel,
          const jif3D::rvec &CovModVec);
      SetupInversion();
      virtual ~SetupInversion();
      };
  /* @} */
  }

#endif /* SETUPINVERSION_H_ */
