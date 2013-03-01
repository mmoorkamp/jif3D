//============================================================================
// Name        : SetupRegularization.h
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPREGULARIZATION_H_
#define SETUPREGULARIZATION_H_

#include "../Inversion/ModelTransforms.h"
#include "../Regularization/MatOpRegularization.h"
#include "../ModelBase/ThreeDModelBase.h"
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the type of regularization used in the inversion
    /*! Choose the type of regularization and the weighting
     * of the different directions depending on command line options.
     */
    class SetupRegularization
      {
    private:
      //here we declare some private variables that are set when
      //parsing the options
      //! The weight for the absolute model vector minimization in the regularization
      double beta;
      //! Do we substract the starting model for calculating the roughness of the current model
      bool substart;
      //! The regularization weight in x-direction
      double xweight;
      //! The regularization weight in y-direction
      double yweight;
      //! The regularization weight in z-direction
      double zweight;
      //! The smallness parameter b for minimum support regularization
      double minsuppb;
    public:
      //! Do we want to substract the starting model from the current model to calculate the roughness
      bool GetSubStart()
        {
          return substart;
        }
      //! Setup the possible command line options
      /*! This function creates an options_description object for boost::program
       * options that describes all possible options for regularization. This can
       * then be used to parse the command line or a configuration file etc.
       * @return The options_description object for the regularization options
       */
      po::options_description SetupOptions();
      //! Setup the regularization objective function depending on the program options
      /*! This function creates a new Regularization objective depending on the options that
       * have been set.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param StartModel The starting model containing a description of the cell geometry
       * @param Transform A transformation object to transform generalized inversion parameters to physical parameters that we want to regularize
       * @param CovModVec The model covariance vector, can be empty if covariance is assumed 1 or has to have one value for each cell in StartModel
       * @return A shared pointer to the configured regularization objective function object
       */
      boost::shared_ptr<jiba::ObjectiveFunction> SetupObjective(
          const po::variables_map &vm, const ThreeDModelBase &StartModel,
          const jiba::rvec &CovModVec);
      SetupRegularization();
      virtual ~SetupRegularization();
      };
  /* @} */
  }

#endif /* SETUPREGULARIZATION_H_ */
