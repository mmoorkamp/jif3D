//============================================================================
// Name        : SetupRegularization.h
// Author      : Mar 2, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPREGULARIZATION_H_
#define SETUPREGULARIZATION_H_

#include "../Global/Jif3DGlobal.h"
#include "GeneralDataSetup.h"
#include "../Inversion/ModelTransforms.h"
#include "../Regularization/RegularizationFunction.h"
#include "../ModelBase/ThreeDModelBase.h"
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

namespace jif3D
  {

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the type of regularization used in the inversion
    /*! Choose the type of regularization and the weighting
     * of the different directions depending on command line options.
     */
    class J3DEXPORT SetupRegularization
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
      //! Number of bins for histogram in entropy regularization
      size_t nbins;
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
       * @return A shared pointer to the configured regularization objective function object
       */
      boost::shared_ptr<jif3D::RegularizationFunction> SetupObjective(
          const po::variables_map &vm, const ThreeDModelBase &StartModel);
      //! Setup the regularization objective function depending on the program options
      /*! This function creates a new Regularization objective depending on the options that
       * have been set. It also returns the tear models that are useful for cross-gradient calculations
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param StartModel The starting model containing a description of the cell geometry
       * @param CovModVec The model covariance vector, can be empty if covariance is assumed 1 or has to have one value for each cell in StartModel
       * @param TearModelX The model describing the tear in the regularization in x-direction
       * @param TearModelY The model describing the tear in the regularization in y-direction
       * @param TearModelZ The model describing the tear in the regularization in z-direction
       * @return A shared pointer to the configured regularization objective function object
       */
      boost::shared_ptr<jif3D::RegularizationFunction> SetupObjective(
          const po::variables_map &vm, const ThreeDModelBase &StartModel,
          const jif3D::rvec &CovModVec, jif3D::ThreeDModelBase &TearModelX,
          jif3D::ThreeDModelBase &TearModelY, jif3D::ThreeDModelBase &TearModelZ);
      SetupRegularization();
      virtual ~SetupRegularization();
      };
  /* @} */
  }

#endif /* SETUPREGULARIZATION_H_ */
