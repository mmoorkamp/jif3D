//============================================================================
// Name        : SetupCoupling.h
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPCOUPLING_H_
#define SETUPCOUPLING_H_

#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/MatOpRegularization.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../MT/X3DModel.h"
#include "../Gravity/ThreeDGravityModel.h"
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the coupling for the joint inversion depending on command line parameters
    /*! This class sets up how the different methods are connected within the joint inversion.
     * Currently there are two possibilities: Direct parameter coupling (the default), or
     * cross gradient coupling. As the coupling has an impact on the parameter transformations,
     * the number of objective functions and the length of the model vector, the member functions
     * need information about most parts of the joint inversion.
     */
    class SetupCoupling
      {
    private:
      //! The transformation between generalized model parameters and slowness including selecting the right range from the model vector
      boost::shared_ptr<jiba::GeneralModelTransform> SlowTrans;
      //! The transformation between generalized model parameters and conductivity including selecting the right range from the model vector
      boost::shared_ptr<jiba::GeneralModelTransform> CondTrans;
      //! The transformation between generalized model parameters and density including selecting the right range from the model vector
      boost::shared_ptr<jiba::GeneralModelTransform> DensTrans;
      //! The transformation between generalized model parameters and slowness without selecting the right range from the model vector
      boost::shared_ptr<jiba::GeneralModelTransform> SlowCrossTrans;
      //! The transformation between generalized model parameters and conductivity without selecting the right range from the model vector
      boost::shared_ptr<jiba::GeneralModelTransform> CondCrossTrans;
      //! The transformation between generalized model parameters and density without selecting the right range from the model vector
      boost::shared_ptr<jiba::GeneralModelTransform> DensCrossTrans;
      //! The minimal conductivity in S/m
      double mincond;
      //! The maximum conductivity in S/m
      double maxcond;
      //! The minimum slowness in s/m
      double minslow;
      //! The maximum slowness in s/m
      double maxslow;
      //! The minimum density in g/cm3
      double mindens;
      //! The maximum density in g/cm3
      double maxdens;
    public:
      //! Return an options descriptions object for boost::program_options that contains information about options for coupling the methods
      po::options_description SetupOptions();
      //! Set the transformations between generalized and physical parameters for the different methods
      /*! Depending on whether we use direct parameter coupling or a cross gradient approach we have
       * to use different transformations between generalized and physical parameters for each method.
       * This function looks at the program options and sets the appropriate transforms for seismics, MT, gravity and
       * regularization.
       * @param vm The variables map containing the options set by the user
       * @param TomoTransform The parameter transform for the seismic objective function
       * @param GravityTransform The parameter transform for the gravity objective function
       * @param MTTransform The parameter transform for the MT objective function
       * @param RegTransform The parameter transform for the regularization
       * @param Wavelet Parametrize the inversion by a wavelet transform of the model parameters
       */
      void SetupTransforms(const po::variables_map &vm, boost::shared_ptr<
          jiba::GeneralModelTransform> &TomoTransform, boost::shared_ptr<
          jiba::GeneralModelTransform> &GravityTransform, boost::shared_ptr<
          jiba::GeneralModelTransform> &MTTransform, boost::shared_ptr<
          jiba::GeneralModelTransform> &RegTransform, bool Wavelet = false);
      //! Set the model vector for the inversion, its length and content depends on the coupling method
      /*! For the direct coupling joint inversion the model parameter vector only contains a slowness value
       * for each cell, while for the cross-gradient it contains slowness, density and conductivity. Also
       * depending on the approach we have to set a different number of regularization and coupling objective functions
       * which is done here as well.
       * @param vm The variables map containing the options set by the user
       * @param InvModel Contains the inversion parameter vector at the end
       * @param SeisMod The seismic starting model
       * @param GravMod The gravity starting model
       * @param MTMod The MT starting model
       * @param Objective The joint objective object we want to add the coupling and regularization terms to
       * @param Regularization An object for regularization (gradient/curvature etc.)
       * @param substart Do we want to substract the starting model for roughness calculations
       */
      void SetupModelVector(const po::variables_map &vm, jiba::rvec &InvModel,
          const jiba::ThreeDSeismicModel &SeisMod,
          const jiba::ThreeDGravityModel GravMod,
          const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
          boost::shared_ptr<jiba::MatOpRegularization> Regularization,
          bool substart);
      SetupCoupling();
      virtual ~SetupCoupling();
      };

  }

#endif /* SETUPCOUPLING_H_ */
