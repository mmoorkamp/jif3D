//============================================================================
// Name        : SetupDCCoupling.h
// Author      : July 3, 2014
// Version     :
// Copyright   : 2014, zhanjieshi
//============================================================================

#ifndef SETUPDCCOUPLING_H_
#define SETUPDCCOUPLING_H_

#include "../Global/Jif3DGlobal.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/RegularizationFunction.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../MT/X3DModel.h"
#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <string>

namespace jif3D
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
    class J3DEXPORT SetupDCCoupling
      {
    private:
      //! The transformation between generalized model parameters and slowness including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> SlowTrans;
      //! The transformation between generalized model parameters and conductivity including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> CondTrans;
      //! The transformation between generalized model parameters and density including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> ResTrans;
      //! The transformation between generalized model parameters and slowness for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> SlowCrossTrans;
      //! The transformation between generalized model parameters and conductivity for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> CondCrossTrans;
      //! The transformation between generalized model parameters and density for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> ResCrossTrans;
      //! The transformation between generalized model parameters and slowness for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> SlowRegTrans;
      //! The transformation between generalized model parameters and conductivity for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> CondRegTrans;
      //! The transformation between generalized model parameters and density for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> ResRegTrans;
      //! The minimal conductivity in S/m
      double mincond;
      //! The maximum conductivity in S/m
      double maxcond;
      //! The minimum slowness in s/m
      double minslow;
      //! The maximum slowness in s/m
      double maxslow;
      //! The minimum density in g/cm3
      double minres;
      //! The maximum density in g/cm3
      double maxres;
      //! The a coefficient for the slowness-density transform
      double res_a;
      //! The b coefficient for the slowness-density transform
      double res_b;
      //! The b coefficient for the slowness-density transform
      double res_c;
      //! The a coefficient for the slowness-conductivity transform
      double cond_a;
      //! The b coefficient for the slowness-conductivity transform
      double cond_b;
      //! The c coefficient for the slowness-conductivity transform
      double cond_c;
      //! The name for model file containing information where the parameter relationship should be applied
      std::string RelModelName;
      //! The replacement value for density where the parameter relationship is not valid
      double ResReplace;
      //! The replacement value for conductivity where the parameter relationship is not valid
      double CondReplace;
      //! Internal function to setup coupling and regularization when using a fixed parameter relationship
      void SetupFixedCouplingModel(jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::ThreeDSeismicModel &SeisMod,
          const jif3D::ThreeDDCResistivityModel &DCMod, const jif3D::ThreeDMTModel &MTMod,
          jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart);
    public:
      //! Return an options descriptions object for boost::program_options that contains information about options for coupling the methods
      po::options_description SetupOptions();
      //! Set the transformations between generalized and physical parameters for the different methods
      /*! Depending on whether we use direct parameter coupling or a cross gradient approach we have
       * to use different transformations between generalized and physical parameters for each method.
       * This function looks at the program options and sets the appropriate transforms for seismics, MT, gravity and
       * regularization.
       * @param vm The variables map containing the options set by the user
       * @param GeometryModel The model containing the geometry of the grid, will read in within this routine
       * @param TomoTransform The parameter transform for the seismic objective function
       * @param GravityTransform The parameter transform for the gravity objective function
       * @param MTTransform The parameter transform for the MT objective function
       * @param Wavelet Parametrize the inversion by a wavelet transform of the model parameters
       */
      void SetupTransforms(const po::variables_map &vm, ThreeDSeismicModel &GeometryModel,
          boost::shared_ptr<jif3D::GeneralModelTransform> &TomoTransform,
          boost::shared_ptr<jif3D::GeneralModelTransform> &DCTransform,
          boost::shared_ptr<jif3D::GeneralModelTransform> &MTTransform, bool Wavelet =
              false);
      //! Set the model vector for the inversion, its length and content depends on the coupling method
      /*! For the direct coupling joint inversion the model parameter vector only contains a slowness value
       * for each cell, while for the cross-gradient it contains slowness, density and conductivity. Also
       * depending on the approach we have to set a different number of regularization and coupling objective functions
       * which is done here as well.
       * @param vm The variables map containing the options set by the user
       * @param ModelGeometry The object carrying information about the grid
       * @param InvModel Contains the inversion parameter vector at the end
       * @param SeisMod The seismic starting model
       * @param GravMod The gravity starting model
       * @param MTMod The MT starting model
       * @param Objective The joint objective object we want to add the coupling and regularization terms to
       * @param Regularization An object for regularization (gradient/curvature etc.)
       * @param substart Do we want to substract the starting model for roughness calculations
       */
      void SetupModelVector(const po::variables_map &vm, jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::ThreeDSeismicModel &SeisMod,
          const jif3D::ThreeDDCResistivityModel &DCMod, const jif3D::ThreeDMTModel &MTMod,
          jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart);

      SetupDCCoupling();
      virtual ~SetupDCCoupling();
      };

  }

#endif /* SETUPDCCOUPLING_H_ */
