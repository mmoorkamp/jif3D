//============================================================================
// Name        : SetupCoupling.h
// Author      : Mar 2, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPCOUPLING_H_
#define SETUPCOUPLING_H_

#include "../Global/Jif3DGlobal.h"
#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/RegularizationFunction.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../SurfaceWaves/SurfaceWaveModel.h"
#include "../MT/X3DModel.h"
#include "../Gravity/ThreeDGravityModel.h"
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
    class J3DEXPORT SetupCoupling
      {
    private:
      //! The transformation between generalized model parameters and slowness including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> SlowTrans;
      //! The transformation between generalized model parameters and conductivity including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> CondTrans;
      //! The transformation between generalized model parameters and density including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> DensTrans;
      //! The transformation between generalized model parameters and slowness for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> SlowCrossTrans;
      //! The transformation between generalized model parameters and conductivity for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> CondCrossTrans;
      //! The transformation between generalized model parameters and density for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> DensCrossTrans;
      //! The transformation between generalized model parameters and slowness for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> SlowRegTrans;
      //! The transformation between generalized model parameters and conductivity for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> CondRegTrans;
      //! The transformation between generalized model parameters and density for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> DensRegTrans;
      //! The minimal conductivity in S/m
      double mincond;
      //! The maximum conductivity in S/m
      double maxcond;
      //! The minimum slowness in s/m
      double minslow;
      //! The maximum slowness in s/m
      double maxslow;
      //! The minimum density in kg/m3
      double mindens;
      //! The maximum density in kg/m3
      double maxdens;
      //! The minimum S-wave velocity in m/s
      double minvel;
      //! The maximum S-wave velocity in m/s
      double maxvel;
      //! The a coefficient for the slowness-density transform
      double density_a;
      //! The b coefficient for the slowness-density transform
      double density_b;
      //! The a coefficient for the slowness-conductivity transform
      double cond_a;
      //! The b coefficient for the slowness-conductivity transform
      double cond_b;
      //! The c coefficient for the slowness-conductivity transform
      double cond_c;
      double vpvs;
      //! The name for model file containing information where the parameter relationship should be applied
      std::string RelModelName;
      //! The replacement value for density where the parameter relationship is not valid
      double DensReplace;
      //! The replacement value for conductivity where the parameter relationship is not valid
      double CondReplace;
      //! Internal function to setup coupling and regularization when using the cross-gradient approach
      void SetupCrossGradModel(jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::SurfaceWaveModel &SeisMod,
          const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
          jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart,
          const jif3D::ThreeDModelBase &TearModelX,
          const jif3D::ThreeDModelBase &TearModelY,
          const jif3D::ThreeDModelBase &TearModelZ, const jif3D::rvec &CovVec);
      void SetupMIModel(jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::SurfaceWaveModel &SeisMod,
          const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
          jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart);
      //! Internal function to setup coupling and regularization when using a fixed parameter relationship
      void SetupFixedCouplingModel(jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::SurfaceWaveModel &SeisMod,
          const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
          jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart);
      //! Internal function to setup coupling and regularization when using a parameter relationship designed for salt (unstable at the moment)
      void SetupSaltModel(const po::variables_map &vm, jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::SurfaceWaveModel &SeisMod,
          const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
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
      void SetupTransforms(const po::variables_map &vm, ThreeDModelBase &GeometryModel,
          boost::shared_ptr<jif3D::GeneralModelTransform> &TomoTransform,
          boost::shared_ptr<jif3D::GeneralModelTransform> &GravityTransform,
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
       * @param TearModelX The model describing the tear in the regularization in x-direction
       * @param TearModelY The model describing the tear in the regularization in y-direction
       * @param TearModelZ The model describing the tear in the regularization in z-direction
       * @param CovVec The vector containing the diagonal elements of the model covariance
       */
      void SetupModelVector(const po::variables_map &vm, jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::SurfaceWaveModel &SeisMod,
          const jif3D::ThreeDGravityModel &GravMod, const jif3D::ThreeDMTModel &MTMod,
          jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart,
          const jif3D::ThreeDModelBase &TearModelX,
          const jif3D::ThreeDModelBase &TearModelY,
          const jif3D::ThreeDModelBase &TearModelZ, const jif3D::rvec CovVec);

      SetupCoupling();
      virtual ~SetupCoupling();
      };

  }

#endif /* SETUPCOUPLING_H_ */
