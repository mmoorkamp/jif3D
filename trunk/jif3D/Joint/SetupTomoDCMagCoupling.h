//============================================================================
// Name        : SetupTomoDCMagCoupling.h
// Author      : August 3, 2014
// Version     : 
// Copyright   : 2014, zhanjie
//============================================================================

#ifndef SETUPTOMODCMAGCOUPLING_H_
#define SETUPTOMODCMAGCOUPLING_H_

#include "../Inversion/ModelTransforms.h"
#include "../Inversion/JointObjective.h"
#include "../Regularization/RegularizationFunction.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include "../Magnetics/ThreeDMagneticModel.h"
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>
#include <string>

namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the coupling for the joint inversion depending on command line parameters
    /*! This class sets up how the tomo, dcresistivity and magnetic methods are connected within the joint inversion.
     * Currently, only cross gradient coupling is used.
     */
    class SetupTomoDCMagCoupling
      {
    private:
      //! The transformation between generalized model parameters and slowness including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> SlowTrans;
      //! The transformation between generalized model parameters and resistivity including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> ResTrans;
      //! The transformation between generalized model parameters and susceptibility including selecting the right range from the model vector
      boost::shared_ptr<jif3D::GeneralModelTransform> SuscepTrans;
      //! The transformation between generalized model parameters and slowness for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> SlowCrossTrans;
      //! The transformation between generalized model parameters and resistivity for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> ResCrossTrans;
      //! The transformation between generalized model parameters and susceptibility for the cross-gradient functional
      boost::shared_ptr<jif3D::GeneralModelTransform> SuscepCrossTrans;
      //! The transformation between generalized model parameters and slowness for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> SlowRegTrans;
      //! The transformation between generalized model parameters and resistivity for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> ResRegTrans;
      //! The transformation between generalized model parameters and susceptibility for the regularization functional
      boost::shared_ptr<jif3D::ChainedTransform> SuscepRegTrans;
      //! The minimum slowness in s/m
      double minslow;
      //! The maximum slowness in s/m
      double maxslow;
      //! The minimal resistivity in ohm.m
      double minres;
      //! The maximum resistivity in ohm.m
      double maxres;
      //! The minimum susceptibility in SI
      double minsuscep;
      //! The maximum susceptibility in SI
      double maxsuscep;
      //! Internal function to setup coupling and regularization when using the cross-gradient approach
      void SetupCrossGradModel(jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::ThreeDSeismicModel &SeisMod,
          const jif3D::ThreeDDCResistivityModel &DCMod, const jif3D::ThreeDMagneticModel &MagMod,
          jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart);

    public:
      //! Return an options descriptions object for boost::program_options that contains information about options for coupling the methods
      po::options_description SetupOptions();
      //! Set the transformations between generalized and physical parameters for the different methods
      /*! This function looks at the program options and sets the appropriate transforms for seismics, DCResistivity, magnetics and
       * regularization with cross gradient coupling.
       * @param vm The variables map containing the options set by the user
       * @param GeometryModel The model containing the geometry of the grid, will read in within this routine
       * @param TomoTransform The parameter transform for the seismic objective function
       * @param DCTransform The parameter transform for the DCResistivity objective function
       * @param MagTransform The parameter transform for the Magnetics objective function
       * @param Wavelet Parametrize the inversion by a wavelet transform of the model parameters
       */
      void SetupTransforms(const po::variables_map &vm, ThreeDSeismicModel &GeometryModel,
          boost::shared_ptr<jif3D::GeneralModelTransform> &TomoTransform,
          boost::shared_ptr<jif3D::GeneralModelTransform> &DCTransform,
          boost::shared_ptr<jif3D::GeneralModelTransform> &MagTransform, bool Wavelet =
              false);
      //! Set the model vector for the inversion
      /*! for the cross-gradient it contains slowness, resistivity and susceptibility.
       * we also have to set a different number of regularization and coupling objective functions
       * which is done here as well.
       * @param vm The variables map containing the options set by the user
       * @param ModelGeometry The object carrying information about the grid
       * @param InvModel Contains the inversion parameter vector at the end
       * @param SeisMod The seismic starting model
       * @param DCMod The resistivity starting model
       * @param MagMod The magnetics starting model
       * @param Objective The joint objective object we want to add the coupling and regularization terms to
       * @param Regularization An object for regularization (gradient/curvature etc.)
       * @param substart Do we want to substract the starting model for roughness calculations
       */
      void SetupModelVector(const po::variables_map &vm, jif3D::rvec &InvModel,
          const jif3D::ThreeDModelBase &ModelGeometry,
          const jif3D::ThreeDSeismicModel &SeisMod,
          const jif3D::ThreeDDCResistivityModel &DCMod, const jif3D::ThreeDMagneticModel &MagMod,
          jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::RegularizationFunction> Regularization, bool substart);

      SetupTomoDCMagCoupling();
      virtual ~SetupTomoDCMagCoupling();
      };

  }

#endif /* SETUPTOMODCMAGCOUPLING_H_ */
