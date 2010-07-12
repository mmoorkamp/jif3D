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
      boost::shared_ptr<jiba::GeneralModelTransform> SlowTrans;
      boost::shared_ptr<jiba::GeneralModelTransform> CondTrans;
      boost::shared_ptr<jiba::GeneralModelTransform> DensTrans;
    public:
      po::options_description SetupOptions();
      void SetupTransforms(const po::variables_map &vm, boost::shared_ptr<
          jiba::GeneralModelTransform> &TomoTransform, boost::shared_ptr<
          jiba::GeneralModelTransform> &GravityTransform, boost::shared_ptr<
          jiba::GeneralModelTransform> &MTTransform, boost::shared_ptr<
          jiba::GeneralModelTransform> &RegTransform);
      void SetupModelVector(const po::variables_map &vm, jiba::rvec &InvModel,
          const jiba::ThreeDSeismicModel &SeisMod,
          const jiba::ThreeDGravityModel GravMod,
          const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective,
          boost::shared_ptr<jiba::MatOpRegularization> Regularization, bool substart);
      SetupCoupling();
      virtual ~SetupCoupling();
      };

  }

#endif /* SETUPCOUPLING_H_ */
