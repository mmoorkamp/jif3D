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
#include "../Tomo/ThreeDSeismicModel.h"
#include "../MT/X3DModel.h"
#include "../Gravity/ThreeDGravityModel.h"
#include <boost/program_options.hpp>

namespace jiba
  {
    namespace po = boost::program_options;
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
          const jiba::ThreeDMTModel &MTMod, jiba::JointObjective &Objective);
      SetupCoupling();
      virtual ~SetupCoupling();
      };

  }

#endif /* SETUPCOUPLING_H_ */
