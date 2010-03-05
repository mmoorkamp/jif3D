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
    class SetupGravity
      {
    private:
      jiba::ThreeDGravityModel GravModel;
    public:
      po::options_description SetupOptions();
      void
          SetupObjective(const po::variables_map &vm,
              jiba::JointObjective &Objective,
              const ThreeDSeismicModel &StartModel, boost::shared_ptr<jiba::GeneralModelTransform> Transform);
      const jiba::ThreeDGravityModel &GetModel() {return GravModel;}
      SetupGravity();
      virtual ~SetupGravity();
      };

  }

#endif /* SETUPGRAVITY_H_ */
