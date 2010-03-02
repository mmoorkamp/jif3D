//============================================================================
// Name        : SetupRegularization.h
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPREGULARIZATION_H_
#define SETUPREGULARIZATION_H_

#include "../Inversion/JointObjective.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include <boost/program_options.hpp>

namespace jiba
  {
    namespace po = boost::program_options;
    class SetupRegularization
      {
    public:
      po::options_description SetupOptions();
      void
          SetupObjective(const po::variables_map &vm,
              jiba::JointObjective &Objective,
              const ThreeDSeismicModel &StartModel, boost::shared_ptr<jiba::GeneralModelTransform> Transform);
      SetupRegularization();
      virtual ~SetupRegularization();
      };

  }

#endif /* SETUPREGULARIZATION_H_ */
