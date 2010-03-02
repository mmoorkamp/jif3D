//============================================================================
// Name        : SetupTomo.h
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPTOMO_H_
#define SETUPTOMO_H_

#include "../Inversion/JointObjective.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include <boost/program_options.hpp>

namespace jiba
  {
    namespace po = boost::program_options;
    class SetupTomo
      {
    public:
      po::options_description SetupOptions();
      void
      SetupObjective(const po::variables_map &vm,
          jiba::JointObjective &Objective, ThreeDSeismicModel &StartModel,
          boost::shared_ptr<jiba::GeneralModelTransform> Transform);
      SetupTomo();
      virtual ~SetupTomo();
      };

  }

#endif /* SETUPTOMO_H_ */
