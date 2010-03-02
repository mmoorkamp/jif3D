//============================================================================
// Name        : SetupInversion.h
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPINVERSION_H_
#define SETUPINVERSION_H_

#include "../Inversion/ObjectiveFunction.h"
#include "../Inversion/GradientBasedOptimization.h"
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    namespace po = boost::program_options;
    class SetupInversion
      {

    public:
      po::options_description SetupOptions();
      void ConfigureInversion(const po::variables_map &vm, boost::shared_ptr<
          jiba::GradientBasedOptimization> &Optimizer, boost::shared_ptr<
          jiba::ObjectiveFunction> ObjFunction);
      SetupInversion();
      virtual ~SetupInversion();
      };

  }

#endif /* SETUPINVERSION_H_ */
