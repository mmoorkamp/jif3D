//============================================================================
// Name        : SetupCoupling.h
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPCOUPLING_H_
#define SETUPCOUPLING_H_

#include "../Inversion/ModelTransforms.h"
#include <boost/program_options.hpp>

namespace jiba
  {
    namespace po = boost::program_options;
    class SetupCoupling
      {
    public:
      po::options_description SetupOptions();
      void ConfigureCoupling(
          boost::shared_ptr<jiba::GeneralModelTransform> &InversionTransform,
          boost::shared_ptr<jiba::GeneralModelTransform> &TomoTransform,
          boost::shared_ptr<jiba::GeneralModelTransform> &MTTransform,
          boost::shared_ptr<jiba::GeneralModelTransform> &GravityTransform,
          boost::shared_ptr<jiba::GeneralModelTransform> &RegTransform);
      SetupCoupling();
      virtual ~SetupCoupling();
      };

  }

#endif /* SETUPCOUPLING_H_ */
