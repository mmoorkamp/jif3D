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

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the core optimization
    /*! Depending on command line options this class sets up different optimization algorithm
     * that can be used for joint inversion.
     */
    class SetupInversion
      {

    public:
      po::options_description SetupOptions();
      void ConfigureInversion(const po::variables_map &vm, boost::shared_ptr<
          jiba::GradientBasedOptimization> &Optimizer, boost::shared_ptr<
          jiba::ObjectiveFunction> ObjFunction, const jiba::rvec &InvModel, const jiba::rvec &CovModVec);
      SetupInversion();
      virtual ~SetupInversion();
      };
    /* @} */
  }

#endif /* SETUPINVERSION_H_ */
