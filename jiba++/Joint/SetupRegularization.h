//============================================================================
// Name        : SetupRegularization.h
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPREGULARIZATION_H_
#define SETUPREGULARIZATION_H_

#include "../Inversion/ModelTransforms.h"
#include "../Regularization/MatOpRegularization.h"
#include "../ModelBase/ThreeDModelBase.h"
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

namespace jiba
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the type of regularization used in the inversion
    /*! Choose the type of regularization and the weighting
     * of the different directions depending on command line options.
     */
    class SetupRegularization
      {
    private:
      double beta;
      bool substart;
    public:
      bool GetSubStart()
        {
          return substart;
        }
      po::options_description SetupOptions();
      boost::shared_ptr<jiba::MatOpRegularization> SetupObjective(
          const po::variables_map &vm, const ThreeDModelBase &StartModel,
          boost::shared_ptr<jiba::GeneralModelTransform> Transform,
          const jiba::rvec &CovModVec);
      SetupRegularization();
      virtual ~SetupRegularization();
      };
  /* @} */
  }

#endif /* SETUPREGULARIZATION_H_ */
