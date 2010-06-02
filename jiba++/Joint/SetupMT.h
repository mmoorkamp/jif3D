//============================================================================
// Name        : SetupMT.h
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPMT_H_
#define SETUPMT_H_

#include "../MT/X3DModel.h"
#include "../Inversion/JointObjective.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include <boost/program_options.hpp>

namespace jiba
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the MT part of the joint inversion
    /*!This class reads the information about grid sizes, measurement configuration, data etc.
     * from the command line and through interactive input. It configures the objective
     * function for MT data and adds it to the joint objective.
     */
    class SetupMT
      {
    private:
      // The relative data error to assume for construction of the data variance
      double relerr;
      std::string FineModelName;
      jiba::X3DModel MTModel;
    public:
      po::options_description SetupOptions();
      void
      SetupObjective(const po::variables_map &vm,
          jiba::JointObjective &Objective,
          const ThreeDSeismicModel &StartModel, boost::shared_ptr<
              jiba::GeneralModelTransform> Transform);
      const jiba::X3DModel &GetModel()
        {
          return MTModel;
        }
      SetupMT();
      virtual ~SetupMT();
      };
  /* @} */
  }

#endif /* SETUPMT_H_ */
