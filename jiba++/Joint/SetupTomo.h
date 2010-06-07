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


    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the tomography part of the joint inversion
    /*!This class reads the information about grid sizes, measurement configuration, data etc.
     * from the command line and through interactive input. It configures the objective
     * function for seismic refraction data and adds it to the joint objective.
     */
    class SetupTomo
      {
    private:
      // The picking  error in ms to assume for construction of the data variance
      double pickerr;
      // Storage for the name of the refinement model, can optionally be set on the command line
      std::string FineModelName;
    public:
      po::options_description SetupOptions();
      void
      SetupObjective(const po::variables_map &vm,
          jiba::JointObjective &Objective, ThreeDSeismicModel &StartModel,
          boost::shared_ptr<jiba::GeneralModelTransform> Transform);
      SetupTomo();
      virtual ~SetupTomo();
      };
    /* @} */
  }

#endif /* SETUPTOMO_H_ */
