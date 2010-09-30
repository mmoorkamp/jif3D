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
      //! The relative data error to assume for construction of the data variance
      double relerr;
      //! The file name for a model with cell refinements
      std::string FineModelName;
      //! The object containing the geometry of the MT inversion model, has to match the geometry of the starting model
      jiba::X3DModel MTModel;
    public:
      //! Setup the program options for the MT part of the inversion
      po::options_description SetupOptions();
      //! Setup the MT objective function
      /*! Setup the objective function for inverting magnetotelluric data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created MT objective is added to
       * @param StartModel The starting model geometry
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param NeedStartModel Do we need to ask the user for a starting model (true) or do we generate the geometry from the parameter StartModel
       */
      void
      SetupObjective(const po::variables_map &vm,
          jiba::JointObjective &Objective, const ThreeDModelBase &StartModel,
          boost::shared_ptr<jiba::GeneralModelTransform> Transform, bool NeedStartModel = true);
      //! Return the MT model that has been set for the inversion
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
