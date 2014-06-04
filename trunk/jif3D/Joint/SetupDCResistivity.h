//============================================================================
// Name        : SetupDCResistivity.h
// Author      : 4 Jun 2014
// Version     : 
// Copyright   : 2014, mm489
//============================================================================

#ifndef SETUPDCRESISTIVITY_H_
#define SETUPDCRESISTIVITY_H_

#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/JointObjective.h"
#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include "../DCResistivity/DCResistivityCalculator.h"
#include <boost/program_options.hpp>

namespace jif3D
  {

    namespace po = boost::program_options;

    class SetupDCResistivity
      {
    private:
      double relerr;
      double minerr;
      boost::shared_ptr<DCResistivityCalculator> Calculator;
      //! Stores the grid for the DC resistivity model and the starting model
      jif3D::ThreeDDCResistivityModel Model;
      //! Possible pointer to the scalar Magnetics objective function, gets assigned below depending on user input
      boost::shared_ptr<jif3D::ThreeDModelObjective<DCResistivityCalculator> > DCObjective;
    public:
      boost::shared_ptr<DCResistivityCalculator> GetCalculator()
        {
          return Calculator;
        }
      //! read-only access to the objective function for DC resistivity
      const jif3D::ThreeDModelObjective<DCResistivityCalculator> &GetObjective()
        {
          return *DCObjective;
        }

      //! Return the DC resistivity starting model and measurement positions  that were read in
      /*! The model object also contains the measurement positions.
       * @return The model object containing the starting model and measurement positions.
       */
      const jif3D::ThreeDDCResistivityModel &GetModel() const
        {
          return Model;
        }

      //! Return an options descriptions object for boost::program_options that contains information about Magnetics options
      po::options_description SetupOptions();
      //! Setup the objective function and add to the joint objective
      /*! Setup the objective function for inverting scalar and tensorial data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created Magnetics objective(s) are added to
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param xorigin The origin for the inversion grid in x-direction
       * @param yorigin The origin for the inversion grid in y-direction
       * @param TempDir A directory to store temporary files with sensitivity information
       * @return True if the weight for the DC resistivity objective is greater zero, i.e. we added an objective function to JointObjective, false otherwise
       */
      bool
      SetupObjective(const po::variables_map &vm, jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::GeneralModelTransform> &Transform,
          double xorigin = 0.0, double yorigin = 0.0);
      SetupDCResistivity();
      virtual ~SetupDCResistivity();
      };

  } /* namespace jif3D */

#endif /* SETUPDCRESISTIVITY_H_ */
