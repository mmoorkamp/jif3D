//============================================================================
// Name        : SetupDCResistivity.h
// Author      : 4 Jun 2014
// Version     :
// Copyright   : 2014, mm489
//============================================================================

#ifndef SETUPDCRESISTIVITY_H_
#define SETUPDCRESISTIVITY_H_

#include "../Global/Jif3DGlobal.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/JointObjective.h"
#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include "../DCResistivity/DCResistivityCalculator.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the dcresistivity part of the joint inversion
    /*!This class reads the information about grid sizes, measurement configuration, data etc.
     * from the command line and through interactive input. It configures the objective
     * function for dc data and adds it to the joint objective.
     * Also for the joint inversion the DC model is always considered the starting model in
     * terms of geometry.
     */
    class J3DEXPORT SetupDCResistivity
      {
    private:
      //! A shared pointer to the objective function object for DC data
      boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::DCResistivityCalculator> > DCObjective;
      //! The data error in ohm.m to assume for construction of the data variance
      double relerr, minerr;
      //! Cell size in m for  the refinement model, can optionally be set on the command line
      double CellSize;
      //! The name of the file with the dc resistivity data
      std::string dcdatafilename;
      //! The name of the starting model
      std::string dcmodelfilename;
      //! The weight for the DC data in the joint inversion
      double dclambda;
      //! Stores the grid for the DC resistivity model and the starting model
      jif3D::ThreeDDCResistivityModel DCModel;

    public:

      //! read-only access to the objective function for DC resistivity
      const jif3D::ThreeDModelObjective<DCResistivityCalculator> &GetDCObjective()
        {
          return *DCObjective;
        }

      //! Return the DC resistivity starting model and measurement positions  that were read in
      /*! The model object also contains the measurement positions.
       * @return The model object containing the starting model and measurement positions.
       */
      const jif3D::ThreeDDCResistivityModel &GetModel() const
        {
          return DCModel;
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
          double xorigin = 0.0, double yorigin = 0.0, boost::filesystem::path TempDir =
                  boost::filesystem::current_path());
      SetupDCResistivity();
      virtual ~SetupDCResistivity();
      };

  } /* namespace jif3D */

#endif /* SETUPDCRESISTIVITY_H_ */
