//============================================================================
// Name        : SetupGravity.h
// Author      : Mar 1, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPGRAVITY_H_
#define SETUPGRAVITY_H_

#include "GeneralDataSetup.h"
#include "../Global/Jif3DGlobal.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/JointObjective.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../Gravity/ThreeDGravityModel.h"
#include "../Gravity/ScalarGravityData.h"
#include "../Gravity/TensorGravityData.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the gravity objective function for joint inversion
    /*! This class reads the information about grid sizes, measurement configuration etc.
     * from the command line and through interactive input.
     */
    class J3DEXPORT SetupGravity : public GeneralDataSetup
      {
    private:
#ifdef GRAVDISK
      typedef typename jif3D::DiskGravMagCalculator<jif3D::ScalarGravityData> ScalarCalculatorType;
      typedef typename jif3D::DiskGravMagCalculator<jif3D::TensorGravityData> TensorCalculatorType;
#endif
#ifdef GRAVMEM
      typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ScalarGravityData> ScalarCalculatorType;
      typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::TensorGravityData> TensorCalculatorType;
#endif
#ifdef GRAVCALC
      typedef typename jif3D::MinMemGravMagCalculator<jif3D::ScalarGravityData> ScalarCalculatorType;
      typedef typename jif3D::MinMemGravMagCalculator<jif3D::TensorGravityData> TensorCalculatorType;
#endif
      //! The minimum density limit in the inversion
      double mindens;
      //! The maximum density limit in the inversion
      double maxdens;
      //! The relative error for the scalar data to assume for construction of the data variance
      double scalrelerr;
      //! The relative error for the ftg data to assume for construction of the data variance
      double ftgrelerr;
      //! The minimum error for the scalar data to assume for construction of the data variance
      double scalminerr;
      //! The minimum error for the ftg data to assume for construction of the data variance
      double ftgminerr;
      //! Stores the grid for the scalar gravity model and the starting model
      jif3D::ThreeDGravityModel ScalGravModel;
      //! Possible pointer to the scalar gravity objective function, gets assigned below depending on user input
      boost::shared_ptr<jif3D::ThreeDModelObjective<ScalarCalculatorType> > ScalGravObjective;
      //! Possible pointer to the tensor gravity objective function, gets assigned below depending on user input
      boost::shared_ptr<jif3D::ThreeDModelObjective<TensorCalculatorType> > FTGObjective;
      //! Does the user want scalar gravity calculations and have we set up everything?
      bool HaveScal;
      //! Does the user want tensor gravity calculations and have we set up everything?
      bool HaveFTG;
    public:

      //! Return an options descriptions object for boost::program_options that contains information about gravity options
      virtual po::options_description SetupOptions() override;
      //! Setup the objective function and add to the joint objective
      /*! Setup the objective function for inverting scalar and tensorial data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created gravity objective(s) are added to
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param xorigin The origin for the inversion grid in x-direction
       * @param yorigin The origin for the inversion grid in y-direction
       * @param TempDir A directory to store temporary files with sensitivity information
       * @return True if the weight for one of the gravity objectives is greater zero, i.e. we added an objective function to JointObjective, false otherwise
       */
      virtual bool
      SetupObjective(const po::variables_map &vm, jif3D::JointObjective &Objective,
          jif3D::ThreeDModelBase &InversionMesh, jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
          std::vector<std::string> &SegmentNames,
          std::vector<parametertype> &SegmentTypes,
          boost::filesystem::path TempDir =
              boost::filesystem::current_path()) override;
      virtual void IterationOutput(const std::string &filename, const jif3D::rvec &ModelVector) override;
      virtual void FinalOutput(const jif3D::rvec &FinalModelVector) override;
      SetupGravity();
      virtual ~SetupGravity();
      };
  /* @} */
  }

#endif /* SETUPGRAVITY_H_ */
