//============================================================================
// Name        : SetupMagnetics.h
// Author      : Mar 1, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPMAGNETICS_H_
#define SETUPMAGNETICS_H_

#include "../Global/Jif3DGlobal.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../Magnetics/ThreeDSusceptibilityModel.h"
#include "../Magnetics/TotalFieldMagneticData.h"
#include "../Inversion/JointObjective.h"
#include "GeneralDataSetup.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>



namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the Magnetics objective function for joint inversion
    /*! This class reads the information about grid sizes, measurement configuration etc.
     * from the command line and through interactive input.
     */
    class J3DEXPORT SetupMagnetics : public GeneralDataSetup
      {
    public:

#ifdef MAGDISK
      typedef typename jif3D::DiskGravMagCalculator<jif3D::TotalFieldMagneticData> MagCalculatorType;
#endif
#ifdef MAGMEM
      typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::TotalFieldMagneticData> MagCalculatorType;
#endif
#ifdef MAGCALC
      typedef typename jif3D::MinMemGravMagCalculator<jif3D::TotalFieldMagneticData> MagCalculatorType;
#endif
    private:
      double minsus;
      double maxsus;
      double maglambda;
      double inclination;
      double declination;
      double fieldstrength;
      boost::shared_ptr<MagCalculatorType> Calculator;
      //! The relative error for the scalar data to assume for construction of the data variance
      double relerr;
      //! The minimum error for the scalar data to assume for construction of the data variance
      double minerr;
      //! Stores the grid for the scalar Magnetics model and the starting model
      jif3D::ThreeDSusceptibilityModel Model;
      //! Possible pointer to the scalar Magnetics objective function, gets assigned below depending on user input
      boost::shared_ptr<jif3D::ThreeDModelObjective<MagCalculatorType> > MagObjective;
    public:
      boost::shared_ptr<MagCalculatorType> GetCalculator()
        {
          return Calculator;
        }
      double GetInclination() const
        {
          return inclination;
        }
      double GetDeclination() const
        {
          return declination;
        }
      double GetFielStrength() const
        {
          return fieldstrength;
        }
      //! read-only access to the objective function for scalar Magnetics data
      const jif3D::ThreeDModelObjective<MagCalculatorType>& GetObjective()
        {
          return *MagObjective;
        }

      //! Return the scalar Magnetics starting model and measurement positions  that were read in
      /*! The model object also contains the measurement positions.
       * @return The model object containing the starting model and measurement positions.
       */
      const jif3D::ThreeDSusceptibilityModel& GetModel() const
        {
          return Model;
        }

      //! Return an options descriptions object for boost::program_options that contains information about Magnetics options
      virtual po::options_description SetupOptions() override;
      //! Setup the objective function and add to the joint objective
      /*! Setup the objective function for inverting scalar and tensorial data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created Magnetics objective(s) are added to
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param xorigin The origin for the inversion grid in x-direction
       * @param yorigin The origin for the inversion grid in y-direction
       * @param TempDir A directory to store temporary files with sensitivity information
       * @return True if the weight for one of the Magnetics objectives is greater zero, i.e. we added an objective function to JointObjective, false otherwise
       */
      virtual bool
      SetupObjective(const po::variables_map &vm, jif3D::JointObjective &Objective,
          jif3D::ThreeDModelBase &InversionMesh, jif3D::rvec &CovModVec,
          std::vector<size_t> &startindices,
          std::vector<std::string> &SegmentNames,
          std::vector<parametertype> &SegmentTypes,
          boost::filesystem::path TempDir =
              boost::filesystem::current_path()) override;
      virtual void IterationOutput(const std::string &filename, const jif3D::rvec &ModelVector) override;
      virtual void FinalOutput(const jif3D::rvec &FinalModelVector) override;

      SetupMagnetics();
      virtual ~SetupMagnetics();
      };
  /* @} */
  }

#endif /* SETUPMAGNETICS_H_ */
