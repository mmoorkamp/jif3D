/*
 * SetupMagnetization.h
 *
 *  Created on: Jun 2, 2022
 *      Author: max
 */

#ifndef JOINT_SETUPMAGNETIZATION_H_
#define JOINT_SETUPMAGNETIZATION_H_

#include "GeneralDataSetup.h"
#include "../Magnetics/ThreeComponentMagneticData.h"
#include "../Magnetics/ThreeDMagnetizationModel.h"
#include "../GravMag/DiskGravMagCalculator.h"
#include "../GravMag/MinMemGravMagCalculator.h"
#include "../GravMag/FullSensitivityGravMagCalculator.h"
#include "../Inversion/ThreeDModelObjective.h"
#include <boost/shared_ptr.hpp>

namespace jif3D
  {

    class SetupMagnetization: public GeneralDataSetup
      {
      private:
#ifdef MAGDISK
      typedef typename jif3D::DiskGravMagCalculator<jif3D::ThreeComponentMagneticData> MagCalculatorType;
#endif
#ifdef MAGMEM
      typedef typename jif3D::FullSensitivityGravMagCalculator<jif3D::ThreeComponentMagneticData> MagCalculatorType;
#endif
#ifdef MAGCALC
      typedef typename jif3D::MinMemGravMagCalculator<jif3D::ThreeComponentMagneticData> MagCalculatorType;
#endif
      double minmag;
      double maxmag;
      double maglambda;

      boost::shared_ptr<MagCalculatorType> Calculator;
      //! The relative error for the scalar data to assume for construction of the data variance
      double relerr;
      //! The minimum error for the scalar data to assume for construction of the data variance
      double minerr;
      //! Stores the grid for the scalar Magnetics model and the starting model
      jif3D::ThreeDMagnetizationModel Model;
      //! Possible pointer to the scalar Magnetics objective function, gets assigned below depending on user input
      boost::shared_ptr<jif3D::ThreeDModelObjective<MagCalculatorType> > MagObjective;

    public:
      virtual po::options_description SetupOptions() override;
      virtual bool
      SetupObjective(const po::variables_map &vm,
          jif3D::JointObjective &Objective, jif3D::ThreeDModelBase &InversionMesh,
          jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
          std::vector<std::string> &SegmentNames,
          std::vector<parametertype> &SegmentTypes,
          boost::filesystem::path TempDir = boost::filesystem::current_path()) override;
      virtual void IterationOutput(const std::string &filename,
          const jif3D::rvec &ModelVector) override;
      virtual void FinalOutput(const std::string &filename ,const jif3D::rvec &FinalModelVector) override;
      SetupMagnetization();
      virtual ~SetupMagnetization();
      };

  } /* namespace jif3D */

#endif /* JOINT_SETUPMAGNETIZATION_H_ */
