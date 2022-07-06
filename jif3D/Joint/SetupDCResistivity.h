//============================================================================
// Name        : SetupDCResistivity.h
// Author      : 12 Wen 2021
// Version     :
// Copyright   : 2021, zhanjie and mm489
//============================================================================

#ifndef SETUPDCRESISTIVITY_H_
#define SETUPDCRESISTIVITY_H_

#include "../Global/Jif3DGlobal.h"
#include "GeneralDataSetup.h"

#include "../Inversion/ThreeDModelObjective.h"
#include "../Inversion/JointObjective.h"
#include "../DCResistivity/ThreeDDCResistivityModel.h"
#include "../DCResistivity/DCResistivityData.h"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "../DCResistivity/DCResistivityCalculator.h"

namespace jif3D
  {
    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the dcresistivity part of the joint inversion
    /*!This class reads the information about grid sizes, measurement configuration, data etc.
     * from the command line and through interactive input. It configures the objective
     * function for dc data and adds it to the joint objective.
     * Also for the joint inversion the DC model is always considered the starting model in
     * terms of geometry.
     */
    class J3DEXPORT SetupDCResistivity  : public GeneralDataSetup
      {
    private:
      //! A shared pointer to the objective function object for DC data
      boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::DCResistivityCalculator> > DCObjective;
      //! The data error in ohm.m to assume for construction of the data variance
      double relerr, minerr;
      double maxres;
      double minres;
      //! Cell size in m for  the refinement model, can optionally be set on the command line
      double CellSize;
      //! The name of the starting model
      std::string dcmodelfilename;
      //! The name of the file with the dc resistivity data
      std::string dcdatafilename;
      //! The weight for the DC data in the joint inversion
      double dclambda;
    	jif3D::DCResistivityData DCData;

      //! The resistivity starting model
      jif3D::ThreeDDCResistivityModel DCModel;
    public:
      //! Read only access to the starting model for DCResistivity
      /*! Read only access to the starting model for DCResistivity
       * @return The object containing the starting model.
       */
      const jif3D::ThreeDDCResistivityModel &GetModel() const
        {
          return DCModel;
        }
      //! read-only access to the objective function for DC resistivity data
      const jif3D::ThreeDModelObjective<DCResistivityCalculator> &GetDCObjective()
        {
          return *DCObjective;
        }

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
      virtual void FinalOutput(const std::string &filename ,const jif3D::rvec &FinalModelVector) override;
      SetupDCResistivity();
      virtual ~SetupDCResistivity();
      };
  /* @} */
  }

#endif /* SETUPDCRESISTIVITY_H_ */
