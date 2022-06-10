//============================================================================
// Name        : SetupTomo.h
// Author      : Mar 2, 2010
// Version     :
// Copyright   : 2010, mmoorkamp
//============================================================================

#ifndef SETUPTOMO_H_
#define SETUPTOMO_H_

#include "GeneralDataSetup.h"
#include "../Global/Jif3DGlobal.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../Tomo/TomographyCalculator.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include <boost/program_options.hpp>

namespace jif3D
  {
    namespace po = boost::program_options;

    /** \addtogroup joint Joint inversion routines */
    /* @{ */

    //! Setup the tomography part of the joint inversion
    /*!This class reads the information about grid sizes, measurement configuration, data etc.
     * from the command line and through interactive input. It configures the objective
     * function for seismic refraction data and adds it to the joint objective.
     * Also for the joint inversion the tomography model is always considered the starting model in
     * terms of geometry.
     */
    class J3DEXPORT SetupTomo: public GeneralDataSetup
      {
    private:
      //! A shared pointer to the objective function object for seismic tomography data
      boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::TomographyCalculator> > TomoObjective;
      //! The picking  error in ms to assume for construction of the data variance
      double pickerr;
      //! Cell size in m for  the refinement model, can optionally be set on the command line
      double CellSize;
      //! The name of the starting model
      std::string modelfilename;
      //! The name of the file with the travel time picks
      std::string datafilename;
      //! The weight for the tomography data in the joint inversion
      double tomolambda;
      //! The tomography starting model
      jif3D::ThreeDSeismicModel TomoModel;
      double minslow;
      double maxslow;
      jif3D::TomographyData TomoData;
      bool WriteRays;
    public:
      //! Read only access to the starting model for seismic tomography
      /*! Read only access to the starting model for seismic tomography
       * @return The object containing the starting model.
       */
      const jif3D::ThreeDSeismicModel& GetModel() const
        {
          return TomoModel;
        }
      //! read only access to the objective function object for seismic tomography data
      const jif3D::ThreeDModelObjective<jif3D::TomographyCalculator>& GetTomoObjective()
        {
          return *TomoObjective;
        }
      //! Setup the program options for the tomography part of the inversion
      po::options_description SetupOptions();
      //! Setup the tomography objective function
      /*! Setup the objective function for inverting seismic tomography data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created tomography objective is added to
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param xorigin The origin for the inversion grid in x-direction
       * @param yorigin The origin for the inversion grid in y-direction
       * @return True if the weight the tomography objective is greater zero, i.e. we added an objective function to JointObjective, false otherwise
       */
      virtual bool
      SetupObjective(const boost::program_options::variables_map &vm,
          jif3D::JointObjective &Objective, jif3D::ThreeDModelBase &InversionMesh,
          jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
          std::vector<std::string> &SegmentNames,
          std::vector<parametertype> &SegmentTypes, boost::filesystem::path TempDir =
              boost::filesystem::current_path()) override;
      virtual void IterationOutput(const std::string &filename,
          const jif3D::rvec &ModelVector) override;
      virtual void FinalOutput(const jif3D::rvec &FinalModelVector) override;
      SetupTomo();
      virtual ~SetupTomo();
      };
  /* @} */
  }

#endif /* SETUPTOMO_H_ */
