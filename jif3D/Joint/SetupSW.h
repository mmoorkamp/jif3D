/*
 * SetupSW.h
 *
 *  Created on: 29 Oct 2019
 *      Author: bweise
 */

#ifndef JOINT_SETUPSW_H_
#define JOINT_SETUPSW_H_

#include <boost/program_options.hpp>
#include "../Global/Jif3DGlobal.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../SurfaceWaves/SurfaceWaveCalculator.h"
#include "../SurfaceWaves/SurfaceWaveModel.h"

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
    class J3DEXPORT SetupSW
      {
    private:
      //! A shared pointer to the objective function object for surface wave tomography data
      boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::SurfaceWaveCalculator> > SurfaceWaveObjective;
      //! Relative error of phase delays for construction of the data variance
      double relerr;
      //! Minimal error of phase delays for construction of the data variance
      double minerr;
      //! The name of the starting model
      std::string modelfilename;
      //! The name of the file with the phase delays
      std::string datafilename;
      //! The weight for the tomography data in the joint inversion
      double swlambda;
      //! The tomography starting model
      jif3D::SurfaceWaveModel SWModel;
    public:
      //! Read only access to the starting model for surface wave tomography
      /*! Read only access to the starting model for surface wave tomography
       * @return The object containing the starting model.
       */
      const jif3D::SurfaceWaveModel &GetModel() const
        {
          return SWModel;
        }
      //! read only access to the objective function object for surface wave tomography data
      const jif3D::ThreeDModelObjective<jif3D::SurfaceWaveCalculator>& GetSurfaceWaveObjective()
        {
          return *SurfaceWaveObjective;
        }
      //! Setup the program options for the tomography part of the inversion
      po::options_description SetupOptions();
      //! Setup the tomography objective function
      /*! Setup the objective function for inverting seismic surface wave data based on
       * the program options.
       * @param vm The variable map from boost::program_options that contains the actually set options
       * @param Objective An existing JointObjective object that the newly created tomography objective is added to
       * @param Transform A transformation object to transform generalized to physical parameters
       * @param xorigin The origin for the inversion grid in x-direction
       * @param yorigin The origin for the inversion grid in y-direction
       * @return True if the weight the tomography objective is greater zero, i.e. we added an objective function to JointObjective, false otherwise
       */
      bool
      SetupObjective(const po::variables_map &vm, jif3D::JointObjective &Objective,
          boost::shared_ptr<jif3D::GeneralModelTransform> Transform, double xorigin = 0.0,
          double yorigin = 0.0);
      SetupSW();
      virtual ~SetupSW();
      };
}
#endif /* JOINT_SETUPSW_H_ */
