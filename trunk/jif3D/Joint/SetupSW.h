//
// Created by wangchao on 2022/4/14.
//

#ifndef JOINT_SETUPDSW_H_
#define JOINT_SETUPDSW_H_

#include <boost/program_options.hpp>
#include "GeneralDataSetup.h"

#include "../Global/Jif3DGlobal.h"
#include "../Inversion/JointObjective.h"
#include "../Inversion/ThreeDModelObjective.h"
#include "../SurfaceWaves/SurfaceWaveCalculator.h"
#include "../SurfaceWaves/SurfaceWaveModel.h"

namespace jif3D
  {
    namespace po = boost::program_options;

    class J3DEXPORT SetupSW : public GeneralDataSetup
      {
    private:
      boost::shared_ptr<jif3D::ThreeDModelObjective<jif3D::SurfaceWaveCalculator> > DSurfaceWaveObjective;
      double relerr;
      double minerr;
      double minvs;
      double maxvs;
      std::string modelfilename;
      std::string datafilename;
      double dswlambda;
      jif3D::SurfaceWaveModel DSurfaceWaveModel;
    public:
      const jif3D::SurfaceWaveModel& GetModel() const
        {
          return DSurfaceWaveModel;
        }
      //! read only access to the objective function object for surface wave tomography data
      const jif3D::ThreeDModelObjective<jif3D::SurfaceWaveCalculator>& GetDSurfaceWaveObjective()
        {
          return *DSurfaceWaveObjective;
        }
      //! Setup the program options for the tomography part of the inversion
      po::options_description SetupOptions();
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
      SetupSW();
      virtual ~SetupSW();

      };

  }

#endif //JOINT_SETUPDSW_H
