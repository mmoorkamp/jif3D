/*
 * GeneralDataSetup.h
 *
 *  Created on: May 30, 2022
 *      Author: max
 */

#ifndef JOINT_GENERALDATASETUP_H_
#define JOINT_GENERALDATASETUP_H_

#include "../Global/Jif3DGlobal.h"
#include "../ModelTransforms/MultiSectionTransform.h"
#include "../ModelBase/ThreeDModelBase.h"
#include "../Inversion/JointObjective.h"
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <string>

#ifdef HAVEHPX
#include <hpx/modules/program_options.hpp>
#endif
#ifdef HAVEOPENMP
#include <boost/program_options.hpp>
#endif

#ifdef HAVEHPX
namespace po = hpx::program_options;
#endif
#ifdef HAVEOPENMP
namespace po = boost::program_options;
#endif

namespace jif3D
  {



    class GeneralDataSetup
      {
    private:
      std::string ParameterName;
    protected:
      jif3D::rvec StartingParameters;
      boost::shared_ptr<jif3D::MultiSectionTransform> Transform;
    public:
      enum  parametertype
        {
        gridparameter, other
        };
      //! Access to the underlying parameter transform
      boost::shared_ptr<jif3D::MultiSectionTransform> GetTransform()
        {
          return Transform;
        }
      std::string GetParameterName()
        {
          return ParameterName;
        }
      jif3D::rvec GetStartingParameters()
        {
          return StartingParameters;
        }
      virtual po::options_description SetupOptions() = 0;
      virtual bool
      SetupObjective(const po::variables_map &vm,
          jif3D::JointObjective &Objective, jif3D::ThreeDModelBase &InversionMesh,
          jif3D::rvec &CovModVec, std::vector<size_t> &startindices,
          std::vector<std::string> &SegmentNames,
          std::vector<parametertype> &SegmentTypes, boost::filesystem::path TempDir =
              boost::filesystem::current_path()) = 0;
      virtual void IterationOutput(const std::string &filename,
          const jif3D::rvec &ModelVector) = 0;
      virtual void FinalOutput(const std::string &filename ,const jif3D::rvec &FinalModelVector) =0;
      GeneralDataSetup(std::string Name);
      virtual ~GeneralDataSetup();
      };

  } /* namespace jif3D */

#endif /* JOINT_GENERALDATASETUP_H_ */
