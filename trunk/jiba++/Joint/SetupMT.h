//============================================================================
// Name        : SetupMT.h
// Author      : Mar 1, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#ifndef SETUPMT_H_
#define SETUPMT_H_

#include "../MT/X3DModel.h"
#include "../Inversion/JointObjective.h"
#include "../Tomo/ThreeDSeismicModel.h"
#include <boost/program_options.hpp>

namespace jiba
  {
    namespace po = boost::program_options;
    class SetupMT
      {
    private:
      jiba::X3DModel MTModel;
    public:
      po::options_description SetupOptions();
      void
      SetupObjective(const po::variables_map &vm,
          jiba::JointObjective &Objective,
          const ThreeDSeismicModel &StartModel, boost::shared_ptr<
              jiba::GeneralModelTransform> Transform);
      const jiba::X3DModel &GetModel() {return MTModel;}
      SetupMT();
      virtual ~SetupMT();
      };

  }

#endif /* SETUPMT_H_ */
