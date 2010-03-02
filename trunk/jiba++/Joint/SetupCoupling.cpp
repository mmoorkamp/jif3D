//============================================================================
// Name        : SetupCoupling.cpp
// Author      : Mar 2, 2010
// Version     : 
// Copyright   : 2010, mmoorkamp
//============================================================================


#include "SetupCoupling.h"

namespace jiba
  {

    SetupCoupling::SetupCoupling()
      {
      }

    SetupCoupling::~SetupCoupling()
      {
      }

    po::options_description SetupCoupling::SetupOptions()
      {
        po::options_description desc("Gravity options");

        return desc;
      }

    void SetupCoupling::ConfigureCoupling(
        boost::shared_ptr<jiba::GeneralModelTransform> &InversionTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &TomoTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &MTTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &GravityTransform,
        boost::shared_ptr<jiba::GeneralModelTransform> &RegTransform)
      {
        const double minslow = 1e-4;
        const double maxslow = 0.005;

        TomoTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
            new jiba::TanhTransform(minslow, maxslow));
        GravityTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
            new jiba::DensityTransform(TomoTransform));
        MTTransform = boost::shared_ptr<jiba::GeneralModelTransform>(
            new jiba::ConductivityTransform(TomoTransform));

        InversionTransform = TomoTransform;
        RegTransform = TomoTransform;
      }
  }
