/*
 * JointMethods.h
 *
 *  Created on: May 14, 2020
 *      Author: max
 */

#ifndef JOINT_JOINTMETHODS_H_
#define JOINT_JOINTMETHODS_H_
#include "../Global/Jif3DGlobal.h"
#include <string>
#include <vector>

namespace jif3D
  {
    const std::string MTMethodString = "MT";
    const std::string SurfaceWaveMethodString = "SW";
    const std::string TravelTimeMethodString = "TT";
    const std::string ScalarGravityMethodString = "ScalGrav";
    const std::string FTGMethodString = "FTG";
    const std::string MagneticsMethodString = "Mag";
    const std::string DCResistivityMethodString = "DC";

    const std::vector<string> MethodStrings =
      { MTMethodString, SurfaceWaveMethodString, TravelTimeMethodString,
          ScalarGravityMethodString, FTGMethodString, MagneticsMethodString,
          DCResistivityMethodString };

    enum class PhysicalProperties
      {
      Conductivity, SVelocity, PVelocity, Density, Susceptibility
      };
  }

#endif /* JOINT_JOINTMETHODS_H_ */
