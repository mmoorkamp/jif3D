/*
 * realinfo.h
 *
 *  Created on: 19 Aug 2014
 *      Author: mm489
 */

#ifndef REALINFO_H_
#define REALINFO_H_

#include "../Global/Serialization.h"
#include <string>
#include "../MT/X3DModel.h"

struct realinfo
  {
  jif3D::X3DModel Model;
  double bg_conductivity;
  double phase1cond;
  double phase2cond;
  double phase1frac;
  std::string tempdir;
  std::string x3dname;

  //! Provide serialization to be able to store objects
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
    {
      ar & Model;
      ar & bg_conductivity;
      ar & phase1cond;
      ar & phase2cond;
      ar & phase1frac;
      ar & tempdir;
      ar & x3dname;
    }
  realinfo(jif3D::X3DModel M, double bgc, double p1c, double p2c, double p1f,
      std::string td, std::string x3d) :
      Model(M), bg_conductivity(bgc), phase1cond(p1c), phase2cond(p2c), phase1frac(p1f), tempdir(
          td), x3dname(x3d)
    {

    }
  realinfo() :
      Model(), bg_conductivity(-1.0), phase1cond(-1.0), phase2cond(-1.0), phase1frac(
          -1.0), tempdir(), x3dname()
    {

    }
  };




#endif /* REALINFO_H_ */
