/*
 * ReadAnyModel.h
 *
 *  Created on: 7 Apr 2015
 *      Author: mm489
 */

#ifndef READANYMODEL_H_
#define READANYMODEL_H_

#include <string>
#include <boost/shared_ptr.hpp>
#include "../ModelBase/ThreeDModelBase.h"

namespace jif3D
  {
    boost::shared_ptr<jif3D::ThreeDModelBase> ReadAnyModel(const std::string &Filename);

  }
#endif /* READANYMODEL_H_ */
