/*
 * Entropy.h
 *
 *  Created on: 20 Oct 2020
 *      Author: moorkamp
 */

#ifndef MI_ENTROPY_H_
#define MI_ENTROPY_H_

#include "../Global/Serialization.h"
#include "../Global/VecMat.h"

namespace jif3D
  {
    double shan_entropy(const jif3D::rvec &x);
    jif3D::rvec diff_shan_entropy(const jif3D::rvec &x);
  }
#endif /* MI_ENTROPY_H_ */
