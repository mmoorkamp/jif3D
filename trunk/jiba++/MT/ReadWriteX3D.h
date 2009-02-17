//============================================================================
// Name        : ReadWriteX3D.h
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef READWRITEX3D_H_
#define READWRITEX3D_H_

#include "../ModelBase/ThreeDModelBase.h"
#include <string>

namespace jiba {

  //! Read an ascii 3D MT model in the format used by x3D, not feature complete !
  /*! Read in a file with a 3D MT model in the format used by x3D
   * @param filename The name of the file
   * @param XCellSizes The size of the cells in x-direction [m]
   * @param YCellSizes The size of the cells in y-direction [m]
   * @param ZCellSizes The size of the cells in z-direction [m]
   * @param Data The conductivities in Siemens
   */
  void Read3DModelFromX3D(const std::string &filename,
          ThreeDModelBase::t3DModelDim &XCellSizes,
          ThreeDModelBase::t3DModelDim &YCellSizes,
          ThreeDModelBase::t3DModelDim &ZCellSizes,
          ThreeDModelBase::t3DModelData &Data);

}


#endif /* READWRITEX3D_H_ */
