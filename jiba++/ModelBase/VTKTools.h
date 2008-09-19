//============================================================================
// Name        : VTKTools.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef VTKTOOLS_H_
#define VTKTOOLS_H_

#include "ThreeDModelBase.h"
#include "../Gravity/ThreeDGravityModel.h"
#include <string>
namespace jiba
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */

    //! Write a 3D model and its geometry in a .vtk file for plotting
    void Write3DModelToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data);

    //! Write scalar data with 3D coordinate information in a .vtk file for plotting
    void Write3DDataToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDGravityModel::tScalarMeasVec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ);

    //! Write \f$ 3 \times 3\f$ tensor data with 3D coordinate information in a .vtk file for plotting
    void Write3DTensorDataToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDGravityModel::tTensorMeasVec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ);
  /* @} */

  }
#endif /* VTKTOOLS_H_ */
