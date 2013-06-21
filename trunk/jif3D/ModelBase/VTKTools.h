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
/*! \file VTKTools.h
 * A collection of function to write 3D models and data with 3D positioning information
 * to files in simple .vtk format for plotting with Paraview or Visit
 */
namespace jif3D
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */

    //! Write a 3D model and its geometry into a .vtk file for plotting
    void Write3DModelToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data);

    //! Write scalar data with 3D coordinate information into a .vtk file for plotting
    void Write3DDataToVTK(const std::string &filename,
        const std::string &DataName,
        const jif3D::rvec &Data,
        const ThreeDModelBase::tMeasPosVec &PosX,
        const ThreeDModelBase::tMeasPosVec &PosY,
        const ThreeDModelBase::tMeasPosVec &PosZ);

    //! Write \f$ 3 \times 3\f$ tensor data with 3D coordinate information into a .vtk file for plotting
    void Write3DTensorDataToVTK(const std::string &filename,
        const std::string &DataName,
        const jif3D::rvec &Data,
        const ThreeDGravityModel::tMeasPosVec &PosX,
        const ThreeDGravityModel::tMeasPosVec &PosY,
        const ThreeDGravityModel::tMeasPosVec &PosZ);
  /* @} */

  }
#endif /* VTKTOOLS_H_ */
