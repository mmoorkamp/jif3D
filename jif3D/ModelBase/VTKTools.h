//============================================================================
// Name        : VTKTools.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef VTKTOOLS_H_
#define VTKTOOLS_H_

#include "ThreeDModelBase.h"
#include <string>
#include "../Global/VecMat.h"

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
        const ThreeDModelBase::t3DModelData &Data,
        double xorigin = 0, double yorigin = 0, double zorigin = 0);

    //! Write a 3D model with vector valued cells and its geometry into a .vtk file for plotting
    void Write3DVectorModelToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &CompX,
        const ThreeDModelBase::t3DModelData &CompY,
        const ThreeDModelBase::t3DModelData &CompZ);

    //! Read a 3D model and its geometry from a .vtk file
    void Read3DModelFromVTK(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes,
        ThreeDModelBase::t3DModelData &Data);

    //! Write scalar data with 3D coordinate information into a .vtk file for plotting
    void Write3DDataToVTK(const std::string &filename,
        const std::string &DataName,
        const jif3D::rvec &Data,
        const ThreeDModelBase::tMeasPosVec &PosX,
        const ThreeDModelBase::tMeasPosVec &PosY,
        const ThreeDModelBase::tMeasPosVec &PosZ);

    //! Write three-component vector data with 3D coordinate information into a .vtk file for plotting
    void Write3DVectorDataToVTK(const std::string &filename,
        const std::string &DataName,
        const jif3D::rvec &Data,
        const ThreeDModelBase::tMeasPosVec &PosX,
        const ThreeDModelBase::tMeasPosVec &PosY,
        const ThreeDModelBase::tMeasPosVec &PosZ);

    //! Write \f$ 3 \times 3\f$ tensor data with 3D coordinate information into a .vtk file for plotting
    void Write3DTensorDataToVTK(const std::string &filename,
        const std::string &DataName,
        const jif3D::rvec &Data,
        const ThreeDModelBase::tMeasPosVec &PosX,
        const ThreeDModelBase::tMeasPosVec &PosY,
        const ThreeDModelBase::tMeasPosVec &PosZ);
  /* @} */

  }
#endif /* VTKTOOLS_H_ */
