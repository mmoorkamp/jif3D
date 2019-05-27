//============================================================================
// Name        : VTKTools.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef VTKTOOLS_H_
#define VTKTOOLS_H_

#include <string>

#include "ThreeDModelBase.h"

#include "../Global/VecMat.h"
#include "../Global/Jif3DGlobal.h"

/*! \file VTKTools.h
 * A collection of function to write 3D models and data with 3D positioning information
 * to files in simple .vtk format for plotting with Paraview or Visit
 */
namespace jif3D
  {
    /** \addtogroup modelbase Basic classes and routines for 3D models */
    /* @{ */

    //! Write a 3D model and its geometry into a .vtk file for plotting
    J3DEXPORT void Write3DModelToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDModelBase::t3DModelDim &XCellCoords,
        const ThreeDModelBase::t3DModelDim &YCellCoords,
        const ThreeDModelBase::t3DModelDim &ZCellCoords,
        const ThreeDModelBase::t3DModelData &Data);

    //! Write a 3D model with vector valued cells and its geometry into a .vtk file for plotting
    J3DEXPORT void Write3DVectorModelToVTK(const std::string &filename,
        const std::string &DataName,
        const ThreeDModelBase::t3DModelDim &XCellCoords,
        const ThreeDModelBase::t3DModelDim &YCellCoords,
        const ThreeDModelBase::t3DModelDim &ZCellCoords,
        const ThreeDModelBase::t3DModelData &CompX,
        const ThreeDModelBase::t3DModelData &CompY,
        const ThreeDModelBase::t3DModelData &CompZ);

    //! Read a 3D model and its geometry from a .vtk file
    J3DEXPORT void Read3DModelFromVTK(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellCoords,
        ThreeDModelBase::t3DModelDim &YCellCoords,
        ThreeDModelBase::t3DModelDim &ZCellCoords,
        ThreeDModelBase::t3DModelData &Data);

    J3DEXPORT void Read3DVectorModelFromVTK(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellCoords,
        ThreeDModelBase::t3DModelDim &YCellCoords,
        ThreeDModelBase::t3DModelDim &ZCellCoords, ThreeDModelBase::t3DModelData &CompX,
         ThreeDModelBase::t3DModelData &CompY, ThreeDModelBase::t3DModelData &CompZ);

    //! Write scalar data with 3D coordinate information into a .vtk file for plotting
    J3DEXPORT void Write3DDataToVTK(const std::string &filename,
        const std::string &DataName,
        const jif3D::rvec &Data,
        const std::vector<double> &PosX,
        const std::vector<double> &PosY,
        const std::vector<double> &PosZ);

    //! Write three-component vector data with 3D coordinate information into a .vtk file for plotting
    J3DEXPORT void Write3DVectorDataToVTK(const std::string &filename,
        const std::string &DataName,
        const jif3D::rvec &Data,
        const std::vector<double> &PosX,
        const std::vector<double> &PosY,
        const std::vector<double> &PosZ);

    //! Write \f$ 3 \times 3\f$ tensor data with 3D coordinate information into a .vtk file for plotting
    J3DEXPORT void Write3DTensorDataToVTK(const std::string &filename,
        const std::string &DataName,
        const jif3D::rvec &Data,
        const std::vector<double> &PosX,
        const std::vector<double> &PosY,
        const std::vector<double> &PosZ);
  /* @} */

  }
#endif /* VTKTOOLS_H_ */
