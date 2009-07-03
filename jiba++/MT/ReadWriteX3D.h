//============================================================================
// Name        : ReadWriteX3D.h
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================


#ifndef READWRITEX3D_H_
#define READWRITEX3D_H_

#include "../ModelBase/ThreeDModelBase.h"
#include "X3DModel.h"
#include <string>

namespace jiba
  {

    //! Read an ascii 3D MT model in the format used by x3D, not feature complete !
    /*! Read in a file with a 3D MT model in the format used by x3D
     * @param filename The name of the file
     * @param XCellSizes The size of the cells in x-direction [m]
     * @param YCellSizes The size of the cells in y-direction [m]
     * @param ZCellSizes The size of the cells in z-direction [m]
     * @param Data The conductivities in Siemens
     * @param bg_conductivities The conductivities in Siemens of the background layers
     * @param bg_thicknesses The thicknesses in m of the background layers
     */
    void Read3DModelFromX3D(const std::string &filename,
        ThreeDModelBase::t3DModelDim &XCellSizes,
        ThreeDModelBase::t3DModelDim &YCellSizes,
        ThreeDModelBase::t3DModelDim &ZCellSizes,
        ThreeDModelBase::t3DModelData &Data,
        std::vector<double> &bg_conductivities,
        std::vector<double> &bg_thicknesses);
    //! Write a 3D model to an ascii file compatible with x3D
    /*! Write a file with a 3D MT model in the format used by x3D
     * @param filename The name of the file
     * @param XCellSizes The size of the cells in x-direction [m]
     * @param YCellSizes The size of the cells in y-direction [m]
     * @param ZCellSizes The size of the cells in z-direction [m]
     * @param Data The conductivities in Siemens
     * @param bg_conductivities The conductivities in Siemens of the background layers
     * @param bg_thicknesses The thicknesses in m of the background layers
     */
    void Write3DModelForX3D(const std::string &filename,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const ThreeDModelBase::t3DModelData &Data,
        const std::vector<double> &bg_conductivities,
        const std::vector<double> &bg_thicknesses);
    //! Write the file a.project that controls the forward calculation parameters for x3D
    void WriteProjectFile(std::vector<double> &Frequencies,
        X3DModel::ProblemType Type, const std::string &ResultFilename,
        const std::string &ModelFilename);
  }

#endif /* READWRITEX3D_H_ */
