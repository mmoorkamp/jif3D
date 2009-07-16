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
#include <vector>
#include <complex>
#include <boost/multi_array.hpp>

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
    void WriteProjectFile(const std::vector<double> &Frequencies,
        X3DModel::ProblemType Type, const std::string &ResultFilename,
        const std::string &ModelFilename);

    //!Write a file containing information about sources (electric or magnetic)
    void WriteSourceFile(const std::string &filename, const double SourceDepth,
        const boost::multi_array<std::complex<double>, 2> &XPolMoments,
        const boost::multi_array<std::complex<double>, 2> &YPolMoments);
    void ReadEMO(const std::string &filename,
        std::vector<std::complex<double> > &Ex, std::vector<
            std::complex<double> > &Ey, std::vector<std::complex<double> > &Hx,
        std::vector<std::complex<double> > &Hy);

    void ReadEMA(const std::string &filename,
        std::vector<std::complex<double> > &Ex, std::vector<
            std::complex<double> > &Ey, std::vector<std::complex<double> > &Ez);

    std::vector<std::complex<double> > ResortFields(const std::vector<
        std::complex<double> > &InField, const size_t nx, const size_t ny,
        const size_t nz);
  }

#endif /* READWRITEX3D_H_ */
