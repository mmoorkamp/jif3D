//============================================================================
// Name        : ReadWriteX3D.h
// Author      : Feb 17, 2009
// Version     :
// Copyright   : 2009, mmoorkamp
//============================================================================

#ifndef READWRITEX3D_H_
#define READWRITEX3D_H_

#include <string>
#include <vector>
#include <complex>
#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>
#include "../ModelBase/ThreeDModelBase.h"
#include "X3DModel.h"

namespace jif3D
  {
    /** \addtogroup mtmodelling Forward modelling of magnetotelluric data */
    /* @{ */
    //! Read an ascii 3D MT model in the format used by x3D, not feature complete !
    /*! Read in a file with a 3D MT model in the format used by x3D. Note that
     * we do not support all possible models that can be described in the model
     * file, but only the ones that are relevant to our joint inversion. For
     * example we assume that there is one modeling domain, all layers are located
     * next to each other and each layer has the same number of cells.
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
        ThreeDModelBase::t3DModelDim &ZCellSizes, ThreeDModelBase::t3DModelData &Data,
        std::vector<double> &bg_conductivities, std::vector<double> &bg_thicknesses);
    //! Write a 3D model to an ascii file compatible with x3D
    /*! Write a file with a 3D MT model in the format used by x3D.
     * @param filename The name of the file
     * @param XCellSizes The size of the cells in x-direction [m]
     * @param YCellSizes The size of the cells in y-direction [m]
     * @param ZCellSizes The size of the cells in z-direction [m]
     * @param ObservationDepths The depth to each of the observation sites in m
     * @param Data The conductivities in Siemens
     * @param bg_conductivities The conductivities in Siemens of the background layers
     * @param bg_thicknesses The thicknesses in m of the background layers
     */
    void Write3DModelForX3D(const std::string &filename,
        const ThreeDModelBase::t3DModelDim &XCellSizes,
        const ThreeDModelBase::t3DModelDim &YCellSizes,
        const ThreeDModelBase::t3DModelDim &ZCellSizes,
        const std::vector<double> &ObservationDepths,
        const ThreeDModelBase::t3DModelData &Data,
        const std::vector<double> &bg_conductivities,
        const std::vector<double> &bg_thicknesses);

    //! Write the file a.project that controls the forward calculation parameters for x3D
    /*! The project file is always called a.project and controls what x3d calculates.
     * @param RootDir The start of the directory name without the frequency index
     * @param Frequencies The vector of frequencies in Hz for which we want to calculate data
     * @param Type The type of data we want to calculate, MT, electric dipole or magnetic dipole fields
     * @param ResultFilename The filename root for the results
     * @param ModelFilename The filename for the model in x3d format
     */
    void WriteProjectFile(const boost::filesystem::path &RootDir,
        const std::vector<double> &Frequencies, X3DModel::ProblemType Type,
        const std::string &ResultFilename, const std::string &ModelFilename);

    //!Write a file containing information about sources (electric or magnetic)
    /*! Write a file that contains the dipole moments of electric or magnetic dipoles
     * for x3d.
     * @param filename The name of the ascii file we want to write
     * @param SourceXIndex The indices in x-direction within the model for each source
     * @param SourceYIndex The indices in y-direction within the model for each source
     * @param SourceDepths The depth of the dipoles in m
     * @param XPolMoments The x-component of the dipole moment for each source position
     * @param YPolMoments The y-component of the dipole moment for each source position
     * @param ZCellBoundaries The cell boundaries of the modelling domain in z-direction
     * @param ZCellSizes The cell sizes of the modelling domain in z-direction
     * @param maxx The number of cells of the modelling domain in x-direction
     * @param maxy The number of cells of the modelling domain in y-direction
     */
    void WriteSourceFile(const std::string &filename,
        const std::vector<size_t> &SourceXIndex, const std::vector<size_t> &SourceYIndex,
        const std::vector<double> &SourceDepths,
        const std::vector<std::complex<double> > &XPolMoments,
        const std::vector<std::complex<double> > &YPolMoments,
        const jif3D::ThreeDModelBase::t3DModelDim &ZCellBoundaries,
        const jif3D::ThreeDModelBase::t3DModelDim &ZCellSizes,
        const size_t maxx,
        const size_t maxy);

    //! Read one .emo file produced by x3d that contains the electric and magnetic fields at the observation sites
    /*! x3d produces files with ending .emo that contain the electric and magnetic fields at the
     * observation points that were specified for the forward calculation. This function reads in
     * the contents of the file.
     * @param filename The complete filename including the ending
     * @param Ex The x-component of the electric field in V/m
     * @param Ey The y-component of the electric field in V/m
     * @param Hx The x-component of the magnetic field in A/m
     * @param Hy The y-component of the magnetic field in A/m
     */
    void ReadEMO(const std::string &filename, std::vector<std::complex<double> > &Ex,
        std::vector<std::complex<double> > &Ey, std::vector<std::complex<double> > &Hx,
        std::vector<std::complex<double> > &Hy);

    //! Read one .ema file produced by x3d that contains the electric field at all model cells
    /*! x3d produces files with ending .ema that contain the electric field in all model cells.
     * This information is needed to calculate the gradient using the adjoint approach. This function reads in
     * the contents of the file.
     *
     * @param filename The complete filename including the ending
     * @param Ex The x-component of the electric field at all cells in V/m
     * @param Ey The y-component of the electric field at all cells in V/m
     * @param Ez The z-component of the electric field at all cells in V/m
     * @param ncellsx The number of cells in x-direction
     * @param ncellsy The number of cells in y-direction
     * @param ncellsz The number of cells in z-direction
     */
    void ReadEMA(const std::string &filename, std::vector<std::complex<double> > &Ex,
        std::vector<std::complex<double> > &Ey, std::vector<std::complex<double> > &Ez,
        const size_t ncellsx, const size_t ncellsy, const size_t ncellsz);

  /* @} */
  }

#endif /* READWRITEX3D_H_ */
