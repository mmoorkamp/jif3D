//============================================================================
// Name        : ThreeDGravityModel.h
// Author      : Max Moorkamp
// Version     :
// Copyright   : 2008, MM
//============================================================================

#ifndef THREEDGRAVITYMODEL_H_
#define THREEDGRAVITYMODEL_H_

#include <string>
#include <boost/numeric/ublas/vector.hpp>
#include "../ModelBase/ThreeDModelBase.h"
#include "../Global/VecMat.h"

/*! \file ThreeDGravityModel.h
 * This file contains the class ThreeDGravityModel and associated helper functions and constants
 *
 */
namespace jiba
  {

    /** \addtogroup gravity Gravity forward modelling, display and inversion
     *
     * This module contains a library of functions and some programs for scalar gravity and FTG
     * forward modelling and inversion.
     *
     * There are three big building blocks for the gravity forward modelling. The class ThreeDGravityModel
     * stores a gridded density model including its geometry, a 1D layered background and the position of measurements.
     * In addition it contains some functionality to read and write models to and from files for storage and plotting.
     *
     * There are two class hierarchies that work together to manage the calculation of gravity data. The classes
     * derived from ThreeDGravityImplementation determine what kind of data is calculated, e.g. scalar or FTG, and
     * how it is calculated, e.g. on the CPU using parallelization or on a GPU using the CUDA interface. These classes
     * are not intended for direct use however, but wrapped by the classes derived from ThreeDGravityCalculator.
     *
     * This hierarchy of classes manages what kind of additional data is stored to accelerate the
     * forward calculation. For example, MinMemGravityCalculator just forwards to its Implementation class
     * and does not store any sensitivities. This means that at every call the data is calculated completely new.
     * The classes derived from CachedGravityCalculator in contrast store the sensitivity matrix or a processed version.
     * At each forward calculation they examine the model and whether it requires new sensitivities. If it does they forward
     * the call to the Implementation object, otherwise they perform a simple matrix-vector product to quickly obtain
     * the data.
     */
    /* @{ */


    //! We store the 3x3 matrix for gravimetric measurements in a ublas matrix with real entries
    typedef rmat GravimetryMatrix;

    //! The class used to store the gravity model and the location of the measurement points
    class ThreeDGravityModel: public ThreeDModelBase
      {
    public:
      typedef std::vector<double> tBackgroundVec;
      typedef std::vector<double> tMeasPosVec;
    private:
      //! Create a dimension for the measurement positions in a netcdf file
      NcDim *WriteDimensionToNetCDF(NcFile &NetCDFFile,
          const std::string &SizeName, const tMeasPosVec &Position) const;
      //! the x-coordinates of the measurement points
      tMeasPosVec MeasPosX;
      //! the y-coordinates of the measurement points
      tMeasPosVec MeasPosY;
      //! the z-coordinates of the measurement points
      tMeasPosVec MeasPosZ;
      //! The densities of the background layers
      tBackgroundVec bg_densities;
      //! The thicknesses of the background layers
      tBackgroundVec bg_thicknesses;
      //! Write out the values for a the measurement to an ascii file
      void
          PlotMeasAscii(const std::string &filename, rvec &Data) const;
    public:
      //! return read only access to the stored density values
      const t3DModelData &GetDensities() const
        {
          return ThreeDModelBase::GetData();
        }
      //! return a reference to stored densities
      t3DModelData &SetDensities()
        {
          return ThreeDModelBase::SetData();
        }
      //! Set the density of the background, it extends to infinity in horizontal directions and to the depth specified by the thicknesses in vertical direction
      void SetBackgroundDensities(const tBackgroundVec value)
        {
          bg_densities.clear();
          copy(value.begin(), value.end(), back_inserter(bg_densities));
        }
      //! Return the densities of the background layers
      const tBackgroundVec &GetBackgroundDensities() const
        {
          return bg_densities;
        }
      //! Set the thicknesses of the background layers, the individual thicknesses are given in m the total thickness of the background layers does not need to coincide with the gridded domain
      void SetBackgroundThicknesses(const tBackgroundVec value)
        {
          bg_thicknesses.clear();
          copy(value.begin(), value.end(), back_inserter(bg_thicknesses));
        }
      //! Return the thicknesses of the background layers in m
      const tBackgroundVec &GetBackgroundThicknesses() const
        {
          return bg_thicknesses;
        }
      //! Add a measurement point to the model
      void AddMeasurementPoint(const double xcoord, const double ycoord,
          const double zcoord)
        {
          MeasPosX.push_back(xcoord);
          MeasPosY.push_back(ycoord);
          MeasPosZ.push_back(zcoord);
        }
      //! remove all information about measurement points
      void ClearMeasurementPoints()
        {
          MeasPosX.clear();
          MeasPosY.clear();
          MeasPosZ.clear();
        }
      //! Return the x-coordinates of all measurement points read-only
      /*! This function provides read-only access to the x-coordinates
       * of the measurement points. The only way to modify the position of
       * the measurements is to delete them with ClearMeasurementPoints and
       * add new ones with AddMeasurementPoint. This ensures that we have all
       * three coordinate values for all points.
       * @return A vector with the x-coordinates of all measurement points in m
       */
      const tMeasPosVec &GetMeasPosX() const
        {
          return MeasPosX;
        }
      //! Return the y-coordinates of all measurement points read-only
      const tMeasPosVec &GetMeasPosY() const
        {
          return MeasPosY;
        }
      //! Return the z-coordinates of all measurement points read-only
      const tMeasPosVec &GetMeasPosZ() const
        {
          return MeasPosZ;
        }
      //! Write the density model and all associated information in a netcdf file
      void WriteNetCDF(const std::string filename) const;
      //! Write the density model in VTK format, at the moment the best format for plotting
      void WriteVTK(const std::string filename)
        {
          ThreeDModelBase::WriteVTK(filename, "Density");
        }
      //! Read the density model and all associated information from a netcdf file
      void ReadNetCDF(const std::string filename);
      //! Read an igmas xyz model file
      void ReadIgmas(const std::string filename);
      //! Save all scalar Synthetic data in a netcdf file
      void SaveScalarMeasurements(const std::string filename);
      //! Save all tensor Synthetic data in a netcdf file
      void SaveTensorMeasurements(const std::string filename);
      //! Write all scalar synthetic data in an ascii file for spatial plotting in gmt
      void PlotScalarMeasurements(const std::string filename);
      //! Write all tensor synthetic data in an ascii file for spatial plotting in gmt
      void PlotTensorMeasurements(const std::string filename_root);
      //! Read the Measurement positions from a netcdf file
      void ReadMeasPosNetCDF(const std::string filename);
      //! Read the Measurement positions from an ascii file
      void ReadMeasPosAscii(const std::string filename);
      ThreeDGravityModel();
      virtual ~ThreeDGravityModel();
      };
  /* @} */
  }

#endif /*THREEDGRAVITYMODEL_H_*/
