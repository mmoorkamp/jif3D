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

/*! \file This file contains the class ThreeDGravityModel and associated helper functions and constants
 *
 */
namespace jiba
  {

    //! The constant of gravity
    static const double Grav_const = 6.67428e-8; // in units cm^3/g s

    //! We store the 3x3 matrix for gravimetric measurements in a ublas matrix with real entries
    typedef rmat GravimetryMatrix;

    //! Calculate a single geometric term for the graviational acceleration due to a rectangular prism
    double CalcGravTerm(const double x, const double y, const double z);
    //! Calculate the geometric term  due to a rectangular prism
    double CalcGravBoxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    //! Calculate the geometic part of the gravimetry matrix for a single rectangular prism
    GravimetryMatrix CalcTensorBoxTerm(const double meas_x,
        const double meas_y, const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    double CalcUxxTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    double CalcUxyTerm(const double meas_x, const double meas_y,
        const double meas_z, const double ul_corner_x,
        const double ul_corner_y, const double ul_corner_z,
        const double x_size, const double y_size, const double z_size);

    //! the gravitational acceleration of an semi-infinite slab
    double CalcGravSemiInfSheet(const double hor_dist, const double ver_dist,
        const double thick, const double density);
    //! Calculate one of the terms of diagonal elements of the gravimetric matrxi
    double CalcFTGDiagonalTerm(const double a, const double b, const double c);
    //! Calculate one of the terms of off-diagonal elements of the gravimetric matrxi
    double CalcFTGOffDiagonalTerm(const double value, const double x,
        const double y, const double z);
    //! the gravitational acceleration of an infinite slab
    inline double CalcInfSheetTerm(const double thick)
      {
        return 2.0 * M_PI * Grav_const * thick;
      }

    void ConstructDepthWeighting(const ThreeDModelBase::t3DModelDim &XSizes, const ThreeDModelBase::t3DModelDim &YSizes,
        const ThreeDModelBase::t3DModelDim &ZSizes, const double z0, rvec &WeightVector);
    //! We store the 3x3 matrix for gravimetric measurements in a ublas matrix with real entries
    typedef rmat GravimetryMatrix;

    //! The class used to store the gravity model and calculate the gravity values at the measurement points
    class ThreeDGravityModel: public ThreeDModelBase
      {
    public:
      typedef std::vector<double> tScalarMeasVec;
      typedef std::vector<GravimetryMatrix> tTensorMeasVec;
      typedef std::vector<double> tMeasPosVec;
    private:
      //! Calculate the response of the 1D background
      double CalcBackground(const double xmeas, const double ymeas,
          const double zmeas, const double xwidth, const double ywidth,
          const double zwidth, const size_t meas_index);
      //! For the given model, calculate the gravity at measurement points x, y and z due to the discretized part of the mesh
      double CalcScalarMeas(const double x_meas, const double y_meas,
          const double z_meas, const size_t meas_index);
      //! For the given model, calculate the gravimetric matrix at points x,y, and z
      GravimetryMatrix CalcTensorMeas(const double x_meas, const double y_meas,
          const double z_meas, const size_t meas_index);
      //! Correct for the fact that the grid does not go to infinity
      GravimetryMatrix AdjustTensorBackground(const double x_meas,
          const double y_meas, const double z_meas, const double xwidth,
          const double ywidth, const double zwidth, const size_t meas_index);
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
      tScalarMeasVec bg_densities;
      //! The thicknesses of the background layers
      tScalarMeasVec bg_thicknesses;
      //! Do we want to store the sensitivity matrix for scalar gravity
      const bool StoreScalarSensitivities;
      //! Do we want to store the sensitivity matrix for tensor gravity
      const bool StoreTensorSensitivities;
      //! Have we already calculated the scalar sensitivities
      bool HaveCalculatedScalarSensitivities;
      //! Have we already calculated the tensor sensitivities
      bool HaveCalculatedTensorSensitivities;
      //! The Matrix for the scalar sensitivities
      rmat ScalarSensitivities;
      //! The Matrix for the tensor sensitivities
      rmat TensorSensitivities;
      //!We store the synthetic tensor data so we can plot or save it without recalculation
      tTensorMeasVec TensorResults;
      //!We store the synthetic scalar data so we can plot or save it without recalculation
      tScalarMeasVec ScalarResults;
      //! We implement this virtual function to make sure sensitivities are recalculated when the cell sizes change
      virtual void SetCellSizesAction(t3DModelDim &sizes)
        {
          HaveCalculatedScalarSensitivities = false;
          HaveCalculatedTensorSensitivities = false;
        }
      //! Write out the values for a the measurement to an ascii file
      void
          PlotMeasAscii(const std::string &filename, tScalarMeasVec &Data) const;
    public:
      //! For the given model, calculate the scalar gravity at all measurement points
      tScalarMeasVec CalcGravity();
      //! For the given model, calculate the FTG matrix at all measurement points
      tTensorMeasVec CalcTensorGravity();
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
      //! Set the density of the background, it extends to infinity in horizontal directions and to the depth of the model in z-direction
      void SetBackgroundDensities(const tScalarMeasVec value)
        {
          bg_densities.clear();
          copy(value.begin(), value.end(), back_inserter(bg_densities));
        }
      //! Set the thicknesses of the background layers, the individual thicknesses are given in m
      void SetBackgroundThicknesses(const tScalarMeasVec value)
        {
          bg_thicknesses.clear();
          copy(value.begin(), value.end(), back_inserter(bg_thicknesses));
        }
      //! Add a measurement point to the model
      void AddMeasurementPoint(const double xcoord, const double ycoord,
          const double zcoord)
        {
          MeasPosX.push_back(xcoord);
          MeasPosY.push_back(ycoord);
          MeasPosZ.push_back(zcoord);
          HaveCalculatedScalarSensitivities = false; // we have to recalculate
          HaveCalculatedTensorSensitivities = false;
        }
      //! remove all information about measurement points
      void ClearMeasurementPoints()
        {
          MeasPosX.clear();
          MeasPosY.clear();
          MeasPosZ.clear();
          //the sensitivities will no longer be valid
          HaveCalculatedScalarSensitivities = false;
          HaveCalculatedTensorSensitivities = false;
        }
      //! Return the x-coordinates of all measurement points read-only
      /*! This function provides read-only access to the x-coordinates
       * of the measurement points. The only way to modify the position of
       * the measurements is to delete them with ClearMeasurementPoints and
       * add new ones with AddMeasurementPoint. This ensures that we have all
       * three coordinate values for all points.
       * @return A vector with the x-coordinates of all measurement points in m
       */
      const tMeasPosVec &GetMeasPosX()
        {
          return MeasPosX;
        }
      //! Return the y-coordinates of all measurement points read-only
      const tMeasPosVec &GetMeasPosY()
        {
          return MeasPosY;
        }
      //! Return the z-coordinates of all measurement points read-only
      const tMeasPosVec &GetMeasPosZ()
        {
          return MeasPosZ;
        }
      //! Get the sensitivity matrix for scalar gravity measurements
      const rmat &GetScalarSensitivities() const
        {
          return ScalarSensitivities;
        }
      //! Get the sensitivity matrix for tensor gravity measurements
      const rmat &GetTensorSensitivities() const
        {
          return TensorSensitivities;
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
      //! When we create the object we have to specify whether we want to store scalar and/or tensor sensitivities
      ThreeDGravityModel(const bool storescalar = false,
          const bool storetensor = false);
      virtual ~ThreeDGravityModel();
      };

  }

#endif /*THREEDGRAVITYMODEL_H_*/
