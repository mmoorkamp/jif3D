/*
 * GeneralData.h
 *
 *  Created on: May 20, 2019
 *      Author: max
 */

#ifndef DATABASE_GENERALDATA_H_
#define DATABASE_GENERALDATA_H_

#include "../Global/Serialization.h"
#include "../Global/FatalException.h"
#include "../Global/Jif3DGlobal.h"
#include <vector>
#include <string>
namespace jif3D
  {

    class GeneralData
      {
    public:
      //! The type of the measurement position vector, this is a std::vector because we want to easily append elements
      typedef std::vector<double> tMeasPosVec;
    private:
      //! the x-coordinates of the measurement points
      tMeasPosVec MeasPosX;
      //! the y-coordinates of the measurement points
      tMeasPosVec MeasPosY;
      //! the z-coordinates of the measurement points
      tMeasPosVec MeasPosZ;
      //! The object we use to store the data, the order and size per measurement position depends on the concrete data type
      std::vector<double> Data;
      //! The object we use to store the errors, always has the same size as the data vector
      std::vector<double> Errors;
    public:
      //! Provide serialization to be able to store objects and, more importantly for simpler  parallelization
      template<class Archive>
      void serialize(Archive & ar, const unsigned int version)
        {
          ar & MeasPosX;
          ar & MeasPosY;
          ar & MeasPosZ;
          ar & Data;
          ar & Errors;
        }
      //! Add a measurement point to the model
      void AddMeasurementPoint(const double xcoord, const double ycoord,
          const double zcoord)
        {
          MeasPosX.push_back(xcoord);
          MeasPosY.push_back(ycoord);
          MeasPosZ.push_back(zcoord);
        }
      void SetMeasurementPoints(const std::vector<double> &xcoords,
          const std::vector<double> &ycoords, const std::vector<double> &zcoords)
        {
          const size_t nmeas = xcoords.size();
          if (ycoords.size() != nmeas)
            {
              throw jif3D::FatalException(
                  "Trying to set coordinates with different number of values", __FILE__,
                  __LINE__);
            }
          if (zcoords.size() != nmeas)
            {
              throw jif3D::FatalException(
                  "Trying to set coordinates with different number of values", __FILE__,
                  __LINE__);
            }
          MeasPosX = xcoords;
          MeasPosY = ycoords;
          MeasPosZ = zcoords;
        }
      //! remove all information about measurement points
      void ClearMeasurementPoints()
        {
          MeasPosX.clear();
          MeasPosY.clear();
          MeasPosZ.clear();
        }
      void CopyMeasurementConfiguration(const jif3D::GeneralData &Old)
        {
          MeasPosX = Old.MeasPosX;
          MeasPosY = Old.MeasPosY;
          MeasPosZ = Old.MeasPosZ;
        }
      //! Return the x-coordinates (Northing) of all measurement points read-only
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
      //! Return the y-coordinates (Easting)of all measurement points read-only
      const tMeasPosVec &GetMeasPosY() const
        {
          return MeasPosY;
        }
      //! Return the z-coordinates (Depth) of all measurement points read-only
      const tMeasPosVec &GetMeasPosZ() const
        {
          return MeasPosZ;
        }
      const std::vector<double> &GetData() const
        {
          return Data;
        }
      const std::vector<double> &GetErrors() const
        {
          return Errors;
        }
      void SetDataAndErrors(const std::vector<double> &D, const std::vector<double> &E)
        {
          if (D.size() != E.size())
            {
              throw jif3D::FatalException("Data and Errors do not have the same size !",
              __FILE__, __LINE__);
            }
          Data = D;
          Errors = E;
        }
      void SetErrors(const std::vector<double> &E)
        {
          if (Data.size() != E.size())
            {
              throw jif3D::FatalException("Data and Errors do not have the same size !",
              __FILE__, __LINE__);
            }
          Errors = E;
        }
      //! Write all information to a netcdf file
      virtual void ReadNetCDF(const std::string &filename) = 0;
      virtual void WriteNetCDF(const std::string &filename) const = 0;
      void WriteMeasurementPoints(const std::string &filename) const;
      GeneralData();
      virtual ~GeneralData();
              };

  } /* namespace jif3D */

#endif /* DATABASE_GENERALDATA_H_ */
